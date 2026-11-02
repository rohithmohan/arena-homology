#!/usr/bin/env python3
"""Fetch chain sequences for every Arena target from RCSB (data.rcsb.org GraphQL).

Inputs : data/derived/targets.json
Outputs: data/raw/rcsb/*.json            (raw GraphQL responses, cached)
         data/derived/chain_sequences.tsv  (target, role, chain, auth/label ids, seq, ...)
         data/derived/all_chains.fasta
         data/derived/antigen_chains.fasta, heavy_vdomain.fasta
         data/derived/fetch_report.txt

Chain IDs in the Arena (from SAbDab) are matched against *auth* asym IDs
(auth_asym_ids); we verify this by checking that the auth-mapped heavy/light chains
look like antibody V domains, and record how often a *label* interpretation would
fail or disagree.
"""
import json, os, re, sys, time, urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DER = os.path.join(ROOT, 'data', 'derived')
RAW = os.path.join(ROOT, 'data', 'raw', 'rcsb')
os.makedirs(RAW, exist_ok=True)

QUERY = """query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id
 polymer_entities{rcsb_id
  rcsb_polymer_entity_container_identifiers{auth_asym_ids asym_ids}
  entity_poly{pdbx_seq_one_letter_code_can type}
  rcsb_polymer_entity{pdbx_description}}}}"""


def gql(ids):
    body = json.dumps({'query': QUERY, 'variables': {'ids': ids}}).encode()
    for attempt in range(5):
        try:
            req = urllib.request.Request('https://data.rcsb.org/graphql', data=body,
                                         headers={'Content-Type': 'application/json'})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception as e:  # network retry
            print('retry', attempt, e, file=sys.stderr)
            time.sleep(2 ** (attempt + 1))
    raise RuntimeError('RCSB query failed')


# Heavy-chain FR4 "WGxG" (J-region) motif; light chain "FGxG". Used to (a) sanity-check
# chain assignment and (b) cut the heavy chain to its variable domain so that
# Fab constant domains (CH1) don't dominate VH clustering.
H_FR4 = re.compile(r'W[GA].G')
L_FR4 = re.compile(r'FG.G')
H_FR1 = re.compile(r'C..[ST]')       # Cys22 region, loose
VH_LIKE = re.compile(r'W[VFYI]RQ')   # FR2 start of VH/VHH


def vdomain(seq, kind):
    """Return (variable-domain sequence, method). Cut 11 residues after FR4 W/F
    (WGQGTxVTVSS / FGxGTKLEIK), searching only the first 150 residues."""
    pat = H_FR4 if kind == 'H' else L_FR4
    best = None
    for m in pat.finditer(seq[:160]):
        if m.start() >= 80:
            best = m
            break
    if best:
        return seq[:best.start() + 11], 'fr4_motif'
    return seq[:125], 'fallback_first125'


def main():
    targets = json.load(open(os.path.join(DER, 'targets.json')))
    pdbs = sorted({t['pdb'] for t in targets})
    entries = {}
    for i in range(0, len(pdbs), 40):
        chunk = pdbs[i:i + 40]
        cache = os.path.join(RAW, f'chunk_{i:04d}.json')
        if os.path.exists(cache):
            data = json.load(open(cache))
        else:
            data = gql(chunk)
            json.dump(data, open(cache, 'w'))
        for e in (data.get('data') or {}).get('entries') or []:
            if e:
                entries[e['rcsb_id'].upper()] = e
    report = []
    missing_entries = [p for p in pdbs if p.upper() not in entries]
    report.append(f'PDB entries requested={len(pdbs)} returned={len(entries)} missing={missing_entries}')

    def chain_maps(e):
        auth, label = {}, {}
        for pe in e['polymer_entities']:
            ids = pe['rcsb_polymer_entity_container_identifiers']
            info = dict(entity=pe['rcsb_id'], seq=(pe['entity_poly']['pdbx_seq_one_letter_code_can'] or '').replace('\n', ''),
                        ptype=pe['entity_poly']['type'], desc=pe['rcsb_polymer_entity']['pdbx_description'])
            for a in ids.get('auth_asym_ids') or []:
                auth[a] = info
            for l in ids.get('asym_ids') or []:
                label[l] = info
        return auth, label

    rows, problems = [], []
    label_disagree = label_missing = 0
    for t in targets:
        key = t['pdb'] + '_' + t['abId'].replace('sabdab2_', '')
        e = entries.get(t['pdb'].upper())
        if e is None:
            problems.append(f'{key}: entry not found at RCSB')
            continue
        auth, label = chain_maps(e)
        roles = [('H', t['Hchain'])] + ([('L', t['Lchain'])] if t['Lchain'] else []) + [('AG', c) for c in t['agChains']]
        for role, ch in roles:
            info = auth.get(ch)
            if info is None:
                problems.append(f'{key}: {role} chain {ch} not found among auth ids {sorted(auth)}')
                continue
            linfo = label.get(ch)
            if linfo is None:
                label_missing += 1
            elif linfo['entity'] != info['entity']:
                label_disagree += 1
            seq = info['seq']
            row = dict(target=key, pdb=t['pdb'], abId=t['abId'], format=t['format'], role=role, chain=ch,
                       entity=info['entity'], ptype=info['ptype'], desc=info['desc'], length=len(seq), seq=seq)
            if role in ('H', 'L'):
                vd, how = vdomain(seq, role)
                row['vdomain'] = vd
                row['vdomain_method'] = how
                row['vh_like'] = bool(VH_LIKE.search(seq[:60])) if role == 'H' else bool(L_FR4.search(seq[:130]))
            rows.append(row)
    report.append(f'chain lookups ok={len(rows)} problems={len(problems)}')
    # Refine V-domain boundaries with ANARCI (IMGT HMMs) when available.
    try:
        from anarci import anarci as run_anarci
        ab = [r for r in rows if r['role'] in ('H', 'L')]
        _, details, _ = run_anarci([(str(i), r['seq']) for i, r in enumerate(ab)], scheme='imgt',
                                   assign_germline=False, ncpu=4)
        for r, d in zip(ab, details):
            if not d:
                continue
            want = ('H',) if r['role'] == 'H' else ('K', 'L')
            dom = next((x for x in d if x['chain_type'] in want), d[0])
            r['vdomain'] = r['seq'][dom['query_start']:dom['query_end'] + 1]
            r['vdomain_method'] = f"anarci_{dom['chain_type']}"
            r['vh_like'] = dom['chain_type'] in want
    except Exception as ex:  # noqa
        report.append(f'ANARCI unavailable ({ex}); using FR4-motif heuristic for V-domain trimming')
    report.append(f'label-ID interpretation: chain id absent as label id in {label_missing} lookups, '
                  f'maps to a different entity in {label_disagree} lookups (auth IDs used throughout)')
    hrows = [r for r in rows if r['role'] == 'H']
    from collections import Counter
    report.append(f'heavy chains recognised as VH (ANARCI H domain, or VH motif if no ANARCI): '
                  f'{sum(r["vh_like"] for r in hrows)}/{len(hrows)}; trim methods: '
                  f'{dict(Counter(r["vdomain_method"] for r in hrows))}; '
                  f'VH length range {min(len(r["vdomain"]) for r in hrows)}-{max(len(r["vdomain"]) for r in hrows)}')
    report.append('heavy chains not recognised as VH: ' + str([f"{r['target']}:{r['desc']}:{r['vdomain_method']}"
                  for r in hrows if not r['vh_like']]))
    lrows = [r for r in rows if r['role'] == 'L']
    report.append(f'light chains recognised as VL: {sum(r["vh_like"] for r in lrows)}/{len(lrows)}')
    nonprot = [f"{r['target']}:{r['role']}:{r['chain']}:{r['ptype']}" for r in rows if r['ptype'] != 'polypeptide(L)']
    report.append(f'non-protein chains: {nonprot}')
    # antigen chains that look like antibodies (e.g. anti-idiotype targets) - informational
    agab = [f"{r['target']}:{r['chain']}:{r['desc']}" for r in rows if r['role'] == 'AG' and VH_LIKE.search(r['seq'][:60])]
    report.append(f'antigen chains with VH-like motif (antibody-as-antigen?): {agab}')
    report += ['PROBLEM ' + p for p in problems]

    cols = ['target', 'pdb', 'abId', 'format', 'role', 'chain', 'entity', 'ptype', 'desc', 'length',
            'vdomain_method', 'vh_like', 'seq', 'vdomain']
    with open(os.path.join(DER, 'chain_sequences.tsv'), 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
    with open(os.path.join(DER, 'all_chains.fasta'), 'w') as f:
        for r in rows:
            f.write(f">{r['target']}|{r['role']}|{r['chain']} {r['desc']}\n{r['seq']}\n")
    # clustering inputs: one record per unique (target, sequence) to avoid homo-oligomer copies
    with open(os.path.join(DER, 'antigen_chains.fasta'), 'w') as f:
        seen = set()
        for r in rows:
            if r['role'] == 'AG' and (r['target'], r['seq']) not in seen and r['seq']:
                seen.add((r['target'], r['seq']))
                f.write(f">{r['target']}|AG|{r['chain']}\n{r['seq']}\n")
    with open(os.path.join(DER, 'heavy_vdomain.fasta'), 'w') as f:
        for r in hrows:
            f.write(f">{r['target']}|H|{r['chain']}\n{r['vdomain']}\n")
    open(os.path.join(DER, 'fetch_report.txt'), 'w').write('\n'.join(report) + '\n')
    print('\n'.join(report))


if __name__ == '__main__':
    main()
