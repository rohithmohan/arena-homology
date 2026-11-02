#!/usr/bin/env python3
"""Search every Arena target against the whole PDB to measure how "seen" it was before a date.

Steps
 1. Download the PDB reference set: wwPDB pdb_seqres.txt.gz and the initial release date of
    every current entry (RCSB holdings + GraphQL, 1,000 ids per request). Cached in data/pdb/.
 2. Mask common fusion partners (MBP, BRIL/cyt b562, T4 lysozyme, GFP, DsRed/mCherry, IgG1 Fc,
    thioredoxin) in the target antigen chains with X, unless the partner covers >=80 % of the
    chain (then the partner *is* the antigen and nothing is masked).
 3. Antigen search: masked antigen chains vs all unique PDB protein sequences with
    `mmseqs easy-search -s 7.5`; peptides (<20 aa) additionally by exact substring match.
    Coverage is measured over the unmasked part of the target antigen chain.
 4. Antibody side: ANARCI (IMGT numbering + germline assignment) on every target heavy chain and
    on every antibody-like chain in PDB entries that carry a >=90 %-identity / >=80 %-coverage
    hit to a target antigen; flag entries whose heavy chain has the same CDR-H3 (equal length,
    >=80 % identity) and the same IGHV family as the target ("seen complex").

Everything is stored date-free per (target, entry) so 07 can apply any cutoff:
  data/derived/novelty/target_entry_hits.tsv
     target, entry, release, ag_id50 (max identity over antigen chains among hits covering
     >=50 % of that chain), ag_id80 (same, >=80 % coverage), complex (0/1), h3_match details
  data/derived/novelty/target_antibody_anarci.tsv  (CDR-H3, v gene, species per target)
  data/derived/novelty/fusion_masking.tsv
  data/derived/novelty/antigen_chains_masked.fasta
"""
import csv, gzip, json, os, re, shutil, subprocess, sys, time, urllib.request
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DER = os.path.join(ROOT, 'data', 'derived')
NOV = os.path.join(DER, 'novelty')
PDB = os.path.join(ROOT, 'data', 'pdb')
MMSEQS = os.environ.get('MMSEQS') or shutil.which('mmseqs') or os.path.join(ROOT, 'tools', 'mmseqs', 'bin', 'mmseqs')
THREADS = str(os.cpu_count() or 4)
FUSION = {  # UniProt accession -> label
    'P0AEX9': 'MBP (E. coli malE)', 'P0ABE7': 'BRIL / cytochrome b562', 'P00720': 'T4 lysozyme',
    'P42212': 'GFP (A. victoria)', 'Q9U6Y8': 'DsRed (mCherry parent)', 'P01857': 'human IgG1 constant (Fc)',
    'P0AA25': 'thioredoxin (E. coli trxA)'}


def get(url, data=None, tries=5):
    for a in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={'Content-Type': 'application/json'} if data else {})
            with urllib.request.urlopen(req, timeout=300) as r:
                return r.read()
        except Exception as e:
            print('retry', url[:60], e, file=sys.stderr)
            time.sleep(2 ** (a + 1))
    raise RuntimeError(url)


def run(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def read_fasta(path):
    name, seq = None, []
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name, seq = line[1:], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, ''.join(seq)


# ---------------------------------------------------------------- 1. PDB reference
def pdb_reference():
    os.makedirs(PDB, exist_ok=True)
    seqres = os.path.join(PDB, 'pdb_seqres.txt.gz')
    if not os.path.exists(seqres):
        open(seqres, 'wb').write(get('https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz'))
    dates_path = os.path.join(PDB, 'release_dates.tsv')
    if not os.path.exists(dates_path):
        ids = json.loads(get('https://data.rcsb.org/rest/v1/holdings/current/entry_ids'))
        q = 'query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id rcsb_accession_info{initial_release_date}}}'
        rows = []
        for i in range(0, len(ids), 1000):
            d = json.loads(get('https://data.rcsb.org/graphql',
                               json.dumps({'query': q, 'variables': {'ids': ids[i:i + 1000]}}).encode()))
            for e in d['data']['entries'] or []:
                if e and e['rcsb_accession_info']:
                    rows.append((e['rcsb_id'].upper(), e['rcsb_accession_info']['initial_release_date'][:10]))
        with open(dates_path + '.tmp', 'w') as f:
            f.write('entry\trelease\n' + ''.join(f'{a}\t{b}\n' for a, b in rows))
        os.replace(dates_path + '.tmp', dates_path)
    dates = {r['entry']: r['release'] for r in csv.DictReader(open(dates_path), delimiter='\t')}
    # unique protein sequences -> chains
    uniq_fa = os.path.join(PDB, 'pdb_protein_unique.fasta')
    map_path = os.path.join(PDB, 'pdb_protein_unique_map.tsv')
    if not os.path.exists(uniq_fa):
        seq2id, members = {}, defaultdict(list)
        for name, seq in read_fasta(seqres):
            if ' mol:protein ' not in name + ' ':
                continue
            if 'mol:protein' not in name:
                continue
            ident = name.split()[0]
            entry, chain = ident.rsplit('_', 1)
            if seq not in seq2id:
                seq2id[seq] = f'U{len(seq2id)}'
            members[seq2id[seq]].append(f'{entry.upper()}_{chain}')
        with open(uniq_fa, 'w') as f:
            for s, u in seq2id.items():
                f.write(f'>{u}\n{s}\n')
        with open(map_path, 'w') as f:
            for u, m in members.items():
                f.write(f'{u}\t{",".join(m)}\n')
    useq = dict(read_fasta(uniq_fa))
    umap = {}
    for line in open(map_path):
        u, m = line.rstrip('\n').split('\t')
        umap[u] = m.split(',')
    return dates, useq, umap


# ---------------------------------------------------------------- 2. fusion masking
def fusion_mask(chains):
    """chains: list of (rec_id, seq). Returns {rec_id: masked_seq}, report rows."""
    fdir = os.path.join(NOV, 'fusion')
    os.makedirs(fdir, exist_ok=True)
    lib = os.path.join(fdir, 'fusion_partners.fasta')
    if not os.path.exists(lib):
        with open(lib, 'w') as f:
            for acc, lab in FUSION.items():
                seq = ''.join(get(f'https://rest.uniprot.org/uniprotkb/{acc}.fasta').decode().splitlines()[1:])
                f.write(f'>{acc}\n{seq}\n')
    q = os.path.join(fdir, 'query.fasta')
    with open(q, 'w') as f:
        for r, s in chains:
            f.write(f'>{r}\n{s}\n')
    out = os.path.join(fdir, 'hits.m8')
    run([MMSEQS, 'easy-search', q, lib, out, os.path.join(fdir, 'tmp'), '-s', '7.5', '-e', '1e-5',
         '--threads', THREADS, '-v', '1', '--format-output', 'query,target,fident,alnlen,qstart,qend,qlen'])
    shutil.rmtree(os.path.join(fdir, 'tmp'), ignore_errors=True)
    regions = defaultdict(list)
    for qid, tid, fid, alen, qs, qe, ql in csv.reader(open(out), delimiter='\t'):
        if float(fid) >= 0.5 and int(alen) >= 50:
            regions[qid].append((int(qs), int(qe), tid, float(fid)))
    masked, rep = {}, []
    for r, s in chains:
        cov = [False] * len(s)
        for qs, qe, _, _ in regions.get(r, []):
            for i in range(qs - 1, qe):
                cov[i] = True
        frac = sum(cov) / len(s) if s else 0
        if regions.get(r) and frac < 0.8:
            masked[r] = ''.join('X' if c else a for a, c in zip(s, cov))
            action = 'masked'
        else:
            masked[r] = s
            action = 'partner_is_antigen_not_masked' if regions.get(r) else 'none'
        if regions.get(r):
            parts = sorted({FUSION[t] for _, _, t, _ in regions[r]})
            rep.append(dict(chain=r, partners='; '.join(parts), masked_fraction=round(frac, 3), action=action,
                            length=len(s)))
    return masked, rep


# ---------------------------------------------------------------- 4. ANARCI helpers
def anarci_heavy(seqs, ncpu, species=('human', 'mouse'), on_chunk=None):
    """seqs: list of (id, seq). Returns {id: dict(h3, vgene, species)} for chains with an H domain.
    species: germline reference set (ANARCI's default is human + mouse). on_chunk(done_ids, out) is
    called after each chunk so callers can checkpoint."""
    from anarci import anarci as run_anarci
    out = {}
    B = 1000
    for i in range(0, len(seqs), B):
        chunk = seqs[i:i + B]
        numbering, details, _ = run_anarci(chunk, scheme='imgt', assign_germline=True, ncpu=ncpu,
                                           allowed_species=list(species))
        for (sid, _), num, det in zip(chunk, numbering, details):
            if not det:
                continue
            for k, d in enumerate(det):
                if d['chain_type'] != 'H':
                    continue
                h3 = ''.join(a for (pos, ins), a in num[k][0] if 105 <= pos <= 117 and a != '-')
                gl = d.get('germlines') or {}
                v = gl.get('v_gene') or [(None, None), 0]
                (gspecies, vgene), vid = v[0], v[1]
                out[sid] = dict(h3=h3, vgene=vgene or '', species=gspecies or '', v_identity=vid or 0)
                break
        if on_chunk:
            on_chunk([sid for sid, _ in chunk], out)
        print(f'  ANARCI {min(i + B, len(seqs))}/{len(seqs)}', file=sys.stderr)
    return out


def vfamily(vgene):
    m = re.match(r'(IGHV\d+)', vgene or '')
    return m.group(1) if m else ''


def h3_identity(a, b):
    if not a or len(a) != len(b):
        return 0.0
    return sum(x == y for x, y in zip(a, b)) / len(a)


# ---------------------------------------------------------------- main
def main():
    os.makedirs(NOV, exist_ok=True)
    targets = json.load(open(os.path.join(DER, 'targets.json')))
    keys = [t['pdb'] + '_' + t['abId'].replace('sabdab2_', '') for t in targets]
    own = {k: t['pdb'].upper() for k, t in zip(keys, targets)}
    rows = list(csv.DictReader(open(os.path.join(DER, 'chain_sequences.tsv')), delimiter='\t'))

    print('PDB reference ...', file=sys.stderr)
    dates, useq, umap = pdb_reference()
    print(f'  {len(dates)} entries with release dates, {len(useq)} unique protein sequences', file=sys.stderr)

    # antigen chains, one record per (target, distinct sequence)
    ag = []
    seen = set()
    for r in rows:
        if r['role'] == 'AG' and (r['target'], r['seq']) not in seen:
            seen.add((r['target'], r['seq']))
            ag.append((f"{r['target']}|{r['chain']}", r['seq']))
    masked, frep = fusion_mask(ag)
    with open(os.path.join(NOV, 'fusion_masking.tsv'), 'w') as f:
        f.write('chain\tpartners\tmasked_fraction\taction\tlength\n')
        for d in frep:
            f.write('\t'.join(str(d[c]) for c in ('chain', 'partners', 'masked_fraction', 'action', 'length')) + '\n')
    agfa = os.path.join(NOV, 'antigen_chains_masked.fasta')
    with open(agfa, 'w') as f:
        for r, s in ag:
            f.write(f'>{r}\n{masked[r]}\n')
    unmasked_len = {r: sum(c != 'X' for c in masked[r]) for r, _ in ag}

    # ---- 3. antigen search
    hits_path = os.path.join(NOV, 'antigen_vs_pdb.m8')
    if not os.path.exists(hits_path):
        print('antigen search ...', file=sys.stderr)
        tmp = os.path.join(NOV, 'tmp_ag')
        run([MMSEQS, 'easy-search', agfa, os.path.join(PDB, 'pdb_protein_unique.fasta'), hits_path + '.tmp', tmp,
             '-s', '7.5', '-e', '1e-3', '--max-seqs', '20000', '--threads', THREADS, '-v', '1',
             '--format-output', 'query,target,fident,alnlen,qstart,qend,qlen,nident'])
        shutil.rmtree(tmp, ignore_errors=True)
        os.replace(hits_path + '.tmp', hits_path)
    # chain-level best (identity, coverage) per unique PDB sequence
    best = defaultdict(dict)  # query -> {U: (fid, cov)}
    for qid, uid, fid, alen, qs, qe, ql, nid in csv.reader(open(hits_path), delimiter='\t'):
        s = masked[qid]
        cov = sum(1 for c in s[int(qs) - 1:int(qe)] if c != 'X') / max(1, unmasked_len[qid])
        fid = float(fid)
        prev = best[qid].get(uid)
        if prev is None or (fid, cov) > prev:
            best[qid][uid] = (fid, cov)
    # peptides: exact substring match (MMseqs k-mer prefilter is unreliable below ~20 aa)
    peptides = [(r, s) for r, s in ag if len(s) < 20]
    if peptides:
        for uid, s in useq.items():
            for r, p in peptides:
                if p in s:
                    best[r][uid] = (1.0, 1.0)

    # per (target, entry): max identity over antigen chains among hits with cov >= .5 / >= .8
    te50, te80 = defaultdict(float), defaultdict(float)
    ag80_entries = defaultdict(set)
    for qid, d in best.items():
        tgt = qid.split('|')[0]
        for uid, (fid, cov) in d.items():
            if cov < 0.5:
                continue
            for ch in umap[uid]:
                e = ch.rsplit('_', 1)[0]
                if e == own[tgt]:
                    continue
                te50[(tgt, e)] = max(te50[(tgt, e)], fid)
                if cov >= 0.8:
                    te80[(tgt, e)] = max(te80[(tgt, e)], fid)
                    if fid >= 0.9:
                        ag80_entries[tgt].add(e)
    print(f'  {len(te50)} (target, entry) antigen hits at >=50 % coverage', file=sys.stderr)

    # ---- 4. antibodies
    tab_path = os.path.join(NOV, 'target_antibody_anarci.tsv')
    hrows = [r for r in rows if r['role'] == 'H']
    tab = anarci_heavy([(r['target'], r['vdomain']) for r in hrows], ncpu=int(THREADS))
    tab_all = anarci_heavy([(r['target'], r['vdomain']) for r in hrows], ncpu=int(THREADS),
                           species=('human', 'mouse', 'rat', 'rabbit', 'rhesus', 'pig', 'alpaca', 'cow'))
    for path, tb in ((tab_path, tab), (tab_path.replace('.tsv', '_allspecies.tsv'), tab_all)):
        with open(path, 'w') as f:
            f.write('target\tformat\th3\th3_len\tvgene\tspecies\tv_identity\n')
            for r in hrows:
                a = tb.get(r['target'], {})
                f.write(f"{r['target']}\t{r['format']}\t{a.get('h3','')}\t{len(a.get('h3',''))}\t{a.get('vgene','')}"
                        f"\t{a.get('species','')}\t{a.get('v_identity','')}\n")

    # candidate entries for complex check: >=90 % / >=80 % antigen hit (own entry excluded)
    cand_entries = set().union(*ag80_entries.values()) if ag80_entries else set()
    print(f'  {len(cand_entries)} candidate entries for the complex check', file=sys.stderr)
    entry_chains = defaultdict(list)
    for u, chs in umap.items():
        for ch in chs:
            e = ch.rsplit('_', 1)[0]
            if e in cand_entries:
                entry_chains[e].append(u)
    cand_u = sorted({u for us in entry_chains.values() for u in us if 90 <= len(useq[u]) <= 600})
    # antibody-like prefilter: MMseqs search against the target V domains (VH, VHH, VL)
    vfa = os.path.join(NOV, 'target_vdomains.fasta')
    with open(vfa, 'w') as f:
        for r in rows:
            if r['role'] in ('H', 'L'):
                f.write(f">{r['target']}|{r['role']}\n{r['vdomain']}\n")
    cfa = os.path.join(NOV, 'candidate_chains.fasta')
    with open(cfa, 'w') as f:
        for u in cand_u:
            f.write(f'>{u}\n{useq[u]}\n')
    abhits = os.path.join(NOV, 'candidate_vs_vdomains.m8')
    run([MMSEQS, 'easy-search', cfa, vfa, abhits, os.path.join(NOV, 'tmp_ab'), '-s', '5', '-e', '1e-3',
         '--threads', THREADS, '-v', '1', '--format-output', 'query,target,fident'])
    shutil.rmtree(os.path.join(NOV, 'tmp_ab'), ignore_errors=True)
    ab_u = sorted({r[0] for r in csv.reader(open(abhits), delimiter='\t')})
    print(f'  {len(ab_u)} antibody-like unique chains -> ANARCI', file=sys.stderr)
    cache = os.path.join(NOV, 'pdb_heavy_anarci.tsv')
    known = {}
    if os.path.exists(cache):
        for r in csv.DictReader(open(cache), delimiter='\t'):
            known[r['uid']] = r
    todo = [(u, useq[u]) for u in ab_u if u not in known]

    def checkpoint(done, res):
        for u in done:
            a = res.get(u)
            known[u] = dict(uid=u, h3=a['h3'] if a else '', vgene=a['vgene'] if a else '',
                            species=a['species'] if a else '', is_heavy='1' if a else '0')
        with open(cache + '.tmp', 'w') as f:
            f.write('uid\th3\tvgene\tspecies\tis_heavy\n')
            for u, r in known.items():
                f.write(f"{u}\t{r['h3']}\t{r['vgene']}\t{r['species']}\t{r['is_heavy']}\n")
        os.replace(cache + '.tmp', cache)
    if todo:
        anarci_heavy(todo, ncpu=int(THREADS), on_chunk=checkpoint)

    complex_hit = {}
    for tgt, ents in ag80_entries.items():
        t = tab.get(tgt)
        if not t or not t['h3']:
            continue
        fam = vfamily(t['vgene'])
        for e in ents:
            bestm = None
            for u in entry_chains.get(e, []):
                k = known.get(u)
                if not k or k.get('is_heavy') != '1' or not k['h3']:
                    continue
                idt = h3_identity(t['h3'], k['h3'])
                if idt >= 0.8 and fam and vfamily(k['vgene']) == fam:
                    if bestm is None or idt > bestm[0]:
                        bestm = (idt, k['h3'], k['vgene'])
            if bestm:
                complex_hit[(tgt, e)] = bestm

    out = os.path.join(NOV, 'target_entry_hits.tsv')
    with open(out, 'w') as f:
        f.write('target\tentry\trelease\tag_id50\tag_id80\tcomplex\th3_identity\th3_match\tvgene_match\n')
        n_missing = 0
        for (tgt, e), v in sorted(te50.items()):
            rel = dates.get(e)
            if rel is None:
                n_missing += 1
                continue
            c = complex_hit.get((tgt, e))
            f.write(f"{tgt}\t{e}\t{rel}\t{v:.4f}\t{te80.get((tgt, e), 0):.4f}\t{1 if c else 0}\t"
                    f"{c[0]:.3f}\t{c[1]}\t{c[2]}\n" if c else
                    f"{tgt}\t{e}\t{rel}\t{v:.4f}\t{te80.get((tgt, e), 0):.4f}\t0\t\t\t\n")
    print(f'wrote {out}; entries without a release date (obsolete/unreleased): {n_missing}', file=sys.stderr)
    print(f'seen-complex (target, entry) pairs, any date: {len(complex_hit)}; '
          f'targets with any: {len({t for t, _ in complex_hit})}', file=sys.stderr)


if __name__ == '__main__':
    main()
