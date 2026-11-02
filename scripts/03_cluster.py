#!/usr/bin/env python3
"""Cluster Arena targets by antigen (and antibody VH) sequence homology with MMseqs2.

Primary (as specified): mmseqs easy-cluster on antigen chains at 30/50/90 % identity,
coverage 0.8 (cov-mode 0 = bidirectional), then union-find over targets: two targets are
in one group if any of their antigen chains share an MMseqs cluster.

Antibody check: VH / VHH variable domains (ANARCI-trimmed) clustered at 90 % / cov 0.8.
Joint group = (antigen-50 group, VH-90 cluster) pair.

Extra sensitivity (not in the spec, reported separately): antigen all-vs-all search with
edges where identity >= 50 % and the *shorter* chain is >= 80 % covered, which links
sub-constructs (e.g. RBD vs full spike ectodomain) that bidirectional coverage separates.

Outputs: data/derived/groups.csv, data/derived/cluster_report.txt, data/derived/mmseqs/*
"""
import csv, json, os, shutil, subprocess, sys
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DER = os.path.join(ROOT, 'data', 'derived')
MM = os.path.join(DER, 'mmseqs')
MMSEQS = os.environ.get('MMSEQS') or shutil.which('mmseqs') or os.path.join(ROOT, 'tools', 'mmseqs', 'bin', 'mmseqs')
COMMON = ['-c', '0.8', '--cov-mode', '0', '-s', '7.5', '--threads', '4', '-v', '1']


def run(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def easy_cluster(fasta, prefix, ident):
    tmp = os.path.join(MM, 'tmp_' + os.path.basename(prefix))
    run([MMSEQS, 'easy-cluster', fasta, prefix, tmp, '--min-seq-id', str(ident)] + COMMON)
    shutil.rmtree(tmp, ignore_errors=True)
    clus = {}
    for rep, mem in csv.reader(open(prefix + '_cluster.tsv'), delimiter='\t'):
        clus[mem] = rep
    return clus


class UF:
    def __init__(self, items):
        self.p = {i: i for i in items}

    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[max(ra, rb)] = min(ra, rb)


def components(targets, links):
    """links: iterable of sets of targets that must be merged."""
    uf = UF(targets)
    for s in links:
        s = list(s)
        for o in s[1:]:
            uf.union(s[0], o)
    return {t: uf.find(t) for t in targets}


def relabel(assign, prefix, order):
    """Name groups PREFIX_001.. by decreasing size, ties broken by first target order."""
    members = defaultdict(list)
    for t in order:
        members[assign[t]].append(t)
    ranked = sorted(members.values(), key=lambda m: (-len(m), order.index(m[0])))
    out = {}
    for i, m in enumerate(ranked, 1):
        for t in m:
            out[t] = f'{prefix}_{i:03d}'
    return out


def main():
    os.makedirs(MM, exist_ok=True)
    targets = json.load(open(os.path.join(DER, 'targets.json')))
    order = [t['pdb'] + '_' + t['abId'].replace('sabdab2_', '') for t in targets]
    meta = {k: t for k, t in zip(order, targets)}
    ag_fa = os.path.join(DER, 'antigen_chains.fasta')
    vh_fa = os.path.join(DER, 'heavy_vdomain.fasta')
    tgt_of = lambda rec: rec.split('|')[0]

    groups = {t: {} for t in order}
    for ident in (0.3, 0.5, 0.9):
        tag = f'ag{int(ident*100)}'
        clus = easy_cluster(ag_fa, os.path.join(MM, tag), ident)
        by_cluster = defaultdict(set)
        for mem, rep in clus.items():
            by_cluster[rep].add(tgt_of(mem))
        missing = set(order) - {tgt_of(m) for m in clus}
        assert not missing, missing
        assign = relabel(components(order, by_cluster.values()), tag.upper(), order)
        for t in order:
            groups[t][tag] = assign[t]
        groups_n_raw = len(set(clus.values()))
        print(f'{tag}: {groups_n_raw} mmseqs clusters -> {len(set(assign.values()))} target groups')

    vh = easy_cluster(vh_fa, os.path.join(MM, 'vh90'), 0.9)
    vh_assign = relabel({tgt_of(m): rep for m, rep in vh.items()}, 'VH90', order)
    for t in order:
        groups[t]['vh90'] = vh_assign[t]
    joint_raw = {t: groups[t]['ag50'] + '|' + groups[t]['vh90'] for t in order}
    joint = relabel(joint_raw, 'JOINT', order)
    for t in order:
        groups[t]['joint'] = joint[t]

    # sensitivity: shorter-chain coverage graph at 50 % identity
    hits = os.path.join(MM, 'ag_allvsall.m8')
    tmp = os.path.join(MM, 'tmp_search')
    run([MMSEQS, 'easy-search', ag_fa, ag_fa, hits, tmp, '-s', '7.5', '--threads', '4', '-v', '1',
         '--format-output', 'query,target,fident,qcov,tcov,evalue', '-e', '1e-3'])
    shutil.rmtree(tmp, ignore_errors=True)
    links = []
    for q, s, fid, qc, tc, ev in csv.reader(open(hits), delimiter='\t'):
        if float(fid) >= 0.5 and max(float(qc), float(tc)) >= 0.8:
            links.append({tgt_of(q), tgt_of(s)})
    sc = relabel(components(order, links), 'AG50S', order)
    for t in order:
        groups[t]['ag50_shortcov'] = sc[t]

    keys = ['ag30', 'ag50', 'ag90', 'vh90', 'joint', 'ag50_shortcov']
    sizes = {k: Counter(groups[t][k] for t in order) for k in keys}
    with open(os.path.join(DER, 'groups.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['target'] + [x for k in keys for x in (k, k + '_size')])
        for t in order:
            w.writerow([t] + [x for k in keys for x in (groups[t][k], sizes[k][groups[t][k]])])

    rep = []
    for k in keys:
        c = sizes[k]
        dist = Counter(c.values())
        n_eff = len(c)
        rep.append(f'== {k}: {n_eff} groups over {len(order)} targets; singletons={dist.get(1,0)}; '
                   f'targets in groups of size>1: {sum(v for v in c.values() if v > 1)}')
        rep.append('   size distribution (size:count): ' + ', '.join(f'{s}:{n}' for s, n in sorted(dist.items())))
        for g, n in c.most_common(15):
            names = Counter(meta[t]['antigen'] for t in order if groups[t][k] == g)
            fmt = Counter(meta[t]['format'] for t in order if groups[t][k] == g)
            nm = '; '.join(f'{a} x{m}' if m > 1 else a for a, m in names.most_common(6))
            more = f' (+{len(names)-6} more names)' if len(names) > 6 else ''
            rep.append(f'   {g} n={n} [{dict(fmt)}] {nm}{more}')
    open(os.path.join(DER, 'cluster_report.txt'), 'w').write('\n'.join(rep) + '\n')
    print('\n'.join(rep))


if __name__ == '__main__':
    main()
