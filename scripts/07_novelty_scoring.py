#!/usr/bin/env python3
"""Training-set novelty analysis, clonotype-based antibody redundancy, and fusion-masked clustering.

Inputs : models_cutoffs.csv, data/derived/novelty/* (from 06), data/derived/groups.csv, scores.
Outputs (repo root):
  novelty_per_target_model.csv   target x model: antigen max identity before cutoff, stratum, seen complex
  novelty_strata.csv             per model x stratum: n, mean DockQ (+CI), acceptable %
  novel_only_leaderboard.csv     per-model novel set (<50 % antigen identity) and common-subset versions
  protenix_natural_experiment.csv
  mol1_consistency_check.csv
  antibody_clonotypes.csv        per target clonotype + joint (ag50 & clonotype) groups
  clonotype_rescoring.csv
  fusion_masking_effect.csv      groupings before/after masking fusion partners
  data/derived/results_novelty.json

Conventions
  * "Before cutoff" = PDB initial release date strictly earlier than the model's cutoff date.
  * The target's own entry is excluded from the PDB search (06). The one exception is a model
    whose cutoff is *after* the target's release (only Protenix-v1-20250630): its training set
    contained the target itself, so the target is counted as seen (identity 1.0, seen complex)
    and flagged ``own_entry_in_training``.
  * Multi-chain antigens: identity is the max over the target's antigen chains.
  * Models with unknown cutoff are excluded from cutoff-specific results.
  * CIs: paired cluster bootstrap over antigen-50 % groups (B=2000), same resample for all models.
"""
import csv, importlib.util, json, os
from collections import Counter, defaultdict
import numpy as np
from scipy.stats import kendalltau

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DER = os.path.join(ROOT, 'data', 'derived')
NOV = os.path.join(DER, 'novelty')
B = 2000
SEED = 20261103
STRATA = [('<30%', 0, 0.3), ('30-50%', 0.3, 0.5), ('50-90%', 0.5, 0.9), ('>=90%', 0.9, 9)]


def load_module(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, 'scripts', fname))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


R4 = load_module('rescore', '04_rescore.py')
C3 = load_module('cluster', '03_cluster.py')


def wcsv(name, header, rows):
    with open(os.path.join(ROOT, name), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def main():
    rng = np.random.default_rng(SEED)
    targets, scores, models, keys, S, groups, ag_desc = R4.load()
    n, M = S.shape
    mids = [m['id'] for m in models]
    names = [m['name'] for m in models]
    tidx = {k: i for i, k in enumerate(keys)}
    trel = np.array([t['release'] for t in targets])
    fmt = np.array([t['format'] for t in targets])
    ag50 = [groups[k]['ag50'] for k in keys]
    gidx, G = R4.group_index(ag50)
    C = rng.multinomial(G, np.full(G, 1.0 / G), size=B).astype(float)  # one resample matrix for everything

    def boot_subset_means(mask):
        """mask: n x M bool (target in model's subset). Returns point (M,), boot (B,M)."""
        mk = mask.astype(float)
        gs = np.zeros((G, M)); gn = np.zeros((G, M))
        np.add.at(gs, gidx, S * mk)
        np.add.at(gn, gidx, mk)
        with np.errstate(invalid='ignore', divide='ignore'):
            point = gs.sum(0) / gn.sum(0)
            boot = (C @ gs) / (C @ gn)
        return point, boot

    def boot_vec(v, mask=None):
        """Cluster-bootstrap mean of a per-target vector v over targets in mask."""
        m = np.ones(n, bool) if mask is None else mask
        gs = np.zeros(G); gn = np.zeros(G)
        np.add.at(gs, gidx[m], v[m]); np.add.at(gn, gidx[m], 1)
        with np.errstate(invalid='ignore', divide='ignore'):
            b = (C @ gs) / (C @ gn)
        return v[m].mean() if m.any() else np.nan, np.nanpercentile(b, 2.5), np.nanpercentile(b, 97.5)

    # ---------------- cutoffs
    cut = {}
    cut_rows = list(csv.DictReader(open(os.path.join(ROOT, 'models_cutoffs.csv'))))
    for r in cut_rows:
        if r['cutoff_date'] and r['cutoff_date'] != 'unknown':
            cut[r['model_id']] = r['cutoff_date']
    known = [j for j in range(M) if mids[j] in cut]
    print('models with known cutoff:', [(names[j], cut[mids[j]]) for j in known])
    print('unknown:', [names[j] for j in range(M) if j not in known])

    # ---------------- per (target, entry) hits
    hits = defaultdict(list)  # target -> list of (release, id50, id80, complex)
    for r in csv.DictReader(open(os.path.join(NOV, 'target_entry_hits.tsv')), delimiter='\t'):
        hits[r['target']].append((r['release'], float(r['ag_id50']), float(r['ag_id80']), r['complex'] == '1'))

    def novelty_at(k, date):
        """max antigen identity (>=50 % cov) and seen-complex flag using entries released before date."""
        ident, cplx = 0.0, False
        for rel, i50, i80, c in hits.get(k, []):
            if rel < date:
                ident = max(ident, i50)
                cplx = cplx or c
        return ident, cplx

    ID = np.full((n, M), np.nan)
    CX = np.zeros((n, M), bool)
    OWN = np.zeros((n, M), bool)
    for j in known:
        for i, k in enumerate(keys):
            ident, cplx = novelty_at(k, cut[mids[j]])
            if trel[i] < cut[mids[j]]:
                OWN[i, j] = True
                ident, cplx = 1.0, True
            ID[i, j], CX[i, j] = ident, cplx

    def stratum(x):
        for lab, lo, hi in STRATA:
            if lo <= x < hi:
                return lab
        return ''

    rows = []
    for i, k in enumerate(keys):
        for j in known:
            rows.append([k, fmt[i], trel[i], groups[k]['ag50'], mids[j], names[j], cut[mids[j]], round(ID[i, j], 4),
                         stratum(ID[i, j]), int(CX[i, j]), int(OWN[i, j]), S[i, j]])
    wcsv('novelty_per_target_model.csv', ['target', 'format', 'target_release', 'ag50_group', 'model_id', 'model',
                                          'cutoff', 'antigen_max_identity_before_cutoff', 'stratum',
                                          'seen_complex', 'own_entry_in_training', 'dockq'], rows)

    res = dict(cutoffs=cut, known=[mids[j] for j in known], n_targets=n, B=B, seed=SEED)
    acc = (S >= 0.23).astype(float)

    # ---------------- strata
    srows, res['strata'] = [], {}
    for j in known:
        d = {}
        for lab, lo, hi in STRATA:
            m = (ID[:, j] >= lo) & (ID[:, j] < hi)
            mean, lo_, hi_ = boot_vec(S[:, j], m) if m.any() else (np.nan, np.nan, np.nan)
            d[lab] = dict(n=int(m.sum()), mean=float(mean), lo=float(lo_), hi=float(hi_),
                          acc=float(acc[m, j].mean()) if m.any() else np.nan)
            srows.append([mids[j], names[j], cut[mids[j]], 'antigen ' + lab, int(m.sum()), round(mean, 4),
                          round(lo_, 4), round(hi_, 4), round(d[lab]['acc'] * 100, 1) if m.any() else ''])
        for lab, m in [('seen complex', CX[:, j]), ('no seen complex', ~CX[:, j])]:
            mean, lo_, hi_ = boot_vec(S[:, j], m) if m.any() else (np.nan, np.nan, np.nan)
            d[lab] = dict(n=int(m.sum()), mean=float(mean), lo=float(lo_), hi=float(hi_),
                          acc=float(acc[m, j].mean()) if m.any() else np.nan)
            srows.append([mids[j], names[j], cut[mids[j]], lab, int(m.sum()), round(mean, 4), round(lo_, 4),
                          round(hi_, 4), round(d[lab]['acc'] * 100, 1) if m.any() else ''])
        d['own_entry_in_training'] = int(OWN[:, j].sum())
        d['all'] = dict(mean=float(S[:, j].mean()), acc=float(acc[:, j].mean()))
        res['strata'][mids[j]] = d
    wcsv('novelty_strata.csv', ['model_id', 'model', 'cutoff', 'stratum', 'n_targets', 'mean_dockq', 'ci_lo', 'ci_hi',
                                'acceptable_pct'], srows)

    # ---------------- novel-only leaderboards
    K = np.array(known)

    def leaderboard(mask, label, K=K):
        K = np.array(K)
        point, boot = boot_subset_means(mask)
        pk, bk = point[K], boot[:, K]
        rk = R4.ranks(pk)
        brk = R4.ranks(np.nan_to_num(bk, nan=-1))
        accp = np.array([acc[mask[:, j], j].mean() if mask[:, j].any() else np.nan for j in range(M)])[K]
        out = []
        for a, j in enumerate(K):
            out.append(dict(model_id=mids[j], model=names[j], cutoff=cut[mids[j]], n=int(mask[:, j].sum()),
                            mean=float(pk[a]), lo=float(np.nanpercentile(bk[:, a], 2.5)),
                            hi=float(np.nanpercentile(bk[:, a], 97.5)), rank=int(rk[a]),
                            rlo=int(np.percentile(brk[:, a], 2.5)), rhi=int(np.percentile(brk[:, a], 97.5)),
                            p1=float((brk[:, a] == 1).mean()), acc=float(accp[a]),
                            full_mean=float(S[:, j].mean()), full_rank=0))
        full_rank = R4.ranks(S.mean(0)[K])
        for a, o in enumerate(out):
            o['full_rank'] = int(full_rank[a])
        tau = kendalltau(S.mean(0)[K], pk)[0]
        return dict(label=label, rows=sorted(out, key=lambda o: o['rank']), tau_vs_full=float(tau))

    novel = ID < 0.5
    lb_own = leaderboard(np.where(np.isnan(ID), False, novel), 'own novel set (<50 % antigen identity before own cutoff)')
    common_all = np.all(novel[:, K], axis=1)
    lb_common = leaderboard(np.repeat(common_all[:, None], M, 1), 'common subset: novel for every model with a known cutoff')
    pre24 = [j for j in K if cut[mids[j]] < '2024-09-09']
    common_pre = np.all(novel[:, pre24], axis=1)
    lb_common_pre = leaderboard(np.repeat(common_pre[:, None], M, 1),
                                'common subset excluding cutoffs after 2024-09-09 (Protenix-v1-20250630)', K=pre24)
    nocx = ~CX
    lb_nocx = leaderboard(np.where(np.isnan(ID), False, nocx), 'own set without a seen complex')
    res['leaderboards'] = dict(own=lb_own, common=lb_common, common_pre2024=lb_common_pre, no_complex=lb_nocx,
                               n_common=int(common_all.sum()), n_common_pre=int(common_pre.sum()))
    lrows = []
    for key, lb in res['leaderboards'].items():
        if not isinstance(lb, dict):
            continue
        for o in lb['rows']:
            lrows.append([key, lb['label'], o['model_id'], o['model'], o['cutoff'], o['n'], round(o['mean'], 4),
                          round(o['lo'], 4), round(o['hi'], 4), o['rank'], f"{o['rlo']}-{o['rhi']}", round(o['p1'], 3),
                          round(o['acc'] * 100, 1), round(o['full_mean'], 4), o['full_rank'], round(lb['tau_vs_full'], 3)])
    wcsv('novel_only_leaderboard.csv', ['subset', 'description', 'model_id', 'model', 'cutoff', 'n_targets',
                                        'mean_dockq', 'ci_lo', 'ci_hi', 'rank_among_known', 'rank95',
                                        'p_rank1', 'acceptable_pct', 'full_mean_dockq', 'full_rank_among_known',
                                        'kendall_tau_vs_full'], lrows)

    # ---------------- Protenix natural experiment
    pe = {}
    if 'protenix1_21' in cut and 'protenix1_25' in cut:
        a, b = mids.index('protenix1_21'), mids.index('protenix1_25')
        d = S[:, b] - S[:, a]
        i21, i25, c21, c25 = ID[:, a], ID[:, b], CX[:, a], CX[:, b]
        own25 = OWN[:, b]
        cats = np.array(['novel at both (<50%)'] * n, dtype=object)
        cats[(i21 >= 0.5)] = 'antigen seen by both (>=50%)'
        cats[(i21 < 0.5) & (i25 >= 0.5)] = 'antigen newly seen 2021->2025'
        cats[(~c21) & c25] = 'complex newly seen 2021->2025'
        cats[own25] = 'target itself in 2025 training set'
        prow = []
        for c in ['novel at both (<50%)', 'antigen seen by both (>=50%)', 'antigen newly seen 2021->2025',
                  'complex newly seen 2021->2025', 'target itself in 2025 training set']:
            m = cats == c
            if not m.any():
                prow.append([c, 0, '', '', '', '', '', '', ''])
                pe[c] = dict(n=0)
                continue
            mean, lo, hi = boot_vec(d, m)
            pe[c] = dict(n=int(m.sum()), d=float(mean), lo=float(lo), hi=float(hi),
                         s21=float(S[m, a].mean()), s25=float(S[m, b].mean()),
                         acc21=float(acc[m, a].mean()), acc25=float(acc[m, b].mean()))
            prow.append([c, int(m.sum()), round(S[m, a].mean(), 4), round(S[m, b].mean(), 4), round(mean, 4),
                         round(lo, 4), round(hi, 4), round(acc[m, a].mean() * 100, 1), round(acc[m, b].mean() * 100, 1)])
        # difference-in-differences vs novel-at-both (paired bootstrap)
        base = cats == 'novel at both (<50%)'
        for c in list(pe):
            m = cats == c
            if c == 'novel at both (<50%)' or not m.any() or not base.any():
                continue
            gs1 = np.zeros(G); gn1 = np.zeros(G); gs0 = np.zeros(G); gn0 = np.zeros(G)
            np.add.at(gs1, gidx[m], d[m]); np.add.at(gn1, gidx[m], 1)
            np.add.at(gs0, gidx[base], d[base]); np.add.at(gn0, gidx[base], 1)
            with np.errstate(invalid='ignore', divide='ignore'):
                bd = (C @ gs1) / (C @ gn1) - (C @ gs0) / (C @ gn0)
            pe[c]['did'] = float(d[m].mean() - d[base].mean())
            pe[c]['did_lo'], pe[c]['did_hi'] = float(np.nanpercentile(bd, 2.5)), float(np.nanpercentile(bd, 97.5))
        m_after = ~own25
        gs1 = np.zeros(G); gn1 = np.zeros(G); gs0 = np.zeros(G); gn0 = np.zeros(G)
        np.add.at(gs1, gidx[own25], d[own25]); np.add.at(gn1, gidx[own25], 1)
        np.add.at(gs0, gidx[m_after], d[m_after]); np.add.at(gn0, gidx[m_after], 1)
        with np.errstate(invalid='ignore', divide='ignore'):
            bd = (C @ gs1) / (C @ gn1) - (C @ gs0) / (C @ gn0)
        pe['_own_vs_after'] = dict(did=float(d[own25].mean() - d[m_after].mean()),
                                   lo=float(np.nanpercentile(bd, 2.5)), hi=float(np.nanpercentile(bd, 97.5)))
        mean, lo, hi = boot_vec(d, m_after)
        pe['_targets_released_after_2025_cutoff'] = dict(n=int(m_after.sum()), d=float(mean), lo=float(lo), hi=float(hi))
        mean, lo, hi = boot_vec(d)
        pe['_all'] = dict(n=n, d=float(mean), lo=float(lo), hi=float(hi))
        wcsv('protenix_natural_experiment.csv', ['stratum', 'n_targets', 'protenix_v1_2021_mean', 'protenix_v1_2025_mean',
                                                 'diff_2025_minus_2021', 'ci_lo', 'ci_hi', 'acc_2021_pct', 'acc_2025_pct'],
             prow + [['all targets released after 2025-06-30', pe['_targets_released_after_2025_cutoff']['n'], '', '',
                      round(pe['_targets_released_after_2025_cutoff']['d'], 4),
                      round(pe['_targets_released_after_2025_cutoff']['lo'], 4),
                      round(pe['_targets_released_after_2025_cutoff']['hi'], 4), '', ''],
                     ['all targets', n, '', '', round(pe['_all']['d'], 4), round(pe['_all']['lo'], 4),
                      round(pe['_all']['hi'], 4), '', '']])
    res['protenix'] = pe

    # ---------------- Mol-1 consistency check
    mc = {}
    if 'mol1' in cut:
        jm = mids.index('mol1')
        c0 = cut['mol1']
        peers = [j for j in known if cut[mids[j]] == c0 and j != jm]
        later = [j for j in known if cut[mids[j]] > c0 and not (mids[j] == 'protenix1_25')]
        # "seen later": not seen at Mol-1's cutoff, but seen before the target's own release
        at_cut = ID[:, jm]
        at_rel = np.array([novelty_at(k, trel[i])[0] for i, k in enumerate(keys)])
        cx_cut = CX[:, jm]
        cx_rel = np.array([novelty_at(k, trel[i])[1] for i, k in enumerate(keys)])
        cat = np.array(['novel until target release'] * n, dtype=object)
        cat[at_cut >= 0.5] = 'antigen seen before 2021-09-30'
        cat[(at_cut < 0.5) & (at_rel >= 0.5)] = 'antigen first seen after 2021-09-30'
        cat[(~cx_cut) & cx_rel] = 'complex first seen after 2021-09-30'
        cat[cx_cut] = 'complex seen before 2021-09-30'
        rel_peer = S[:, jm] - S[:, peers].mean(1) if peers else np.full(n, np.nan)
        mrow = []
        base = cat == 'novel until target release'
        order = ['novel until target release', 'antigen first seen after 2021-09-30', 'complex first seen after 2021-09-30',
                 'antigen seen before 2021-09-30', 'complex seen before 2021-09-30']
        for c in order:
            m = cat == c
            entry = dict(n=int(m.sum()))
            if m.any():
                entry['mol1'] = float(S[m, jm].mean())
                entry['peers'] = float(S[m][:, peers].mean()) if peers else np.nan
                mean, lo, hi = boot_vec(rel_peer, m)
                entry.update(rel=float(mean), lo=float(lo), hi=float(hi))
                if later:
                    rl = S[:, later].mean(1) - S[:, peers].mean(1)
                    mean2, lo2, hi2 = boot_vec(rl, m)
                    entry.update(later_rel=float(mean2), later_lo=float(lo2), later_hi=float(hi2))
                if c != order[0] and base.any():
                    gs1 = np.zeros(G); gn1 = np.zeros(G); gs0 = np.zeros(G); gn0 = np.zeros(G)
                    np.add.at(gs1, gidx[m], rel_peer[m]); np.add.at(gn1, gidx[m], 1)
                    np.add.at(gs0, gidx[base], rel_peer[base]); np.add.at(gn0, gidx[base], 1)
                    with np.errstate(invalid='ignore', divide='ignore'):
                        bd = (C @ gs1) / (C @ gn1) - (C @ gs0) / (C @ gn0)
                    entry['did'] = float(rel_peer[m].mean() - rel_peer[base].mean())
                    entry['did_lo'], entry['did_hi'] = float(np.nanpercentile(bd, 2.5)), float(np.nanpercentile(bd, 97.5))
            mc[c] = entry
            mrow.append([c, entry['n'], round(entry.get('mol1', np.nan), 4), round(entry.get('peers', np.nan), 4),
                         round(entry.get('rel', np.nan), 4), round(entry.get('lo', np.nan), 4), round(entry.get('hi', np.nan), 4),
                         round(entry.get('did', np.nan), 4), round(entry.get('did_lo', np.nan), 4),
                         round(entry.get('did_hi', np.nan), 4), round(entry.get('later_rel', np.nan), 4)])
        # unique wins: a top model beating the best of the other top models by >= 0.3 DockQ
        from scipy.stats import fisher_exact
        top = [mids.index(x) for x in ('mol1', 'opendde', 'opendde_abag', 'protenix2') if x in mids]
        cats4 = np.array(['novel until target release' if c == 'novel until target release' else
                          'antigen seen' if c.startswith('antigen') else c for c in cat], dtype=object)
        uw = {}
        for j in top:
            others = [o for o in top if o != j]
            win = S[:, j] - S[:, others].max(1) >= 0.3
            uw[names[j]] = {c: [int(win[cats4 == c].sum()), int((cats4 == c).sum())] for c in
                            ['novel until target release', 'antigen seen', 'complex seen before 2021-09-30',
                             'complex first seen after 2021-09-30']}
            if mids[j] == 'mol1':
                inside = cats4 == 'complex first seen after 2021-09-30'
                tab = [[int(win[inside].sum()), int((~win[inside]).sum())],
                       [int(win[~inside].sum()), int((~win[~inside]).sum())]]
                mc['_unique_win_fisher_p'] = float(fisher_exact(tab, alternative='greater')[1])
        mc['_unique_wins'] = uw
        lst = []
        for i, k in enumerate(keys):
            if cat[i] != 'complex first seen after 2021-09-30':
                continue
            first = min((h for h in hits[k] if h[3] and c0 <= h[0] < trel[i]), key=lambda h: h[0])
            ents = [r for r in csv.DictReader(open(os.path.join(NOV, 'target_entry_hits.tsv')), delimiter='\t')
                    if r['target'] == k and r['complex'] == '1' and c0 <= r['release'] < trel[i]]
            fe = min(ents, key=lambda r: r['release'])
            row = dict(target=k, antigen='; '.join(ag_desc.get(k, []))[:60], first_entry=fe['entry'],
                       first_release=fe['release'], n_entries=len(ents),
                       is_arena_target=fe['entry'] in {t['pdb'].upper() for t in targets},
                       mol1=float(S[i, jm]), best_other_top=float(S[i, [o for o in top if o != jm]].max()),
                       peers=float(S[i, peers].mean()))
            for x in ('boltz2', 'rf3_25', 'promera', 'protenix1_25'):
                if x in mids:
                    row[x] = float(S[i, mids.index(x)])
                    row[x + '_cutoff_after_first'] = cut.get(x, '9999') > fe['release']
            lst.append(row)
        mc['_complex_after_targets'] = sorted(lst, key=lambda r: r['first_release'])
        wcsv('mol1_complexes_first_seen_after_cutoff.csv', list(lst[0].keys()) if lst else ['target'],
             [[r[c] for c in lst[0].keys()] for r in lst] if lst else [])
        mc['_peers'] = [names[j] for j in peers]
        mc['_later'] = [names[j] for j in later]
        wcsv('mol1_consistency_check.csv', ['stratum', 'n_targets', 'mol1_mean', 'same_cutoff_peers_mean',
                                            'mol1_minus_peers', 'ci_lo', 'ci_hi', 'diff_vs_novel_stratum', 'did_ci_lo',
                                            'did_ci_hi', 'later_cutoff_models_minus_peers'], mrow)
    res['mol1'] = mc

    # ---------------- B. clonotypes
    ab = {r['target']: r for r in csv.DictReader(open(os.path.join(NOV, 'target_antibody_anarci_allspecies.tsv')), delimiter='\t')}

    def gene(r):
        return (r['species'], r['vgene'].split('*')[0]) if r['vgene'] else None

    def same_clonotype(a, b):
        ra, rb = ab.get(a), ab.get(b)
        if not ra or not rb or not ra['h3'] or not rb['h3'] or gene(ra) is None or gene(ra) != gene(rb):
            return False
        return len(ra['h3']) == len(rb['h3']) and \
            sum(x == y for x, y in zip(ra['h3'], rb['h3'])) / len(ra['h3']) >= 0.8

    clono_links, joint_links = [], []
    for x in range(n):
        for y in range(x + 1, n):
            if same_clonotype(keys[x], keys[y]):
                clono_links.append({keys[x], keys[y]})
                if ag50[x] == ag50[y]:
                    joint_links.append({keys[x], keys[y]})
    clono = C3.relabel(C3.components(keys, clono_links), 'CLONO', keys)
    joint = C3.relabel(C3.components(keys, joint_links), 'AGCL', keys)
    vh90 = [groups[k]['vh90'] for k in keys]
    crow = []
    for k in keys:
        r = ab.get(k, {})
        crow.append([k, groups[k]['ag50'], groups[k]['vh90'], groups[k]['vh90_size'], clono[k],
                     Counter(clono.values())[clono[k]], joint[k], Counter(joint.values())[joint[k]],
                     r.get('vgene', ''), r.get('species', ''), r.get('h3', ''), r.get('h3_len', '')])
    wcsv('antibody_clonotypes.csv', ['target', 'ag50_group', 'vh90_group', 'vh90_size', 'clonotype', 'clonotype_size',
                                     'ag50_clonotype_group', 'ag50_clonotype_size', 'ighv_gene', 'species', 'cdr_h3',
                                     'cdr_h3_len'], crow)
    # VH90 multi-member clusters: how many distinct clonotypes / antigens do they contain?
    vh_multi = defaultdict(list)
    for k in keys:
        vh_multi[groups[k]['vh90']].append(k)
    vh_report = []
    for g, mem in sorted(vh_multi.items(), key=lambda kv: -len(kv[1])):
        if len(mem) < 2:
            continue
        vh_report.append(dict(group=g, n=len(mem), n_clonotypes=len({clono[k] for k in mem}),
                              n_ag50=len({groups[k]['ag50'] for k in mem}),
                              antigens='; '.join(sorted({d for k in mem for d in ag_desc.get(k, [])[:1]}))[:160],
                              h3s=', '.join(ab.get(k, {}).get('h3', '?') for k in mem),
                              vgenes=', '.join(sorted({ab.get(k, {}).get('vgene', '?') for k in mem}))))

    raw_mean = S.mean(0)
    cres = {}
    crows = []
    for lab, labels in [('clonotype', [clono[k] for k in keys]), ('ag50+clonotype', [joint[k] for k in keys]),
                        ('vh90 (old)', vh90), ('ag50+vh90 (old joint)', [groups[k]['joint'] for k in keys])]:

        gi, Gc = R4.group_index(labels)
        wm, _, _ = R4.weighted(S, gi, Gc)
        bw = R4.cluster_boot(S, gi, Gc, rng, 'weighted')
        lo, hi, rlo, rhi, p1, _ = R4.summarize_boot(bw)
        tau = kendalltau(raw_mean, wm)[0]
        taus = []
        for _ in range(1000):
            w2, _, _ = R4.weighted(S, rng.permutation(gi), Gc)
            taus.append(kendalltau(raw_mean, w2)[0])
        sizes = Counter(labels)
        cres[lab] = dict(n_groups=Gc, n_in_multi=int(sum(v for v in sizes.values() if v > 1)),
                         largest=int(max(sizes.values())), tau=float(tau), tau_null_median=float(np.median(taus)),
                         p_tau=float((np.array(taus) <= tau).mean()), dockq=wm.tolist(), rank=R4.ranks(wm).tolist(),
                         lo=lo.tolist(), hi=hi.tolist(), rlo=rlo.tolist(), rhi=rhi.tolist(), p1=p1.tolist())
        for j in np.argsort(-wm):
            crows.append([lab, Gc, names[j], round(wm[j], 4), round(wm[j] - raw_mean[j], 4), cres[lab]['rank'][j],
                          f'{rlo[j]}-{rhi[j]}', round(p1[j], 3), round(tau, 3), round(cres[lab]['p_tau'], 3)])
    wcsv('clonotype_rescoring.csv', ['grouping', 'n_groups', 'model', 'weighted_dockq', 'delta_vs_raw', 'rank',
                                     'rank95', 'p_rank1', 'kendall_tau_vs_raw', 'p_perm_tau'], crows)
    res['clonotype'] = dict(results=cres, vh90_multi=vh_report,
                            n_with_h3=sum(1 for k in keys if ab.get(k, {}).get('h3')),
                            n_with_gene=sum(1 for k in keys if ab.get(k, {}).get('vgene')))

    # ---------------- C. fusion-partner masking: re-cluster on masked antigen chains
    fm = list(csv.DictReader(open(os.path.join(NOV, 'fusion_masking.tsv')), delimiter='\t'))
    masked_fa = os.path.join(NOV, 'antigen_chains_masked.fasta')
    mdir = os.path.join(NOV, 'masked_clusters')
    os.makedirs(mdir, exist_ok=True)
    C3.MM = mdir
    tgt_of = lambda rec: rec.split('|')[0]
    newg = {}
    for ident in (0.3, 0.5, 0.9):
        tag = f'ag{int(ident * 100)}'
        clus = C3.easy_cluster(masked_fa, os.path.join(mdir, tag), ident)
        byc = defaultdict(set)
        for mem, rep in clus.items():
            byc[rep].add(tgt_of(mem))
        newg[tag] = C3.relabel(C3.components(keys, byc.values()), tag.upper() + 'M', keys)
    # shorter-chain coverage on masked sequences
    import subprocess, shutil
    hits_m = os.path.join(mdir, 'allvsall.m8')
    subprocess.run([C3.MMSEQS, 'easy-search', masked_fa, masked_fa, hits_m, os.path.join(mdir, 'tmp'), '-s', '7.5',
                    '--threads', '4', '-v', '1', '--format-output', 'query,target,fident,qcov,tcov,evalue', '-e', '1e-3'],
                   check=True, stdout=subprocess.DEVNULL)
    shutil.rmtree(os.path.join(mdir, 'tmp'), ignore_errors=True)
    links = [{tgt_of(q), tgt_of(s)} for q, s, fid, qc, tc, ev in csv.reader(open(hits_m), delimiter='\t')
             if float(fid) >= 0.5 and max(float(qc), float(tc)) >= 0.8]
    newg['ag50_shortcov'] = C3.relabel(C3.components(keys, links), 'AG50SM', keys)
    frow, fres = [], {}
    for tag in ['ag30', 'ag50', 'ag90', 'ag50_shortcov']:
        old = {k: groups[k][tag] for k in keys}
        new = newg[tag]
        # pairs together before but not after (split by masking), and vice versa
        def pairs(assign):
            by = defaultdict(list)
            for k in keys:
                by[assign[k]].append(k)
            return {(a, b) for mem in by.values() for x, a in enumerate(mem) for b in mem[x + 1:]}
        po, pn = pairs(old), pairs(new)
        split, merged = sorted(po - pn), sorted(pn - po)
        gi, Gm = R4.group_index([new[k] for k in keys])
        wm, _, _ = R4.weighted(S, gi, Gm)
        go, Go = R4.group_index([old[k] for k in keys])
        wo, _, _ = R4.weighted(S, go, Go)
        fres[tag] = dict(groups_before=Go, groups_after=Gm, split_pairs=split[:50], n_split=len(split),
                         n_merged=len(merged), tau_before_after=float(kendalltau(wo, wm)[0]),
                         max_abs_score_change=float(np.abs(wm - wo).max()),
                         rank_changes=[(names[j], int(R4.ranks(wo)[j]), int(R4.ranks(wm)[j])) for j in range(M)
                                       if R4.ranks(wo)[j] != R4.ranks(wm)[j]])
        for a, b in split:
            frow.append([tag, 'split by masking', a, b, '; '.join(ag_desc.get(a, [])), '; '.join(ag_desc.get(b, []))])
        for a, b in merged:
            frow.append([tag, 'merged after masking', a, b, '; '.join(ag_desc.get(a, [])), '; '.join(ag_desc.get(b, []))])
    wcsv('fusion_masking_effect.csv', ['grouping', 'change', 'target_a', 'target_b', 'antigen_a', 'antigen_b'], frow)
    res['fusion'] = dict(masking=fm, effect=fres)

    json.dump(res, open(os.path.join(DER, 'results_novelty.json'), 'w'), default=float)

    # ---------------- console summary
    print('\n== novel-only leaderboard (own set) ==')
    for o in lb_own['rows']:
        print(f"{o['model']:22s} n={o['n']:3d} {o['mean']:.3f} [{o['lo']:.3f},{o['hi']:.3f}] #{o['rank']} "
              f"[{o['rlo']}-{o['rhi']}] p1={o['p1']:.2f} | full {o['full_mean']:.3f} #{o['full_rank']}")
    print('tau', round(lb_own['tau_vs_full'], 3))
    for key in ['common', 'common_pre2024']:
        lb = res['leaderboards'][key]
        print(f"\n== {lb['label']} n={res['leaderboards']['n_' + key.replace('common', 'common').split('_')[0]] if False else ''}")
        for o in lb['rows'][:8]:
            print(f"{o['model']:22s} n={o['n']:3d} {o['mean']:.3f} [{o['lo']:.3f},{o['hi']:.3f}] #{o['rank']} [{o['rlo']}-{o['rhi']}]")
    print('\n== strata (n / mean) ==')
    for j in known:
        d = res['strata'][mids[j]]
        print(f"{names[j]:22s} " + ' '.join(f"{lab}:{d[lab]['n']}/{d[lab]['mean']:.3f}" for lab, _, _ in STRATA)
              + f"  cx:{d['seen complex']['n']}/{d['seen complex']['mean']:.3f}")
    print('\n== protenix ==')
    for c, v in pe.items():
        print(c, v)
    print('\n== mol1 ==')
    for c, v in mc.items():
        print(c, v)
    print('\n== clonotype ==')
    for lab, v in cres.items():
        print(lab, {k: v[k] for k in ('n_groups', 'n_in_multi', 'largest', 'tau', 'tau_null_median', 'p_tau')})
    for v in vh_report[:6]:
        print(v)
    print('\n== fusion ==')
    for tag, v in fres.items():
        print(tag, {k: v[k] for k in ('groups_before', 'groups_after', 'n_split', 'n_merged', 'tau_before_after',
                                      'max_abs_score_change', 'rank_changes')})


if __name__ == '__main__':
    main()
