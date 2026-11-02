#!/usr/bin/env python3
"""Re-score the CoFold Arena Ab-Ag leaderboard with homology-cluster weighting.

For every grouping (raw = each target its own group; antigen 30/50/90; joint Ag50+VH90;
plus sensitivity groupings ag50_shortcov and vh90) we compute, per model:
  * weighted mean DockQ and weighted acceptable/medium/high rates, weight = 1/group size
  * the mean of group means (must equal the weighted mean - asserted)
  * cluster bootstrap (resample whole groups, B=2000, same resample for every model):
    95 % percentile CI of the score, 95 % rank interval, P(rank 1)
  * for raw scores, both an ordinary per-target bootstrap and a cluster bootstrap
  * Kendall tau vs raw ranking, per-model rank / score change
  * leave-one-group-out (10 largest groups): change of each model's raw and weighted score
  * Fab / VHH breakdown (group sizes recomputed within each subset)

Outputs (in arena_homology/):
  targets_with_clusters.csv, leaderboard_comparison.csv, leaderboard_long.csv,
  logo_largest_groups.csv, format_breakdown.csv, data/derived/results.json
"""
import csv, json, os
import numpy as np
from scipy.stats import kendalltau

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DER = os.path.join(ROOT, 'data', 'derived')
B = 2000
SEED = 20261102
THR = {'acc': 0.23, 'med': 0.49, 'high': 0.80}
PRIMARY = ['raw', 'ag30', 'ag50', 'ag90', 'joint']
SENS = ['ag50_shortcov', 'vh90']
GROUPINGS = PRIMARY + SENS
EXCLUDE = {'esmfold2_fast'}  # not on the leaderboard (has no MSA-mode scores)


def load():
    targets = json.load(open(os.path.join(DER, 'targets.json')))
    scores = json.load(open(os.path.join(DER, 'scores.json')))
    models = [m for m in json.load(open(os.path.join(DER, 'models.json'))) if m['id'] not in EXCLUDE]
    keys = [t['pdb'] + '_' + t['abId'].replace('sabdab2_', '') for t in targets]
    S = np.array([[scores[k][m['id']]['msa'] for m in models] for k in keys], dtype=float)
    groups = {r['target']: r for r in csv.DictReader(open(os.path.join(DER, 'groups.csv')))}
    ag_desc = {}
    for r in csv.DictReader(open(os.path.join(DER, 'chain_sequences.tsv')), delimiter='\t'):
        if r['role'] == 'AG':
            ag_desc.setdefault(r['target'], [])
            if r['desc'] not in ag_desc[r['target']]:
                ag_desc[r['target']].append(r['desc'])
    return targets, scores, models, keys, S, groups, ag_desc


def group_index(labels):
    uniq = {g: i for i, g in enumerate(dict.fromkeys(labels))}
    return np.array([uniq[g] for g in labels]), len(uniq)


def metric_matrix(S):
    return {'dockq': S, **{k: (S >= v).astype(float) for k, v in THR.items()}}


def weighted(X, gidx, G):
    size = np.bincount(gidx, minlength=G)
    w = 1.0 / size[gidx]
    wm = (w[:, None] * X).sum(0) / w.sum()
    gsum = np.zeros((G, X.shape[1]))
    np.add.at(gsum, gidx, X)
    mgm = (gsum / size[:, None]).mean(0)
    assert np.allclose(wm, mgm, atol=1e-12), 'weighted mean != mean of group means'
    return wm, gsum, size


def ranks(v):
    """rank 1 = best; ties broken by order (scores rarely tie)."""
    order = np.argsort(-v, axis=-1, kind='stable')
    r = np.empty_like(order)
    np.put_along_axis(r, order, np.arange(1, v.shape[-1] + 1)[None, :].repeat(v.shape[0], 0) if v.ndim == 2
                      else np.arange(1, v.shape[-1] + 1), axis=-1)
    return r


def cluster_boot(X, gidx, G, rng, estimator):
    """Resample groups with replacement. estimator 'weighted' = mean of group means;
    'raw' = plain mean over all targets in resampled groups."""
    _, gsum, size = weighted(X, gidx, G)
    C = rng.multinomial(G, np.full(G, 1.0 / G), size=B).astype(float)  # B x G counts
    if estimator == 'weighted':
        return C @ (gsum / size[:, None]) / G
    return (C @ gsum) / (C @ size)[:, None]


def summarize_boot(boot):
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
    rk = ranks(boot)
    rlo, rhi = np.percentile(rk, [2.5, 97.5], axis=0)
    p1 = (rk == 1).mean(0)
    return lo, hi, rlo.round().astype(int), rhi.round().astype(int), p1, boot.std(0)


def main():
    targets, scores, models, keys, S, groups, ag_desc = load()
    n, M = S.shape
    names = [m['name'] for m in models]
    mids = [m['id'] for m in models]
    X = metric_matrix(S)
    rng = np.random.default_rng(SEED)
    labels = {'raw': keys}
    for g in GROUPINGS[1:]:
        labels[g] = [groups[k][g] for k in keys]

    # sanity: published numbers
    raw_mean = S.mean(0)
    pub = {'mol1': 0.359, 'opendde': 0.355, 'protenix2': 0.327, 'boltz2': 0.175, 'afmult_np': 0.106}
    for mid, v in pub.items():
        assert round(raw_mean[mids.index(mid)], 3) == v, (mid, raw_mean[mids.index(mid)])
    assert round(X['acc'][:, mids.index('mol1')].mean() * 100, 1) == 48.8

    res = {'models': names, 'model_ids': mids, 'n_targets': n, 'B': B, 'seed': SEED, 'groupings': {}}
    raw_rank = ranks(raw_mean)
    for g in GROUPINGS:
        gidx, G = group_index(labels[g])
        out = {'n_groups': G}
        for met, Xm in X.items():
            wm, _, size = weighted(Xm, gidx, G)
            out[met] = wm.tolist()
        wmean = np.array(out['dockq'])
        out['rank'] = ranks(wmean).tolist()
        bw = cluster_boot(S, gidx, G, rng, 'weighted')
        lo, hi, rlo, rhi, p1, sd = summarize_boot(bw)
        out.update(ci_lo=lo.tolist(), ci_hi=hi.tolist(), rank_lo=rlo.tolist(), rank_hi=rhi.tolist(),
                   p_rank1=p1.tolist(), boot_sd=sd.tolist())
        ba = cluster_boot(X['acc'], gidx, G, rng, 'weighted')
        out['acc_ci_lo'], out['acc_ci_hi'] = np.percentile(ba, [2.5, 97.5], axis=0).tolist()
        # raw (unweighted) score but with this grouping's cluster bootstrap -> honest CI for the published mean
        br = cluster_boot(S, gidx, G, rng, 'raw')
        lo2, hi2, rlo2, rhi2, p12, sd2 = summarize_boot(br)
        out.update(rawscore_clusterboot_ci_lo=lo2.tolist(), rawscore_clusterboot_ci_hi=hi2.tolist(),
                   rawscore_clusterboot_rank_lo=rlo2.tolist(), rawscore_clusterboot_rank_hi=rhi2.tolist(),
                   rawscore_clusterboot_p_rank1=p12.tolist(), rawscore_clusterboot_sd=sd2.tolist())
        tau, p = kendalltau(raw_mean, wmean)
        out['kendall_tau_vs_raw'] = float(tau)
        out['kendall_p'] = float(p)
        out['effective_n'] = float(G)
        out['kish_ess'] = float((1 / np.bincount(gidx)[gidx]).sum() ** 2 / ((1 / np.bincount(gidx)[gidx]) ** 2).sum())
        res['groupings'][g] = out

    # permutation null: shuffle target->group assignment keeping the group-size distribution.
    # Tells whether the observed shift is specific to homology structure or what any
    # reweighting of this shape would do.
    res['perm_null'] = {}
    for g in ['ag30', 'ag50', 'ag90']:
        gidx, G = group_index(labels[g])
        obs = np.array(res['groupings'][g]['dockq']) - raw_mean
        taus, deltas = [], []
        for _ in range(1000):
            wm, _, _ = weighted(S, rng.permutation(gidx), G)
            taus.append(kendalltau(raw_mean, wm)[0])
            deltas.append(wm - raw_mean)
        deltas = np.array(deltas)
        # two-sided p for each model's delta relative to the mean shift across models (relative gain)
        rel_obs = obs - obs.mean()
        rel_null = deltas - deltas.mean(1, keepdims=True)
        p_rel = ((np.abs(rel_null) >= np.abs(rel_obs)[None, :]).mean(0))
        res['perm_null'][g] = dict(tau_null_median=float(np.median(taus)), tau_null_lo=float(np.percentile(taus, 2.5)),
                                   p_tau_le_obs=float((np.array(taus) <= res['groupings'][g]['kendall_tau_vs_raw']).mean()),
                                   delta_null_sd=deltas.std(0).tolist(), rel_delta_p=p_rel.tolist(),
                                   mean_shift_obs=float(obs.mean()), mean_shift_null=float(deltas.mean()))
        print(g, 'perm null tau median', round(np.median(taus), 3), 'P(tau<=obs)', res['perm_null'][g]['p_tau_le_obs'],
              'mean shift obs', round(obs.mean(), 4), 'null', round(deltas.mean(), 4),
              'models with rel-delta p<0.05:', [names[j] for j in range(M) if p_rel[j] < 0.05])

    # pairwise: top-model differences, paired via same resample
    top = np.argsort(-raw_mean)[:5]
    pair = {}
    for g in ['raw', 'ag30', 'ag50', 'ag90', 'joint', 'ag50_shortcov']:
        gidx, G = group_index(labels[g])
        bw = cluster_boot(S, gidx, G, rng, 'weighted')
        d = {}
        for a in top:
            for b in top:
                if b == a:
                    continue
                diff = bw[:, a] - bw[:, b]
                d[f'{mids[a]}-{mids[b]}'] = dict(lo=float(np.percentile(diff, 2.5)), hi=float(np.percentile(diff, 97.5)),
                                                 p_gt0=float((diff > 0).mean()))
        pair[g] = d
    # raw score, cluster bootstrap on ag50 for the pairwise
    gidx, G = group_index(labels['ag50'])
    br = cluster_boot(S, gidx, G, rng, 'raw')
    pair['raw_score_ag50_clusterboot'] = {f'{mids[a]}-{mids[b]}': dict(
        lo=float(np.percentile(br[:, a] - br[:, b], 2.5)), hi=float(np.percentile(br[:, a] - br[:, b], 97.5)),
        p_gt0=float((br[:, a] > br[:, b]).mean())) for a in top for b in top if a != b}
    res['pairwise_top5'] = pair
    res['top5'] = [mids[i] for i in top]

    # leave-one-group-out for the 10 largest groups of each antigen grouping
    logo_rows = []
    res['logo'] = {}
    for g in ['ag30', 'ag50', 'ag90', 'joint', 'ag50_shortcov']:
        lab = np.array(labels[g])
        uniq, cnt = np.unique(lab, return_counts=True)
        largest = [u for u, _ in sorted(zip(uniq, cnt), key=lambda x: (-x[1], x[0]))[:10]]
        gidx, G = group_index(labels[g])
        wfull = np.array(res['groupings'][g]['dockq'])
        res['logo'][g] = []
        for grp in largest:
            keep = lab != grp
            raw_drop = S[keep].mean(0) - raw_mean
            gi2, G2 = group_index(list(lab[keep]))
            w_drop, _, _ = weighted(S[keep], gi2, G2)
            w_drop = w_drop - wfull
            rk_raw = ranks(S[keep].mean(0))
            members = [k for k in np.array(keys)[~keep]]
            names_in = sorted({d for k in members for d in ag_desc.get(k, [])})
            label = '; '.join(names_in[:3]) + (f' (+{len(names_in)-3})' if len(names_in) > 3 else '')
            entry = dict(group=grp, size=int((~keep).sum()), antigens=label,
                         d_raw=raw_drop.tolist(), d_weighted=w_drop.tolist(),
                         rank_raw_after=rk_raw.tolist(),
                         group_mean=S[~keep].mean(0).tolist())
            res['logo'][g].append(entry)
            for j in range(M):
                logo_rows.append([g, grp, entry['size'], label, names[j], round(S[~keep, j].mean(), 4),
                                  round(raw_drop[j], 4), round(w_drop[j], 4), int(raw_rank[j]), int(rk_raw[j])])

    # format breakdown (group sizes recomputed within subset)
    fmt_rows = []
    res['format'] = {}
    fmts = np.array([t['format'] for t in targets])
    for f in ['Fab', 'VHH']:
        sel = fmts == f
        res['format'][f] = {'n': int(sel.sum())}
        base = S[sel].mean(0)
        for g in GROUPINGS:
            lab = list(np.array(labels[g])[sel])
            gidx, G = group_index(lab)
            wm, _, _ = weighted(S[sel], gidx, G)
            wa, _, _ = weighted(X['acc'][sel], gidx, G)
            bw = cluster_boot(S[sel], gidx, G, rng, 'weighted')
            lo, hi, rlo, rhi, p1, _ = summarize_boot(bw)
            rk = ranks(wm)
            tau = kendalltau(base, wm)[0] if g != 'raw' else 1.0
            res['format'][f][g] = dict(n_groups=G, dockq=wm.tolist(), acc=wa.tolist(), rank=rk.tolist(),
                                       ci_lo=lo.tolist(), ci_hi=hi.tolist(), rank_lo=rlo.tolist(),
                                       rank_hi=rhi.tolist(), p_rank1=p1.tolist(), tau=float(tau))
            for j in range(M):
                fmt_rows.append([f, g, int(sel.sum()), G, names[j], round(wm[j], 4), round(lo[j], 4), round(hi[j], 4),
                                 int(rk[j]), int(rlo[j]), int(rhi[j]), round(p1[j], 3), round(wa[j] * 100, 1),
                                 round(tau, 3)])

    # ---------- write CSVs ----------
    def w(path, header, rows):
        with open(os.path.join(ROOT, path), 'w', newline='') as fh:
            cw = csv.writer(fh)
            cw.writerow(header)
            cw.writerows(rows)

    w('logo_largest_groups.csv', ['grouping', 'group', 'size', 'antigens', 'model', 'model_mean_in_group',
                                  'delta_raw_score_if_dropped', 'delta_weighted_score_if_dropped',
                                  'raw_rank', 'raw_rank_if_dropped'], logo_rows)
    w('format_breakdown.csv', ['format', 'grouping', 'n_targets', 'n_groups', 'model', 'score', 'ci_lo', 'ci_hi',
                               'rank', 'rank_lo95', 'rank_hi95', 'p_rank1', 'acceptable_pct', 'kendall_tau_vs_raw'],
      fmt_rows)

    raw = res['groupings']['raw']
    hdr = ['model_id', 'model', 'org', 'raw_dockq', 'raw_rank', 'raw_ci_lo_targetboot', 'raw_ci_hi_targetboot',
           'raw_rank95_targetboot', 'raw_p_rank1_targetboot',
           'raw_ci_lo_ag50clusterboot', 'raw_ci_hi_ag50clusterboot', 'raw_rank95_ag50clusterboot',
           'raw_acc_pct', 'raw_med_pct', 'raw_high_pct']
    for g in GROUPINGS[1:]:
        hdr += [f'{g}_dockq', f'{g}_delta', f'{g}_rank', f'{g}_rank_change', f'{g}_ci_lo', f'{g}_ci_hi',
                f'{g}_rank95', f'{g}_p_rank1', f'{g}_acc_pct', f'{g}_med_pct', f'{g}_high_pct']
    rows = []
    a50 = res['groupings']['ag50']
    for j in np.argsort(raw_rank):
        r = [mids[j], names[j], models[j]['org'], round(raw['dockq'][j], 4), raw['rank'][j],
             round(raw['ci_lo'][j], 4), round(raw['ci_hi'][j], 4), f"{raw['rank_lo'][j]}-{raw['rank_hi'][j]}",
             round(raw['p_rank1'][j], 3),
             round(a50['rawscore_clusterboot_ci_lo'][j], 4), round(a50['rawscore_clusterboot_ci_hi'][j], 4),
             f"{a50['rawscore_clusterboot_rank_lo'][j]}-{a50['rawscore_clusterboot_rank_hi'][j]}",
             round(raw['acc'][j] * 100, 1), round(raw['med'][j] * 100, 1), round(raw['high'][j] * 100, 1)]
        for g in GROUPINGS[1:]:
            o = res['groupings'][g]
            r += [round(o['dockq'][j], 4), round(o['dockq'][j] - raw['dockq'][j], 4), o['rank'][j],
                  raw['rank'][j] - o['rank'][j], round(o['ci_lo'][j], 4), round(o['ci_hi'][j], 4),
                  f"{o['rank_lo'][j]}-{o['rank_hi'][j]}", round(o['p_rank1'][j], 3),
                  round(o['acc'][j] * 100, 1), round(o['med'][j] * 100, 1), round(o['high'][j] * 100, 1)]
        rows.append(r)
    w('leaderboard_comparison.csv', hdr, rows)

    long = []
    for g in GROUPINGS:
        o = res['groupings'][g]
        for j in range(M):
            long.append([g, o['n_groups'], mids[j], names[j], round(o['dockq'][j], 4), round(o['ci_lo'][j], 4),
                         round(o['ci_hi'][j], 4), o['rank'][j], o['rank_lo'][j], o['rank_hi'][j],
                         round(o['p_rank1'][j], 3), round(o['acc'][j] * 100, 1), round(o['acc_ci_lo'][j] * 100, 1),
                         round(o['acc_ci_hi'][j] * 100, 1), round(o['med'][j] * 100, 1),
                         round(o['high'][j] * 100, 1), round(o['kendall_tau_vs_raw'], 3)])
    w('leaderboard_long.csv', ['grouping', 'n_groups', 'model_id', 'model', 'weighted_dockq', 'ci_lo', 'ci_hi',
                               'rank', 'rank_lo95', 'rank_hi95', 'p_rank1', 'acceptable_pct', 'acc_ci_lo',
                               'acc_ci_hi', 'medium_pct', 'high_pct', 'kendall_tau_vs_raw'], long)

    # per-target table
    allm = json.load(open(os.path.join(DER, 'models.json')))
    hdr = ['target', 'pdb', 'abId', 'antigen_arena', 'antigen_rcsb', 'agUniProt', 'format', 'Hchain', 'Lchain',
           'agChains', 'release', 'resolution', 'method']
    for g in GROUPINGS[1:]:
        hdr += [f'{g}_group', f'{g}_size', f'{g}_weight']
    hdr += [f'msa_{m["id"]}' for m in models] + [f'nomsa_{m["id"]}' for m in allm]
    trows = []
    for i, (k, t) in enumerate(zip(keys, targets)):
        r = [k, t['pdb'], t['abId'], t['antigen'], ' | '.join(ag_desc.get(k, [])), t['agUniProt'], t['format'],
             t['Hchain'], t['Lchain'] or '', ';'.join(t['agChains']), t['release'], t['resolution'], t['method']]
        for g in GROUPINGS[1:]:
            sz = int(groups[k][g + '_size'])
            r += [groups[k][g], sz, round(1 / sz, 6)]
        r += [scores[k][m['id']]['msa'] for m in models]
        r += [scores[k][m['id']].get('nomsa', '') for m in allm]
        trows.append(r)
    w('targets_with_clusters.csv', hdr, trows)
    json.dump(res, open(os.path.join(DER, 'results.json'), 'w'))

    # console summary
    print(f"{'model':22s} {'raw':>6s} " + ' '.join(f'{g:>16s}' for g in GROUPINGS[1:]))
    for j in np.argsort(raw_rank):
        line = f"{names[j]:22s} {raw['dockq'][j]:.3f}#{raw['rank'][j]:<2d}"
        for g in GROUPINGS[1:]:
            o = res['groupings'][g]
            line += f"   {o['dockq'][j]:.3f}#{o['rank'][j]:<2d}[{o['rank_lo'][j]}-{o['rank_hi'][j]}]"
        print(line)
    for g in GROUPINGS:
        o = res['groupings'][g]
        print(g, 'groups', o['n_groups'], 'kish', round(o['kish_ess'], 1), 'tau', round(o['kendall_tau_vs_raw'], 3),
              'P(rank1):', {names[j]: round(o['p_rank1'][j], 2) for j in range(M) if o['p_rank1'][j] > 0.01})
    print('raw score, ag50 cluster-boot P(rank1):',
          {names[j]: round(a50['rawscore_clusterboot_p_rank1'][j], 2) for j in range(M)
           if a50['rawscore_clusterboot_p_rank1'][j] > 0.01})


if __name__ == '__main__':
    main()
