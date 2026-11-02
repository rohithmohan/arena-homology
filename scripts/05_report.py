#!/usr/bin/env python3
"""Build arena_homology/report.html from the analysis outputs (self-contained, no external assets)."""
import csv, json, os, html
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DER = os.path.join(ROOT, 'data', 'derived')
R = json.load(open(os.path.join(DER, 'results.json')))
T = list(csv.DictReader(open(os.path.join(ROOT, 'targets_with_clusters.csv'))))
FETCH = open(os.path.join(DER, 'fetch_report.txt')).read()
names, mids = R['models'], R['model_ids']
G = R['groupings']
GLABEL = {'raw': 'Raw (published)', 'ag30': 'Antigen 30%', 'ag50': 'Antigen 50%', 'ag90': 'Antigen 90%',
          'joint': 'Antigen 50% + VH 90%', 'ag50_shortcov': 'Antigen 50%, shorter-chain cov.', 'vh90': 'VH 90% only'}
M = len(names)
raw_order = sorted(range(M), key=lambda j: G['raw']['rank'][j])
esc = html.escape


def f3(x):
    return f'{x:.3f}'


# ---------- data for charts ----------
models = []
for j in range(M):
    d = {'id': mids[j], 'name': names[j]}
    for g, o in G.items():
        d[g] = dict(s=o['dockq'][j], lo=o['ci_lo'][j], hi=o['ci_hi'][j], r=o['rank'][j], rlo=o['rank_lo'][j],
                    rhi=o['rank_hi'][j], p1=o['p_rank1'][j], acc=o['acc'][j])
    models.append(d)

buckets = [(1, 1, '1'), (2, 2, '2'), (3, 3, '3'), (4, 5, '4–5'), (6, 9, '6–9'), (10, 14, '10–14'), (15, 99, '15+')]
hist = {}
for g in ['ag30', 'ag50', 'ag90']:
    sizes = Counter(t[g + '_group'] for t in T)
    hist[g] = []
    for lo, hi, lab in buckets:
        ng = sum(1 for s in sizes.values() if lo <= s <= hi)
        nt = sum(s for s in sizes.values() if lo <= s <= hi)
        hist[g].append(dict(label=lab, groups=ng, targets=nt))

top6 = raw_order[:6]
logo = {}
for g in ['ag30', 'ag50', 'ag90']:
    logo[g] = [dict(group=e['group'], size=e['size'], antigens=e['antigens'],
                    d=[e['d_raw'][j] for j in top6], gm=[e['group_mean'][j] for j in top6],
                    after=[e['rank_raw_after'][j] for j in top6]) for e in R['logo'][g]]

# largest groups tables
largest = {}
for g in ['ag30', 'ag50', 'ag90', 'ag50_shortcov', 'joint']:
    by = defaultdict(list)
    for t in T:
        by[t[g + '_group']].append(t)
    rows = []
    for grp, mem in sorted(by.items(), key=lambda kv: (-len(kv[1]), kv[0]))[:15]:
        descs = Counter(d for t in mem for d in set(t['antigen_rcsb'].split(' | ')))
        fmt = Counter(t['format'] for t in mem)
        mean_by_model = [sum(float(t['msa_' + mids[j]]) for t in mem) / len(mem) for j in top6[:3]]
        rows.append(dict(group=grp, n=len(mem), fab=fmt.get('Fab', 0), vhh=fmt.get('VHH', 0),
                         antigens='; '.join(f'{d} ×{c}' if c > 1 else d for d, c in descs.most_common(5))
                                  + (f' (+{len(descs) - 5} more)' if len(descs) > 5 else ''),
                         pdbs=', '.join(t['pdb'] for t in mem[:8]) + (' …' if len(mem) > 8 else ''),
                         means=mean_by_model))
    largest[g] = rows

fmt = {}
for f in ['Fab', 'VHH']:
    o = R['format'][f]
    order = sorted(range(M), key=lambda j: o['raw']['rank'][j])[:8]
    fmt[f] = dict(n=o['n'], rows=[dict(name=names[j], **{g: dict(s=o[g]['dockq'][j], r=o[g]['rank'][j],
                                                                 rlo=o[g]['rank_lo'][j], rhi=o[g]['rank_hi'][j],
                                                                 p1=o[g]['p_rank1'][j])
                                                          for g in ['raw', 'ag30', 'ag50', 'ag90', 'joint']})
                                  for j in order],
                  groups={g: o[g]['n_groups'] for g in ['raw', 'ag30', 'ag50', 'ag90', 'joint']},
                  tau={g: o[g]['tau'] for g in ['ag30', 'ag50', 'ag90', 'joint']})

DATA = dict(models=models, raw_order=raw_order, hist=hist, logo=logo, top6=[names[j] for j in top6],
            glabel=GLABEL, largest=largest, fmt=fmt,
            pair=R['pairwise_top5'], top5=R['top5'])

# ---------- headline numbers ----------
mol = mids.index('mol1')
sd_t = sum(G['raw']['boot_sd']) / M
ratio = {g: sum(G[g]['rawscore_clusterboot_sd'][j] / G['raw']['boot_sd'][j] for j in range(M)) / M for g in ['ag30', 'ag50', 'ag90']}
wratio = {g: sum(G[g]['boot_sd'][j] / G['raw']['boot_sd'][j] for j in range(M)) / M for g in ['ag30', 'ag50', 'ag90', 'joint']}
shift = {g: sum(G[g]['dockq'][j] - G['raw']['dockq'][j] for j in range(M)) / M for g in ['ag30', 'ag50', 'ag90', 'joint']}
PN = R['perm_null']
pw = R['pairwise_top5']


def pair(g, a, b):
    d = pw[g][f'{a}-{b}']
    return f"{d['lo']:+.3f} to {d['hi']:+.3f}"


def rankchg(g):
    out = []
    for j in raw_order:
        dr = G['raw']['rank'][j] - G[g]['rank'][j]
        if dr:
            out.append(f"{esc(names[j])} {G['raw']['rank'][j]}→{G[g]['rank'][j]}")
    return ', '.join(out)


stats = [
    ('Targets', '422', '422 unique antigen UniProt IDs'),
    ('Groups at 30% / 50% / 90%', f"{G['ag30']['n_groups']} / {G['ag50']['n_groups']} / {G['ag90']['n_groups']}",
     f"{sum(1 for t in T if int(t['ag30_size'])>1)} targets share a 30% group with another"),
    ('Kendall τ vs raw (30/50/90)', f"{G['ag30']['kendall_tau_vs_raw']:.2f} / {G['ag50']['kendall_tau_vs_raw']:.2f} / {G['ag90']['kendall_tau_vs_raw']:.2f}",
     'all 18 models'),
    ('Mol-1 P(rank 1)', f"{G['raw']['p_rank1'][mol]:.2f} → {G['ag50']['p_rank1'][mol]:.2f}",
     f"raw → antigen-50% weighting ({G['ag30']['p_rank1'][mol]:.2f} at 30%)"),
]

# ---------- leaderboard table ----------
def cell(g, j):
    o = G[g]
    d = o['dockq'][j] - G['raw']['dockq'][j]
    dr = G['raw']['rank'][j] - o['rank'][j]
    arrow = f'<span class="up">▲{dr}</span>' if dr > 0 else (f'<span class="down">▼{-dr}</span>' if dr < 0 else '')
    return (f'<td class="num">{f3(o["dockq"][j])} <span class="mut">({d:+.3f})</span><br>'
            f'<span class="mut">#{o["rank"][j]} [{o["rank_lo"][j]}–{o["rank_hi"][j]}]</span> {arrow}</td>')


lb_rows = []
a50 = G['ag50']
for j in raw_order:
    r = G['raw']
    lb_rows.append(
        f'<tr><td>{esc(names[j])}</td>'
        f'<td class="num">{f3(r["dockq"][j])}<br><span class="mut">#{r["rank"][j]} [{r["rank_lo"][j]}–{r["rank_hi"][j]}]</span></td>'
        f'<td class="num">{f3(r["ci_lo"][j])}–{f3(r["ci_hi"][j])}<br><span class="mut">{f3(a50["rawscore_clusterboot_ci_lo"][j])}–{f3(a50["rawscore_clusterboot_ci_hi"][j])}</span></td>'
        f'<td class="num">{r["acc"][j]*100:.1f} → {G["ag50"]["acc"][j]*100:.1f}</td>'
        + ''.join(cell(g, j) for g in ['ag30', 'ag50', 'ag90', 'joint', 'ag50_shortcov']) + '</tr>')

tau_rows = ''.join(
    f'<tr><td>{esc(GLABEL[g])}</td><td class="num">{G[g]["n_groups"]}</td><td class="num">{G[g]["kish_ess"]:.0f}</td>'
    f'<td class="num">{G[g]["kendall_tau_vs_raw"]:.3f}</td>'
    f'<td class="num">{(PN[g]["tau_null_median"] if g in PN else float("nan")):.3f}</td>'
    f'<td class="num">{(PN[g]["p_tau_le_obs"] if g in PN else float("nan")):.3f}</td>'
    f'<td class="num">{sum(G[g]["dockq"][j]-G["raw"]["dockq"][j] for j in range(M))/M:+.4f}</td>'
    f'<td class="num">{sum(G[g]["boot_sd"][j]/G["raw"]["boot_sd"][j] for j in range(M))/M:.2f}×</td></tr>'
    for g in ['ag30', 'ag50', 'ag90', 'joint', 'ag50_shortcov', 'vh90']).replace('nan', '–')

fetch_lines = '\n'.join(esc(l) for l in FETCH.strip().splitlines())

PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Arena Homology Reweighting</title>
<meta name="claude:viewport-fit-sides" content="cover">
<style>
:root{color-scheme:light;--page:#f9f9f7;--surface:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--muted:#898781;--grid:#e1e0d9;--axis:#c3c2b7;
--border:rgba(11,11,11,.10);--s1:#2a78d6;--s2:#eb6834;--s3:#1baf7a;--neg:#e34948;--pos:#2a78d6;--mid:#f0efec;--good:#006300;--bad:#b42c2c}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){color-scheme:dark;--page:#0d0d0d;--surface:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--muted:#898781;--grid:#2c2c2a;--axis:#383835;
--border:rgba(255,255,255,.10);--s1:#3987e5;--s2:#d95926;--s3:#199e70;--neg:#e66767;--pos:#3987e5;--mid:#383835;--good:#0ca30c;--bad:#e66767}}
:root[data-theme="dark"]{color-scheme:dark;--page:#0d0d0d;--surface:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--muted:#898781;--grid:#2c2c2a;--axis:#383835;
--border:rgba(255,255,255,.10);--s1:#3987e5;--s2:#d95926;--s3:#199e70;--neg:#e66767;--pos:#3987e5;--mid:#383835;--good:#0ca30c;--bad:#e66767}
*{box-sizing:border-box}
body{margin:0;background:var(--page);color:var(--ink);font:15px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1080px;margin:0 auto;padding-block:32px 64px;padding-left:max(16px,env(safe-area-inset-left,0px));padding-right:max(16px,env(safe-area-inset-right,0px))}
h1{font-size:28px;line-height:1.2;margin:0 0 6px;font-weight:650}
h2{font-size:19px;margin:44px 0 6px;font-weight:620}
h3{font-size:15px;margin:22px 0 6px;font-weight:620}
p,li{color:var(--ink2);max-width:76ch}
p strong,li strong{color:var(--ink)}
.lede{font-size:16px;color:var(--ink2);margin:0 0 20px}
.card{background:var(--surface);border:1px solid var(--border);border-radius:12px;padding:18px 18px 12px;margin:14px 0}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:12px;margin:18px 0}
.tile{background:var(--surface);border:1px solid var(--border);border-radius:12px;padding:14px 16px}
.tile .l{font-size:13px;color:var(--ink2)}.tile .v{font-size:24px;font-weight:620;margin:2px 0}.tile .s{font-size:12.5px;color:var(--muted)}
.answer{border-left:3px solid var(--s1);padding:2px 0 2px 16px;margin:18px 0}
.answer p{margin:6px 0}
.controls{display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin:4px 0 10px;font-size:13.5px;color:var(--ink2)}
select{font:inherit;color:var(--ink);background:var(--surface);border:1px solid var(--axis);border-radius:8px;padding:4px 8px}
.legend{display:flex;flex-wrap:wrap;gap:16px;font-size:13px;color:var(--ink2);margin:2px 0 8px}
.legend span{display:inline-flex;align-items:center;gap:6px}
.key-dot{width:10px;height:10px;border-radius:50%;display:inline-block}
.key-rect{width:12px;height:10px;border-radius:2px;display:inline-block}
.key-line{width:16px;height:2px;display:inline-block}
.scroll{overflow-x:auto;-webkit-overflow-scrolling:touch}
svg text{font-family:inherit}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{padding:7px 8px;border-bottom:1px solid var(--grid);text-align:left;vertical-align:top}
th{color:var(--ink2);font-weight:600;font-size:12.5px;white-space:nowrap}
tbody td:first-child{white-space:nowrap}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
.mut{color:var(--muted)}.up{color:var(--good);font-weight:600}.down{color:var(--bad);font-weight:600}
.cap{font-size:12.5px;color:var(--muted);margin:6px 0 2px}
details{margin:10px 0}summary{cursor:pointer;color:var(--ink2);font-size:13.5px}
pre{white-space:pre-wrap;font-size:12px;color:var(--ink2);background:var(--surface);border:1px solid var(--border);border-radius:8px;padding:10px;overflow-x:auto}
code{font-size:12.5px}
#tip{position:fixed;pointer-events:none;background:var(--surface);color:var(--ink);border:1px solid var(--border);border-radius:8px;
padding:8px 10px;font-size:12.5px;box-shadow:0 4px 16px rgba(0,0,0,.14);display:none;z-index:10;max-width:300px}
#tip .tv{font-weight:650;font-size:13.5px}#tip .tl{color:var(--ink2)}
.heat td.hc{text-align:center;font-variant-numeric:tabular-nums;border:2px solid var(--surface);border-radius:4px}
.hit{cursor:default}
.hit:focus{outline:none}
.hit:focus-visible>.m{stroke:var(--ink);stroke-width:2}
</style></head><body><main>
<h1>Does near-duplicate antigen homology move the CoFold Arena Ab–Ag leaderboard?</h1>
<p class="lede">422 antibody–antigen targets, 18 leaderboard models, MSA-mode DockQ. Targets are grouped by antigen sequence homology (MMseqs2), each group gets equal weight, and the results are compared with the published, unweighted leaderboard. Data snapshot: bundle <code>index-D1t1x7qL.js</code>, fetched 2026-11-02.</p>

<div class="tiles">__TILES__</div>

<div class="answer">
<p><strong>Short answer: duplication doesn't change who's on top, but it does make the top more decisive.</strong> Mol-1 stays #1 under every weighting, and the top four (Mol-1, OpenDDE, OpenDDE-AbAg, Protenix-v2) stay the top four. Under antigen-50% weighting Mol-1's lead over OpenDDE grows from __LEAD_RAW__ to __LEAD_50__ DockQ. Its bootstrap probability of being #1 rises from __P1_RAW__ to __P1_50__ (__P1_30__ at 30%), because the largest homology groups (influenza HA, HLA class I) are ones where Mol-1 does relatively poorly. Even so, the Mol-1 vs OpenDDE gap is not resolved: its 95% CI is __PAIR_50__ at 50% and __PAIR_30__ at 30%.</p>
<p>Most of the effect sits mid-table, where models are 0.01 DockQ apart. IntelliFold-2-Flash climbs 3–4 places and Boltz-2 drops three; the other moves are one or two places among models whose rank intervals overlap. Weighting raises nearly every model's score (mean +__SHIFT30__ at 30%), because the over-represented antigen families are harder than average. Scores move by at most ~0.02 DockQ. The honest error bars widen by about __WR50__× (cluster bootstrap with weighting vs the naive per-target bootstrap).</p>
</div>

<h2>1 · Raw vs homology-weighted scores</h2>
<p>Each row is a leaderboard model, ordered by published rank. Dots are mean DockQ; whiskers are 95% bootstrap CIs. The raw CI comes from the ordinary per-target bootstrap; the weighted CI comes from a cluster bootstrap that resamples whole homology groups.</p>
<div class="card">
<div class="controls"><label for="selA">Weighting</label><select id="selA"></select></div>
<div class="legend" id="legA"></div>
<div id="chartA"></div>
<p class="cap">Hover a dot for the value, CI, rank interval and P(rank 1). All values are also in the leaderboard table below.</p>
</div>

<h2>2 · Rank changes</h2>
<p>Published rank (left) against rank under each weighting. Lines are coloured by the direction of the move under antigen-30% weighting, the strongest weighting. Grey lines don't move at 30%.</p>
<div class="card">
<div class="legend"><span><i class="key-line" style="background:var(--pos)"></i>moves up at 30%</span><span><i class="key-line" style="background:var(--neg)"></i>moves down at 30%</span><span><i class="key-line" style="background:var(--axis)"></i>unchanged at 30%</span></div>
<div class="scroll"><div id="chartB" style="min-width:640px"></div></div>
</div>

<h2>3 · How much duplication is there?</h2>
<p>The number of targets that sit in a homology group of a given size. Bars at size 1 are targets with no homologue in the set. Everything to the right is down-weighted.</p>
<div class="card">
<div class="legend"><span><i class="key-rect" style="background:var(--s1)"></i>Antigen 30%</span><span><i class="key-rect" style="background:var(--s2)"></i>Antigen 50%</span><span><i class="key-rect" style="background:var(--s3)"></i>Antigen 90%</span></div>
<div id="chartC"></div>
<p class="cap">Hover a bar for the number of groups and targets. The VH-90% clustering is much flatter: __VH_SINGLE__ of 422 heavy chains are singletons. The joint antigen + VH grouping has only __JOINT_MULTI__ targets in groups of size &gt;1.</p>
</div>

<h3>Largest groups</h3>
<div class="controls"><label for="selL">Grouping</label><select id="selL"></select></div>
<div class="card scroll"><table id="tblL"></table></div>
<p class="cap">Antigen names are RCSB entity descriptions of the Arena's antigen chains. The Arena's own <code>antigen</code> field sometimes holds a ligand name (e.g. “2-acetamido-2-deoxy-β-D-glucopyranose”), a SAbDab parsing quirk.</p>

<h2>4 · Which groups move the scores</h2>
<p>Leave-one-group-out on the published (unweighted) mean. Each cell shows how much a top-6 model's raw score changes, ×1000, when that group is removed. Blue means the model would gain, i.e. the group pulls it down. Red means the group props it up.</p>
<div class="card">
<div class="controls"><label for="selD">Grouping</label><select id="selD"></select></div>
<div class="legend"><span><i class="key-rect" style="background:var(--pos)"></i>score rises if group dropped</span><span><i class="key-rect" style="background:var(--neg)"></i>score falls if group dropped</span></div>
<div class="scroll"><table class="heat" id="tblD"></table></div>
<p class="cap">No single group shifts a top-6 model by more than ~0.006 DockQ, or any model by more than ~0.007 (0.010 for the merged 26-target HA group in the shorter-chain-coverage grouping). Dropping one group never changes who is #1. Dropping the neuraminidase group (6 Fabs, OpenDDE 0.42 vs OpenDDE-AbAg 0.21) swaps #2 and #3. The SARS-CoV-2 XBB.1.5 RBD/S1 group (Mol-1 0.06 vs OpenDDE-AbAg 0.73) and the HA and HLA groups move the top models in opposite directions. The full table for all 18 models is in <code>logo_largest_groups.csv</code>.</p>
</div>

<h2>5 · Fab vs VHH</h2>
<p>Group sizes are recomputed within each subset. Each cell shows score · rank [95% rank interval] · P(rank 1).</p>
<div class="controls"><label for="selF">Format</label><select id="selF"><option>Fab</option><option>VHH</option></select></div>
<div class="card scroll"><table id="tblF"></table></div>
<p class="cap" id="capF"></p>

<h2>6 · Full leaderboard comparison</h2>
<div class="card scroll"><table>
<thead><tr><th>Model</th><th class="num">Raw<br><span class="mut">#rank [95%]</span></th><th class="num">Raw 95% CI<br><span class="mut">target / ag50-cluster boot</span></th><th class="num">Accept. % <br><span class="mut">raw → ag50</span></th>
<th class="num">Antigen 30%</th><th class="num">Antigen 50%</th><th class="num">Antigen 90%</th><th class="num">Ag50 + VH90</th><th class="num">Ag50 short-cov</th></tr></thead>
<tbody>__LB_ROWS__</tbody></table></div>
<p class="cap">Weighted cells show weighted mean DockQ (Δ vs raw) · rank [cluster-bootstrap 95% rank interval], with ▲/▼ for rank change. Medium and high rates are in <code>leaderboard_comparison.csv</code>.</p>

<h3>Ranking agreement and the permutation null</h3>
<div class="card scroll"><table>
<thead><tr><th>Grouping</th><th class="num">Groups</th><th class="num">Kish ESS</th><th class="num">Kendall τ</th><th class="num">τ under null</th><th class="num">P(τ<sub>null</sub> ≤ τ)</th><th class="num">Mean Δ score</th><th class="num">CI width vs naive</th></tr></thead>
<tbody>__TAU_ROWS__</tbody></table></div>
<p class="cap">Null: 1,000 random reassignments of targets to groups with the same size distribution. It asks whether the reshuffle comes from homology itself or would happen under any reweighting this uneven. At 30% the observed τ is lower than almost every null draw. At 50% and 90% it is within the range of noise. No single model's relative gain under weighting is individually significant against the null, Mol-1's included. "CI width" is the mean bootstrap SD of the weighted score over the per-target bootstrap SD of the raw score. The raw score's own cluster-bootstrap SD is only __RATIO50__× the naive one at 50%, so the published CIs are only slightly overconfident.</p>

<h2>Method and assumptions</h2>
<ul>
<li><strong>Data.</strong> The <code>Cd</code> (targets), <code>He</code> (scores) and <code>Pd</code> (models) literals were brace-matched out of the bundle and evaluated in an empty <code>vm</code> context. Leaderboard = <code>msa</code> scores, <code>esmfold2_fast</code> excluded (it has no <code>msa</code> scores). Every published check reproduces exactly: Mol-1 0.359 / 48.8% acceptable, OpenDDE 0.355, Protenix-v2 0.327, Boltz-2 0.175, AF-Multimer 0.106.</li>
<li><strong>Sequences.</strong> Taken from the RCSB GraphQL API as SEQRES canonical sequences. All 420 PDB entries resolved, including 12PG, 46GH, 10ND and 13FL, which are ordinary 4-character IDs (new digit-leading issue) and not extended IDs. <strong>The chain IDs are auth IDs:</strong> read as label IDs, 309 of 1,223 lookups don't exist and 257 hit a different entity. Two PDB entries contribute two antibody targets each. No target was dropped.</li>
<li><strong>Antigen clustering.</strong> <code>mmseqs easy-cluster --min-seq-id {0.3,0.5,0.9} -c 0.8 --cov-mode 0 -s 7.5</code> on the unique antigen chains of each target (462 sequences), then union-find over targets. Only the chains the Arena lists as antigen are used. For pMHC targets that means the MHC heavy chain and peptide (sometimes β2m), so different peptides on HLA-A2 group together through the HLA chain.</li>
<li><strong>Coverage caveat.</strong> Bidirectional 80% coverage keeps HA1-only constructs and full-length HA0 apart, and likewise RBD vs ectodomain. The spec'd thresholds therefore <em>under</em>-merge sub-construct duplicates. The extra “shorter-chain coverage” grouping links two targets when identity ≥50% and the shorter chain is ≥80% covered (all-vs-all <code>easy-search</code>). That grouping merges 26 HA targets (HA1-only and full-length) into one group and gives results close to antigen-30%.</li>
<li><strong>Antibody clustering.</strong> V domains were trimmed with ANARCI (IMGT) before clustering at 90% identity, so that Fab CH1 domains don't merge unrelated Fabs. The two shark VNARs (not modelled by ANARCI) fall back to their first 125 residues. Joint group = identical (antigen-50% group, VH-90% cluster) pair. The VH clusters with mixed antigens (e.g. a 5-member cluster of synthetic-library Fabs) are why the joint grouping also requires a shared antigen group.</li>
<li><strong>Scoring.</strong> Weight = 1/group size. The weighted mean was asserted equal to the mean of group means for every grouping, metric and subset. The bootstrap uses B=2,000 and seed 20261102; every model is scored on the same resample, and rank intervals are 2.5–97.5 percentiles of the per-resample rank.</li>
<li><strong>Not done.</strong> No epitope-level clustering: two antibodies against different epitopes of the same HA still share a group. No training-set overlap analysis (which models saw homologues before their cutoffs). This analysis measures within-benchmark redundancy only.</li>
</ul>
<details><summary>Sequence fetch log</summary><pre>__FETCH__</pre></details>
<p class="cap">Files: <code>targets_with_clusters.csv</code>, <code>leaderboard_comparison.csv</code>, <code>leaderboard_long.csv</code>, <code>logo_largest_groups.csv</code>, <code>format_breakdown.csv</code>; scripts in <code>scripts/</code>, run with <code>run_all.sh</code>.</p>
</main>
<div id="tip" role="tooltip"></div>
<script>
const D = __DATA__;
const GL = D.glabel;
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const NS = 'http://www.w3.org/2000/svg';
function el(tag, attrs, parent){const e=document.createElementNS(NS,tag);for(const k in attrs)e.setAttribute(k,attrs[k]);if(parent)parent.appendChild(e);return e;}
function txt(parent, x, y, s, attrs={}){const t=el('text',Object.assign({x,y,fill:css('--ink2'),'font-size':12},attrs),parent);t.textContent=s;return t;}
const tip=document.getElementById('tip');
function showTip(ev, rows){tip.replaceChildren();rows.forEach(([v,l],i)=>{const d=document.createElement('div');const b=document.createElement('span');b.className=i===0?'tv':'';b.textContent=v;d.appendChild(b);if(l){const s=document.createElement('span');s.className='tl';s.textContent=' '+l;d.appendChild(s);}tip.appendChild(d);});
 tip.style.display='block';const r=ev.target.getBoundingClientRect?ev.target.getBoundingClientRect():{left:ev.clientX,top:ev.clientY};
 const x=(ev.clientX??r.left)+14, y=(ev.clientY??r.top)+14;tip.style.left=Math.min(x, innerWidth-tip.offsetWidth-8)+'px';tip.style.top=Math.min(y, innerHeight-tip.offsetHeight-8)+'px';}
function hideTip(){tip.style.display='none';}
function bindTip(node, rows){node.setAttribute('tabindex','0');node.classList.add('hit');node.addEventListener('pointermove',e=>showTip(e,rows));node.addEventListener('pointerleave',hideTip);
 node.addEventListener('focus',e=>{const r=node.getBoundingClientRect();showTip({clientX:r.right,clientY:r.top,target:node},rows);});node.addEventListener('blur',hideTip);}
function fill(sel, opts, def){opts.forEach(([v,l])=>{const o=document.createElement('option');o.value=v;o.textContent=l;if(v===def)o.selected=true;sel.appendChild(o);});}
const f3=x=>x.toFixed(3);

// ---------- A: dot + CI ----------
const selA=document.getElementById('selA');
fill(selA, ['ag30','ag50','ag90','joint','ag50_shortcov','vh90'].map(g=>[g,GL[g]]),'ag50');
function drawA(){
 const g=selA.value, host=document.getElementById('chartA'); host.replaceChildren();
 const leg=document.getElementById('legA'); leg.replaceChildren();
 [['--s1','Raw (per-target bootstrap CI)'],['--s2',GL[g]+' weighted (cluster bootstrap CI)']].forEach(([c,l])=>{const s=document.createElement('span');const k=document.createElement('i');k.className='key-dot';k.style.background=`var(${c})`;s.appendChild(k);s.appendChild(document.createTextNode(l));leg.appendChild(s);});
 const W=Math.max(320, host.clientWidth), narrow=W<560, left=narrow?118:160, right=narrow?14:70, rowH=30, top=8, H=top+D.raw_order.length*rowH+30;
 const svg=el('svg',{width:W,height:H,viewBox:`0 0 ${W} ${H}`,role:'img','aria-label':'Raw vs weighted DockQ by model'},host);
 const x0=0.08,x1=0.42, X=v=>left+(v-x0)/(x1-x0)*(W-left-right);
 for(let t=0.1;t<=0.401;t+=(narrow?0.1:0.05)){const x=X(t);el('line',{x1:x,x2:x,y1:top,y2:H-24,stroke:css('--grid'),'stroke-width':1},svg);txt(svg,x,H-8,t.toFixed(2),{'text-anchor':'middle',fill:css('--muted'),'font-size':11});}
 D.raw_order.forEach((j,i)=>{const m=D.models[j], y=top+i*rowH+rowH/2;
  txt(svg,left-10,y+4,(narrow&&m.name.length>16?m.name.slice(0,15)+'…':m.name),{'text-anchor':'end',fill:css('--ink'),'font-size':narrow?11:12.5});
  [['raw','--s1',-5],[g,'--s2',5]].forEach(([gg,c,dy])=>{const v=m[gg];const yy=y+dy;
   el('line',{x1:X(v.lo),x2:X(v.hi),y1:yy,y2:yy,stroke:css(c),'stroke-width':2,'stroke-linecap':'round'},svg);
   const grp=el('g',{},svg); el('rect',{x:X(v.lo)-4,y:yy-8,width:X(v.hi)-X(v.lo)+8,height:16,fill:'transparent'},grp);
   el('circle',{class:'m',cx:X(v.s),cy:yy,r:4.5,fill:css(c),stroke:css('--surface'),'stroke-width':2},grp);
   bindTip(grp,[[f3(v.s)+' DockQ',`${m.name} · ${gg==='raw'?'Raw':GL[gg]}`],[`95% CI ${f3(v.lo)}–${f3(v.hi)}`],[`rank #${v.r} (95% ${v.rlo}–${v.rhi}) · P(#1) ${v.p1.toFixed(2)}`],[`acceptable ${(v.acc*100).toFixed(1)}%`]]);});
  if(!narrow){const d=m[g].s-m.raw.s;txt(svg,W-right+8,y+4,(d>=0?'+':'')+d.toFixed(3),{fill:css('--muted'),'font-size':11.5});}
 });
 if(!narrow) txt(svg,W-right+8,H-8,'Δ',{fill:css('--muted'),'font-size':11});
}

// ---------- B: slope ----------
function drawB(){
 const host=document.getElementById('chartB'); host.replaceChildren();
 const cols=['raw','ag90','ag50','ag30','joint','ag50_shortcov'], short={raw:'Raw',ag90:'Ag 90%',ag50:'Ag 50%',ag30:'Ag 30%',joint:'Ag50+VH90',ag50_shortcov:'Ag50 short-cov'};
 const W=Math.max(640, host.clientWidth), left=172, right=172, top=30, rowH=24, n=D.models.length, H=top+n*rowH+10;
 const svg=el('svg',{width:W,height:H,viewBox:`0 0 ${W} ${H}`,role:'img','aria-label':'Rank change slope chart'},host);
 const X=i=>left+i*(W-left-right)/(cols.length-1), Y=r=>top+(r-1)*rowH+rowH/2;
 cols.forEach((c,i)=>{el('line',{x1:X(i),x2:X(i),y1:top-6,y2:H-10,stroke:css('--grid'),'stroke-width':1},svg);txt(svg,X(i),16,short[c],{'text-anchor':'middle',fill:css('--ink2'),'font-size':12,'font-weight':600});});
 const order=D.models.map((m,j)=>j).sort((a,b)=>Math.abs(D.models[a].raw.r-D.models[a].ag30.r)-Math.abs(D.models[b].raw.r-D.models[b].ag30.r));
 order.forEach(j=>{const m=D.models[j], d=m.raw.r-m.ag30.r, c=d>0?css('--pos'):d<0?css('--neg'):css('--axis');
  const g=el('g',{},svg);
  const pts=cols.map((cc,i)=>`${X(i)},${Y(m[cc].r)}`).join(' ');
  el('polyline',{points:pts,fill:'none',stroke:'transparent','stroke-width':12},g);
  el('polyline',{class:'m',points:pts,fill:'none',stroke:c,'stroke-width':2,'stroke-linejoin':'round','stroke-linecap':'round'},g);
  cols.forEach((cc,i)=>el('circle',{cx:X(i),cy:Y(m[cc].r),r:4,fill:c,stroke:css('--surface'),'stroke-width':2},g));
  txt(svg,X(0)-10,Y(m.raw.r)+4,`${m.raw.r}. ${m.name}`,{'text-anchor':'end',fill:css('--ink'),'font-size':12});
  txt(svg,X(cols.length-1)+10,Y(m.ag50_shortcov.r)+4,`${m.ag50_shortcov.r}. ${m.name}`,{fill:css('--ink2'),'font-size':12});
  bindTip(g,[[m.name,''],...cols.map(cc=>[`#${m[cc].r}`,`${short[cc]} · ${f3(m[cc].s)}`])]);
  g.addEventListener('pointerenter',()=>g.parentNode.appendChild(g));
 });
}

// ---------- C: histogram ----------
function drawC(){
 const host=document.getElementById('chartC'); host.replaceChildren();
 const keys=['ag30','ag50','ag90'], cols={ag30:'--s1',ag50:'--s2',ag90:'--s3'}, B=D.hist.ag30.map(b=>b.label);
 const W=Math.max(320,host.clientWidth), H=260, left=44, right=10, top=14, bottom=40;
 const svg=el('svg',{width:W,height:H,viewBox:`0 0 ${W} ${H}`,role:'img','aria-label':'Targets by group size'},host);
 const ymax=350, Y=v=>top+(1-v/ymax)*(H-top-bottom), band=(W-left-right)/B.length, bw=Math.min(18,(band-16)/3);
 for(let t=0;t<=ymax;t+=50){el('line',{x1:left,x2:W-right,y1:Y(t),y2:Y(t),stroke:css(t===0?'--axis':'--grid'),'stroke-width':1},svg);txt(svg,left-8,Y(t)+4,String(t),{'text-anchor':'end',fill:css('--muted'),'font-size':11});}
 txt(svg,left,H-6,'group size (targets per group)',{fill:css('--muted'),'font-size':11});
 B.forEach((lab,i)=>{const cx=left+band*i+band/2; txt(svg,cx,H-bottom+16,lab,{'text-anchor':'middle',fill:css('--ink2'),'font-size':11.5});
  keys.forEach((k,ki)=>{const b=D.hist[k][i], x=cx+(ki-1)*(bw+2)-bw/2, y=Y(b.targets), h=Y(0)-y;
   const g=el('g',{},svg); el('rect',{x:x-2,y:top,width:bw+4,height:Y(0)-top,fill:'transparent'},g);
   if(h>0){const r=Math.min(4,h);el('path',{class:'m',d:`M${x},${Y(0)} V${y+r} Q${x},${y} ${x+r},${y} H${x+bw-r} Q${x+bw},${y} ${x+bw},${y+r} V${Y(0)} Z`,fill:css(cols[k])},g);}
   bindTip(g,[[`${b.targets} targets`,`in ${b.groups} group${b.groups===1?'':'s'}`],[GL[k]+' · size '+lab]]);
   if(b.targets>0&&i>0) txt(svg,x+bw/2,y-4,String(b.targets),{'text-anchor':'middle',fill:css('--muted'),'font-size':9.5});
  });});
}

// ---------- largest groups table ----------
const selL=document.getElementById('selL');
fill(selL,['ag30','ag50','ag90','ag50_shortcov','joint'].map(g=>[g,GL[g]]),'ag50');
function drawL(){const t=document.getElementById('tblL');t.replaceChildren();const g=selL.value;
 const th=document.createElement('thead');const tr=document.createElement('tr');['Group','n','Fab/VHH','Antigens (RCSB)',...D.top6.slice(0,3).map(n=>n+' mean'),'PDB'].forEach((h,i)=>{const c=document.createElement('th');c.textContent=h;if(i===1||i===2||(i>=4&&i<=6))c.className='num';tr.appendChild(c);});th.appendChild(tr);t.appendChild(th);
 const tb=document.createElement('tbody');D.largest[g].forEach(r=>{const tr=document.createElement('tr');
  [[r.group],[String(r.n),'num'],[`${r.fab}/${r.vhh}`,'num'],[r.antigens],...r.means.map(v=>[f3(v),'num']),[r.pdbs,'mut']].forEach(([v,c])=>{const td=document.createElement('td');td.textContent=v;if(c)td.className=c;tr.appendChild(td);});tb.appendChild(tr);});t.appendChild(tb);}

// ---------- D: LOGO heat table ----------
const selD=document.getElementById('selD');
fill(selD,['ag30','ag50','ag90'].map(g=>[g,GL[g]]),'ag50');
function drawD(){const t=document.getElementById('tblD');t.replaceChildren();const g=selD.value;const rows=D.logo[g];
 const vmax=Math.max(...rows.flatMap(r=>r.d.map(Math.abs)));
 const th=document.createElement('thead');const tr=document.createElement('tr');['Group','n','Antigens',...D.top6].forEach((h,i)=>{const c=document.createElement('th');c.textContent=h;if(i===1)c.className='num';tr.appendChild(c);});th.appendChild(tr);t.appendChild(th);
 const tb=document.createElement('tbody');rows.forEach(r=>{const tr=document.createElement('tr');
  [[r.group],[String(r.size),'num'],[r.antigens]].forEach(([v,c])=>{const td=document.createElement('td');td.textContent=v;if(c)td.className=c;tr.appendChild(td);});
  r.d.forEach((v,k)=>{const td=document.createElement('td');td.className='hc';const a=Math.min(1,Math.abs(v)/vmax);
   td.style.background=`color-mix(in oklab, var(${v>=0?'--pos':'--neg'}) ${Math.round(a*85)}%, var(--mid))`;
   td.style.color=a>0.55?'#ffffff':css('--ink');td.textContent=(v>=0?'+':'')+(v*1000).toFixed(1);
   bindTip(td,[[`${(v>=0?'+':'')}${v.toFixed(4)} DockQ`,`if ${r.group} dropped`],[`${D.top6[k]} mean in group ${f3(r.gm[k])}`],[`raw rank after drop #${r.after[k]}`]]);tr.appendChild(td);});
  tb.appendChild(tr);});t.appendChild(tb);}

// ---------- F: format table ----------
const selF=document.getElementById('selF');
function drawF(){const f=selF.value, F=D.fmt[f], t=document.getElementById('tblF');t.replaceChildren();const gs=['raw','ag30','ag50','ag90','joint'];
 const th=document.createElement('thead');const tr=document.createElement('tr');['Model',...gs.map(g=>`${GL[g]} (${F.groups[g]} groups)`)].forEach((h,i)=>{const c=document.createElement('th');c.textContent=h;if(i)c.className='num';tr.appendChild(c);});th.appendChild(tr);t.appendChild(th);
 const tb=document.createElement('tbody');F.rows.forEach(r=>{const tr=document.createElement('tr');const td0=document.createElement('td');td0.textContent=r.name;tr.appendChild(td0);
  gs.forEach(g=>{const v=r[g];const td=document.createElement('td');td.className='num';td.textContent=`${f3(v.s)} · #${v.r} [${v.rlo}–${v.rhi}] · ${v.p1.toFixed(2)}`;tr.appendChild(td);});tb.appendChild(tr);});t.appendChild(tb);
 document.getElementById('capF').textContent=`${f}: ${F.n} targets. Kendall τ vs raw (all 18 models): `+['ag30','ag50','ag90','joint'].map(g=>`${GL[g]} ${F.tau[g].toFixed(3)}`).join(' · ')+'. Top 8 by raw rank shown; all 18 are in format_breakdown.csv.';}

function drawAll(){drawA();drawB();drawC();drawL();drawD();drawF();}
[selA,selL,selD,selF].forEach(s=>s.addEventListener('change',drawAll));
let rt;addEventListener('resize',()=>{clearTimeout(rt);rt=setTimeout(drawAll,150);});
matchMedia('(prefers-color-scheme: dark)').addEventListener('change',drawAll);
new MutationObserver(drawAll).observe(document.documentElement,{attributes:true,attributeFilter:['data-theme']});
drawAll();
</script></body></html>"""

op = mids.index('opendde')
tiles = ''.join(f'<div class="tile"><div class="l">{esc(l)}</div><div class="v">{v}</div><div class="s">{esc(s)}</div></div>'
                for l, v, s in stats)
vh_single = sum(1 for t in T if t['vh90_size'] == '1')
joint_multi = sum(1 for t in T if int(t['joint_size']) > 1)
rep = {
    '__TILES__': tiles,
    '__LEAD_RAW__': f"{G['raw']['dockq'][mol]-G['raw']['dockq'][op]:.3f}",
    '__LEAD_50__': f"{G['ag50']['dockq'][mol]-G['ag50']['dockq'][op]:.3f}",
    '__P1_RAW__': f"{G['raw']['p_rank1'][mol]:.2f}", '__P1_50__': f"{G['ag50']['p_rank1'][mol]:.2f}",
    '__P1_30__': f"{G['ag30']['p_rank1'][mol]:.2f}",
    '__PAIR_50__': pair('ag50', 'mol1', 'opendde'), '__PAIR_30__': pair('ag30', 'mol1', 'opendde'),
    '__SHIFT30__': f"{shift['ag30']:.3f}", '__WR50__': f"{wratio['ag50']:.2f}", '__RATIO50__': f"{ratio['ag50']:.2f}",
    '__VH_SINGLE__': str(vh_single), '__JOINT_MULTI__': str(joint_multi),
    '__LB_ROWS__': '\n'.join(lb_rows), '__TAU_ROWS__': tau_rows, '__FETCH__': fetch_lines,
    '__DATA__': json.dumps(DATA).replace('</', '<\\/'),
}
for k, v in rep.items():
    PAGE = PAGE.replace(k, v)
assert '__' not in PAGE.replace('__proto__', ''), [x for x in rep if x in PAGE]
open(os.path.join(ROOT, 'report.html'), 'w').write(PAGE)
art = PAGE.split('<meta name="viewport" content="width=device-width,initial-scale=1">', 1)[1]
art = art.replace('</head><body>', '', 1).replace('</body></html>', '')
out = os.environ.get('ARTIFACT_OUT')
if out:
    open(out, 'w').write(art.lstrip())
    print('wrote artifact variant', out)
print('wrote report.html', len(PAGE) // 1024, 'KB')
