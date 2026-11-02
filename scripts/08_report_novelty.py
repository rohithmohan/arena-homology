#!/usr/bin/env python3
"""Build report_novelty.html (training-set novelty, clonotypes, fusion masking) from 07's outputs.

Set ARTIFACT_OUT=<path> to also write a variant without the document skeleton (for claude.ai artifacts).
"""
import csv, html, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DER = os.path.join(ROOT, 'data', 'derived')
R = json.load(open(os.path.join(DER, 'results_novelty.json')))
R0 = json.load(open(os.path.join(DER, 'results.json')))
CUT = list(csv.DictReader(open(os.path.join(ROOT, 'models_cutoffs.csv'))))
NARR = json.load(open(os.path.join(ROOT, 'scripts', 'novelty_narrative.json')))
esc = html.escape
names = R0['models']
mids = R0['model_ids']
name_of = dict(zip(mids, names))
STR = ['<30%', '30-50%', '50-90%', '>=90%']


def f3(x):
    return '–' if x is None or x != x else f'{x:.3f}'


def pct(x):
    return '–' if x is None or x != x else f'{x*100:.1f}'


# ---------- cutoff table
cut_rows = []
for r in CUT:
    src = r['source_url']
    link = f'<a href="{esc(src)}" target="_blank" rel="noopener">source</a>' if src.startswith('http') else esc(src)
    cut_rows.append(f"<tr><td>{esc(r['model'])}</td><td class='num'>{esc(r['cutoff_date'])}</td><td>{esc(r['cutoff_type'])}</td>"
                    f"<td>{esc(r['confidence'])}</td><td>{link}</td><td class='ev'>{esc(r['evidence'])}</td></tr>")

# ---------- strata
known = R['known']
order = sorted(known, key=lambda m: -R['strata'][m]['all']['mean'])
strata_rows = []
comp = []
for m in order:
    d = R['strata'][m]
    cells = ''.join(f"<td class='num'>{f3(d[s]['mean'])}<br><span class='mut'>n={d[s]['n']} · {pct(d[s]['acc'])}%</span></td>"
                    for s in STR)
    cx = d['seen complex']
    ncx = d['no seen complex']
    own = f" <span class='flag'>{d['own_entry_in_training']} targets in training</span>" if d['own_entry_in_training'] else ''
    strata_rows.append(f"<tr><td>{esc(name_of[m])}{own}</td><td class='num'>{esc(R['cutoffs'][m])}</td>"
                       f"<td class='num'>{f3(d['all']['mean'])}</td>{cells}"
                       f"<td class='num'>{f3(cx['mean'])}<br><span class='mut'>n={cx['n']}</span></td>"
                       f"<td class='num'>{f3(ncx['mean'])}<br><span class='mut'>n={ncx['n']}</span></td></tr>")
    comp.append(dict(model=name_of[m], n=[d[s]['n'] for s in STR], mean=[d[s]['mean'] for s in STR],
                     cx=cx['n']))

# ---------- leaderboards
def lb_rows(lb):
    out = []
    for o in lb['rows']:
        dr = o['full_rank'] - o['rank']
        arrow = f'<span class="up">▲{dr}</span>' if dr > 0 else (f'<span class="down">▼{-dr}</span>' if dr < 0 else '')
        out.append(f"<tr><td>{esc(o['model'])}</td><td class='num'>{o['n']}</td><td class='num'>{f3(o['mean'])}</td>"
                   f"<td class='num'>{f3(o['lo'])}–{f3(o['hi'])}</td><td class='num'>#{o['rank']} [{o['rlo']}–{o['rhi']}] {arrow}</td>"
                   f"<td class='num'>{o['p1']:.2f}</td><td class='num'>{pct(o['acc'])}</td>"
                   f"<td class='num'>{f3(o['full_mean'])} <span class='mut'>#{o['full_rank']}</span></td></tr>")
    return '\n'.join(out)


L = R['leaderboards']
chart_lb = {k: [dict(model=o['model'], mean=o['mean'], lo=o['lo'], hi=o['hi'], full=o['full_mean'], n=o['n'],
                     rank=o['rank'], rlo=o['rlo'], rhi=o['rhi'], p1=o['p1']) for o in L[k]['rows']]
            for k in ('own', 'common', 'common_pre2024', 'no_complex')}

# ---------- protenix
PE = R.get('protenix', {})
pe_order = ['novel at both (<50%)', 'antigen seen by both (>=50%)', 'antigen newly seen 2021->2025',
            'complex newly seen 2021->2025', 'target itself in 2025 training set']
pe_rows, pe_chart = [], []
for c in pe_order:
    v = PE.get(c, {})
    if not v or not v.get('n'):
        pe_rows.append(f"<tr><td>{esc(c)}</td><td class='num'>0</td><td colspan='5' class='mut'>no targets</td></tr>")
        continue
    did = f"{v['did']:+.3f} <span class='mut'>({v['did_lo']:+.3f} to {v['did_hi']:+.3f})</span>" if 'did' in v else '–'
    pe_rows.append(f"<tr><td>{esc(c)}</td><td class='num'>{v['n']}</td><td class='num'>{f3(v['s21'])}</td><td class='num'>{f3(v['s25'])}</td>"
                   f"<td class='num'>{v['d']:+.3f} <span class='mut'>({v['lo']:+.3f} to {v['hi']:+.3f})</span></td>"
                   f"<td class='num'>{pct(v['acc21'])} → {pct(v['acc25'])}</td><td class='num'>{did}</td></tr>")
    pe_chart.append(dict(label=c, n=v['n'], d=v['d'], lo=v['lo'], hi=v['hi']))

# ---------- mol1
MC = R.get('mol1', {})
mc_order = ['novel until target release', 'antigen first seen after 2021-09-30', 'complex first seen after 2021-09-30',
            'antigen seen before 2021-09-30', 'complex seen before 2021-09-30']
mc_rows, mc_chart = [], []
for c in mc_order:
    v = MC.get(c, {})
    if not v or not v.get('n'):
        mc_rows.append(f"<tr><td>{esc(c)}</td><td class='num'>0</td><td colspan='5' class='mut'>no targets</td></tr>")
        continue
    did = f"{v['did']:+.3f} <span class='mut'>({v['did_lo']:+.3f} to {v['did_hi']:+.3f})</span>" if 'did' in v else '–'
    lat = f"{v['later_rel']:+.3f} <span class='mut'>({v['later_lo']:+.3f} to {v['later_hi']:+.3f})</span>" if 'later_rel' in v else '–'
    mc_rows.append(f"<tr><td>{esc(c)}</td><td class='num'>{v['n']}</td><td class='num'>{f3(v['mol1'])}</td><td class='num'>{f3(v['peers'])}</td>"
                   f"<td class='num'>{v['rel']:+.3f} <span class='mut'>({v['lo']:+.3f} to {v['hi']:+.3f})</span></td>"
                   f"<td class='num'>{did}</td><td class='num'>{lat}</td></tr>")
    mc_chart.append(dict(label=c, n=v['n'], rel=v['rel'], lo=v['lo'], hi=v['hi'],
                         lrel=v.get('later_rel'), llo=v.get('later_lo'), lhi=v.get('later_hi')))

# ---------- clonotypes
CL = R['clonotype']
cl_rows = []
for lab in ['vh90 (old)', 'clonotype', 'ag50+vh90 (old joint)', 'ag50+clonotype']:
    v = CL['results'][lab]
    top = sorted(range(len(names)), key=lambda j: v['rank'][j])[:4]
    top_txt = ', '.join(f"{names[j]} {v['dockq'][j]:.3f}" for j in top)
    cl_rows.append(f"<tr><td>{esc(lab)}</td><td class='num'>{v['n_groups']}</td><td class='num'>{v['n_in_multi']}</td>"
                   f"<td class='num'>{v['largest']}</td><td class='num'>{v['tau']:.3f}</td><td class='num'>{v['tau_null_median']:.3f}</td>"
                   f"<td class='num'>{v['p_tau']:.2f}</td><td>{esc(top_txt)}</td></tr>")
vh_rows = ''.join(f"<tr><td>{esc(x['group'])}</td><td class='num'>{x['n']}</td><td class='num'>{x['n_clonotypes']}</td>"
                  f"<td class='num'>{x['n_ag50']}</td><td>{esc(x['vgenes'])}</td><td class='mono'>{esc(x['h3s'])}</td></tr>"
                  for x in CL['vh90_multi'])

# ---------- fusion
FU = R['fusion']
fm_rows = ''.join(f"<tr><td>{esc(x['chain'].split('|')[0])}</td><td>{esc(x['partners'])}</td>"
                  f"<td class='num'>{float(x['masked_fraction'])*100:.0f}%</td><td>{esc(x['action'].replace('_', ' '))}</td></tr>"
                  for x in FU['masking'])
fe_rows = ''.join(f"<tr><td>{esc(t)}</td><td class='num'>{v['groups_before']} → {v['groups_after']}</td>"
                  f"<td class='num'>{v['n_split']}</td><td class='num'>{v['n_merged']}</td><td class='num'>{v['tau_before_after']:.3f}</td>"
                  f"<td class='num'>{v['max_abs_score_change']:.4f}</td><td>{esc('; '.join(f'{a} {b}→{c}' for a, b, c in v['rank_changes']) or 'none')}</td></tr>"
                  for t, v in FU['effect'].items())
split_list = '; '.join(f"{a.split('_')[0]}–{b.split('_')[0]}" for v in FU['effect'].values() for a, b in v['split_pairs'])

DATA = dict(comp=comp, lb=chart_lb, pe=pe_chart, mc=mc_chart, strata=STR)

PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Arena Training-Set Novelty</title>
<meta name="claude:viewport-fit-sides" content="cover">
<style>
:root{color-scheme:light;--page:#f9f9f7;--surface:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--muted:#898781;--grid:#e1e0d9;--axis:#c3c2b7;
--border:rgba(11,11,11,.10);--s1:#2a78d6;--s2:#eb6834;--s3:#1baf7a;--neg:#e34948;--pos:#2a78d6;--mid:#f0efec;--good:#006300;--bad:#b42c2c;
--q1:#cde2fb;--q2:#86b6ef;--q3:#3987e5;--q4:#1c5cab;--flagbg:#fdf0e3;--flagink:#8a4a12}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){color-scheme:dark;--page:#0d0d0d;--surface:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--muted:#898781;--grid:#2c2c2a;--axis:#383835;
--border:rgba(255,255,255,.10);--s1:#3987e5;--s2:#d95926;--s3:#199e70;--neg:#e66767;--pos:#3987e5;--mid:#383835;--good:#0ca30c;--bad:#e66767;
--q1:#184f95;--q2:#256abf;--q3:#5598e7;--q4:#9ec5f4;--flagbg:#3a2a1a;--flagink:#f2b98a}}
:root[data-theme="dark"]{color-scheme:dark;--page:#0d0d0d;--surface:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--muted:#898781;--grid:#2c2c2a;--axis:#383835;
--border:rgba(255,255,255,.10);--s1:#3987e5;--s2:#d95926;--s3:#199e70;--neg:#e66767;--pos:#3987e5;--mid:#383835;--good:#0ca30c;--bad:#e66767;
--q1:#184f95;--q2:#256abf;--q3:#5598e7;--q4:#9ec5f4;--flagbg:#3a2a1a;--flagink:#f2b98a}
*{box-sizing:border-box}
body{margin:0;background:var(--page);color:var(--ink);font:15px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1080px;margin:0 auto;padding-block:32px 64px;padding-left:max(16px,env(safe-area-inset-left,0px));padding-right:max(16px,env(safe-area-inset-right,0px))}
h1{font-size:28px;line-height:1.2;margin:0 0 6px;font-weight:650;text-wrap:balance}
h2{font-size:19px;margin:44px 0 6px;font-weight:620;text-wrap:balance}
h3{font-size:15px;margin:22px 0 6px;font-weight:620}
p,li{color:var(--ink2);max-width:76ch}
p strong,li strong{color:var(--ink)}
a{color:var(--s1)}
.lede{font-size:16px;margin:0 0 20px}
.card{background:var(--surface);border:1px solid var(--border);border-radius:12px;padding:18px 18px 12px;margin:14px 0}
.answer{border-left:3px solid var(--s1);padding:2px 0 2px 16px;margin:18px 0}
.answer p{margin:8px 0}
.controls{display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin:4px 0 10px;font-size:13.5px;color:var(--ink2)}
select{max-width:100%;min-width:0;font:inherit;color:var(--ink);background:var(--surface);border:1px solid var(--axis);border-radius:8px;padding:4px 8px}
.legend{display:flex;flex-wrap:wrap;gap:16px;font-size:13px;color:var(--ink2);margin:2px 0 8px}
.legend span{display:inline-flex;align-items:center;gap:6px}
.key-dot{width:10px;height:10px;border-radius:50%;display:inline-block}
.key-rect{width:12px;height:10px;border-radius:2px;display:inline-block}
.scroll{overflow-x:auto;-webkit-overflow-scrolling:touch}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{padding:7px 8px;border-bottom:1px solid var(--grid);text-align:left;vertical-align:top}
th{color:var(--ink2);font-weight:600;font-size:12.5px;white-space:nowrap}
tbody td:first-child{white-space:nowrap}
td.num,th.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
td.ev{min-width:260px;color:var(--ink2);font-size:12.5px}
td.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px}
.mut{color:var(--muted)}.up{color:var(--good);font-weight:600}.down{color:var(--bad);font-weight:600}
.flag{display:inline-block;font-size:11px;font-weight:600;color:var(--flagink);background:var(--flagbg);border-radius:6px;padding:1px 6px;margin-left:6px}
.cap{font-size:12.5px;color:var(--muted);margin:6px 0 2px;max-width:90ch}
code{font-size:12.5px}
#tip{position:fixed;pointer-events:none;background:var(--surface);color:var(--ink);border:1px solid var(--border);border-radius:8px;
padding:8px 10px;font-size:12.5px;box-shadow:0 4px 16px rgba(0,0,0,.14);display:none;z-index:10;max-width:300px}
#tip .tv{font-weight:650;font-size:13.5px}#tip .tl{color:var(--ink2)}
.hit{cursor:default}.hit:focus{outline:none}.hit:focus-visible>.m{stroke:var(--ink);stroke-width:2}
</style></head><body><main>
<h1>How much of the CoFold Arena Ab–Ag leaderboard comes from targets models have effectively seen?</h1>
<p class="lede">For each model with a documented training cutoff, every one of the 422 targets is searched against the whole PDB as it stood before that cutoff. Two things are checked: whether a homologous antigen was already deposited, and whether the same antibody–antigen complex was (same antigen at ≥90% identity, plus a heavy chain with the same CDR-H3 and IGHV family). Companion to <code>report.html</code> (homology reweighting).</p>

<div class="answer">__SUMMARY__</div>

<h2>1 · Training cutoffs</h2>
<p>Models marked unknown are excluded from every cutoff-specific result below. Dates are PDB release-date cutoffs as stated (or inferred, where marked) by each model's authors.</p>
<div class="card scroll"><table><thead><tr><th>Model</th><th class="num">Cutoff</th><th>Type</th><th>Confidence</th><th>Source</th><th>Evidence</th></tr></thead>
<tbody>__CUT_ROWS__</tbody></table></div>

<h2>2 · Novelty strata per model</h2>
<p>Each target's antigen is placed in a stratum by its highest identity to any PDB entry released before the model's cutoff. Hits must cover ≥50% of the target antigen chain, with the target's own entry excluded and fusion partners masked. The bar shows how a model's 422 targets split across strata; the table gives mean DockQ and acceptable rate per stratum.</p>
<div class="card">
<div class="legend"><span><i class="key-rect" style="background:var(--q1)"></i>&lt;30%</span><span><i class="key-rect" style="background:var(--q2)"></i>30–50%</span><span><i class="key-rect" style="background:var(--q3)"></i>50–90%</span><span><i class="key-rect" style="background:var(--q4)"></i>≥90% antigen identity before cutoff</span></div>
<div id="chartS"></div>
<p class="cap">Hover a segment for its count and the model's mean DockQ in that stratum.</p>
</div>
<div class="card scroll"><table><thead><tr><th>Model</th><th class="num">Cutoff</th><th class="num">All</th><th class="num">&lt;30%</th><th class="num">30–50%</th><th class="num">50–90%</th><th class="num">≥90%</th><th class="num">Seen complex</th><th class="num">No seen complex</th></tr></thead>
<tbody>__STRATA_ROWS__</tbody></table></div>
<p class="cap">Cells: mean DockQ, then n and acceptable rate (DockQ ≥ 0.23). Rows are sorted by full-set mean. 95% CIs per stratum are in <code>novelty_strata.csv</code>.</p>

<h2>3 · Novel-only leaderboard</h2>
<p>Each model is scored only on targets whose antigen had &lt;50% identity to anything in the PDB before that model's cutoff. Because those sets differ by model, the common-subset views score every model on the same targets. CIs and rank intervals come from the paired cluster bootstrap over antigen-50% groups (B = 2,000).</p>
<div class="card">
<div class="controls"><label for="selL">Target set</label><select id="selL">
<option value="own">Each model's own novel set</option><option value="common">Common subset: novel for every model</option>
<option value="common_pre2024">Common subset, excluding Protenix-v1-20250630</option><option value="no_complex">Each model's set without a seen complex</option></select></div>
<div class="legend"><span><i class="key-dot" style="background:var(--s1)"></i>Full set (published mean)</span><span><i class="key-dot" style="background:var(--s2)"></i>Novel-only mean, 95% CI</span></div>
<div id="chartL"></div>
</div>
<h3>Each model's own novel set</h3>
<div class="card scroll"><table><thead><tr><th>Model</th><th class="num">n</th><th class="num">Mean DockQ</th><th class="num">95% CI</th><th class="num">Rank [95%]</th><th class="num">P(#1)</th><th class="num">Accept. %</th><th class="num">Full set</th></tr></thead>
<tbody>__LB_OWN__</tbody></table></div>
<p class="cap">__LB_OWN_CAP__</p>
<h3>Common subset (__N_COMMON__ targets novel for every model with a known cutoff)</h3>
<div class="card scroll"><table><thead><tr><th>Model</th><th class="num">n</th><th class="num">Mean DockQ</th><th class="num">95% CI</th><th class="num">Rank [95%]</th><th class="num">P(#1)</th><th class="num">Accept. %</th><th class="num">Full set</th></tr></thead>
<tbody>__LB_COMMON__</tbody></table></div>
<h3>Common subset without Protenix-v1-20250630 (__N_COMMON_PRE__ targets; that model is left out because it trained on some of them)</h3>
<div class="card scroll"><table><thead><tr><th>Model</th><th class="num">n</th><th class="num">Mean DockQ</th><th class="num">95% CI</th><th class="num">Rank [95%]</th><th class="num">P(#1)</th><th class="num">Accept. %</th><th class="num">Full set</th></tr></thead>
<tbody>__LB_COMMON_PRE__</tbody></table></div>
<p class="cap">Ranks are among models with a known cutoff only. ▲/▼ compare with the same models' full-set ranks.</p>

<h2>4 · Protenix-v1: same model, two cutoffs</h2>
<p>Protenix-v1 (2021-09-30 cutoff) and Protenix-v1-20250630 share an architecture and differ in training data, so their per-target difference isolates what later data adds. 163 targets were released before 2025-06-30. For those, the 2025 checkpoint's training set contained the target structure itself.</p>
<div class="card">
<div id="chartP"></div>
<p class="cap">Dots: mean per-target DockQ difference (2025 checkpoint − 2021 checkpoint) in each stratum, with 95% cluster-bootstrap CI. Strata are exclusive, assigned in priority order from the bottom row up.</p>
</div>
<div class="card scroll"><table><thead><tr><th>Stratum</th><th class="num">n</th><th class="num">v1 (2021)</th><th class="num">v1 (2025)</th><th class="num">Difference (95% CI)</th><th class="num">Accept. %</th><th class="num">vs novel stratum</th></tr></thead>
<tbody>__PE_ROWS__</tbody></table></div>
<p class="cap">__PE_CAP__</p>

<h2>5 · Mol-1 consistency check</h2>
<p>Mol-1's stated cutoff is 2021-09-30. If its training data match that date, it should have no particular advantage on targets whose antigen or complex first appeared in the PDB <em>after</em> 2021-09-30. This compares Mol-1 with the mean of the other models that share the 2021-09-30 cutoff (__PEERS__). As a positive control it also shows how much models with later cutoffs (__LATER__) gain on the same strata, which tells us whether an effect of this kind is detectable at all.</p>
<div class="card">
<div class="legend"><span><i class="key-dot" style="background:var(--s1)"></i>Mol-1 − same-cutoff peers</span><span><i class="key-dot" style="background:var(--s2)"></i>Later-cutoff models − same-cutoff peers (control)</span></div>
<div id="chartM"></div>
</div>
<div class="card scroll"><table><thead><tr><th>Stratum</th><th class="num">n</th><th class="num">Mol-1</th><th class="num">Peers</th><th class="num">Mol-1 − peers (95% CI)</th><th class="num">vs novel stratum</th><th class="num">Later-cutoff − peers</th></tr></thead>
<tbody>__MC_ROWS__</tbody></table></div>
<p class="cap">__MC_CAP__</p>

<h2>6 · Antibody redundancy, redone with clonotypes</h2>
<p>Whole-V-domain clustering at 90% identity is replaced with clonotype-style grouping: same IGHV gene (species-aware, from ANARCI germline assignment), the same CDR-H3 length, and ≥80% CDR-H3 identity, linked by single linkage. VHHs are included; the two shark VNARs have no IGHV assignment and stay singletons.</p>
<div class="card scroll"><table><thead><tr><th>Grouping</th><th class="num">Groups</th><th class="num">Targets in shared groups</th><th class="num">Largest</th><th class="num">Kendall τ vs raw</th><th class="num">τ null median</th><th class="num">P(τ<sub>null</sub> ≤ τ)</th><th>Top 4 (weighted DockQ)</th></tr></thead>
<tbody>__CL_ROWS__</tbody></table></div>
<p class="cap">__CL_CAP__</p>
<h3>What the old VH-90% clusters contained</h3>
<div class="card scroll"><table><thead><tr><th>VH-90% group</th><th class="num">n</th><th class="num">Clonotypes</th><th class="num">Antigen groups</th><th>IGHV genes</th><th>CDR-H3s</th></tr></thead>
<tbody>__VH_ROWS__</tbody></table></div>

<h2>7 · Fusion-partner masking</h2>
<p>Antigen chains that carry an MBP, BRIL/cyt b562, T4 lysozyme, GFP, DsRed/mCherry, IgG1 Fc or thioredoxin segment (≥50% identity over ≥50 residues) had that segment masked before re-clustering. A chain is left unmasked when the partner covers ≥80% of it, because then the partner is the antigen. The novelty search above also uses the masked chains.</p>
<div class="card scroll"><table><thead><tr><th>Target</th><th>Partner</th><th class="num">Chain covered</th><th>Action</th></tr></thead><tbody>__FM_ROWS__</tbody></table></div>
<div class="card scroll"><table><thead><tr><th>Grouping</th><th class="num">Groups before → after</th><th class="num">Pairs split</th><th class="num">Pairs merged</th><th class="num">τ (weighted, before vs after)</th><th class="num">Max |Δ score|</th><th>Rank changes</th></tr></thead><tbody>__FE_ROWS__</tbody></table></div>
<p class="cap">Pairs split by masking: __SPLIT__. The primary groupings in <code>report.html</code> are unchanged; this is a sensitivity check.</p>

<h2>Method and assumptions</h2>
<ul>__METHOD__</ul>
<p class="cap">Files: <code>models_cutoffs.csv</code>, <code>novelty_per_target_model.csv</code>, <code>novelty_strata.csv</code>, <code>novel_only_leaderboard.csv</code>, <code>protenix_natural_experiment.csv</code>, <code>mol1_consistency_check.csv</code>, <code>antibody_clonotypes.csv</code>, <code>clonotype_rescoring.csv</code>, <code>fusion_masking_effect.csv</code>. Scripts: <code>scripts/06_pdb_novelty.py</code>, <code>07_novelty_scoring.py</code>, <code>08_report_novelty.py</code>.</p>
</main>
<div id="tip" role="tooltip"></div>
<script>
const D = __DATA__;
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const NS = 'http://www.w3.org/2000/svg';
function el(tag, attrs, parent){const e=document.createElementNS(NS,tag);for(const k in attrs)e.setAttribute(k,attrs[k]);if(parent)parent.appendChild(e);return e;}
function txt(parent, x, y, s, attrs={}){const t=el('text',Object.assign({x,y,fill:css('--ink2'),'font-size':12},attrs),parent);t.textContent=s;return t;}
const tip=document.getElementById('tip');
function showTip(ev, rows){tip.replaceChildren();rows.forEach(([v,l],i)=>{const d=document.createElement('div');const b=document.createElement('span');b.className=i===0?'tv':'';b.textContent=v;d.appendChild(b);if(l){const s=document.createElement('span');s.className='tl';s.textContent=' '+l;d.appendChild(s);}tip.appendChild(d);});
 tip.style.display='block';const x=ev.clientX+14, y=ev.clientY+14;tip.style.left=Math.min(x, innerWidth-tip.offsetWidth-8)+'px';tip.style.top=Math.min(y, innerHeight-tip.offsetHeight-8)+'px';}
function hideTip(){tip.style.display='none';}
function bindTip(node, rows){node.setAttribute('tabindex','0');node.classList.add('hit');node.addEventListener('pointermove',e=>showTip(e,rows));node.addEventListener('pointerleave',hideTip);
 node.addEventListener('focus',()=>{const r=node.getBoundingClientRect();showTip({clientX:r.right,clientY:r.top},rows);});node.addEventListener('blur',hideTip);}
const f3=x=>(x==null||isNaN(x))?'–':x.toFixed(3), sg=x=>(x>=0?'+':'')+x.toFixed(3);

// stacked composition bars
function drawS(){const host=document.getElementById('chartS');host.replaceChildren();
 const W=Math.max(320,host.clientWidth), narrow=W<560, left=narrow?118:170, right=10, rowH=26, top=4, H=top+D.comp.length*rowH+26;
 const svg=el('svg',{width:W,height:H,viewBox:`0 0 ${W} ${H}`,role:'img','aria-label':'Targets per novelty stratum by model'},host);
 const X=v=>left+v/422*(W-left-right), cols=['--q1','--q2','--q3','--q4'];
 [0,100,200,300,422].forEach(t=>{el('line',{x1:X(t),x2:X(t),y1:top,y2:H-22,stroke:css('--grid'),'stroke-width':1},svg);txt(svg,X(t),H-6,String(t),{'text-anchor':'middle',fill:css('--muted'),'font-size':11});});
 D.comp.forEach((m,i)=>{const y=top+i*rowH+5, h=rowH-10;txt(svg,left-10,y+h/2+4,narrow&&m.model.length>16?m.model.slice(0,15)+'…':m.model,{'text-anchor':'end',fill:css('--ink'),'font-size':narrow?11:12.5});
  let acc=0;m.n.forEach((c,k)=>{if(!c)return;const x=X(acc), w=X(acc+c)-X(acc);const g=el('g',{},svg);
   el('rect',{class:'m',x:x,y:y,width:Math.max(0,w-2),height:h,rx:k===3||acc+c===422?3:0,fill:css(cols[k])},g);
   bindTip(g,[[`${c} targets`,`${D.strata[k]} antigen identity`],[`${m.model} mean DockQ ${f3(m.mean[k])}`]]);acc+=c;});});
}
// dot + CI generic
function dotCI(hostId, rows, opt){const host=document.getElementById(hostId);host.replaceChildren();
 const W=Math.max(320,host.clientWidth), narrow=W<560, left=narrow?150:(opt.left||220), right=narrow?12:20, rowH=opt.rowH||30, top=6, H=top+rows.length*rowH+30;
 const svg=el('svg',{width:W,height:H,viewBox:`0 0 ${W} ${H}`,role:'img','aria-label':opt.label},host);
 const [x0,x1]=opt.domain, X=v=>left+(v-x0)/(x1-x0)*(W-left-right);
 opt.ticks.forEach(t=>{const x=X(t);el('line',{x1:x,x2:x,y1:top,y2:H-24,stroke:css(t===0&&opt.zero?'--axis':'--grid'),'stroke-width':1},svg);txt(svg,x,H-8,opt.fmt(t),{'text-anchor':'middle',fill:css('--muted'),'font-size':11});});
 rows.forEach((r,i)=>{const y=top+i*rowH+rowH/2;let lab=r.label;if(narrow&&lab.length>24)lab=lab.slice(0,23)+'…';txt(svg,left-10,y+4,lab,{'text-anchor':'end',fill:css('--ink'),'font-size':narrow?11:12.5});
  r.series.forEach(s=>{if(s.v==null||isNaN(s.v))return;const yy=y+(s.dy||0);
   if(s.lo!=null&&!isNaN(s.lo))el('line',{x1:X(s.lo),x2:X(s.hi),y1:yy,y2:yy,stroke:css(s.c),'stroke-width':2,'stroke-linecap':'round'},svg);
   const g=el('g',{},svg);el('rect',{x:Math.min(X(s.lo??s.v),X(s.v))-6,y:yy-9,width:Math.abs(X(s.hi??s.v)-X(s.lo??s.v))+12,height:18,fill:'transparent'},g);
   el('circle',{class:'m',cx:X(s.v),cy:yy,r:4.5,fill:css(s.c),stroke:css('--surface'),'stroke-width':2},g);bindTip(g,s.tip);});});
}
function drawL(){const k=document.getElementById('selL').value;
 const rows=D.lb[k].map(o=>({label:o.model,series:[{v:o.full,c:'--s1',dy:-5,tip:[[f3(o.full),'full-set mean'],[o.model]]},
  {v:o.mean,lo:o.lo,hi:o.hi,c:'--s2',dy:5,tip:[[f3(o.mean),`novel-only, n=${o.n}`],[`95% CI ${f3(o.lo)}–${f3(o.hi)}`],[`rank #${o.rank} (95% ${o.rlo}–${o.rhi}) · P(#1) ${o.p1.toFixed(2)}`]]}]}));
 const lo=Math.min(...D.lb[k].map(o=>o.lo).filter(x=>!isNaN(x)),...D.lb[k].map(o=>o.full));const hi=Math.max(...D.lb[k].map(o=>o.hi).filter(x=>!isNaN(x)),...D.lb[k].map(o=>o.full));
 const a=Math.floor(lo*20)/20,b=Math.ceil(hi*20)/20;const ticks=[];for(let t=a;t<=b+1e-9;t+=0.05)ticks.push(+t.toFixed(2));
 dotCI('chartL',rows,{domain:[a,b],ticks,fmt:t=>t.toFixed(2),label:'Novel-only leaderboard',left:170});}
function sym(vals){const m=Math.max(0.05,...vals.filter(v=>v!=null&&!isNaN(v)).map(Math.abs));const step=[0.05,0.1,0.2,0.25,0.5].find(st=>Math.ceil(m/st)<=3)||0.5;const k=Math.ceil(m/step);const ticks=[];for(let i=-k;i<=k;i++)ticks.push(+(i*step).toFixed(3));return {d:[-k*step,k*step],ticks};}
function drawP(){const vals=D.pe.flatMap(r=>[r.lo,r.hi]);const {d,ticks}=sym(vals);
 dotCI('chartP',D.pe.map(r=>({label:`${r.label} (n=${r.n})`,series:[{v:r.d,lo:r.lo,hi:r.hi,c:'--s1',tip:[[sg(r.d),'DockQ, 2025 − 2021'],[`95% CI ${sg(r.lo)} to ${sg(r.hi)}`],[r.label]]}]})),
  {domain:d,ticks,fmt:t=>(t>0?'+':'')+t.toFixed(2),zero:true,label:'Protenix-v1 cutoff comparison',left:300});}
function drawM(){const vals=D.mc.flatMap(r=>[r.lo,r.hi,r.llo??0,r.lhi??0]);const {d,ticks}=sym(vals);
 dotCI('chartM',D.mc.map(r=>({label:`${r.label} (n=${r.n})`,series:[{v:r.rel,lo:r.lo,hi:r.hi,c:'--s1',dy:-5,tip:[[sg(r.rel),'Mol-1 − peers'],[`95% CI ${sg(r.lo)} to ${sg(r.hi)}`],[r.label]]},
  {v:r.lrel,lo:r.llo,hi:r.lhi,c:'--s2',dy:5,tip:[[r.lrel==null?'–':sg(r.lrel),'later-cutoff models − peers'],[r.llo==null?'':`95% CI ${sg(r.llo)} to ${sg(r.lhi)}`]]}]})),
  {domain:d,ticks,fmt:t=>(t>0?'+':'')+t.toFixed(2),zero:true,label:'Mol-1 consistency check',left:300});}
function drawAll(){drawS();drawL();drawP();drawM();}
document.getElementById('selL').addEventListener('change',drawL);
let rt;addEventListener('resize',()=>{clearTimeout(rt);rt=setTimeout(drawAll,150);});
matchMedia('(prefers-color-scheme: dark)').addEventListener('change',drawAll);
new MutationObserver(drawAll).observe(document.documentElement,{attributes:true,attributeFilter:['data-theme']});
drawAll();
</script></body></html>"""

uw = MC.get('_unique_wins', {}).get('Mol-1', {})
uw_txt = ', '.join(f"{v[0]}/{v[1]} {k}" for k, v in uw.items())
others_uw = '; '.join(f"{m}: " + ', '.join(f"{v[0]}/{v[1]}" for v in d.values()) for m, d in MC.get('_unique_wins', {}).items() if m != 'Mol-1')
fill = {'__TAU__': f"{L['own']['tau_vs_full']:.2f}", '__TAU_NOCX__': f"{L['no_complex']['tau_vs_full']:.2f}",
        '__OWNDID__': (f"{PE['_own_vs_after']['did']:+.3f} (95% CI {PE['_own_vs_after']['lo']:+.3f} to {PE['_own_vs_after']['hi']:+.3f})" if '_own_vs_after' in PE else '–'),
        '__ALLD__': (f"{PE['_all']['d']:+.3f} (95% CI {PE['_all']['lo']:+.3f} to {PE['_all']['hi']:+.3f})" if '_all' in PE else '–'),
        '__UW__': esc(uw_txt) + ' (the other top models, same strata: ' + esc(others_uw) + ')',
        '__FP__': f"{MC.get('_unique_win_fisher_p', float('nan')):.2f}"}
for k in ('lb_own_cap', 'pe_cap', 'mc_cap', 'cl_cap'):
    for a, b in fill.items():
        NARR[k] = NARR[k].replace(a, b)
rep = {
    '__SUMMARY__': NARR['summary'], '__CUT_ROWS__': '\n'.join(cut_rows), '__STRATA_ROWS__': '\n'.join(strata_rows),
    '__LB_OWN__': lb_rows(L['own']), '__LB_COMMON__': lb_rows(L['common']), '__LB_COMMON_PRE__': lb_rows(L['common_pre2024']),
    '__N_COMMON__': str(L['n_common']), '__N_COMMON_PRE__': str(L['n_common_pre']),
    '__LB_OWN_CAP__': NARR['lb_own_cap'], '__PE_ROWS__': '\n'.join(pe_rows), '__PE_CAP__': NARR['pe_cap'],
    '__PEERS__': esc(', '.join(MC.get('_peers', []))), '__LATER__': esc(', '.join(MC.get('_later', []))),
    '__MC_ROWS__': '\n'.join(mc_rows), '__MC_CAP__': NARR['mc_cap'], '__CL_ROWS__': '\n'.join(cl_rows),
    '__CL_CAP__': NARR['cl_cap'], '__VH_ROWS__': vh_rows, '__FM_ROWS__': fm_rows, '__FE_ROWS__': fe_rows,
    '__SPLIT__': esc(split_list or 'none'), '__METHOD__': NARR['method'],
    '__DATA__': json.dumps(DATA).replace('</', '<\\/'),
}
for k, v in rep.items():
    PAGE = PAGE.replace(k, v)
leftover = [k for k in rep if k in PAGE]
assert not leftover, leftover
open(os.path.join(ROOT, 'report_novelty.html'), 'w').write(PAGE)
out = os.environ.get('ARTIFACT_OUT')
if out:
    art = PAGE.split('<meta name="viewport" content="width=device-width,initial-scale=1">', 1)[1]
    open(out, 'w').write(art.replace('</head><body>', '', 1).replace('</body></html>', '').lstrip())
    print('wrote artifact variant', out)
print('wrote report_novelty.html', len(PAGE) // 1024, 'KB')
