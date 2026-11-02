# arena-homology

How much do near-duplicate targets affect the [CoFold Arena](https://cofoldarena.ai) antibody–antigen leaderboard?

The Arena already deduplicates by exact antigen UniProt ID. Related antigens still appear many times, though: influenza hemagglutinin strains, HLA class I alleles, SARS-CoV-2 RBD/S1 variants, HIV Env SOSIPs. This repo clusters the 422 targets by sequence homology with MMseqs2. It then re-scores all 18 leaderboard models with each cluster weighted equally and compares the result with the published leaderboard.

Open `report.html` in a browser for the full write-up with charts.

## Result in brief

- **Group counts:** the 422 targets fall into 297 / 315 / 355 antigen groups at 30% / 50% / 90% identity (coverage 0.8). The largest groups are flu HA (17), HLA class I (14), neuraminidase (11), SFTSV Gn (9) and HIV Env (5).
- **Who's on top:** this doesn't change. Mol-1 is #1 under every weighting, and the top four (Mol-1, OpenDDE, OpenDDE-AbAg, Protenix-v2) stay the top four.
- **Mol-1's lead:** it widens from 0.004 to 0.013 DockQ (antigen 50%). Mol-1's bootstrap P(rank 1) rises from 0.51 to 0.70 (0.82 at 30%). The Mol-1 vs OpenDDE gap is still not resolved (95% CI −0.020 to +0.047).
- **Where the movement is:** mid-table. IntelliFold-2-Flash goes up 3–4 places and Boltz-2 down 3. Kendall τ against the published ranking is 0.87 / 0.91 / 0.91 at 30 / 50 / 90%.
- **Score shift:** weighting raises almost every model's score (+0.010 on average at 30%), because the duplicated antigen families are harder than average.
- **Antibody-side duplication:** this hardly matters. Only 9 targets share both an antigen group and a 90% VH cluster.

Data snapshot: bundle `assets/index-D1t1x7qL.js`, fetched 2026-11-02.

## Running it

```bash
pip install -r requirements.txt       # numpy, scipy, anarci
# also needs: node, and HMMER's hmmscan on PATH (apt-get install hmmer / conda install -c bioconda hmmer)
./run_all.sh
```

`run_all.sh` re-fetches the Arena bundle (the asset hash is read from the homepage). It downloads the MMseqs2 release 18 static binary into `tools/` if `mmseqs` is not on PATH, or uses `$MMSEQS` if set. The scripts then run in order:

| Script | What it does |
|---|---|
| `scripts/01_extract.js` | Brace-matches the `Cd` (targets), `He` (per-target DockQ) and `Pd` (models) literals out of the bundle and evaluates them in an empty `vm` context. |
| `scripts/02_fetch_sequences.py` | Fetches SEQRES sequences from the RCSB GraphQL API (cached in `data/raw/rcsb/`), maps chains by **auth** asym ID, trims VH/VL to the V domain with ANARCI. |
| `scripts/03_cluster.py` | MMseqs2 `easy-cluster` on antigen chains (30/50/90%) and VH domains (90%), then union-find over targets. |
| `scripts/04_rescore.py` | Weighted scores, cluster bootstrap (B=2000), Kendall τ, permutation null, leave-one-group-out, Fab/VHH split. It first asserts that the published leaderboard numbers reproduce. |
| `scripts/05_report.py` | Builds `report.html`. Set `ARTIFACT_OUT=<path>` to also write a version without the document skeleton, for publishing as a claude.ai artifact. |

A re-run fetches the live site, so the results will change if the Arena adds targets or models. `data/derived/` holds the snapshot used for the committed results.

## Outputs

- `targets_with_clusters.csv`: one row per target. It holds the metadata, the RCSB antigen description, group ID, size and weight for each grouping, and every model's MSA and no-MSA DockQ.
- `leaderboard_comparison.csv`: one row per model. Raw score and rank with per-target and cluster-bootstrap CIs, then for each grouping: weighted DockQ, Δ, rank, rank change, CI, 95% rank interval, P(rank 1), and acceptable/medium/high rates.
- `leaderboard_long.csv`: the same data in long format.
- `logo_largest_groups.csv`: leave-one-group-out for the 10 largest groups of each grouping.
- `format_breakdown.csv`: Fab vs VHH, with group sizes recomputed within each subset.
- `data/derived/`: extracted Arena data, sequences (FASTA + `chain_sequences.tsv`), MMseqs2 cluster tables, `cluster_report.txt`, `fetch_report.txt` and `results.json`.

## Assumptions

- **Leaderboard definition:** `msa` scores, with `esmfold2_fast` excluded (it has no `msa` scores), giving 18 models. Acceptable / medium / high thresholds are DockQ ≥ 0.23 / 0.49 / 0.80.
- **PDB IDs:** 12PG, 46GH, 10ND and 13FL are ordinary 4-character IDs, not extended IDs. All 420 entries resolved and no target was dropped.
- **Chain IDs:** SAbDab chain IDs are auth IDs. Read as label IDs, 566 of 1,223 lookups fail or hit the wrong entity.
- **Which chains are clustered:** only the chains the Arena lists as antigen. For pMHC targets these include the MHC heavy chain, so different peptides on HLA-A2 group together.
- **Coverage 0.8 (bidirectional) under-merges sub-constructs.** HA1-only constructs separate from full-length HA0, and RBD from full spike. The extra `ag50_shortcov` grouping fixes this: it links two targets at ≥50% identity when the shorter chain is ≥80% covered.
- **Heavy chains** are trimmed to the V domain (ANARCI, IMGT) before clustering, so that Fab CH1 domains don't merge unrelated antibodies. Two shark VNARs fall back to their first 125 residues.
- **Joint group:** the (antigen-50% group, VH-90% cluster) pair.
- **Not done:** no epitope-level grouping, and no training-set overlap analysis.
