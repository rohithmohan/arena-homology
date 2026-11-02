# arena-homology

How much do near-duplicate targets affect the [CoFold Arena](https://cofoldarena.ai) antibody–antigen leaderboard?

The Arena already deduplicates by exact antigen UniProt ID. Related antigens still appear many times, though: influenza hemagglutinin strains, HLA class I alleles, SARS-CoV-2 RBD/S1 variants, HIV Env SOSIPs. This repo clusters the 422 targets by sequence homology with MMseqs2. It then re-scores all 18 leaderboard models with each cluster weighted equally and compares the result with the published leaderboard.

Open `report.html` in a browser for the full write-up with charts. A follow-up, `report_novelty.html`, asks how much of each model's score comes from targets it could have seen in training. It also redoes antibody redundancy with clonotypes and masks fusion partners.

## Result in brief

- **Group counts:** the 422 targets fall into 297 / 315 / 355 antigen groups at 30% / 50% / 90% identity (coverage 0.8). The largest groups are flu HA (17), HLA class I (14), neuraminidase (11), SFTSV Gn (9) and HIV Env (5).
- **Who's on top:** this doesn't change. Mol-1 is #1 under every weighting, and the top four (Mol-1, OpenDDE, OpenDDE-AbAg, Protenix-v2) stay the top four.
- **Mol-1's lead:** it widens from 0.004 to 0.013 DockQ (antigen 50%). Mol-1's bootstrap P(rank 1) rises from 0.51 to 0.70 (0.82 at 30%). The Mol-1 vs OpenDDE gap is still not resolved (95% CI −0.020 to +0.047).
- **Where the movement is:** mid-table. IntelliFold-2-Flash goes up 3–4 places and Boltz-2 down 3. Kendall τ against the published ranking is 0.87 / 0.91 / 0.91 at 30 / 50 / 90%.
- **Score shift:** weighting raises almost every model's score (+0.010 on average at 30%), because the duplicated antigen families are harder than average.
- **Antibody-side duplication:** this hardly matters. Only 9 targets share both an antigen group and a 90% VH cluster.

Data snapshot: bundle `assets/index-D1t1x7qL.js`, fetched 2026-11-02.

## Follow-up: training-set novelty (scripts 06–08)

- **Antigen overlap:** 76% of targets (320/422) had an antigen at ≥90% identity in the PDB before 2021-09-30, the cutoff 12 of the 18 models share. Only 17% were below 50%. Scores on familiar antigens are no higher than on novel ones.
- **Novel-only leaderboard:** re-scoring each model only on targets whose antigen was <50% identical to anything before its own cutoff keeps the same top four. Mol-1 stays #1 (rank interval 1–4, n = 71); the mid-table reshuffles (τ 0.75).
- **Seen complexes:** same antigen plus the same CDR-H3 and IGHV family. These are 17 targets (4%) for 2021-cutoff models, including rituximab–CD20, Nb35–Gs and the GFP enhancer nanobody. Top models score far higher on them, but removing them lowers any 2021 model's mean by ≤0.014 DockQ.
- **Protenix natural experiment:** Protenix-v1-20250630 vs Protenix-v1 (2021). Training on the exact target structure is worth +0.051 DockQ (CI +0.003 to +0.100) relative to targets released after the 2025 cutoff. The extra four years of other data add ≈0.
- **Mol-1 consistency check:** results are consistent with the stated 2021-09-30 cutoff in aggregate. Two individual targets are worth a look: KRas 8VGQ, and the RSV F Fab 5B11 in 8WZ4/8WZ5. There Mol-1 alone is accurate, and the complex was first deposited after 2021. Nine targets can't separate later exposure from capability.
- **Clonotypes:** grouping by same IGHV gene, CDR-H3 length and ≥80% CDR-H3 identity links only 13 targets, and the weighted ranking equals the published one. The old VH-90% clusters were mostly synthetic-library Fabs sharing a framework.
- **Fusion masking:** masking MBP, BRIL, T4L, GFP/RFP, Fc and thioredoxin splits only the two known fusion merges and changes no rank.

`models_cutoffs.csv` lists every model's PDB cutoff with a source URL and the supporting quote. Four are inferred: AF-Multimer (assumed v2.3), RF3-09/21 (no exact day given), OpenFold3-p2 and OpenDDE-AbAg.

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
| `scripts/06_pdb_novelty.py` | Downloads wwPDB `pdb_seqres` and all release dates, masks fusion partners, searches target antigens against the whole PDB (MMseqs2 `-s 7.5`), and runs ANARCI on target and PDB heavy chains for the seen-complex check. |
| `scripts/07_novelty_scoring.py` | Novelty strata, novel-only and common-subset leaderboards, the Protenix natural experiment, the Mol-1 consistency check, clonotype re-weighting, and fusion-masked re-clustering. |
| `scripts/08_report_novelty.py` | Builds `report_novelty.html`; narrative text lives in `scripts/novelty_narrative.json`. |
| `scripts/05_report.py` | Builds `report.html`. Set `ARTIFACT_OUT=<path>` to also write a version without the document skeleton, for publishing as a claude.ai artifact. |

A re-run fetches the live site, so the results will change if the Arena adds targets or models. `data/derived/` holds the snapshot used for the committed results.

## Outputs

- `targets_with_clusters.csv`: one row per target. It holds the metadata, the RCSB antigen description, group ID, size and weight for each grouping, and every model's MSA and no-MSA DockQ.
- `leaderboard_comparison.csv`: one row per model. Raw score and rank with per-target and cluster-bootstrap CIs, then for each grouping: weighted DockQ, Δ, rank, rank change, CI, 95% rank interval, P(rank 1), and acceptable/medium/high rates.
- `leaderboard_long.csv`: the same data in long format.
- `logo_largest_groups.csv`: leave-one-group-out for the 10 largest groups of each grouping.
- `format_breakdown.csv`: Fab vs VHH, with group sizes recomputed within each subset.
- Follow-up outputs: `models_cutoffs.csv`, `novelty_per_target_model.csv`, `novelty_strata.csv`, `novel_only_leaderboard.csv`, `protenix_natural_experiment.csv`, `mol1_consistency_check.csv`, `mol1_complexes_first_seen_after_cutoff.csv`, `antibody_clonotypes.csv`, `clonotype_rescoring.csv`, `fusion_masking_effect.csv`.
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
