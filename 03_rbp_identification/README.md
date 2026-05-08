# 03 — RBP Identification

Identify receptor-binding proteins (RBPs) — phage tail spikes / tail fibers — from annotated phage proteomes. Feeds Module 05 (structure prediction) and Module 04 (embedding).

## Method

**PhageRBPdetect** (Boeckaerts et al. 2022, *Viruses* 14:1329) — dual-track HMM + XGBoost classifier:
1. Pfam HMM scan against known RBP-related domains (tail fiber, tail spike, carbohydrate-binding, etc.)
2. ESM-2 embedding + XGBoost on Pfam-missed ORFs (catches RBPs that are too divergent for HMM)

Reported performance: precision-recall AUC 93.8%, F1 84.0%.

GitHub: https://github.com/dimiboeckaerts/PhageRBPdetection

## Inputs

`inputs/` reads from upstream `02_annotation/outputs/phage_proteins/<acc>.faa`.

## Outputs

`outputs/`
- `<phage_acc>_rbp_candidates.csv` — ranked RBP candidates per phage with `(orf_id, score, evidence_track)`
- `<phage_acc>_rbp_sequences.faa` — RBP candidate FASTA for downstream use

## Status

⏳ Not started — to be set up during dry lab sprint (5/7–5/17). First target: phiL7 (NCBI EU717894).
