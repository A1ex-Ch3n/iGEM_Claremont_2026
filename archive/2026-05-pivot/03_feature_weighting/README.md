# Step 3 — Feature Weighting / 特徵加權

## Overview / 概述

This step extracts six quantitative features for each phage–host pair and combines them into a normalized feature matrix. The matrix feeds both a linear regression baseline (Step 5) and the final XGBoost classifier.

本步驟為每對噬菌體–宿主配對提取六個量化特徵，並合併為標準化特徵矩陣。該矩陣用於線性回歸基準模型（步驟五）及最終 XGBoost 分類器。

---

## Inputs / 輸入

| Source | Description |
|--------|-------------|
| `02_annotation/outputs/phage_proteins/<acc>.faa` | Phage predicted protein sequences (PHANOTATE) |
| `02_annotation/outputs/host_proteins/<acc>.faa` | Host predicted protein sequences (Prodigal) |
| `01_data_ground_truth/outputs/phage_genomes/` | Raw phage DNA (for GC content, CAI) |
| `01_data_ground_truth/outputs/interaction_matrix/final_interaction_matrix.csv` | Ground-truth affinity labels (y) |

---

## Factor Definitions / 因子定義

| # | Owner | Factor | Data Source | Biological Meaning | x Computation |
|---|-------|--------|------------|-------------------|---------------|
| 1 | **Alex** | GC Content | DNA genome | Genomic composition match; similar GC suggests co-evolution | Absolute difference: `\|GC_phage − GC_host\|` |
| 2 | **Olivia** | pI & Acidity | `.faa` (protein) | Charge distribution & acidic AA ratio; affects protein stability in host cytoplasm | Acidity diff: `\|Acid%_phage − Acid%_host\|` (D+E residue fraction) |
| 3 | **Weitao** | Protein Length (Size) | `.faa` (protein) | Average protein length; reflects genome compactness & replication strategy | Mean length ratio: `MeanLen_phage / MeanLen_host` |
| 4 | **Sarah** | GRAVY Hydrophobicity | `.faa` (protein) | Average hydrophobicity; reflects environmental adaptation (temp, membrane) | Absolute difference: `\|GRAVY_phage − GRAVY_host\|` |
| 5 | **Angela** | Sequence Similarity % | `.faa` (protein) | Sequence homology; detects HGT or auxiliary metabolic genes (AMGs) | Normalized ratio: `SharedGenes / TotalPhageGenes` |
| 6 | **Carol** | CAI (Codon Adaptation Index) | DNA + `.faa` | Measures how efficiently phage hijacks host translation machinery | Mean CAI of phage genes against host high-expression codon table |

---

## Factor Computation Details / 特徵計算細節

### Factor 1 — GC Content (Alex)
Compute GC percentage of each phage and host DNA genome. `x = |GC_phage − GC_host|`.

### Factor 2 — pI & Acidity (Olivia)
- **pI:** median isoelectric point of all proteins, computed via `BioPython.ProteinAnalysis` or `R Peptides`.
- **Acidity:** fraction of D (Aspartate) + E (Glutamate) residues across the proteome.
- Use `03_feature_weighting/processes/factor2_pI_acidity/compute_pI_acidity.py`.

### Factor 3 — Protein Length (Weitao)
Mean protein length from `.faa`. `x = MeanLen_phage / MeanLen_host`.

### Factor 4 — GRAVY (Sarah)
`GRAVY = Σ(residue_hydrophobicity) / sequence_length` (Kyte–Doolittle scale, via BioPython).
Plot phage vs host density to verify distribution shift before computing differences.

### Factor 5 — Similarity % (Angela)
All-vs-all BLAST (or MMseqs2) of phage proteins against host proteins.
`x = number_of_phage_proteins_with_hit / total_phage_proteins`.

### Factor 6 — CAI (Carol — key factor)
1. Extract host ribosomal protein genes as the high-expression reference set.
2. Build codon usage table from this set.
3. Compute CAI for each phage gene against the table; take the mean.
Tool: `cai` Python library or BioPython `CodonAdaptationIndex`.

---

## Regression Model / 迴歸分析

**Goal / 目標:** Find `f(x)` such that `y = ax + b`, where `y` is the infection probability (0–1) from the interaction matrix.

**Workflow / 建議流程:**
1. **Normalization / 標準化:** Z-score standardize all six `x` features before regression (GC is ~0.1-scale, protein length is ~100-scale).
2. **Scatter check / 散佈圖:** Each factor owner plots their `x` vs `y` to confirm linearity.
3. **Multiple linear regression** on the normalized feature matrix using `sklearn.linear_model.LinearRegression` or `statsmodels.OLS`.
4. Output `proteins_weights.csv` ranking factor importances for Step 4 (ESM-2 selection) and Step 5 (classifier).

---

## Outputs / 輸出

| File | Description |
|------|-------------|
| `outputs/per_factor/factor1/f01_gc_content.csv` | Per-pair GC difference |
| `outputs/per_factor/factor2/f02_pI_acidity_per_genome.csv` | Per-genome pI & acidity (Olivia's output, already exists) |
| `outputs/per_factor/factor3/f03_length_ratio.csv` | Per-pair mean-length ratio |
| `outputs/per_factor/factor4/f04_gravy.csv` | Per-pair GRAVY difference |
| `outputs/per_factor/factor5/f05_similarity.csv` | Per-pair similarity fraction |
| `outputs/per_factor/factor6/f06_cai.csv` | Per-pair mean CAI |
| `outputs/normalized_factor_matrix.csv` | All six factors, Z-score normalized |
| `outputs/proteins_weights.csv` | Regression coefficients / feature importances |

---

## Current Status / 目前進度

| Factor | Owner | Status |
|--------|-------|--------|
| 2 — pI & Acidity | Olivia | Partial — `compute_pI_acidity.py` exists; per-genome output computed; per-pair output pending |
| 1, 3, 4, 5, 6 | Alex/Weitao/Sarah/Angela/Carol | Not started — scripts need to be written |
| Weighting model / regression | Alex (core) | Not started |
