# iGEM Claremont 2026 — Phage–Host Infectivity Prediction

## Project overview

We are building a machine-learning pipeline to predict which *Xanthomonas*-infecting bacteriophages will lyse which bacterial host strains. The pipeline outputs a **digital phagogram** — a ranked recommendation table — that guides wet-lab plaque-assay validation. This is Claremont Colleges' iGEM 2026 dry-lab contribution.

**Core engineer:** Alex Chen. **Team:** Sarah, Weitao, Olivia, Angela, Carol, Ryan, Leah.

Pipeline diagram: `docs/workflow_chart.jpeg`

---

## Pipeline at a glance

| Step | Folder | Description | Owner(s) | Status |
|------|--------|-------------|----------|--------|
| 1 | `01_data_ground_truth/` | Download phage & host genomes; build interaction matrix | Sarah, Weitao | Partial ✅ |
| 2 | `02_annotation/` | PHANOTATE (phages) + Prodigal (hosts) + pharokka | Weitao, Olivia | Partial 🟡 |
| 3 | `03_feature_weighting/` | 6 host–phage biophysical features + regression | All | Partial 🟡 |
| 4 | `04_protein_embedding/` | ESM-2 embeddings of high-priority proteins | TBD | Not started ⬜ |
| 5 | `05_predictive_modeling/` | XGBoost/CNN/RF classifier; KNN baseline already working | Sarah, TBD | Partial 🟡 |
| 6 | `06_in_silico_knockdown/` | fastISM sequence perturbation → sensitivity heatmaps | TBD | Not started ⬜ |
| 7 | `07_final_output/` | Assemble digital phagogram CSV | TBD | Not started ⬜ |

---

## Repository layout

```
iGEM_Claremont_2026/
├── CLAUDE.md                        ← this file
├── docs/                            ← bilingual guides, diagrams
├── 01_data_ground_truth/            ← NCBI download, interaction matrix
│   ├── inputs/                      ← manual seed lists
│   ├── processes/                   ← scripts (fetch_positive_pairs.py, fetch_negative_pairs.py, download_genomes.py)
│   └── outputs/                     ← canonical CSVs + FASTA pools
├── 02_annotation/                   ← PHANOTATE + Prodigal + pharokka
│   ├── processes/phage_phanotate/   ← PHANOTATE binary + batch_phanotate.py
│   ├── processes/host_prodigal/     ← batch_prodigal.py
│   └── outputs/                     ← per-accession .faa files, pharokka bundles
├── 03_feature_weighting/            ← 6-factor feature engineering + regression
│   ├── processes/factor{1-6}_*/    ← one sub-folder per factor (see README.md)
│   ├── processes/weighting_model/   ← XGBoost/linear regression aggregator
│   └── outputs/                     ← per-factor CSVs, normalized matrix, weights
├── 04_protein_embedding/            ← ESM-2 embeddings (planned)
├── 05_predictive_modeling/          ← KNN baseline (working) + future classifier
├── 06_in_silico_knockdown/          ← fastISM (planned)
├── 07_final_output/                 ← digital phagogram assembler (planned)
├── shared/                          ← cross-step utilities, conda env
│   └── env/environment.yml
├── archive/                         ← superseded scripts + historical task specs
└── docs/                            ← bilingual pipeline guides
```

Each step folder has its own `README.md` describing inputs, processes, outputs, and current gaps.

---

## Data flow contract

1. **`inputs/`** — read-only references to upstream step `outputs/` (or external seeds). Never write generated data here.
2. **`processes/`** — the only place with code. All scripts read from `inputs/` pointers or upstream `outputs/`; write to their step's `outputs/`.
3. **`outputs/`** — canonical products consumed by the next step. Large artifact trees are gitignored (see `.gitignore`); commit `MANIFEST.csv` files instead.

**Canonical files:**
- Interaction matrix: `01_data_ground_truth/outputs/interaction_matrix/final_interaction_matrix.csv`
- Phage accession list: `01_data_ground_truth/outputs/interaction_matrix/xanthomonas_phages_accession_list.csv`
- Feature weights: `03_feature_weighting/outputs/proteins_weights.csv` (once generated)

---

## Team ownership map

| Step | Factor | Owner | Key deliverable |
|------|--------|-------|----------------|
| 1 — Data | Interaction matrix | Sarah | `final_interaction_matrix.csv` |
| 1 — Data | Genome download | Weitao | `phage_genomes/xanthomonas_pool/` |
| 2 — Annotation (phage) | PHANOTATE | Weitao | `02_annotation/outputs/phage_proteins/<acc>.faa` |
| 2 — Annotation (host) | Prodigal | Olivia | `02_annotation/outputs/host_proteins/<acc>.faa` |
| 3 — Factor 1 | GC Content | Alex | `outputs/per_factor/factor1/f01_gc_content.csv` |
| 3 — Factor 2 | pI & Acidity | Olivia | `outputs/per_factor/factor2/f02_pI_acidity*.csv` |
| 3 — Factor 3 | Protein Length | Weitao | `outputs/per_factor/factor3/f03_length_ratio.csv` |
| 3 — Factor 4 | GRAVY | Sarah | `outputs/per_factor/factor4/f04_gravy.csv` |
| 3 — Factor 5 | Sequence Similarity | Angela | `outputs/per_factor/factor5/f05_similarity.csv` |
| 3 — Factor 6 | CAI | Carol | `outputs/per_factor/factor6/f06_cai.csv` |
| 3 — Aggregation | Regression model | Alex | `outputs/proteins_weights.csv` |
| 4 — Embedding | ESM-2 | TBD | `04_protein_embedding/outputs/*.npy` |
| 5 — Modeling (baseline) | KNN imputer | Sarah | `05_predictive_modeling/outputs/baseline_taxonomy_knn/` |
| 5 — Modeling (final) | XGBoost classifier | TBD | Probability matrix |
| 6 — Knockdown | fastISM | TBD | Sensitivity heatmaps |
| 7 — Phagogram | Assembler | TBD | `07_final_output/outputs/digital_phagogram.csv` |

---

## Conventions

- **No hard-coded absolute paths.** Use `pathlib.Path(__file__).parents[N]` to anchor at the repo root or step folder. Scripts already patched to follow this rule.
- **Output filenames:** use `f<NN>_<description>` namespace for factor outputs (e.g. `f02_pI_acidity_per_genome.csv`).
- **Code comments:** English only.
- **Docs:** bilingual EN + ZH (Traditional Chinese), following `docs/prodigal_manual.md` style.
- **Large artifacts:** gitignored. Add accession/filename to a `MANIFEST.csv` in the same `outputs/` folder so the set is reproducible.
- **Tool split:** PHANOTATE for phage ORF calling; Prodigal/pyrodigal for bacterial hosts. Never swap.

---

## Quick start

```bash
# 1. Create environment
conda env create -f shared/env/environment.yml
conda activate igem2026

# 2. Step 1 — build interaction matrix (positive pairs)
python 01_data_ground_truth/processes/fetch_positive_pairs.py

# 3. Step 1 — add negative pairs
python 01_data_ground_truth/processes/fetch_negative_pairs.py

# 4. Step 1 — download phage genomes
python 01_data_ground_truth/processes/download_genomes.py

# 5. Step 2 — annotate phage proteins (PHANOTATE must be installed)
python 02_annotation/processes/phage_phanotate/batch_phanotate.py

# 6. Step 2 — annotate host proteins (Prodigal)
python 02_annotation/processes/host_prodigal/batch_prodigal.py

# 7. Step 3 — compute pI & acidity (Factor 2, example)
python 03_feature_weighting/processes/factor2_pI_acidity/compute_pI_acidity.py \
    --phage-dir 02_annotation/outputs/phage_proteins \
    --out 03_feature_weighting/outputs/per_factor/factor2

# 8. Step 5 — run KNN baseline predictor
python 05_predictive_modeling/processes/baseline_taxonomy_knn/predictor.py
```

---

## Key files to read first (for new contributors)

| File | Why |
|------|-----|
| `03_feature_weighting/README.md` | Full 6-factor spec with biological rationale and computation formulas |
| `docs/integrated_pipeline_guide.md` | Bilingual step-by-step pipeline walkthrough |
| `01_data_ground_truth/outputs/interaction_matrix/final_interaction_matrix.csv` | The ground-truth data everything is built on |
| `05_predictive_modeling/processes/baseline_taxonomy_knn/predictor.py` | Working end-to-end example of the modeling approach |
