# Archive

Contains files that have been superseded or reorganized but are preserved for reference.

## 2026-05 Pivot — pre-pivot artifacts

In May 2026 the project pivoted from a 6-factor weighting / KNN-baseline approach to a closed-loop active-learning pipeline (see `docs/iGEM_2026_Project_Plan.md` and `docs/project_pivot_summary_for_team.md` for the full rationale).

Everything in `2026-05-pivot/` was retired during the pivot. Git history preserves the full state — these copies live here so the new tree stays focused.

| Folder | Why archived |
|---|---|
| `2026-05-pivot/03_feature_weighting/` | 6 hand-crafted features (GC, pI, length, GRAVY, similarity, CAI) — replaced by ESM-2 embeddings + deep ensemble (modules `04_protein_embedding/` + `06_uncertainty_model/`) with substantially stronger literature support. |
| `2026-05-pivot/05_predictive_modeling/` | KNN-imputer baseline — replaced by deep ensemble in `06_uncertainty_model/`. |
| `2026-05-pivot/06_in_silico_knockdown/` | fastISM placeholder — replaced by **wet-lab** receptor knockout (pK18mobsacB) for causal validation. |
| `2026-05-pivot/07_final_output/` | Digital phagogram CSV placeholder — replaced by quantitative motif-level binding atlas (output of cycle iterations, see `08_cycle_data/`). |
| `2026-05-pivot/weitao/` | Weitao's working directory for host genome fetching + protein-size linear regression. Canonical genome data now lives in `00_raw_data/` and `01_data_ground_truth/`; the protein-length feature analysis is part of the retired 6-factor approach. |
| `2026-05-pivot/TASK/` | April 2026 task scratch directory — contained `0428_linear_regression/` (factor 3 spec, superseded) and `0428_phage_genomes/` (phage genome dump, superseded by `00_raw_data/phage_genomes/`). |

## Earlier archive entries

| File / Folder | Original location | Why archived |
|--------------|-------------------|-------------|
| `legacy_master_pipeline.py` | `olivia/04_dry_lab/pipeline/master_pipeline.py` | Mixed Step 1 (NCBI download) and Step 2 (Prodigal annotation on phage genomes). Phage annotation is now done by PHANOTATE (`02_annotation/processes/phage_phanotate/`); host annotation by Prodigal (`02_annotation/processes/host_prodigal/`). Download logic lives in `01_data_ground_truth/processes/download_genomes.py`. |
| `tasks/0428_linear_regression_cn_plan.md` | `TASK/0428_linear_regression/cn_plan.md` | Canonical 6-factor regression spec. Content was translated and expanded into the now-archived `2026-05-pivot/03_feature_weighting/README.md`. Kept here in its original Mandarin form. |
