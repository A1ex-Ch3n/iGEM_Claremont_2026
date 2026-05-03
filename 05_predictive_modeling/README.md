# Step 5 — Predictive Modeling

## Purpose
Predict phage–host infection probability. Three tracks run in parallel: (1) taxonomy-KNN baseline (already working), (2) multi-factor linear regression, (3) XGBoost/CNN/RF classifier on concatenated ESM-2 embeddings.

## Inputs
- `01_data_ground_truth/outputs/interaction_matrix/phage_host_matrix_with_ids.csv` — ground-truth matrix
- `03_feature_weighting/outputs/normalized_factor_matrix.csv` — for linear regression
- `04_protein_embedding/outputs/` — for XGBoost classifier

## Processes (`processes/`)
| Sub-folder | What it does | Status |
|------------|-------------|--------|
| `baseline_taxonomy_knn/predictor.py` | NCBI Taxonomy-weighted KNN imputer; fills NaN cells with predicted probability | Working |
| `linear_regression/` | 6-factor MLR baseline; uses `03_feature_weighting/` outputs | Not started |
| `xgboost_classifier/` | Final classifier on concatenated phage + host ESM-2 embeddings | Not started |

## Outputs (`outputs/`)
| Path | Description |
|------|-------------|
| `baseline_taxonomy_knn/original_matrix.csv` | Known interaction matrix |
| `baseline_taxonomy_knn/predicted_probability.csv` | KNN-imputed probability matrix |
| `baseline_taxonomy_knn/combined_matrix.csv` | Known + predicted combined view |
| `baseline_taxonomy_knn/confidence.csv` | Per-phage confidence (data coverage) |
| `baseline_taxonomy_knn/phage_host_results.xlsx` | Multi-sheet Excel summary |
| `baseline_taxonomy_knn/heatmap_*.png` | Original / predicted / combined heatmaps |

## Owners
Sarah (KNN baseline), TBD (linear regression, classifier)
