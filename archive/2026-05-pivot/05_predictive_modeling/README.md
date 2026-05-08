# Step 5 — Predictive Modeling / 預測模型

## Purpose / 目的

Predict phage–host infection probability using three parallel tracks: (1) taxonomy-KNN baseline (already working), (2) multi-factor linear regression, (3) XGBoost/CNN/RF classifier on concatenated ESM-2 embeddings.

透過三條並行路線預測噬菌體–宿主感染機率：（1）分類學 KNN 基準模型（已完成），（2）多因子線性回歸，（3）基於 ESM-2 嵌入向量的 XGBoost/CNN/RF 分類器。

---

## Inputs / 輸入

| Source | Description / 說明 |
|--------|--------------------|
| `01_data_ground_truth/outputs/interaction_matrix/phage_host_matrix_with_ids.csv` | Ground-truth matrix / 基準互動矩陣 |
| `03_feature_weighting/outputs/normalized_factor_matrix.csv` | For linear regression / 線性回歸用標準化特徵矩陣 |
| `04_protein_embedding/outputs/` | For XGBoost classifier / XGBoost 分類器用嵌入向量 |

---

## Processes / 流程 (`processes/`)

| Sub-folder / 子資料夾 | What it does / 功能 | Status / 狀態 |
|-----------------------|---------------------|---------------|
| `baseline_taxonomy_knn/predictor.py` | NCBI 分類學加權 KNN 插補器；填補 NaN 格為預測機率 / Taxonomy-weighted KNN imputer | ✅ Working |
| `linear_regression/` | 六因子多元線性回歸 / 6-factor MLR baseline | ⬜ Not started |
| `xgboost_classifier/` | 基於 ESM-2 嵌入向量拼接的最終分類器 / Final classifier on ESM-2 embeddings | ⬜ Not started |

---

## Outputs / 輸出 (`outputs/`)

| Path | Description / 說明 |
|------|--------------------|
| `baseline_taxonomy_knn/original_matrix.csv` | Known interaction matrix / 已知互動矩陣 |
| `baseline_taxonomy_knn/predicted_probability.csv` | KNN 預測機率矩陣 / KNN-imputed probability matrix |
| `baseline_taxonomy_knn/combined_matrix.csv` | Known + predicted combined / 已知與預測合併矩陣 |
| `baseline_taxonomy_knn/confidence.csv` | Per-phage confidence (data coverage) / 每個噬菌體的置信度 |
| `baseline_taxonomy_knn/phage_host_results.xlsx` | Multi-sheet Excel summary / 多工作表 Excel 彙整 |
| `baseline_taxonomy_knn/heatmap_*.png` | Original / predicted / combined heatmaps / 熱力圖 |

---

## Owners / 負責人

Sarah（KNN 基準模型），待定 / TBD（線性回歸、分類器）
