# Step 7 — Final Output: Digital Phagogram / 最終輸出：數位噬菌體圖

## Purpose / 目的

Aggregate results from Steps 5 and 6 into a final digital phagogram CSV — a ranked recommendation table of phage–host pairs with prediction score, ML confidence, and binding sensitivity index. This table drives wet-lab prioritization for plaque-assay validation.

將步驟五與步驟六的結果整合為最終數位噬菌體圖 CSV——一個按優先級排序的噬菌體–宿主配對建議表，包含預測分數、機器學習置信度及結合敏感度指數。此表格用於指導濕實驗室噬菌斑實驗的優先安排。

---

## Inputs / 輸入

- `05_predictive_modeling/outputs/` — 預測機率矩陣與置信度分數 / predicted probability matrix, confidence scores
- `06_in_silico_knockdown/outputs/sensitivity_scores.csv` — 結合敏感度分數 / binding sensitivity scores

---

## Processes / 流程 (`processes/`)

- `build_phagogram.py` — 待撰寫 / to be written
- 合併預測結果與敏感度分析，輸出排序建議表 / Merges prediction and sensitivity into a ranked recommendation table

---

## Outputs / 輸出 (`outputs/`)

| File | Columns / 欄位 |
|------|----------------|
| `digital_phagogram.csv` | Phage_ID, Target_Strain, Prediction_Score, Predicted_MIC, Recommendation（High/Medium/Low Priority / Non-Infective）, HRM_Sensitivity_Index |

---

## Owner / 負責人

待定 / TBD

## Status / 進度

⬜ 尚未開始 / Not started.
