# Step 7 — Final Output: Digital Phagogram

## Purpose
Aggregate results from Steps 5 and 6 into a final digital phagogram CSV — a ranked recommendation table of phage–host pairs with prediction score, ML confidence, and binding sensitivity index. This table drives wet-lab prioritization for plaque-assay validation.

## Inputs
- `05_predictive_modeling/outputs/` — predicted probability matrix, confidence scores
- `06_in_silico_knockdown/outputs/sensitivity_scores.csv`

## Processes
- `processes/build_phagogram.py` — to be written; merges prediction and sensitivity into ranked table

## Outputs
| File | Columns |
|------|---------|
| `outputs/digital_phagogram.csv` | Phage_ID, Target_Strain, Prediction_Score, Predicted_MIC, Recommendation, HRM_Sensitivity_Index |

## Owner
TBD

## Status
Not started.
