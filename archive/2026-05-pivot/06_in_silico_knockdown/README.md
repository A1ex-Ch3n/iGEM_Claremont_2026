# Step 6 — In-Silico Knockdown / 電腦模擬基因敲除（fastISM）

## Purpose / 目的

Use fastISM (In-Silico Mutagenesis) to perform sequence perturbation on the highest-probability phage–host pairs identified in Step 5. Generates a binding sensitivity heatmap showing which nucleotide positions most influence predicted infectivity.

對步驟五中預測機率最高的噬菌體–宿主配對，使用 fastISM（電腦模擬突變分析）進行序列擾動。輸出熱力圖，顯示哪些核苷酸位置對感染力預測影響最大。

---

## Inputs / 輸入

- 步驟五篩選出的高機率配對 / Top-scoring phage–host pairs from `05_predictive_modeling/outputs/`
- 對應的噬菌體蛋白質序列 / `02_annotation/outputs/phage_proteins/<acc>.faa`

---

## Processes / 流程 (`processes/`)

- `fastism_knockdown.py` — 待撰寫 / to be written
- 對每個高優先配對的關鍵蛋白質位置計算敏感度分數 / Computes per-position sensitivity score for high-priority pairs

---

## Outputs / 輸出 (`outputs/`)

| File | Description / 說明 |
|------|--------------------|
| `binding_sensitivity_heatmap_<phage>_<host>.png` | Heatmap of positional sensitivity / 位置敏感度熱力圖 |
| `sensitivity_scores.csv` | Per-position sensitivity scores / 每位置敏感度分數表 |

---

## Owner / 負責人

待定 / TBD

## Status / 進度

⬜ 尚未開始。需等待步驟五完成（取得最終機率矩陣以篩選高優先配對）。

Not started. Blocked on Step 5 (final probability matrix).
