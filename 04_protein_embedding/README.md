# Step 4 — Protein Embedding / 蛋白質嵌入向量（ESM-2）

## Purpose / 目的

Embed high-priority phage and host protein sequences using the ESM-2 transformer model. Embeddings replace raw sequence as input to the Step 5 classifier, capturing structural and evolutionary context.

使用 ESM-2 Transformer 模型將高優先級的噬菌體與宿主蛋白質序列轉換為數值嵌入向量。嵌入向量取代原始序列作為步驟五分類器的輸入，捕捉蛋白質結構與演化資訊。

---

## Inputs / 輸入

- `03_feature_weighting/outputs/proteins_weights.csv` — CRITICAL 與 HIGH 優先級蛋白質列表 / ranked protein priority list
- `02_annotation/outputs/phage_proteins/<acc>.faa` 與 `host_proteins/<acc>.faa`

---

## Processes / 流程 (`processes/`)

- ESM-2 包裝腳本（待撰寫）/ ESM-2 wrapper script (to be written)
- 建議模型 / Recommended model: `esm2_t33_650M_UR50D`（標準）或 `esm2_t36_3B_UR50D`（高精度，需 GPU）

---

## Outputs / 輸出 (`outputs/`)

| File | Description / 說明 |
|------|--------------------|
| `<protein_id>.npy` | Per-protein embedding vector (1280-dim or 2560-dim, gitignored) / 每個蛋白質的嵌入向量（已 gitignore） |
| `embedding_index.csv` | Maps accession + protein_id → `.npy` filename / 登錄號與蛋白質 ID 對應表 |

---

## Owner / 負責人

待定 / TBD

## Status / 進度

⬜ 尚未開始。需等待步驟三（蛋白質優先級排名）及步驟二（完整 `.faa` 資料集）完成後啟動。

Not started. Blocked on Step 3 (ranked protein list) and Step 2 (full `.faa` pool).
