# Step 2 — Standardized Annotation / 標準化基因組註釋

## Purpose / 目的

Predict open reading frames and protein sequences for all phage genomes (PHANOTATE) and bacterial host genomes (Prodigal). Also runs pharokka for richer phage functional annotation.

對所有噬菌體基因組（使用 PHANOTATE）及細菌宿主基因組（使用 Prodigal）進行開放閱讀框預測與蛋白質序列提取。同時透過 pharokka 進行更豐富的噬菌體功能性註釋。

---

## Inputs / 輸入

指向 `01_data_ground_truth/outputs/phage_genomes/`（噬菌體）與 `01_data_ground_truth/outputs/host_genomes/`（宿主）。

Points to `01_data_ground_truth/outputs/phage_genomes/` and `01_data_ground_truth/outputs/host_genomes/`.

---

## Processes / 流程 (`processes/`)

| Sub-folder / 子資料夾 | Tool | What it does / 功能 |
|-----------------------|------|---------------------|
| `phage_phanotate/batch_phanotate.py` | PHANOTATE | ORF prediction on phage FASTAs → coordinates + `.faa` / 噬菌體 FASTA 的 ORF 預測 → 座標 + 蛋白質序列 |
| `host_prodigal/batch_prodigal.py` | Prodigal (pyrodigal) | ORF prediction on host FASTAs → `proteins.faa` + `genes.gff` / 宿主 FASTA 的 ORF 預測 |
| `pharokka/` | pharokka | Full annotation bundle: PHROG hits, tRNA, CRISPR / 完整功能性註釋：PHROG 比對、tRNA、CRISPR |

**工具分工說明 / Tool split note:** PHANOTATE 專為噬菌體基因組設計（不假設終止密碼子偏好）；Prodigal/pyrodigal 用於細菌宿主。請勿互換使用。

PHANOTATE is tuned for phage genomes (no stop-codon bias assumptions); Prodigal/pyrodigal is for bacterial hosts. Do not swap them.

---

## Outputs / 輸出 (`outputs/`)

| Path | Description / 說明 |
|------|--------------------|
| `phage_proteins/<acc>.faa` | Translated phage protein sequences (gitignored) / 噬菌體蛋白質序列（已 gitignore） |
| `host_proteins/<acc>.faa` | Translated host protein sequences (gitignored) / 宿主蛋白質序列（已 gitignore） |
| `pharokka_runs/<acc>/` | Full pharokka output bundle per phage (gitignored) / 每個噬菌體的完整 pharokka 輸出（已 gitignore） |
| `phage_phanotate_coords/<acc>.txt` | Raw PHANOTATE coordinate files / 原始 PHANOTATE 座標檔 |

5 個樣本的 pharokka 結果（AB720063.2、AB720064.1、AP008979.1、EU717894.1、JN882298.1）已存於 `pharokka_runs/`（2025-04-28 註釋測試）。

---

## Owners / 負責人

Weitao（PHANOTATE / 噬菌體），Olivia（Prodigal / 宿主）

## Gaps / Next steps / 缺口與下一步

- `batch_phanotate.py` 的 `.faa` 翻譯迴圈未完成（截斷於第 44 行）。需補全：解析 PHANOTATE tab 格式座標 → Biopython 翻譯 → 寫出 `.faa`。
- `batch_phanotate.py` is incomplete — `.faa` translation loop truncated at line 44. Needs: parse PHANOTATE tab-delimited coords, translate via Biopython, write `.faa`.
- Scale pharokka to the full 777-genome pool / 將 pharokka 擴展至全部 777 個基因組。
- `archive/legacy_master_pipeline.py` 使用 Prodigal 處理噬菌體基因組 — 請勿使用，僅供參考。
