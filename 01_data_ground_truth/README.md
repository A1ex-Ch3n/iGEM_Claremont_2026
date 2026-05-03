# Step 1 — Data & Ground Truth / 數據與基準事實

## Purpose / 目的

Download phage and host genomes from NCBI and construct the phage–host interaction matrix (ground truth labels for the classifier).

從 NCBI 下載噬菌體與宿主基因組，並建立噬菌體–宿主感染互動矩陣（作為分類器的基準標籤）。

---

## Inputs / 輸入 (`inputs/`)

| File | Description / 說明 |
|------|-------------------|
| `complete_phage_data_Xanthomonas.csv` | Manual seed list of *Xanthomonas*-infecting phages (from Olivia) / 手動整理之*黃單胞菌*噬菌體種子列表（Olivia 提供） |
| `complete_phage_data_Pseudomonas.csv` | Pseudomonas phage seed list / 假單胞菌噬菌體種子列表 |

---

## Processes / 流程 (`processes/`)

| Script | What it does / 功能 |
|--------|---------------------|
| `fetch_positive_pairs.py` | NCBI Entrez search → positive affinity (=1) phage–host pairs for *Xanthomonas* phages / NCBI Entrez 搜尋 → 陽性感染配對（親和力=1） |
| `fetch_negative_pairs.py` | Generates Affinity=0 negatives via 3 modules: cross-genus, GenBank comment mining, pathovar inference / 透過三個模組生成陰性樣本：跨屬、GenBank 備注挖掘、致病型推斷 |
| `download_genomes.py` | Batch-downloads nucleotide FASTAs for all phage accessions / 批次下載噬菌體核苷酸 FASTA 序列 |

---

## Outputs / 輸出 (`outputs/`)

| Path | Description / 說明 |
|------|--------------------|
| `interaction_matrix/final_interaction_matrix.csv` | **Canonical** 2-D phage × host affinity matrix (1/0/NaN) / **標準** 二維噬菌體×宿主親和力矩陣 |
| `interaction_matrix/phage_host_matrix_with_ids.csv` | Long-format pairs with host RefSeq accessions / 長格式配對表（含宿主 RefSeq 登錄號） |
| `interaction_matrix/phage_host_matrix_raw.csv` | Raw long-format pairs (host names only) / 原始長格式配對（僅含宿主名稱） |
| `interaction_matrix/xanthomonas_phages_accession_list.csv` | Phage NCBI accession list / 噬菌體 NCBI 登錄號列表 |
| `negative_samples/negative_cross_genus.csv` | Module A cross-genus negatives / 模組 A 跨屬陰性樣本 |
| `negative_samples/negative_pv_inference.csv` | Module C pathovar-inference negatives / 模組 C 致病型推斷陰性樣本 |
| `negative_samples/negative_data_combined.csv` | All negatives combined / 所有陰性樣本合併 |
| `phage_genomes/xanthomonas_pool/` | 777 phage genome FASTAs (gitignored; track via MANIFEST) / 777 個噬菌體基因組 FASTA（已 gitignore；以 MANIFEST.csv 追蹤） |
| `host_genomes/ncbi_dataset/` | *Xanthomonas* host genome FASTAs (gitignored) / *黃單胞菌*宿主基因組（已 gitignore） |
| `download_summary.csv` | Download success/error log / 下載成功/失敗紀錄 |

---

## Owners / 負責人

Sarah（互動矩陣、陽性與陰性數據提取），Weitao（基因組下載）

## Gaps / Next steps / 缺口與下一步

- `fetch_negative_pairs.py` should write negatives directly into `final_interaction_matrix.csv` / 應直接將陰性樣本寫入最終互動矩陣
- Add `MANIFEST.csv` listing all 777 accessions in `phage_genomes/xanthomonas_pool/` for git tracking / 建立 MANIFEST.csv 追蹤 777 個登錄號
