# Alex 的 Project 自學指南
# 20 分鐘讀懂整個 pipeline 的 tech stack 和檔案地圖

> 這份文件是寫給你自己的。假設你有基本 Python 能力，但沒有正式學過 ML 或生物資訊學。
> 每個 section 約 2-3 分鐘閱讀。

---

## 先讀這三個檔案，其他的都可以之後再看

| 檔案 | 花多少時間 | 為什麼重要 |
|------|-----------|-----------|
| `CLAUDE.md` | 5 分鐘 | 整個 project 的規則書：資料夾結構、命名規則、notebook 規範 |
| `INTERFACE.md` | 10 分鐘 | 各 module 之間交換資料的格式合約——每個 module 生產什麼、消耗什麼 |
| `docs/pipeline_build_report_2026-05-08.md` | 20 分鐘 | 昨晚建出什麼、發現什麼問題、下一步是什麼（你剛剛生成的那份） |

---

## 整體架構再說一次（用檔案路徑說）

```
你的科學問題
    ↓
00_raw_data/                   777 個噬菌體 + 34 個細菌的 DNA 序列
    ↓
01_data_ground_truth/          哪些噬菌體已知能感染哪些細菌
    ↓
02_annotation/                 DNA → 蛋白質序列（gene calling）
    ↓
03_rbp_identification/         哪些蛋白質是 RBP（噬菌體的「鑰匙」）
    ↓
04_protein_embedding/          蛋白質序列 → 數字向量（讓 ML 能讀）
    ↓
05_structure_prediction/       RBP + 受體蛋白 → 3D 結構 + 結合力估計
    ↓
06_uncertainty_model/          訓練模型，讓它知道自己不確定的地方
    ↓
07_acquisition_function/       （下一個 sprint）決定做哪個實驗
    ↓
08_cycle_data/                 每個 cycle 的濕實驗室數據存放處
```

---

## Module 00 — Raw Data 驗證

### Tech Stack
| 工具 | 是什麼 | 為什麼用 |
|------|--------|---------|
| `pathlib` | Python 標準庫，處理檔案路徑 | 跨平台不用擔心斜線方向 |
| `hashlib` | Python 標準庫，算 SHA-256 | 驗證檔案完整性，確保沒有下載損壞 |
| `Bio.SeqIO` | BioPython 的 FASTA 解析器 | 讀取基因組序列，計算記錄數量 |
| `pytest` | Python 測試框架 | 自動化驗證所有斷言 |

### 要看的檔案
```
✅ 00_raw_data/AGENT_REPORT.md          結果摘要：哪些檔案有問題
✅ 00_raw_data/MANIFEST.csv             所有檔案的清單 + SHA-256（2398 行）
✅ 00_raw_data/processes/01_verify_dataset.ipynb   怎麼做的
⏭️ 00_raw_data/processes/tests/        可跳過，測試邏輯
```

### 你需要知道的一件事
`MANIFEST.csv` 是整個 pipeline 的「收據」——記錄每個原始輸入檔案的指紋。
未來如果有人問「你用的是哪個版本的 phiL7？」，答案就在這裡。

---

## Module 01 — Ground Truth 數據

### Tech Stack
| 工具 | 是什麼 | 為什麼用 |
|------|--------|---------|
| `pandas` | Python 資料分析庫 | 讀寫 CSV，做 join/merge |
| `NCBI datasets CLI` | NCBI 官方命令行工具 | 下載基因組（Bio.Entrez 因 SSL 失敗，改用這個）|
| `BioPython` | 生物資訊學 Python 庫 | 解析 FASTA 格式、計算基因組大小 |

### 要看的檔案
```
✅ 01_data_ground_truth/outputs/interaction_matrix.csv   最重要的輸出
   （2236 對：315 正例 + 1920 負例 + 1 個 phiL7×Xcc 已知對）
✅ 01_data_ground_truth/AGENT_REPORT.md                 為什麼 Bio.Entrez 失敗
⏭️ 01_data_ground_truth/processes/                      notebooks，之後再看
```

### 你需要知道的一件事
`interaction_matrix.csv` 裡的數據來自文獻——哪些噬菌體「已知」能感染哪些細菌。
這是 ML 模型的地基。它現在大部分是負例（不感染）因為負例數據更容易取得。

---

## Module 02 — 基因組注釋

### Tech Stack
| 工具 | 是什麼 | 為什麼用 |
|------|--------|---------|
| `PHANOTATE 1.6.7` | 噬菌體專用 ORF 預測工具 | 能找重疊的基因（噬菌體特有，一般工具會漏掉） |
| `pyrodigal 3.7.1` | Prodigal 的 Python 綁定 | 細菌基因預測，比直接呼叫 CLI 快 |
| `GFF3 格式` | 基因座標的標準格式 | 記錄每個基因的位置（染色體、起點、終點、正負鏈）|
| `BioPython` | | 翻譯 DNA → 蛋白質序列 |

### 要看的檔案
```
✅ 02_annotation/outputs/phage_proteins/EU717894.1.faa     phiL7 的 80 個蛋白質
✅ 02_annotation/outputs/host_proteins/NZ_CP155948.1.faa   Xcc 的 4344 個蛋白質
✅ 02_annotation/outputs/annotation_summary.csv            一行一個基因組的統計
✅ 02_annotation/processes/01_run_phanotate.ipynb          怎麼跑的（雙語 notebook）
⏭️ 02_annotation/outputs/phage_orfs/                      GFF3 格式，暫時不用看
```

### 為什麼要兩個工具？
PHANOTATE 懂噬菌體的「重疊讀碼框」（overlapping ORFs）——噬菌體為了把更多基因塞進小基因組，
基因會彼此重疊。一般工具（Prodigal）假設基因不重疊，所以噬菌體上會漏掉約 15% 的基因。
這是一個常見錯誤，CLAUDE.md 裡明確禁止混用這兩個工具。

---

## Module 03 — RBP 鑑定

### Tech Stack
| 工具 | 是什麼 | 為什麼用 |
|------|--------|---------|
| `PhageRBPdetect v2` | 噬菌體 RBP 偵測工具（Boeckaerts 2022）| 專門找「受體結合蛋白」 |
| `HMMER 3.4 (hmmscan)` | 蛋白質 domain 搜尋引擎 | 用 profile HMM 比對已知 RBP domain |
| `RBPdetect_phageRBPs.hmm` | 92 個 RBP domain 的 profile 資料庫 | PhageRBPdetect 自帶的 HMM profiles |
| `XGBoost` | 梯度提升樹分類器 | ML track：沒有 domain 的 RBP 用這個找（目前被阻塞）|

### 要看的檔案
```
✅ 03_rbp_identification/outputs/EU717894.1_rbp_candidates.csv   58 個蛋白質的評分
✅ 03_rbp_identification/outputs/EU717894.1_rbp_sequences.faa    前 3 名的序列
✅ 03_rbp_identification/AGENT_REPORT.md                         ML track 為什麼失敗
✅ 03_rbp_identification/processes/01_run_phagerbpdetect.ipynb   雙語 notebook
⏭️ 03_rbp_identification/inputs/phagerbpdetect_data/            HMM 資料庫，不用看
```

### 你需要知道的一件事
**rbp_01（712 aa）幾乎確定是尾刺蛋白（tail spike）**——這是 phiL7 用來認識 Xcc 的蛋白質。
它命中了 `Tail_spike_N` 和 `unknown_C54` 兩個 domain。這個蛋白是整個 wet lab 實驗的主角。

---

## Module 04 — 蛋白質 Embedding

### Tech Stack
| 工具 | 是什麼 | 為什麼用 |
|------|--------|---------|
| `fair-esm` | Meta 的 ESM 官方 Python 庫 | 載入 ESM-2 模型 |
| `ESM-2 8M (esm2_t6_8M_UR50D)` | 蛋白質語言模型，8M 參數版 | CPU 上幾秒跑完，smoke test 用 |
| `PyTorch` | 深度學習框架 | ESM-2 的底層，也用來做 tensor 操作 |
| `numpy` | 數值計算庫 | 儲存 embedding 陣列（.npz 格式）|
| `BioPython` | | 讀取 FASTA 輸入序列 |

### 要看的檔案
```
✅ 04_protein_embedding/outputs/embeddings_esm2_t6_8M_phiL7_rbps.npz     RBP 嵌入
✅ 04_protein_embedding/outputs/embeddings_esm2_t6_8M_xcc_receptors.npz  受體嵌入
✅ 04_protein_embedding/outputs/embedding_index.csv                       索引
✅ 04_protein_embedding/processes/01_embed_esm2.ipynb                     主 notebook
⚠️ 04_protein_embedding/inputs/mock_phiL7_rbp_sequences.faa              ← 假序列！
   需要換成 03 的真實輸出再重跑
```

### Embedding 是什麼意思？
把蛋白質序列（一串英文字母：MKLLILT...）轉成一排數字（320 個浮點數）。
這排數字捕捉了蛋白質的生物學意義——功能相似的蛋白質，數字也接近。
就像 word2vec 把「king」和「queen」放在向量空間裡的相近位置一樣。

### ⚠️ 現在需要做的事
Module 04 當時用的是假序列（mock）。現在 Module 03 的真實 RBP 序列已經在：
`03_rbp_identification/outputs/EU717894.1_rbp_sequences.faa`
要打開 `01_embed_esm2.ipynb` 的 Cell 8，把路徑換掉，重跑一次。

---

## Module 05 — 結構預測

### Tech Stack
| 工具 | 是什麼 | 為什麼用 |
|------|--------|---------|
| `Boltz-2 2.0.3` | 蛋白質複合體結構預測（Passaro 2025）| 開源，能預測結構 + 嘗試估計結合力 |
| `pytorch-lightning` | PyTorch 的高層封裝 | Boltz-2 的依賴 |
| `MMSeqs2`（線上）| 多序列比對工具 | 找進化上相關的蛋白質，幫助結構預測 |
| `AlphaFold 3`（未安裝）| Google 的結構預測（Abramson 2024）| 品質更高，但需要 Google 審核才能用 |
| `einops`, `dm-tree` | 矩陣操作工具 | Boltz-2 的依賴 |

### 要看的檔案
```
✅ 05_structure_prediction/AGENT_REPORT.md                      關鍵發現
✅ 05_structure_prediction/outputs/affinity_priors.csv           結合力表（目前只有一行，失敗的）
✅ 05_structure_prediction/inputs/boltz_input_*.fasta            已準備好的輸入
✅ 05_structure_prediction/processes/01_run_boltz2.ipynb         主 notebook
⏭️ 05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__*/msa/   MSA 中間產物，不用看
```

### 最重要的一個發現
**Boltz-2 的親和力預測頭（affinity head）只支援「蛋白質 + 小分子」，不支援「蛋白質 + 蛋白質」。**
我們的 RBP + 受體是兩個蛋白質，所以 `predicted_dG` 永遠是 NaN。
替代方案：用 `confidence`（ipTM 分數，0-1）作為結合品質的代理指標。
這個差異很重要，Module 06 已經考慮進去了。

---

## Module 06 — 不確定性模型

### Tech Stack
| 工具 | 是什麼 | 為什麼用 |
|------|--------|---------|
| `PyTorch` | 深度學習框架 | 建構 MLP、訓練、儲存模型 |
| `numpy` | 數值計算 | 資料矩陣操作 |
| `pandas` | | 讀寫 CSV 輸入輸出 |
| `scikit-learn` | ML 工具庫 | train/val split、校準指標（ECE）|
| `matplotlib + seaborn` | 繪圖庫 | 生成 calibration.png（可靠性圖）|
| `scipy` | 科學計算 | temperature scaling 優化（如果 ECE 太高）|

### 要看的檔案
```
✅ 06_uncertainty_model/AGENT_REPORT.md                                   最重要
✅ 06_uncertainty_model/processes/ensemble.py                             模型架構
✅ 06_uncertainty_model/processes/01_train_deep_ensemble_synthetic.ipynb  主 notebook
✅ 06_uncertainty_model/outputs/cycle_0/model_meta.json                   訓練後的參數記錄
✅ 06_uncertainty_model/outputs/cycle_0/calibration.png                   校準圖
⏭️ 06_uncertainty_model/outputs/cycle_0/ensemble_member_*.pt            5 個模型權重，不用打開
⏭️ 06_uncertainty_model/processes/tests/                                 測試，可跳過
```

### Deep Ensemble 是什麼？
訓練 5 個一模一樣的神經網路，但用不同的隨機起始點（random seed 0, 1, 2, 3, 4）。
然後讓它們各自預測同一個輸入——如果 5 個都說「結合力大概是 7.5」，那就很確定。
如果 5 個的預測差很多（2.3, 5.1, 8.0, 4.7, 6.2），那就非常不確定。
不確定性 = 這 5 個預測的標準差（std）。

### ELISA 到了之後改哪裡？
按照 `AGENT_REPORT.md` 的 "When ELISA Arrives" section，只要改兩個地方：
1. `01_train_deep_ensemble_synthetic.ipynb` 的 **Cell 3**：把假資料換成 ELISA CSV 路徑
2. `model_meta.json` 的 `data_source` 欄位

---

## 哪些檔案可以完全忽略（不用讀）

### 已 Archive 的（和這個版本完全無關）
```
archive/2026-05-pivot/                    舊版 6-factor 迴歸 pipeline
archive/2026-05-pivot/02_annotation_batch_outputs/   批次注釋的副產品
  ├── host_proteins_non_xcc/             Bacillus, E.coli 等無關細菌的蛋白質組
  ├── pharokka_runs/                     pharokka 不在目前 pipeline 裡
  └── phage_phanotate_coords_misc/       非 phiL7 噬菌體的舊坐標
docs/archive/                            舊版 workflow chart 和整合指南
```

### 還在 repo 但暫時不用管的
```
07_acquisition_function/        下一個 sprint 才要建
08_cycle_data/                  6 月 1 日 ELISA 結果才會有東西
scripts/setup_worktrees.sh      隔夜平行構建用的腳本，已完成使命
00_raw_data/data_needs.md       Sarah/Weitao 的待辦事項
```

---

## 看每個 Module 的建議順序

**如果你只有 20 分鐘：**
1. `CLAUDE.md` → 了解規則
2. `docs/pipeline_build_report_2026-05-08.md` → 了解現況
3. `03_rbp_identification/outputs/EU717894.1_rbp_candidates.csv` → 看看找到了什麼

**如果你要跑代碼：**
每個 module 的 `processes/` 裡都有編號的 notebook（01_, 02_...），按順序跑。
每個 module 都有 `requirements.txt` 告訴你要裝什麼。

**如果你要解釋給別人聽：**
`docs/planning/iGEM_2026_Project_Plan.md` → 給 PI 的英文版
`docs/planning/iGEM_2026_项目大纲_中文版.md` → 給團隊的中文版
`docs/pipeline_build_report_2026-05-08.md` → 給自己（也可以上傳 NotebookLM）

---

## 最重要的 3 件事（如果只能記一件事）

1. **資料流是單向的**：`inputs/` → `processes/` → `outputs/`。不要在 `inputs/` 裡寫東西。

2. **Module 04 需要重跑**：embedding 現在用的是假序列。
   把 `01_embed_esm2.ipynb` Cell 8 的路徑改成 Module 03 的真實輸出，跑一次，就好了。

3. **ELISA 到了只改兩行**：Module 06 的架構完全不用動。

---

*更新日期：2026-05-09*
*對應 branch：`active-learning-pipeline`*
