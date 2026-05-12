# Paper Reading Notes — iGEM Claremont 2026
# 論文閱讀筆記

**讀完時間：** 2026-05-11  
**讀了哪些：** 19 篇核心論文（docs/reference/papers/ 全部）  
**目的：** 交叉驗證所有 project documents 中的 assumptions 和引用數字是否正確

---

## ⚠️ 錯誤清單（按嚴重度排序）

### 🔴 CRITICAL — 事實錯誤，必須修正

---

#### Error 1 — ExbD2 不是必需的（Hung 2003 直接推翻）

**我們說：** 「phiL7 通過 **TonB-ExbB-ExbD1-ExbD2** 受體系統感染 Xcc」（所有文件普遍）

**論文實際說：**
> "TonB, ExbB, and ExbD1 are essential for penetration of phage ϕL7"  
> CH620（exbD2::pUD2190）mutant **retained** sensitivity to ϕL7

ExbD2 knockout **不影響** phiL7 感染。論文結論明確只有三個基因 essential：tonB、exbB、exbD1。

**影響：**
- 所有文件裡的「TonB-ExbB-ExbD1-ExbD2」要改成「TonB-ExbB-ExbD1（ExbD2 不是必需的）」
- Wet lab knockout 設計中，exbD2 的優先度要下調
- Module 05 Boltz-2 input 裡 exbD2 的角色需要重新說明

---

#### Error 2 — Boltz-2 的 affinity head 是針對小分子-蛋白，不是蛋白-蛋白（Passaro 2025）

**我們說：** 「Boltz-2 的 affinity head 能直接輸出 binding affinity 估計」用於 RBP × receptor

**論文實際說：**
> "its key distinctive feature is its ability to predict binding affinity, **which measures how tightly small molecules attach to proteins**"

Training data：PubChem、ChEMBL、BindingDB——全是**小分子資料庫**。benchmark 也是 protein-**ligand** benchmark（CDK2, TYK2, JNK1, P38）。

**事實：** 對 protein-protein pair（RBP × TonB），affinity head 輸出 NaN——build report 其實已經正確記錄了這一點，但 planning documents 仍然聲稱我們能得到 zero-shot affinity prior，這是誤導性的。

**我們實際能得到的：** ipTM（interface confidence score），不是 affinity。

**影響：**
- Planning docs 中「Boltz-2 zero-shot 亲和力当作 synthetic prior」需要修正
- 應改為：「Boltz-2 的 ipTM 作為結構信心分數 proxy，不是定量 affinity」

---

#### Error 3 — Greenman 2025 期刊錯誤，且結論被誤讀（Greenman 2025）

**我們說：** 「Greenman et al. 2025, *NAR Genom Bioinform*」，且說 deep ensemble 是最佳 UQ 方法

**論文實際說：**
- 期刊：**PLoS Computational Biology** 21(1): e1012639（不是 NAR）
- 結論：「**there is no single best UQ method** excels across all scenarios」
- 具體：ensemble 在 accuracy 上常是最高的，但 calibration 是最差的之一
- 「uncertainty-based strategies for property optimization **often did not outperform simpler methods**」

**影響：**
- 期刊名稱需要在所有文件中修正
- 「deep ensemble 被認定為最佳方法之一」這個 claim 需要改成更精確的說法

---

#### Error 4 — Hie 2024 用的是 ESM-1b/1v，不是 ESM-2（Hie 2024）

**我們說：** 「Hie et al. 2024 *Nat Biotechnol* 用 ESM 在 ~50 個 antibody 數據點上」

**論文實際說：**
- 用了 **ESM-1b 和 ESM-1v**（不是 ESM-2）
- 訓練在 UniRef50 和 UniRef90
- 每個抗體篩選 **20 個或更少 variant**（不是 50）
- 兩輪實驗，不是主動學習閉環

**額外注意：** Hie 2024 的方法是「語言模型 likelihood 過濾突變」，不是「以測量結果訓練模型後推薦下一輪實驗」的 BALD-style active learning。機制上有根本差異。

---

#### Error 5 — Lee 2009 從未明確指出任何蛋白是 tail spike（Lee 2009）

**Build report 說：** 「Lee et al. 2009 identified gp25 as phiL7's tail spike」

**論文實際說：**
- p20（tail protein III，**1105 aa**）——「sequence analysis suggested that **p20 instead has the potential to play this role**」（host range determination）
- p25（85 aa）——Table 1 標記「No similarity」
- 論文明確指出 phiL7 **lacks** OP1 的 ORF25 tail fiber homologue

**我們的 rbp_01（712 aa）：** 不對應 p20（1105 aa）或 p25（85 aa）。這是 PHANOTATE 找到的 orf_00001，命中 Tail_spike_N HMM domain（PhageRBPdetect 的 HMM 結果），跟 Lee 2009 的注釋是平行的獨立證據。

**結論：** rbp_01 是 tail spike 候選的依據是 PhageRBPdetect HMM（Tail_spike_N domain hit），不是 Lee 2009 的注釋。Build report 的說法需要修正。

---

### 🟡 MODERATE — 不夠精確的 claims

---

#### Error 6 — Boeckaerts 2024 AUC 的條件依賴性

**我們說：** 「AUC 0.82」（固定值）

**論文實際說：** cross-validated ROC AUC **up to 81.8%**，在 100% identity threshold 下。隨著 threshold 降低（更嚴格的 cross-species 測試），AUC 下降到 0.698（75% identity threshold）。

**修正：** 應說「在菌株內驗證下 ROC AUC up to 81.8%，跨菌株應用時下降」

---

#### Error 7 — BALD 原始論文是針對 Gaussian Process Classifier，不是 deep ensemble

**Houlsby 2011 實際：** BALD 在論文中應用於 Gaussian Process Classifier（GPC）。  
**我們的用法：** 把 BALD 的 information-theoretic objective 套用到 deep ensemble 的 variance output 上——這是合理的延伸，但不在原始論文中。  
**另外：** Houlsby 2011 是 classification；我們的問題是 regression（binding score）。

**修正：** 說明我們用的是 BALD 的 objective（equation 2），延伸應用到 deep ensemble regression。

---

#### Error 8 — ALDE 用的是 Thompson sampling，不是 BALD

**Yang 2025 實際：** 用了 DNN ensemble + **Thompson sampling** acquisition function，以及 one-hot encoding（不是 ESM-2 embedding）。  
**我們的說法：** 把 ALDE 作為「BALD 被驗證有效」的依據——但 ALDE 用的是 Thompson sampling。

**修正：** ALDE 驗證了「active learning + uncertainty quantification 在蛋白質工程有效」，不是特別驗證 BALD。BALD 的依據應另外引用。

**另外：** Yang 2025 的 wet lab 結果是「yield 從 12% 提升到 93%」（3 輪），不是「比 random 快 2-5 倍」——那個數字來自 Hie et al. 2022 Cell 的說法。

---

#### Error 9 — PLM-interact 的 AUPR 提升數字需要更謹慎

**我們說：** 「AUPR 提升 16-28%」（對應遷移到 mouse/fly/worm/yeast/*E. coli*）

**論文實際：**
- 測試物種：mouse、fly、worm、yeast、*E. coli*（不包含噬菌體）
- 相對 TUnA（最接近的 baseline）的 AUPR 改善各有差異
- Mouse：TUnA 0.779 → PLM-interact 0.904（+16%）
- Fly：TUnA 0.757 → PLM-interact 0.913（+21%）
- 「virus-host interaction」的說法在 abstract 提到，但指的是人類病毒-宿主，不是噬菌體-細菌

**修正：** 16-28% 數字在特定比較下成立，但需要說清楚是「相對 TUnA baseline，在 5 個物種的跨物種 benchmark 上的 AUPR 提升」，範圍不是固定的。

---

#### Error 10 — Latka 2021 的適用性有限制

**Latka 2021 實際：** Klebsiella 噬菌體的 RBP 是**多醣解聚酶（depolymerase）**，降解莢膜多醣（CPS）。N 端結構錨定 + C 端特異性決定區是這個系統的特性。

**phiL7 × TonB 系統：** TonB 是**內膜蛋白**（iron uptake），phiL7 RBP 不降解多醣，是完全不同的受體結合機制。

**修正：** N/C 端模塊化的概念可以作為一般原則引用，但 Latka 2021 的具體生化機制（depolymerase）不適用於 phiL7 系統，應注明。

---

### ⚪ MINOR — 措辭問題（不影響核心論述）

---

#### Error 11 — Boeckaerts 2022 ML track 的 ESM 版本

Boeckaerts 2022 發表於 2022 年 6 月，ESM-2 論文（Lin et al.）是 2022 年 10 月 bioRxiv preprint。因此 Boeckaerts 2022 用的 ESM 不是 ESM-2，而是舊版 ESM（ESM-1b 或更早）。這個細節不影響我們自己的 pipeline（我們用 ESM-2），但不應說 Boeckaerts 2022「用了 ESM-2」。

---

#### Error 12 — BALD 的嚴格數學說法

我們說「BALD 等價於 mutual information 最大化——預測 entropy 減去預測的 expected entropy」。  
Houlsby 2011 的準確說法：最大化 **I[θ; y | x, D]**，即模型參數 θ 和 label y 之間的條件互信息。  
Equation (2)：argmax H[y|x,D] − E_{θ~p(θ|D)}[H[y|x,θ]]  
我們的描述數學上是對的，只是省略了「參數不確定性」這一層的精確性。

---

## ✅ 確認正確的 claims

| Claim | 來源 | 結論 |
|-------|------|------|
| PhageRBPdetect precision-recall AUC 93.8%，F1 84.0% | Boeckaerts 2022 | ✅ |
| PHANOTATE 用 dynamic programming 處理重疊 ORF | McNair 2019 | ✅ |
| Prodigal 假設 ORF 不重疊（GeneMark-style）| Hyatt 2010 | ✅ |
| ESM-2 訓練時見到 ~65M unique sequences | Lin 2023 | ✅ |
| ESM-2 用 masked language modeling | Lin 2023 | ✅ |
| ESM-2 650M = 1280 維 embedding | Lin 2023 | ✅ |
| Lakshminarayanan 2017：ensemble 輸出 μ(x) + σ²(x)，用 Gaussian NLL loss | 2017 paper | ✅ |
| AF3：MSA + diffusion module，handles protein-protein | Abramson 2024 | ✅ |
| Hung 2003：PMID 12646254，BBRC 302(4):878-884 | Verified | ✅ |
| Hung 2003：tonB、exbB、exbD1 三個 essential | Hung 2003 | ✅ |
| BALD equation 2：argmax H[y|x,D] − E_θ[H[y|x,θ]] | Houlsby 2011 | ✅ |
| PhageHostLearn AUC up to 81.8%（at 100% identity threshold）| Boeckaerts 2024 | ✅ |
| Latka 2021：C 端決定 serotype specificity | Latka 2021 | ✅ |
| PLM-interact trained on human PPI，tested on mouse/fly/worm/yeast/E.coli | Liu 2025 | ✅ |
| Boltz-2 對蛋白-蛋白 pair：affinity head 輸出 NaN | Passaro 2025 + build report | ✅ |

---

## 每篇論文關鍵事實摘要

### Lee 2009 (*Appl Environ Microbiol* 75:7828)
- phiL7：44,080 bp，56% G+C，lytic Siphoviridae
- **59 個 ORF**（PHANOTATE 1.6.7 找到 80 個，因為更新版模型識別更多重疊 ORF）
- 10 個 virion 蛋白（SDS-PAGE + LC-MS/MS）
- **p20（tail protein III，1105 aa）**——建議可能與 host range determination 有關
- phiL7 **lacks** OP1 的 ORF25 tail fiber homologue
- p25（85 aa）——"No similarity"
- **從未說任何蛋白是 tail spike**

### Hung 2003 (*BBRC* 302:878, PMID 12646254)
- ϕL7 感染 *Xcc* P20H（不是 ATCC 33913）
- Tn5 mutagenesis 找到 tonB 突變導致抗性
- **TonB、ExbB、ExbD1 essential；ExbD2 NOT essential**
- CH620（exbD2 mutant）保留對 ϕL7 的敏感性
- 四個基因 co-transcribed（operon）

### Boeckaerts 2022 (*Viruses* 14:1329)
- PhageRBPdetect：HMM track + XGBoost ML track
- precision-recall AUC **93.8%**，F1 **84.0%**（vs PhANNs 69.8%）
- RBP 長度篩選：200 < length < 1500 aa
- N 端結構域 + C 端結合/特異性域的雙區框架
- rbp_01（712 aa）在 200-1500 範圍內 ✅

### Boeckaerts 2024 (*Nat Commun* 15:4355)
- PhageHostLearn：PHANOTATE + PhageRBPdetect + ESM-2（t33_650M_UR50D）+ XGBoost
- cross-validated ROC AUC **up to 81.8%**（100% identity），下降到 0.698（75% identity）
- Klebsiella 特定系統，CPS（capsular polysaccharide）as receptor
- In vitro validation on 28 clinical *K. pneumoniae* isolates

### Yang 2025 (*Nat Commun* 16:714)
- ALDE：active learning + DNN ensemble + **Thompson sampling**（不是 BALD）
- one-hot encoding（不是 ESM-2 embedding）
- 應用：ParPgb 酶的 cyclopropanation，5 個活性位點
- 結果：yield 從 12% → 93%（2 輪 wet lab，總計 216+90+90 個序列）
- 比較的是 AL vs 直接 DE，不是 AL vs random
- Frances Arnold 組（Caltech）

### Lakshminarayanan 2017 (*NeurIPS* 2017)
- Deep Ensembles for predictive uncertainty
- 輸出：μ(x)（預測值）和 σ²(x)（不確定性），用 Gaussian NLL loss
- 不用 bagging；random initialization 就足夠多樣性
- 比 MC Dropout 的 calibration 更好（或相當）
- Google DeepMind 作者

### Houlsby 2011 (arXiv:1112.5745)
- BALD：Bayesian Active Learning by Disagreement
- 應用於 Gaussian Process Classifier（GPC）
- Objective：argmax H[y|x,D] − E_{θ~p(θ|D)}[H[y|x,θ]] (Eqn 2)
- 等價於 I[θ; y|x, D]（model parameters 和 label 的互信息）
- 原始論文是分類問題；延伸到 regression + ensemble 是後人的工作

### Hie 2024 (*Nat Biotechnol* 42:275)
- 語言模型導向的抗體親和力成熟
- 用 **ESM-1b 和 ESM-1v**（不是 ESM-2）
- UniRef50（ESM-1b）和 UniRef90（ESM-1v，~98M total sequences）
- 每個抗體篩選 **20 個或更少 variant**，兩輪
- **機制不同**：語言模型 likelihood 過濾突變，不是 BALD closed-loop
- 7 個抗體，親和力改善最高 160 倍

### Greenman 2025 (*PLoS Comput Biol* 21(1): e1012639)
- Benchmarked on FLIP (Fitness Landscape Inference for Proteins) benchmark
- 三個 datasets：GB1（binding/epistasis）、AAV（viral viability）、Meltome（thermostability）
- 結論：**no single best UQ method**
- Ensemble 常常 accuracy 最高但 calibration 最差
- Uncertainty-based sampling 在 BO 中 often 不比 greedy 好

### Lin 2023 / ESM-2 (*Science* 379:1123, bioRxiv preprint 2022)
- 8M to 15B parameters
- Training：UniRef50/90，**~65M unique sequences** seen during training
- masked language modeling objective
- ESM-2 650M：1280 維；ESM-2 3B：2560 維 ✅
- 比 ESM-1b 650M 更好（相同參數規模下）

### Liu 2025 (*Nat Commun* 16:9012)
- PLM-interact：fine-tune ESM-2（650M）on human PPI data
- Training：474,517 human PPI pairs
- Testing：mouse、fly、worm、yeast、*E. coli*
- AUPR improvements vs TUnA：mouse 16%、fly 21%、等
- **未測試 phage-bacteria** ✅（我們的 claim 正確）
- 「virus-host interaction」在 abstract 提到，指人類病毒（不是噬菌體）

### Passaro 2025 / Boltz-2 (bioRxiv, June 2025)
- 針對**小分子-蛋白 binding affinity**（不是蛋白-蛋白）
- Training data：PubChem、ChEMBL、BindingDB（小分子資料庫）
- Benchmark：protein-**ligand** benchmark（CDK2, TYK2, JNK1, P38）
- Protein-protein pair：affinity head 輸出 NaN（confirmed by build report）
- ipTM = structural confidence，不是 affinity
- Weights 開源（Github: jwohlwend/boltz）

### Latka 2021 (*mBio* 12:e00455-21)
- Klebsiella 噬菌體 RBP 的 N 端結構域（anchor）+ C 端酶活性域（specificity）
- C 端的 depolymerase domain 決定 capsular serotype specificity
- 構建了 16 個 chimeras，通過 VersaTile 方法
- 注意：這個系統是 **CPS depolymerase**，與 phiL7/TonB 的機制完全不同

### Abramson 2024 / AF3 (*Nature* 630:493)
- DOI: 10.1038/s41586-024-07487-w ✅
- Pairformer 替代 evoformer（AF2）
- MSA 的 role 大幅減少
- Diffusion module 直接預測原子坐標
- 能預測 protein-protein、protein-ligand、protein-nucleic acid 等
- 優於 AF-Multimer v2.3 在 protein-antibody 上

### Houlsby 2011 (BALD) — 補充
- arXiv:1112.5745v1，Cambridge 大學
- 作者：Houlsby, Huszár, Ghahramani, Lengyel
- 應用到 GPC，不是 ensemble
- 原始論文不包含 regression 應用

---

## 修正優先序

| 優先度 | 錯誤 | 影響範圍 |
|--------|------|---------|
| P1 🔴 | ExbD2 not essential | 6 個文件，wet lab 設計 |
| P1 🔴 | Boltz-2 affinity = NaN for protein-protein | planning docs 的 synthetic prior claim |
| P1 🔴 | Greenman 2025 期刊錯誤 | 所有 reference lists |
| P2 🟡 | Hie 2024 用 ESM-1b/1v 不是 ESM-2 | papers.md、planning docs |
| P2 🟡 | Lee 2009 不說 gp25 是 tail spike | build report |
| P2 🟡 | ALDE 用 Thompson sampling 不是 BALD | planning docs |
| P3 ⚪ | AUC "0.82" → "up to 81.8%" | planning docs |
| P3 ⚪ | BALD 原始是 GPC，我們是 extension | descriptions |
| P3 ⚪ | Latka 2021 適用性限制 | papers.md |
