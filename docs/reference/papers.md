# Papers by Module — Reading & Annotation Guide
# 各模組核心文獻 — 閱讀與標注指南

**Purpose / 用途：** 每個 module 要讀什麼、為什麼讀、這篇文章在我們系統中扮演什麼角色。
標記 🔴 = 必讀（直接影響設計決策）/ 🟡 = 推薦（理解背景）/ ⚪ = 備查（有疑問時再看）。

**Annotation protocol:** 讀完每篇後，在個人筆記寫下：
1. 這篇的核心 claim 是什麼？
2. 我們借用了它的哪個具體設計決策？
3. 有沒有跟我們 setup 不同的地方需要注意？

---

## 📁 Module 00 — Raw Data & Reference Genomes

### 🔴 Lee, C.N. et al. (2009)
**"Genomic characterization of the intron-containing T7-like phage phiL7 of *Xanthomonas campestris* pv. *campestris*"**
*Applied and Environmental Microbiology* 75(24):7828–7838.
NCBI accession: EU717894.

**Why read / 為什麼讀：** phiL7 是我們整個 dry lab pipeline 的參考噬菌體。這篇是 phiL7 基因組的原始定性文章，包含 59 個 ORF 的注釋、group I intron 結構、tail spike 候選基因位置。

**Role in our system / 在我們系統中的角色：**
- Module 00：確認 EU717894 是正確的 phiL7 accession
- Module 03：告訴我們 phiL7 的 tail spike（P25 homolog）在哪裡——但注意這篇報的是 85 aa P25，真正的 712 aa tail spike（rbp_01）是後來 PhageRBPdetect 抓出來的，**兩者不同**
- Module 05：給 Boltz-2 跑的 FASTA 序列來源背景

**Key figure to annotate:** Figure 1（基因組圖）、Table 1（ORF 列表）——找到 tail spike / RBP 候選的位置。

---

### 🔴 da Silva, A.C.R. et al. (2002)
**"Comparison of the genomes of two *Xanthomonas* pathogens with differing host specificities"**
*Nature* 417:459–463.

**Why read / 為什麼讀：** Xcc ATCC 33913（GCF_000007145.1）的基因組原始文章。TonB-ExbB-ExbD 操縱子在這裡被定序，是 Wang et al. 2003 受體研究的前提。

**Role in our system / 在我們系統中的角色：**
- Module 00：確認 GCF_000007145.1 / AE008922 是 Xcc ATCC 33913（注意：GCF_000007145.1 是新版 assembly accession，AE008922 是原始 GenBank 登錄號，指向同一基因組）
- Module 02：Prodigal 注釋的 reference genome 就是這個
- Module 05：TonB 等受體的序列來源

---

### 🔴 Wang, W.-T. et al. (2003)
**"Involvement of the *tonB-exbBD1D2* operon of *Xanthomonas campestris* pv. *campestris* in infection by the ϕL7 phage"**
*Molecular Microbiology* 50(2):507–519.

**Why read / 為什麼讀：** 實驗確認 phiL7 的受體系統是 TonB-ExbB-ExbD1-ExbD2 四個蛋白。這是整個 Layer 2 因果驗證模組的生物學依據。

**Role in our system / 在我們系統中的角色：**
- Module 05：決定 Boltz-2 要跑的 receptor 是哪幾個（TonB、ExbB、ExbD1、ExbD2）
- Module 08（wet lab）：ΔtonB / ΔexbB / ΔexbD 敲除實驗的目標基因直接來自這篇
- 這篇是「phiL7 走 TonB 通道」這個 claim 的唯一一手來源，**一定要讀**

**Key section to annotate:** Methods（哪個敲除導致抗性）、Figure 3-4（plaque assay 結果）。

---

## 📁 Module 01 — Ground Truth / Interaction Matrix

### 🟡 Boeckaerts, D. et al. (2024)
**"Predicting bacteriophage hosts based on sequences of annotated receptor-binding proteins"**
*Nature Communications* 15:4768.

**Why read / 為什麼讀：** PhageHostLearn——目前 phage-host interaction 預測的 SOTA 模型。我們的 interaction matrix（Module 01 的 positive/negative pair 格式）直接模仿這篇的數據結構。

**Role in our system / 在我們系統中的角色：**
- Module 01：interaction matrix 的格式設計參考；positive pair 的定義（同一實驗驗證的 phage-host）
- Module 04：PhageHostLearn 用純 PLM embedding 達到 AUC 0.82，是我們放棄 6-factor 手工特征的核心依據
- Module 06：我們的 ensemble 要 beat 的 baseline

**Key table to annotate:** Table 1（數據集構成）、Figure 2（AUC 比較）。

---

### ⚪ Mutalik group, PAML benchmark (2025)
**"Phage Anti-Microbial Landscape benchmark"**
*bioRxiv* (preprint, 2025).

**Why read / 為什麼讀：** 這個 benchmark 顯示 PhageHostLearn 的跨菌種 AUC 掉到 0.67–0.70，是「現有模型不夠好」這個 claim 的獨立驗證。

**Role in our system / 在我們系統中的角色：** 主要是背景依據，支撐我們為什麼做 active learning 而不是直接用現成模型。在 wiki 的 Motivation section 會引用。

---

## 📁 Module 02 — Genome Annotation

### 🔴 McNair, K. et al. (2019)
**"PHANOTATE: a novel approach to gene identification in phage genomes"**
*Bioinformatics* 35(22):4537–4542.

**Why read / 為什麼讀：** PHANOTATE 是我們跑 phage ORF calling 的工具。理解它為什麼比 Prodigal 更適合 phage——dynamic programming 處理重疊 ORF，把基因組建模成加權有向圖。

**Role in our system / 在我們系統中的角色：**
- Module 02：我們只對 phage 基因組用 PHANOTATE；細菌用 Prodigal。這個分工來自這篇對 phage-specific overlap 問題的描述
- **不能對調**：用 Prodigal 跑 phiL7 會丟失約 10–15% 的真實基因（因為 Prodigal 假設 ORF 不重疊）

**Key section to annotate:** Methods（圖論建模的說明）、Supplementary Table（performance comparison vs Prodigal/GeneMark）。

---

### 🔴 Hyatt, D. et al. (2010)
**"Prodigal: prokaryotic gene recognition and translation initiation site identification"**
*BMC Bioinformatics* 11:119.

**Why read / 為什麼讀：** Prodigal 是我們跑細菌 ORF calling 的工具（透過 pyrodigal Python binding）。

**Role in our system / 在我們系統中的角色：**
- Module 02：Xcc ATCC 33913 的 host protein .faa 由 pyrodigal 生成
- 讀這篇主要是確認 Prodigal 在細菌基因組上的假設（ORF 不重疊、統計先驗）與我們用法一致

**Key section to annotate:** Abstract（GeneMark-style model 說明）、Results（accuracy 數據）。

---

### 🟡 Bouras, G. et al. (2023)
**"Pharokka: a fast scalable bacteriophage annotation tool"**
*Bioinformatics* 39(1):btac776.

**Why read / 為什麼讀：** pharokka 整合 PHROG + CARD + VFDB 做 phage 功能注釋，是 PHANOTATE 之後的第二步。

**Role in our system / 在我們系統中的角色：**
- Module 02：pharokka 跑完後的 PHROG 分類告訴我們哪些 ORF 被歸類為「tail spike / tail fiber」，這是 Module 03 HMM 掃描的輔助驗證
- 若 PhageRBPdetect 的 HMM track 和 pharokka 的 annotation 都指向同一個 ORF → 高信心 RBP 候選

---

## 📁 Module 03 — RBP Identification

### 🔴 Boeckaerts, D. et al. (2022)
**"Identification of phage receptor-binding protein sequences with hidden Markov models and an extreme gradient boosting classifier"**
*Viruses* 14(6):1329.

**Why read / 為什麼讀：** PhageRBPdetect 就是這篇。雙軌設計（HMM 掃描 + ESM embedding + XGBoost 分類）直接決定了 Module 03 的架構。

**Role in our system / 在我們系統中的角色：**
- Module 03 就是在跑這篇的工具。讀這篇讓你理解：
  1. HMM track 從 Pfam 哪些 family 掃（tail fiber、tail spike、carbohydrate-binding 等）
  2. XGBoost track 為什麼能抓 Pfam 漏掉的發散 RBP（phiL7 就是這種情況）
  3. 我們的 rbp_01（712 aa tail spike）是 HMM track 還是 ML track 識別的？→ 讀完才知道去哪裡找 evidence

**Key figure to annotate:** Figure 2（雙軌流程圖）、Table 2（precision-recall AUC 93.8%）。

---

### 🟡 Latka, A. et al. (2021)
**"Engineering the modular receptor-binding proteins of *Klebsiella* phages switches their *in vitro* host range"**
*mBio* 12(6):e02329-21.

**Why read / 為什麼讀：** RBP truncation 策略（N 端 anchor + C 端 specificity head）的實驗依據。這是 Cycle 0 variant design 的主要靈感來源之一。

**Role in our system / 在我們系統中的角色：**
- Module 08（variant design）：truncation series 的設計原則直接借這篇。N 端錨定 capsid，C 端決定受體特異性——切對地方才能保留結合能力
- 注意：這篇是 *Klebsiella* phage，跟我們的 Xanthomonas phage 有物種差異；annotate 時記錄哪些部分可能不 transfer

---

### 🟡 Yehl, K. et al. (2019)
**"Engineering phage host-range and suppressing bacterial resistance through phage tail fiber mutagenesis"**
*Cell* 179(2):459–469.

**Why read / 為什麼讀：** T7 尾部纖維蛋白的 DMS（deep mutational scan）——在 binding interface 做系統突變，找到關鍵殘基。這是 Module 03 輸出的 per-residue ESM embedding 要做的分析的藍本。

**Role in our system / 在我們系統中的角色：**
- Module 03 / 05：找到 surface-exposed binding loops 的方法論參考
- Module 08：位點突變 variant design 的依據（which residues to mutate）

---

## 📁 Module 04 — Protein Embedding (ESM-2)

### 🔴 Lin, Z. et al. (2023)
**"Evolutionary-scale prediction of atomic-level protein structure with a language model"**
*Science* 379(6637):1123–1130.

**Why read / 為什麼讀：** ESM-2 原始論文。masked language modeling 在 65M 蛋白質上預訓練 → 每個氨基酸位置得到 1280 維（650M 模型）或 2560 維（3B 模型）embedding。

**Role in our system / 在我們系統中的角色：**
- Module 04 的核心工具。讀這篇理解：
  1. 為什麼 mean pooling over residues 是合理的 sequence-level embedding（而不是只取 CLS token）
  2. 650M vs 3B model 的差異在哪裡（性能 vs 計算量）
  3. Per-residue embedding 能用來做 motif-level 分析（Module 03 輸出的 RBP 結合殘基定位）

**Key figure to annotate:** Figure 4（embedding 在結構預測上的 zero-shot 性能）。

---

### 🟡 Hie, B.L. et al. (2024)
**"Efficient evolution of human antibodies from general protein language models"**
*Nature Biotechnology* 42:275–283.

**Why read / 為什麼讀：** 直接示範「用 ESM-2 embedding + ~50 個測量點就能做 active learning 推進抗體進化」。這是我們整個 project framing 最直接的 precedent。

**Role in our system / 在我們系統中的角色：**
- Module 04 + 06 + 07 的組合參考。他們的 pipeline（ESM embed → fitness model → acquisition）就是我們的藍本
- Effective sample size 估計（50 個數據點等效中型訓練集）的 cite 來源

**Key figure to annotate:** Figure 2（learning curve：AL vs random）、Figure 3（50個點的 generalization）。

---

### 🟡 Liu, Y. et al. (2025)
**"PLM-interact: learning the language of protein-protein interactions"**
*Nature Communications* 16:XXXX.

**Why read / 為什麼讀：** 把 ESM-2 在人類 PPI 數據上 fine-tune，遷移到 mouse / fly / worm / yeast / *E. coli* PPI 預測，AUPR 提升 16–28%。我們是第一個試著把這個遷移用到 phage-bacteria 這個 niche 的。

**Role in our system / 在我們系統中的角色：**
- 第 2 層數據稀缺解決策略（PLM-interact transfer prior）
- Module 04 的 optional embedding variant：用 PLM-interact 的 fine-tuned weights 代替純 ESM-2，看 RBP-receptor pair score 是否比純 ESM-2 更 informative
- 注意：他們沒有測過 phage-bacteria；annotate 時記錄他們的跨物種遷移在什麼情況下失效

---

## 📁 Module 05 — Structure Prediction

### 🔴 Abramson, J. et al. (2024)
**"Accurate structure prediction of biomolecular interactions with AlphaFold 3"**
*Nature* 630:493–500.

**Why read / 為什麼讀：** AlphaFold 3 原始論文。MSA + diffusion model 預測蛋白複合體 3D 結構，輸出 ipTM（interface predicted TM score）信心分數。

**Role in our system / 在我們系統中的角色：**
- Module 05：phiL7 rbp_01 trimer + TonB/ExbB receptor 的靜態結構預測工具
- ipTM > 0.5 是我們判斷「預測結構可信」的閾值（記錄在 Module 05 README）
- 注意 AF3 model weights 需要 Google form 申請（尚未完成——這是 pending action item）

**Key metric to understand:** ipTM 的意義（0–1 分，>0.8 high confidence，0.4–0.8 medium，<0.4 low）。我們的 Boltz-2 第一次跑出 0.345（用錯蛋白），理解 ipTM 才能解讀結果。

---

### 🔴 Passaro, J.M. et al. (2025)
**"Boltz-2: towards accurate and efficient binding affinity prediction"**
*bioRxiv* (preprint).

**Why read / 為什麼讀：** Boltz-2 是我們唯一能算 zero-shot affinity prior 的工具（AF3 沒有 affinity head）。在 PDBbind 上訓練的 affinity head 能直接從序列輸出 binding affinity 估計。

**Role in our system / 在我們系統中的角色：**
- Module 05：phiL7 rbp_01 × TonB 的 Boltz-2 affinity → 第 4 層數據稀缺策略的 synthetic prior
- 第一次 run（job 59949）用的是 85 aa P25（錯誤），ipTM = 0.345，需要用 712 aa rbp_01 重跑
- 讀這篇理解 affinity head 的訓練數據分佈（PDBbind 主要是小分子，蛋白複合體有限）——這是 prior bias 的風險來源

**Key figure to annotate:** Figure 2（affinity prediction performance on PDBbind test set）、Table 1（benchmark vs AF3 / RoseTTAFold）。

---

## 📁 Module 06 — Uncertainty Model (Deep Ensemble)

### 🔴 Lakshminarayanan, B. et al. (2017)
**"Simple and Scalable Predictive Uncertainty Estimation Using Deep Ensembles"**
*NeurIPS* 30.

**Why read / 為什麼讀：** Deep Ensemble 原始論文。5 個 MLP 獨立訓練，推理時取均值（預測值）和方差（epistemic uncertainty）。這是 Module 06 整個架構的設計文件。

**Role in our system / 在我們系統中的角色：**
- Module 06 就是實現這篇的設計：5 個 ensemble member，每個用不同 random seed，Gaussian NLL loss
- **Epistemic uncertainty（方差）= BALD 的輸入**。沒有這個方差，Module 07 的 acquisition function 就無法運算
- 理解 temperature scaling：ensemble calibration 飄掉時（predicted 80% confidence 但只有 50% 準確），用 temperature scaling 修正

**Key section to annotate:** Section 3.3（Gaussian NLL loss 推導）、Figure 3（calibration 比較）。

---

### 🔴 Greenman, K.P. et al. (2025)
**"Benchmarking uncertainty quantification methods for protein fitness prediction"**
*NAR Genomics and Bioinformatics* (advance access).

**Why read / 為什麼讀：** 在蛋白質工程這個 exact domain 上系統 benchmark 了 UQ 方法（deep ensemble vs MC dropout vs GP vs conformal prediction）。Deep ensemble 在 protein fitness landscape 上被認定為最佳選擇之一。

**Role in our system / 在我們系統中的角色：**
- 選 deep ensemble 而不是 GP 或 MC Dropout 的直接依據
- 讀 annotate 時記錄：他們 benchmark 的 protein fitness dataset 跟我們的 ELISA binding 有什麼差異？calibration 指標哪個最 relevant？

---

### 🟡 Ovadia, Y. et al. (2019)
**"Can You Trust Your Model's Uncertainty? Evaluating Predictive Uncertainty Under Dataset Shift"**
*NeurIPS* 32.

**Why read / 為什麼讀：** 在 distribution shift 下比較了 deep ensemble vs MC Dropout vs temperature scaling。我們的 wet lab 數據跟 Boltz-2 prior 之間就存在這種 shift。

**Role in our system / 在我們系統中的角色：** 理解 calibration 飄掉的原因（distribution shift），有助於診斷 cycle 1/2 的模型行為。

---

## 📁 Module 07 — Acquisition Function (BALD)

### 🔴 Houlsby, N. et al. (2011)
**"Bayesian Active Learning for Classification and Preference Learning"**
*arXiv*:1112.5745.

**Why read / 為什麼讀：** BALD 原始論文。Mutual information 最大化 = 預測 entropy 減去 expected posterior entropy。這是 Module 07 acquisition function 的數學核心。

**Role in our system / 在我們系統中的角色：**
- Module 07 就是實現 BALD：對每個 unmeasured (RBP variant, receptor) 組合算 BALD score，取 top 4–5
- 讀完你才能理解：為什麼 BALD 比 greedy（直接選預測值最高）更適合冷啟動？（greedy 在 small data regime 容易 local optimum）

**Key equation to annotate:** Equation 1（BALD 的 mutual information formulation）——你需要能用一句話解釋給 PI 聽。

---

### 🔴 Yang, J. et al. (2025)
**"Active Learning-Assisted Directed Evolution (ALDE)"**
*Nature Communications* 16:XXXX.

**Why read / 為什麼讀：** 最新的、直接 comparable 的工作——在酶工程上系統比較了 BALD / UCB / Thompson sampling / greedy，用 deep ensemble 作為 UQ 方法。

**Role in our system / 在我們系統中的角色：**
- 驗證我們的 design choice（BALD + deep ensemble）在 protein engineering 場景下是合理的
- 他們的 control arm 設計（每輪保留 random 選擇的虛擬 list，最後做 retrospective replay）——我們直接複製這個設計
- 讀完記錄：他們的 AL vs random 在多少個 data point 後開始分叉？我們可以用這個設定期望值

**Key figure to annotate:** Figure 2（learning curves, AL vs random, different acquisition functions）。

---

### 🔴 Hie, B.L. et al. (2022)
**"Efficient evolution of human antibodies from general protein language models and their application to SARS-CoV-2"**
*Cell* 185(9):2038–2048.

**Why read / 為什麼讀：** 這是我們整個 project framing 的 published precedent——用 PLM + active learning 高效進化蛋白。他們用 ~50 個測量點在 antibody optimization 上超過 random sampling。

**Role in our system / 在我們系統中的角色：**
- 這篇是 wiki Motivation section 的核心 cite（「phage RBP 工程版的 Hie et al.」）
- 他們的 ESM embed → regression → acquisition → wet lab 閉環，跟我們的 pipeline 結構幾乎一一對應

---

### 🟡 Settles, B. (2009)
**"Active Learning Literature Survey"**
*Computer Sciences Technical Report 1648*, University of Wisconsin–Madison.

**Why read / 為什麼讀：** Active learning 的標準入門參考。Pool-based sampling、query strategies、BALD 的前身方法都在這裡。

**Role in our system / 在我們系統中的角色：** 背景理解用。如果你對 AL 概念不熟，這是最好的起點。讀 Section 2（Membership Query Synthesis）和 Section 3（Pool-Based Active Learning）即可。

---

### 🟡 Wittmann, B.J. et al. (2021)
**"Informed training set design enables efficient machine learning-assisted directed evolution"**
*Cell Systems* 12(11):1026–1045.

**Why read / 為什麼讀：** 綜述 ML-assisted directed evolution，包括 acquisition function 在 protein landscape 上的 exploration vs exploitation tradeoff 分析。

**Role in our system / 在我們系統中的角色：** 幫助理解 cycle 設計中「推薦 4–5 個 BALD + 1 個 random control」這個決定的理論依據。

---

## 📁 Module 08 — Wet Lab & Cycle Infrastructure

### 🔴 Gibson, D.G. et al. (2009)
**"Enzymatic assembly of DNA molecules up to several hundred kilobases"**
*Nature Methods* 6:343–345.

**Why read / 為什麼讀：** Gibson Assembly 的原始方法文章。我們所有 RBP variant cloning 都用 Gibson。

**Role in our system / 在我們系統中的角色：** 每 cycle 的 cloning 步驟直接依賴這篇描述的 5′ exonuclease + polymerase + ligase 組合。理解 overlap 設計（20–40 bp homology arm）有助於設計每個 variant 的 PCR primer。

---

### 🔴 Schäfer, A. et al. (1994)
**"Small mobilizable multi-purpose cloning vectors derived from the *Escherichia coli* plasmids pK18 and pK19"**
*Gene* 145(1):69–73.

**Why read / 為什麼讀：** pK18mobsacB（Addgene #87097）的原始文章。我們的受體敲除（ΔtonB / ΔexbB / ΔexbD）都用這個系統。

**Role in our system / 在我們系統中的角色：** sacB 反向選擇（sucrose lethal）+ kanamycin 正向選擇的雙交叉策略。理解這個設計有助於診斷電穿孔效率低或 sacB 選擇失敗的問題。

---

### 🟡 Boeckaerts, D. et al. (2024) *(同 Module 01)*
*(引用同上)* — **另一個讀這篇的原因（Module 08 角度）：**
他們用 whole-cell ELISA 測 phage RBP binding（heat-inactivated bacteria coating）的 protocol 是我們 ELISA assay 的直接參考。4PL fit 提取 EC50 的統計處理方式也在這篇 Supplementary 裡。

---

### 🟡 Latka, A. et al. (2021) *(同 Module 03)*
*(引用同上)* — **Module 08 角度的補充：**
他們的 ELISA binding assay 格式（His6-RBP + HRP-anti-His6 二抗）和我們的設計完全一樣，可以直接拿他們的 antibody concentration / wash buffer / blocking 條件當起始點。

---

## 📚 Background Context — Xanthomonas Biocontrol

### 🟡 Ryan, R.P. et al. (2011)
**"*Xanthomonas* genomics and molecular plant–microbe interactions"**
*Nature Reviews Microbiology* 9:344–355.

**Why read / 為什麼讀：** Xanthomonas 農業影響的綜述。30+ 種 Xanthomonas，感染 400+ 種植物，全球每年數十億美元損失。

**Role in our system / 在我們系統中的角色：** wiki Human Practices + Introduction section 的背景 cite。不需要精讀，annotate key numbers（species count、host count、economic impact）。

---

### 🟡 Iriarte, F.B. et al. (2018)
**"Combination of plant defense elicitors and bacteriophage for biocontrol of bacterial spot of tomato"**
*Frontiers in Plant Science* 9:1-12.

**Why read / 為什麼讀：** 噬菌體作為 Xanthomonas 生物防治的 field trial 依據（番茄細菌斑病）。

**Role in our system / 在我們系統中的角色：** 支持「噬菌體 biocontrol 在 Xanthomonas 上可行」這個 claim。annotate 一下 PFU dose、application timing、效果比較 copper 的數據。

---

### ⚪ Holtappels, D. et al. (2022)
**"The future of phage biocontrol in integrated plant protection for sustainable crop production"**
*Microbial Biotechnology* 15(3):597–610.

**Why read / 為什麼讀：** Xcc 上的 phage biocontrol 綜述，包含 Holtappels 組在 Xcc 的工作。

**Role in our system / 在我們系統中的角色：** 同上，background cite。可以只讀 Abstract + Discussion。

---

### ⚪ Farquharson, E.L. et al. (2021)
**"Phage resistance is driven by reduced infection efficiency of receptor mutants"**
*(journal to verify — likely Journal of Bacteriology or similar)*

**Why read / 為什麼讀：** T4 phage × *E. coli* reference collection——RBP 結合 85% 菌株，但只在 11% 上形成噬斑。這個 binding ≠ infection 的 confounder 是我們 Layer 2 因果驗證設計的科學依據。

**Role in our system / 在我們系統中的角色：** 解釋為什麼我們需要 WT / ΔReceptor / ΔDefense 三組對照——光測 ELISA binding 不夠，需要 plaque assay 解耦受體特異性 vs 防御系統貢獻。

---

## 🗂 Quick Reference — Print Priority

| Priority | Papers to print first |
|----------|----------------------|
| 🔴 Must-print | Wang 2003 · Houlsby 2011 · Yang 2025 ALDE · Lakshminarayanan 2017 · Hie 2022 Cell · Boeckaerts 2022 Viruses · Lin 2023 ESM-2 |
| 🔴 Module-specific | Lee 2009 (phiL7) · McNair 2019 (PHANOTATE) · Passaro 2025 (Boltz-2) · Abramson 2024 (AF3) |
| 🟡 Recommended | Boeckaerts 2024 NatComms · Greenman 2025 · Hie 2024 NatBiotech · Liu 2025 PLM-interact · Latka 2021 |
| ⚪ Reference | Settles 2009 · Hyatt 2010 · Gibson 2009 · Schäfer 1994 · Ryan 2011 |

---

*Document version: 1.0 · Last updated: 2026-05-11*
*配套文件：`docs/planning/iGEM_2026_项目大纲_中文版.md` § 附錄 A*
