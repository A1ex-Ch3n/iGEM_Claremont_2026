# Phage-Host Interaction Prediction 项目综合文档

> **目的**:这份文档彙整目前所有讨论的内容,作为后续 Claude Code 对话的 context。
>
> **状态**:Pre-experimental planning phase
>
> **目标 audience**:Project team members + Claude Code 后续对话

---

## 目录

1. [项目核心问题与初始想法](#1-项目核心问题与初始想法)
2. [iGEM Toronto 2025 (PHORAGER) 分析](#2-igem-toronto-2025-phorager-分析)
3. [现有学术界研究全景](#3-现有学术界研究全景)
4. [三个 Prediction Levels 的 Frontier 与 Gaps](#4-三个-prediction-levels-的-frontier-与-gaps)
5. [Wet Lab 可行性分析](#5-wet-lab-可行性分析)
6. [Motif-Level 详细方案](#6-motif-level-详细方案)
7. [Motif → Strain Prediction 与 KD 验证策略](#7-motif--strain-prediction-与-kd-验证策略)
8. [Data 需求与 ML Model 配置](#8-data-需求与-ml-model-配置)
9. [推荐方案与下一步](#9-推荐方案与下一步)

---

## 1. 项目核心问题与初始想法

### 我们的目标
建立一个 ML model 来 predict **phage-host interaction**。

### 核心痛点
现有 phage-host interaction 已知的 dataset 太少,无法支撑 ML 训练。

### 初始想法
能不能用 **protein-protein interaction (PPI) data** 来辅助/transfer到 phage-host interaction prediction?

### 关键 insight
PPI data 量远大于 phage-host interaction data。如果能 transfer learning,可解 data scarcity。

---

## 2. iGEM Toronto 2025 (PHORAGER) 分析

### Toronto 团队做了什么
PHORAGER (Phage Host-Optimized Receptor-Activated Generative Engineering Repository):
- **Task**: phage *engineering* (生成新的 RBP),不是 prediction
- **INPUT**: bacterial receptor (target)
- **OUTPUT**: 全新设计的 RBP sequence 能 bind 该 receptor
- **核心技术栈**:
  - ESM3 做 sequence generation
  - Boltz-2 做 structure prediction + affinity scoring
  - MCMC + simulated annealing 做 optimization
  - HADDOCK3 做 physics-based docking 验证
  - CodonTransformer 做 codon optimization

### 他们的实验
- 设计 4 个 E. coli phage 的 RBP:
  - Mu, P2 (target K12 → 改为 truncated R1 LPS)
  - Lambda, HK97 (target LamB → 改为 OmpC)
- Batch 1: ESM3 + Boltz-2 prototype, 实验合成 2 个 (P2_gpH_2, P2_gpH_8)
- Batch 2: 加上 MCMC + simulated annealing, 选了 9 个 cluster representatives

### Toronto 自己承认的限制
1. **Data scarcity**: RBP family 多样性导致 PLM generalize 不好
2. **Glycan modeling 差**: AlphaFold3 / Boltz-2 对 glycan 处理不好,4.4% chirality violation
3. **没有 wet lab quantitative validation**: 只做了 binary spot test,没量 binding affinity
4. **Resistance 没考虑**: pipeline 没 model 抗药性演化

### 我们 vs Toronto 的差异

| 维度 | Toronto PHORAGER | 我们 (proposed) |
|---|---|---|
| Task | 生成新 RBP | 预测 phage 是否 match host |
| Output | 全新 sequences | Probability / score |
| Data scarcity 处理 | 用 Boltz-2 (general biomol model) 当 oracle | 用 general PPI 做 transfer learning |
| Glycan 处理 | CHARMM-GUI + Boltz-2 (admittedly weak) | (open question - our opportunity) |

**关键观察**: 两个 project 互补,不冲突。我们的 prediction model 可以 validate Toronto 的 generated sequences。

---

## 3. 现有学术界研究全景

### A) RBP sequence + PLM 方法 (species/genus level)
- **Gonzales et al. (2023, PLOS One)**: ESM/ProtTrans embeddings → XGBoost → host genus, F1 ~63%
- **PHIStruct (2025)**: SaProt (structure-aware PLM) embeddings,low sequence similarity 时较好
- **PhageHostLearn (Boeckaerts, Nature Comms 2024)**: Klebsiella strain-level, ROC AUC 81.8%, **wet lab validated**
- **GE-PHI (2024)**: knowledge graph + ESM-2

### B) PPI feature 直接当 input (跟我们想法最接近!)
- **Coelho et al. (Sci Rep 2025)**: "A machine learning approach to predict strain-specific phage-host interactions"
  - 直接用 phage-bacteria PPI scores (从 PPIDM database) 当 ML feature
  - Salmonella + E. coli phages, accuracy 78–94%
  - **最接近我们想法的工作**
- **PBIP (2025)**: UniRep embeddings + CNN/Bi-GRU + attention, strain-level Klebsiella

### C) PPI knowledge transfer to virus-host
- **PLM-interact (Liu et al., Nature Comms 2025)** ⭐ **最重要**
  - Fine-tune ESM-2 用 *human* PPI data
  - Transfer 到 mouse/fly/worm/yeast/E. coli → AUPR 提升 16–28%
  - Transfer 到 virus-human PPI prediction 效果好
  - **但他们没有测 phage-bacteria!** ← 这是我们的 niche
- **Dong et al. (BMC Bioinformatics 2021)**: multitask transfer learning, virus-human PPI

### D) Docking / structure-based
- AlphaFold-Multimer (ipTM 当 PPI score)
- Boltz-2 (Toronto 用的, 自带 affinity module)
- HADDOCK3 (physics-based, 比 ML 慢但更精确)

---

## 4. 三个 Prediction Levels 的 Frontier 与 Gaps

### Mental Model

| Level | Question | Output |
|---|---|---|
| **Species/Genus** | Phage 属于 infect 哪个 bacterial taxa? | Multi-class classification |
| **Strain** | Phage X 会不会 infect 这个 *特定* clinical isolate? | Binary per (phage, strain) pair |
| **Motif** | RBP 哪一段决定 specificity? Mutate 哪个 residue 改变 host range? | Sequence/residue-level scoring |

### Level 1: Species/Genus

**代表方法**: HostPhinder (2016), WIsH (2017), VirHostMatcher (2019), vHULK (2022), RaFAH (2021), PHIStruct (2025), MoEPH (2025)

**Frontier**: Bjamalcz benchmarking (2025, Briefings in Bioinformatics) 系统重新评估,SaProt + 严格 benchmark

**Gaps**:
1. Long-tail distribution (rare taxa 没人能 predict)
2. Viral "dark matter" (RefSeq 只有 ~87% 有 host annotation)
3. Annotation 颗粒度太粗 (E. coli K-12 vs O157:H7 天差地别)
4. Performance plateau (~80-90% on well-curated data, polyvalent phages 测不准)
5. Productive infection vs adsorption 的 confound

### Level 2: Strain Level

**代表方法**:
- Beamud et al. (Cell Reports 2023) - Klebsiella, capsule-receptor
- **PhageHostLearn (Boeckaerts, Nature Comms 2024)** - Klebsiella, AUC 81.8%, lab validated
- Gaborieau et al. (Nature Microbiology 2024) - Escherichia 跨 strain
- Coelho et al. (Sci Rep 2025) - PPI feature
- PBIP (2025) - Klebsiella

**Frontier (2025-2026)**:
- **PAML / Phylogeny-Agnostic ML (Mutalik group, Nov 2025, bioRxiv)** - 跨 5 个 datasets, 128,357 interactions, 1058 strains, 560 phages, AUROC 0.67-0.94
- **PhageMind (Jan 2026, arXiv)** - meta-learning (MAML),少 data 就能 adapt
- **Briefings in Bioinformatics review (Feb 2026)** - 系统整理整个 strain-level field

**Gaps**:
1. Sparse + imbalanced PHI matrices
2. Assay-dependent labels (spot test vs plaque vs liquid kill curve)
3. **Most models are species-bound** (PhageHostLearn 只做 Klebsiella, 换种就要重练)
4. Cross-dataset evaluation 揭露泡沫 (within-dataset 0.9+, cross-dataset 掉到 0.6-0.7)
5. Negative data scarcity (没人写 paper 说"phage X 不能 infect strain Y")
6. Bacterial defense systems (CRISPR, R-M, anti-phage systems) 没被 model
7. Receptor diversity 假设单一类型 (PhageHostLearn 假设 K-locus)
8. Productive vs failed infection 区分

### Level 3: Motif Level

**代表方法**:
- PhageRBPdetect (Boeckaerts 2022) - HMM + Pfam, 只到 RBP level
- **SpikeHunter (Yang et al., GigaScience 2024)** - ESM-2 识别 tail spike, 发现 **C-terminal domain swap 是 key specificity determinant**, 230k+ tailspikes 跨 5 pathogens
- Schwarzer et al. / **Latka et al. (mBio 2021)** - 实验证明 Klebsiella RBPs modular: **N-term anchor + C-term specificity**, swap 可换 host range
- Garcia-Doval & van Raaij (PNAS 2012) - T7 phage gp17 tip BC/DE/FG/HI loops 决定 LPS recognition
- **Yehl et al. (Cell 2019)** - T7 deep mutational scanning, 1660 RBP variants × hosts → residue-level

**Frontier (2026 最新)**:
- **GenoPHI (Moriniere et al., April 2026, bioRxiv)** ⭐ - 1,050 genome-wide screens × 255 phages → 19 receptor classes, AlphaFold3 解到 **individual residue level**, 三个 modularity scales:
  - Gene-scale (整个 RBP gene 换)
  - Domain-scale (allelic variation)
  - Residue-scale (point mutations at binding interface)

**Gaps (这层最大空间!)**:
1. **Glycan binding motif 完全没好 model** (所有 frontier 集中在 protein receptor)
2. **Convergent evolution → no homology signal** (打破 alignment-based 方法)
3. Loop region 是 specificity hotspot 但 disordered, AlphaFold confidence 低
4. Lack of paired interface data (PDB 里 phage RBP-receptor complex 可能 <100 个)
5. DMS 太贵 (Yehl 只做了 1 个 phage 1 个 domain)
6. 没有 standardized motif annotation
7. Interpretability gap (PhageHostLearn 给 score 但不告诉哪段在做事)
8. Motif modularity cross-talk 没 model 过

### 三个 Levels 横向对照

| 维度 | Species/Genus | Strain | Motif |
|---|---|---|---|
| Data availability | 高 (10k+ annotations) | 中 (~100k pairs in 5 datasets) | 极低 (~100 cocrystals, 1 saturation scan) |
| SOTA accuracy | 80-90% (within-distribution) | 0.7-0.9 AUROC, drops cross-dataset | 没清楚 benchmark |
| Clinical utility | 太粗,无用 | Critical | 重要 for engineering |
| Glycan vs protein | 不分 | 大致分 | **完全没解** |
| Generalization gap | rare taxa 失效 | 跨 species 失效 | 跨 phage family 失效 |
| Frontier 2025-26 | SaProt + benchmark | PAML, PhageMind | GenoPHI (genome→residue) |
| **Open problem** | Dark matter coverage | Phylogeny-agnostic generalization | **Glycan motifs + interpretability** |

---

## 5. Wet Lab 可行性分析

### Wet Lab × ML Contribution 矩阵

| 方向 | Wet lab 难度 | 设备 | 时间 | ML contribution | 整体 |
|---|---|---|---|---|---|
| Species/genus prediction | ⭐ 现成 dataset | 无 (in silico) | 短 | 低 (saturated) | ⚠️ |
| Strain-level single species | ⭐⭐⭐ strain library + spot test | 中 (BSL-1) | 3-6 月 | 中 | ✅ |
| Strain-level cross-genus | ⭐⭐⭐⭐⭐ 多种 BSL-2 病原 | 高 (BSL-2) | >6 月 | 高 | ❌ 太硬 |
| Motif - DMS | ⭐⭐⭐⭐⭐ NGS | 高 | >1 年 | 很高 | ❌ 太硬 |
| Motif - chimera/domain swap | ⭐⭐⭐ cloning + spot test | 低-中 | 2-4 月 | 高 | ✅✅ |
| **Motif - purified RBP + ELISA** | ⭐⭐ 蛋白表达 + ELISA | 低-中 | 2-3 月 | 中-高 | ✅✅✅ |

### 推荐方向: Purified RBP binding assay + Motif-level prediction

**理由**:
- **Wet lab**: 全 BSL-1, 标准分子生物 + 蛋白纯化 + 96-well ELISA, 6 月内完成
- **Data**: 产生 ~150-500 quantitative binding pairs (这在 phage field 算大)
- **ML**: 直接 attack 三个 open gaps (quantitative motif data, PPI transfer, glycan binding)
- **Differentiation**: 跟 Toronto 互补 (他们 generate, 我们 predict + validate)

---

## 6. Motif-Level 详细方案

### 三种 Motif Granularity

| Granularity | 指什么 | 例子 |
|---|---|---|
| Domain-scale | 整个 functional 模组 (~100-300 aa) | RBP 的 C-term specificity domain vs N-term anchor |
| Loop-scale | functional region (~10-30 aa) | T7 gp17 的 BC/DE/FG/HI loops |
| Residue-scale | 单一 aa 或几个 aa 组合 | T7 gp17 Ala518, Asp520, Val544 |

### Motif Level 能做的四种 Contribution

#### Contribution 1: Motif → host specificity 的 mapping function
- 现有 ML 都是 *whole RBP → host* black box
- 我们做 quantitative motif-level functional mapping (3-5 RBP scaffolds)
- Output: (RBP_variant, receptor, binding_score) quantitative table

#### Contribution 2: Validate (or break) 现有 generative model
- PHORAGER, RFdiffusion, Boltz-2 affinity 都没系统 wet-lab validate motif-level prediction
- 我们的 wet lab dataset 当 ground truth, benchmark:
  - ESM-2 / SaProt embeddings + regression
  - Boltz-2 zero-shot affinity prediction
  - AlphaFold3 ipTM
  - PLM-interact (PPI prior fine-tune)

#### Contribution 3: Motif-level interpretability baseline
- ML attention/saliency 预测重要 residue → wet lab mutate → 看 binding loss
- 教科书级别的 interpretability validation

#### Contribution 4: 第一个 quantitative RBP-glycan binding dataset
- 用 commercial purified LPS (Sigma E. coli O111:B4, O55:B5, K-12) coat plate
- 测 RBP variants 的 binding (Kd / EC50)
- 这个 dataset 学界直接缺的

### 四种 Wet Lab Strategy (按可控性排序)

#### Strategy A: Truncation series ⭐⭐⭐⭐⭐
- 例子: P2 gpH (1-300, 1-400, ..., 100-700, 200-700)
- 每个 = 1 PCR + cloning + expression + ELISA
- **优点**: 单变因 (length); 文献有支持 (N-term anchor, C-term specificity)
- **风险**: trimerization 失败 (对策: 加 GCN4 leucine zipper); inclusion body (对策: SUMO/MBP tag)
- **ML 用途**: minimal binding domain boundary

#### Strategy B: Chimera / domain swap ⭐⭐⭐⭐
- 例子: gpH N-term + gp17 C-term, gp17 N-term + gpH C-term
- **优点**: 直接验证 modularity hypothesis; Latka 2021 已证明 work
- **风险**: junction site 选择 (对策: AlphaFold predict junction, 选 flexible loop); trimerization 跨来源 incompatible (对策: 选 closely related phages)
- **ML 用途**: modularity prior

#### Strategy C: Targeted point mutation ⭐⭐⭐⭐
- 例子: T7 gp17 A518S, D520N, V544I 或 alanine scan loop
- **优点**: 单一变因,可量化 (ΔΔG); literature 有 reference
- **风险**: hotspot 选错没 effect; single mutation 效应小; trimer cooperativity
- **ML 用途**: validate ML attention; mini-DMS

#### Strategy D: Generated novel sequences ⭐⭐
- ESM-3 generate (像 Toronto), 选小 perturbation 版本
- **风险**: 多变因同时改; success rate 低
- **ML 用途**: 跟 Toronto 互补 (predict before generate)

### 强烈建议: A + B + C 混合做

| Strategy | 给 ML 的资讯 |
|---|---|
| A (Truncation) | Domain boundary |
| B (Chimera) | Modularity prior |
| C (Point mutation) | Residue-level resolution |

### 推荐 RBP Candidates

1. **T7 gp17 tip domain (470-553, ~85 aa)** ⭐⭐⭐⭐⭐
   - Garcia-Doval & van Raaij 2012 published expression protocol
   - 结构已解 (PDB: 4A0T)
   - Receptor: E. coli K-12 LPS
   - Monomeric tip 可表达 in BL21
   - **首选**

2. **P22 Gp9 (full RBP, ~666 aa)** ⭐⭐⭐⭐
   - Steinbacher 1996 结构 (PDB: 1TYU)
   - Receptor: Salmonella O-antigen
   - Sigma 商品 antibody → ELISA detection 方便
   - 缺点: trimer

3. **P2 gpH or HK620 RBP** ⭐⭐⭐
   - 跟 Toronto 同 phage scaffold,可直接对比
   - Receptor: E. coli K12 (gpH) / O18 (HK620)

4. **Klebsiella phage K11 RBP (KP32-like)** ⭐⭐⭐
   - Latka mBio 2021 已做过 chimera
   - Receptor: K. pneumoniae capsule
   - 缺点: BSL-2 host (但 RBP 本身可在 E. coli 表达)

---

## 7. Motif → Strain Prediction 与 KD 验证策略

### 关键问题
有没有人做 motif → strain 的 hierarchical prediction?

### 答案: 几乎没人

| Paper | 做了什么 | 距离 motif→strain 多远 |
|---|---|---|
| PhageHostLearn (Boeckaerts) | 提出 layered concept,只做 layer 1 (whole RBP) | 提 layered 但没做 motif level |
| Gaborieau (Nature Microbiol 2024) | Linear mixed model on Escherichia, 用 whole RBP | 没到 motif |
| **Yehl/Voigt 2024 multiobjective** | T7 RBP DMS × multiple strains | 最接近,但只 T7 |
| Toronto PHORAGER | 完全没做 motif → strain | - |

**结论**: 这是个 niche,但 niche 存在有原因 (见下方 caveat)

### KD 实验为什么是强招

KD 给的是 **causal evidence**, 不是 correlation。

```
INPUT:   (RBP_motif_sequence, strain_genome)
OUTPUT:  predicted infection_score
GROUND TRUTH from KD: 
   - WT strain has receptor → infection = X
   - KD strain (receptor knocked down) → infection = Y
   - ΔX = X - Y = receptor-mediated contribution
```

**为什么 causal data 对 ML 重要**:
1. 大部分 phage-host dataset 是 observational
2. Strain 之间有上百基因差异,confounding 严重
3. KD 把 confounding 拿掉,**同一 strain 只差一个 gene**
4. 给 ML 真正的 negative control + receptor-specific signal + epistasis 数据

### 文献中的 KD precedent
- **Mutalik lab 2020 (PLOS Biology)**: CRISPRi 在 E. coli K-12 vs BL21 KD, 14 phages
- **PHAGEPACK (2024 bioRxiv)**: dCas9 + sgRNA library 大规模 KD
- **K. pneumoniae 2025 paper**: knockout wbaP, wbaZ, wzc → adsorption rate 99% → 6-12%
- **U136B + tolC**: knockout tolC → EOP 1 → <10⁻⁷ → complement 回 tolC → 恢复

### ⚠️ CRITICAL CAVEAT (必须面对!)

引用 Doud Frontiers 2025 review:
> "phage T4 RBP **bound to 85% of the 72 strains** in an E. coli reference collection, yet T4 phage **only formed plaques on 11% of the collection** (Farquharson et al., 2021). Receptor-RBP binding may be necessary for infection, but **it is not sufficient**."

```
binding affinity ≠ adsorption ≠ DNA injection ≠ productive infection ≠ lysis
```

Multi-layer decoupling:
1. Binding ≠ adsorption (spatial geometry, ion conditions)
2. Adsorption ≠ DNA injection (secondary receptor, membrane potential)
3. Injection ≠ replication (CRISPR-Cas, R-M, Abi)
4. Replication ≠ lysis (lysogeny vs lytic)

**Motif-level binding 对 strain-level outcome 的 predictive power 可能只解释 20-40% variance**。

**危险**: naïve "motif binding → infection" model accuracy 不会高
**机会**: KD 实验可以 **解构** 这个 multi-layer

### 推荐实验设计: Two-axis Perturbation

```
                Strain WT    Strain ΔReceptor   Strain ΔDefense
RBP WT          [A]          [B]                [C]
RBP motif var1  [D]          [E]                [F]
RBP motif var2  [G]          [H]                [I]
```

每 cell 量两个 readout:
- **In vitro binding** (ELISA)
- **In vivo infection** (spot test / liquid kill curve)

| Comparison | 测量什么 |
|---|---|
| A vs B | Receptor 对 binding/infection 的贡献 |
| A vs D | Motif variation 的贡献 |
| (A vs B) vs (D vs E) | Motif × Receptor epistasis (ML 该学的) |
| A vs C | Defense system 影响 (binding 不变, infection 变) |
| B vs C | 纯粹 downstream effect |

### Layered Model 架构

```
Layer 1 (我们 train + validate):
   f_binding(motif_seq, receptor_seq) → binding affinity score
   ground truth: ELISA

Layer 2 (我们 calibrate):
   f_infection(binding_score, defense_features, strain_metadata) → infection score
   ground truth: KD vs WT spot test
```

- **Layer 1 是 main contribution** (motif-level affinity prediction with PPI prior)
- **Layer 2 是 calibration** (binding 不够,还要 +defense factors)

### KD Wet Lab Feasibility

#### KD 在 E. coli 是 trivial
- **CRISPRi (dCas9 + sgRNA)**: standard tool, Mutalik lab 已 systematically 用过
- **Keio collection**: ~3985 个 E. coli K-12 single gene KO, 直接买不用自己做!
- Targets well-defined:
  - LPS biosynthesis: waaO, waaR, waaJ, waaY, wbbL
  - OMPs: ompC, ompF, ompA, lamB, btuB, fhuA, tolC, fadL (都在 Keio)
  - Capsule: E. coli K12 没有, 用 Klebsiella (BSL-2)

#### 推荐 KO Panel (5-6 个 E. coli K-12 KO)
- WT (BW25113 background)
- ΔlamB (lambda receptor)
- ΔompC (Mu, OmpC-binding phage receptor)
- ΔompF
- ΔtolC (U136B receptor)
- ΔfhuA (T1, T5 receptor)
- ΔwaaO 或 ΔwaaJ (LPS R-core truncation)

对应 phages: lambda, T7, T4, T5, P2, Mu

#### 实验时间
- 取得 Keio: ~$500-1000 + 跨国 distributor 时间
- 验证 KO efficiency: colony PCR / RT-qPCR (1 周)
- Spot test panel: 96-well 一天可测 6 strain × 8 phage × 3 dilution
- **整体 KO 部分 1 个月就能跑完**

### 可测试 Hypothesis

> **H**: 对 motif-level RBP variants, in vitro binding affinity (ELISA) 跟 in vivo infection (spot test on Keio KO panel) 的 correlation 随 receptor 而异。
> - 对 protein receptor (OmpC, LamB, FhuA): correlation 高 (>0.7)
> - 对 glycan receptor (LPS): correlation 中等 (~0.4-0.6)
> - 对 capsule receptor: correlation 低 (需要 depolymerase activity)

**这个 hypothesis 验证后直接告诉 ML field**: 什么 receptor type 适合 binding-based prediction, 什么不适合 → publication-worthy finding。

---

## 8. Data 需求与 ML Model 配置

### Input/Output 定义

#### Layer 1 (主要 contribution)
```
INPUT: (RBP_motif_sequence, receptor_sequence_or_structure)
OUTPUT: binding affinity score (regression)
GROUND TRUTH: ELISA quantitative readout
```

#### Layer 2 (calibration)
```
INPUT: (Layer 1 binding score, strain features)
OUTPUT: infection score (binary or regression)
GROUND TRUTH: spot test on KD strains
```

### Input Representation 详细

| 元素 | Representation |
|---|---|
| RBP motif sequence | Amino acid string → ESM-2 embedding (1280-dim per residue, 或 mean-pooled) |
| Receptor (protein) | Amino acid string → ESM-2 embedding |
| Receptor (glycan) | SMILES → GIN/ChemBERTa, 或 CHARMM-GUI 3D + graph |
| Binding score | Float (Kd, EC50, normalized ELISA OD) |

最简单 baseline: `concat(RBP_embed, receptor_embed) → MLP → score`

### Data Scale 需求

#### 参考 PhageHostLearn (GitHub)
> "Our dataset of around 100 phages and 200 bacteria took 5-6 hours on an 8-core Apple M1"
> "no specialized GPU hardware is strictly needed"
> AUC 0.82 on Klebsiella

#### 我们的 case 因为有 *structure* + *transfer learning*, 需要更少

| 资源 | 数量 | 理由 |
|---|---|---|
| RBP scaffolds | 4-5 个 (P2 gpH, T7 gp17, lambda gpJ, Mu gp49) | Cross-family generalization |
| Variants per scaffold | 8-12 (truncations + chimeras + mutations) | Motif × receptor epistasis |
| **Total RBP variants** | **32-60 个** | |
| Receptor sources | 6-8 个 (LPS chemotypes + OmpC + LamB + FhuA + cells from Keio panel) | |
| **Total (RBP, receptor) pairs** | **200-500 数据点** (3 replicate ELISA each) | |
| KD strain panel for layer 2 | 6 个 (WT + 5 KO from Keio) | |
| Phage × strain spot test | 4-5 phage × 6 strain = 24-30 pairs | Layer 2 calibration |

**Bottom line**: ~300 quantitative pairs Layer 1 够; ~30 binary spot test Layer 2 够。

#### 这够吗?
- **不够** train 全新 deep model from scratch
- **够** fine-tune *pretrained model* (ESM-2 + PLM-interact)
- 只 fine-tune output layer 通常 100-500 examples 就够 → **transfer learning 的核心优势**

#### 与其他 paper 比较

| Paper | Training pairs | Result |
|---|---|---|
| PhageHostLearn | ~100 phage × ~200 bact (sparse) | AUC 0.82 |
| Coelho 2025 | 13 phage × ~50 strains | accuracy 78-94% |
| Yehl 2019 (T7 DMS) | 1660 RBP variants × 3 hosts | regression, R² varies |
| Toronto PHORAGER | 实验只验证 ~10 generated variants | 没真 train ML |
| **我们 (target)** | **~300 quantitative pairs** | **proof-of-concept 可行** |

### 电脑配置需求

#### Layer 1 model (ESM-2 + regression head)

| Component | 规格 |
|---|---|
| GPU | NVIDIA RTX 3060/3090, 或 Colab/Kaggle 免费 GPU |
| GPU 记忆体 | 12 GB 跑 ESM-2 650M; 24 GB 跑 ESM-2 3B |
| RAM | 16-32 GB |
| 硬盘 | 100 GB |
| 运算时间 | Fine-tune: 30 分钟 - 几小时; Pretrain: 1-2 天 |

#### ESM-2 Model Sizes

| Size | Hardware | 适用 |
|---|---|---|
| ESM-2 8M | CPU | toy experiments |
| ESM-2 35M | Laptop | quick prototyping |
| ESM-2 150M | 一般 GPU | reasonable baseline |
| **ESM-2 650M** | **12GB GPU** | **推荐起手** |
| ESM-2 3B | 24GB GPU | 进阶 |
| ESM-2 15B | 多 GPU | 不需要 |

#### Boltz-2 / AlphaFold3
- 只做 inference 不 train
- Boltz-2: ~30 秒-几分钟 per pair on 1 GPU
- 300 pairs × 30 秒 = ~3 小时 (一晚跑完)

#### Layer 2 model
- 完全不用 GPU
- XGBoost / logistic regression on CPU 几分钟跑完

### 两阶段 Training 策略

#### Phase 0 (开始就做): Pretrain on 公开 data
- ~5,000 phage RBP sequences from NCBI (PhageRBPdetect)
- ~150 protein-protein complex from PDB
- ESM-2 / PLM-interact 当 backbone
- Contrastive pretraining 或 masked language modeling
- **完全不需要 wet lab data**

#### Phase 1 (wet lab 开始后): Fine-tune
- 拿 Phase 0 model
- Fine-tune last layer 用 ELISA 数据 (~300 pairs)
- **30 分钟到几小时**搞定

### Critical Data 设计原则

1. **Quantitative,不是 binary** - 保留 ELISA OD 完整 numerical readout
2. **必须有 negative pairs** - WT P2 gpH × R1 LPS (negative), WT P2 gpH × K12 LPS (positive)
3. **Replicate 跟 control** - 至少 3 replicate, BSA 当 negative, lectin 当 positive
4. **Train/test split 要小心** - split by RBP scaffold 或 receptor (不要随便切!)

### MVP First (强烈建议)

#### MVP 阶段 (1-2 个月)
- 1 个 RBP scaffold (P2 gpH, 因为 Toronto 有 Boltz-2 prediction 可对比)
- 6 个 variants (WT + 2 truncations + 2 point mutants + 1 chimera)
- 4 个 receptors (K12 LPS + R1 LPS + OmpC + BSA control)
- → **24 pairs × 3 replicate = 72 data points**
- Train simple ESM-2 + MLP, 看 baseline performance

#### MVP 跑完后再 scale up
- 如果 MVP 看到 signal → 扩到 4-5 RBP × 8-10 variants × 6-8 receptor
- 如果 MVP 完全没 signal → **debug 不要硬上更多 data**

### 最常踩的坑

1. **ELISA 数值变异大** → 同 batch 同天测同 RBP 的所有 variants; 设 internal standard
2. **所有 RBP variants 都不 express 或 mis-fold** → 先验证 WT 能 express 再做 variants; 准备 backup scaffold
3. **ML 学到 "RBP scaffold identity" 而非 motif signal** → 设计 chimera (P2 N + T7 C) 强迫区分
4. **Binding 跟 infection 混淆 (Farquharson)** → 严格 separate Layer 1 跟 Layer 2

---

## 9. 推荐方案与下一步

### 最终推荐方向

**Motif-level RBP-receptor binding prediction with PPI knowledge transfer + KD-based strain-level validation**

### 为什么这个方向赢

| 角度 | 理由 |
|---|---|
| Wet lab 可操控性 | Strategy A/B/C 都是单变因 + standardized cloning + ELISA + BSL-1 + 6 月内完成 |
| Data 量 | 25-50 RBP variants × 5-10 receptors = 125-500 quantitative pairs (在 phage field 算大) |
| ML 跟 SOTA 直接比 | ESM-2, Boltz-2, AlphaFold3, PLM-interact 的 zero-shot 都可直接 benchmark |
| 填的 gap | Quantitative motif data 缺 + glycan binding 缺 + ML interpretability 缺 + PPI transfer to phage 没人做 |
| 跟 Toronto/iGEM scale 相容 | 不需 BSL-2, 不需 NGS, 不需 cryo-EM |
| Differentiation | 跟 Toronto 互补 (他们 generate, 我们 predict), 跟 PhageHostLearn 互补 (他们 strain binary, 我们 motif quantitative) |
| KD 验证给 causal claim | 这个方向 unique 的 advantage |

### Project Hook

**"Motif-level quantitative phage RBP binding atlas + ML prediction benchmark + KD-validated strain-level inference"**

### 三层 Contribution

1. **Wet lab dataset**: 第一个 quantitative motif-level RBP-receptor binding atlas (含 glycan)
2. **ML benchmark**: 第一个系统 evaluate ESM-2/Boltz-2/AlphaFold3/PLM-interact 在 phage RBP context
3. **Causal validation**: 第一个用 KD 把 binding signal 跟 downstream defense 解构,做出 layered interpretable model

### 6-Month Plan

#### Month 1-2: Setup + In silico phase 0
- 选 4-5 个 RBP, design 25 个 constructs (truncation + chimera + mutation)
- Cloning 同时进行
- **同时开始** in silico baseline:
  - 跑 ESM-2 / Boltz-2 zero-shot prediction on all designed constructs
  - Pretrain PLM-interact on phage RBP corpus
  - Set up benchmark pipeline

#### Month 3-4: Wet lab 主体
- Express + purify ~25 RBP variants
- ELISA against 5-10 receptor sources (purified LPS + OMPs + cells)
- 取得 Keio KO collection, 验证 KO efficiency
- → ~150-250 quantitative binding data points

#### Month 5: ML model + KD validation
- Build Layer 1 model: train pipeline, compare 4 strategies
  - ESM-2 + MLP (baseline)
  - PLM-interact pretrain → fine-tune (PPI transfer)
  - Boltz-2 zero-shot
  - Your dataset only (no pretrain)
- Spot test panel on Keio KO strains (Layer 2 ground truth)
- Build Layer 2 calibration model

#### Month 6: Validation + write up
- Hold out RBP variants 做 blind test
- 预测 + 验证一个全新 RBP-receptor pair
- Two-axis perturbation analysis (motif × receptor KO)
- 写 paper / iGEM wiki

### 后续可深入讨论的问题

(供 Claude Code 后续对话用)

- [ ] 详细 ESM-2 fine-tuning code template
- [ ] MVP 阶段具体实验 layout (哪 6 个 variants, 哪 4 个 receptors, 怎么 split)
- [ ] 具体 ELISA optimization protocol (positive/negative control, signal/noise)
- [ ] Keio collection 获取 protocol + 替代方案
- [ ] Layered ML architecture 详细设计 (具体 components, training schedule)
- [ ] PLM-interact 架构 adaptation 细节
- [ ] 哪些 phage-host datasets 公开可用
- [ ] Glycan-protein interaction prediction 最新方法
- [ ] Off-target KD effect 的 control 实验
- [ ] Codon optimization for RBP expression
- [ ] Trimerization 失败的 backup plan (GCN4 leucine zipper 等)
- [ ] BLAST + clustering 流程做 sequence novelty validation
- [ ] Cross-validation strategy (split by scaffold vs by receptor)
- [ ] Negative control 设计 (BSA, irrelevant lectin, non-binding RBP)
- [ ] Boltz-2 inference 的 batch processing 优化

---

## 附录: 关键 References

### Phage-host prediction (ML side)
- Boeckaerts et al. (2024) PhageHostLearn. *Nature Communications*. https://www.nature.com/articles/s41467-024-48675-6
- Coelho et al. (2025) PPI-based phage-host prediction. *Scientific Reports*.
- Liu et al. (2025) PLM-interact. *Nature Communications*.
- Mutalik et al. (2020) CRISPRi phage host range. *PLOS Biology*.

### Phage RBP biology
- Yehl et al. (2019) T7 DMS. *Cell*.
- Latka et al. (2021) Klebsiella RBP modularity. *mBio*.
- Garcia-Doval & van Raaij (2012) T7 gp17 tip. *PNAS*.
- Schwarzer et al. RBP swap.

### Generative protein design
- Abramson et al. (2024) AlphaFold3. *Nature*.
- Passaro et al. (2025) Boltz-2. *bioRxiv*.
- Hayes et al. (ESM3).

### iGEM Toronto 2025
- Wiki: https://2025.igem.wiki/toronto/model/

### Frontier 2025-2026
- PAML (Mutalik group, Nov 2025, bioRxiv) - Phylogeny-Agnostic ML
- PhageMind (Jan 2026, arXiv) - meta-learning
- GenoPHI (Moriniere et al., April 2026, bioRxiv) - genome-wide → residue level
- Briefings in Bioinformatics review (Feb 2026)
- Doud et al. (2025) Frontiers in Cellular and Infection Microbiology - phage therapy AI review

---

**文档版本**: v1.0
**最后更新**: 2026/05/07
**作者状态**: Pre-experimental planning,待 PI/team 讨论后定案
