# iGEM Claremont 2026 — Progress Report
# 項目進度報告

**Date / 日期：** 2026-05-11
**Branch：** active-learning-pipeline
**Author：** Alex Chen

---

# ENGLISH VERSION

## What We Built This Week

### May 7–8: Project Pivot + Overnight Parallel Build

The project pivoted from a 6-factor biophysical scoring pipeline to a **closed-loop active learning system**. Seven AI agents worked in parallel overnight (May 7→8), each building one pipeline module from scratch. By morning, Modules 00–06 all had working code, test suites, and READMEs.

Post-build fixes applied on May 10:
- Re-downloaded phiL7 genome from NCBI (contaminated headers fixed)
- Replaced two invalid bacteria accessions (KY000037, PY746849) with valid assemblies
- Gitignored 630 MB of genome binaries; re-downloadable via fetch scripts
- Set up CARC Laguna HPC environment (boltz2 conda env, 8 debugging iterations to resolve torch/CUDA/trifast version conflicts)

### May 10: First Boltz-2 GPU Run

Successfully ran Boltz-2 on CARC Laguna (NVIDIA L40S, job 59949, 47 seconds). **However:** the run used phiL7 P25 (85 aa) instead of rbp_01 (712 aa) — the wrong protein. ipTM = 0.345 was discarded.

### May 11: Literature Audit + Correct Boltz-2 Run

Read all 19 core papers and cross-checked every quantitative claim in the project documents. Found and corrected 5 factual errors:

1. **ExbD2 NOT essential** (Hung 2003): only TonB, ExbB, ExbD1 are required for phiL7 infection — corrected across all planning docs
2. **Boltz-2 affinity head = small molecule only** (Passaro 2025): protein-protein pairs output NaN; we use ipTM as proxy
3. **Greenman 2025 journal = PLoS Comput Biol** (not NAR Genomics); conclusion is "no single best UQ method"
4. **Hie 2024 used ESM-1b/1v** (not ESM-2); ~20 variants per antibody (not ~50)
5. **Lee 2009 does not name a tail spike** — rbp_01 identified by PhageRBPdetect HMM (Tail_spike_N domain), not Lee annotation

Re-ran Boltz-2 correctly (job 59986) after fixing the FASTA and clearing the processed cache. Results: **ipTM = 0.365, chain A ptm = 0.808**.

Module 04 re-run with real RBP sequences (3 RBPs, not 5 with mocks). Notes file for all papers written: `docs/reference/paper_reading_notes.md`.

---

## What Each Module Does

| Module | Role in the Pipeline | Key Output |
|--------|---------------------|------------|
| **00 Raw Data** | Genome library for training data | 777 phage + 34 bacteria genomes |
| **01 Ground Truth** | Labeled interaction pairs (known phage-host pairs) — the only "ground truth" in the system | interaction_matrix.csv: 2,236 pairs |
| **02 Annotation** | Translates raw DNA into protein sequences | phiL7: 80 ORFs; Xcc: 4,344 ORFs |
| **03 RBP ID** | Finds the "key" protein (RBP) on each phage that determines which host it can infect | rbp_01 (712 aa): primary tail spike candidate |
| **04 Embedding** | Converts protein sequences into vectors that a neural network can process | 1280-dim ESM-2 embeddings (on Laguna) |
| **05 Structure** | Predicts 3D structure of RBP bound to receptor; gives a "structural confidence prior" for the model | PDB file + ipTM score |
| **06 Ensemble** | Trains 5 independent neural networks; their disagreement = uncertainty = what the model doesn't know | (predicted score, uncertainty) per variant |
| **07 BALD** | Uses uncertainty to pick the next experiment: select the variant the model is most uncertain about | Ranked variant list for wet lab |
| **08 Cycle Data** | Ingests wet lab ELISA results, triggers retraining, closes the loop | Per-cycle ELISA data + model checkpoints |

---

## Validation Strategy Analysis

There are three practical validation tiers. Here is an honest comparison:

### Tier 1 — ELISA Only

**What you do:** Express His6-RBP variants, measure binding to heat-inactivated Xanthomonas cells (10⁸ CFU/well), fit 4PL curve → EC50.

**Dry lab required:** Module 06 (ensemble prediction of binding score) + Module 07 (BALD selects which variants to test).

**Biological question answered:** "Which rbp_01 variants bind TonB-expressing Xcc cells more or less strongly than wild-type?"

**What it cannot answer:** Whether binding leads to productive infection. Whether the binding signal is receptor-specific or non-specific. T4 phage RBPs bind 85% of *E. coli* strains but only infect 11% (Farquharson 2021) — ELISA alone cannot distinguish.

**Time / Cost:** 1–2 weeks per cycle, ~$400 reagents per plate batch. **This is what Cycle 0–2 actually runs.**

**Story completeness:** ⭐⭐⭐ — "We found variants that bind better." Clean quantitative story but leaves the mechanism open.

---

### Tier 2 — ELISA + Plaque Assay (WT strain only)

**What you do:** Same ELISA as Tier 1, plus spot assay of phiL7 on the Xanthomonas lawn to confirm lytic infection.

**Dry lab required:** Same as Tier 1; adds a binary "infects / doesn't infect" data layer.

**Biological question answered:** "Does higher ELISA binding correlate with productive infection?" Begins to address the binding ≠ infection confound.

**Time / Cost:** Adds 1–2 days per batch (plaque assay is simple, protocol already exists on Benchling). Near-zero incremental cost.

**Story completeness:** ⭐⭐⭐⭐ — "We found variants that bind better AND confirmed infection." Reviewer-grade validation. **This is the minimum for a credible iGEM story.**

---

### Tier 3 — ELISA + Plaque on WT + ΔtonB + ΔexbB + ΔexbD1

**What you do:** Generate markerless deletion strains (pK18mobsacB system, 4–6 weeks). Then run ELISA binding and plaque assay on each knockout strain alongside wild-type.

**Dry lab required:** Tier 1 + Layer 2 causal decomposition. Comparing: WT vs ΔReceptor binding and infection rates.

**Biological question answered:**
- X − Y = receptor-specific contribution to binding (WT binding minus ΔtonB binding)
- Whether the model's binding score correctly captures *receptor-specific* affinity rather than non-specific surface sticking
- Quantitative decomposition: "X% of the binding signal comes from TonB specificity"

**ExbD2 note:** Hung 2003 confirms ΔexbD2 does NOT affect phiL7 infection — making it a built-in negative control that validates the whole experiment.

**Time / Cost:** 4–6 weeks to generate knockouts (electroporation + double-crossover selection + PCR verification). ~$95 for pK18mobsacB plasmid. Medium technical difficulty — Xanthomonas electroporation protocol already exists on Benchling.

**Story completeness:** ⭐⭐⭐⭐⭐ — "We not only predicted and measured binding, but causally showed that the binding signal is receptor-specific." This is a **Nature-quality mechanistic claim** that most iGEM teams cannot make.

---

### Tier 4 (Optional) — Adding ΔDefense knockouts

**What you do:** Also delete the anti-phage defense systems (CRISPR, restriction-modification).

**Biological question answered:** Full decomposition: receptor contribution vs defense contribution. Allows the claim "X% comes from receptor recognition, Y% from defense evasion."

**Reality check:** CRISPR and RM loci in Xcc are not as well-characterized as tonB/exbB. This adds significant molecular biology work with uncertain return. **Not recommended for iGEM 2026 timeline** — save for a follow-up paper.

---

### Recommendation

| Tier | Feasibility for iGEM 2026 | Story | Verdict |
|------|---------------------------|-------|---------|
| ELISA only | ✅ Easy | ⭐⭐⭐ | Minimum viable |
| + Plaque assay | ✅ Easy (+2 days) | ⭐⭐⭐⭐ | **Do this** |
| + ΔtonB/ΔexbB/ΔexbD1 | 🟡 Medium (4–6 weeks) | ⭐⭐⭐⭐⭐ | **Do this if wet lab starts May 17** |
| + ΔDefense | ❌ Hard, uncertain | — | Post-paper only |

**Bottom line:** Start the receptor knockouts on May 17. If they succeed by early July, you have a five-star story. If they fail, you fall back to Tier 2 and still have a strong submission. The ELISA + plaque + receptor knockout combination is exactly what Hung 2003 did — and that paper was published in *Biochemical and Biophysical Research Communications*, which validates the experimental feasibility.


---

## Discussion-Ready Outputs (as of 2026-05-11)

> For use in PI meetings, team reviews, or paper planning.
> All files listed are committed on branch `active-learning-pipeline`.

### 1. rbp_01 × TonB 3D Complex Structure
**File:** `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/boltz_results_*/predictions/*.pdb`
**View with:** PyMOL or ChimeraX (free)
**What to discuss:**
- Chain A (rbp_01, 712 aa): Boltz-2 predicts a well-folded monomer (chain A ptm = **0.808**)
- Interface confidence (ipTM = **0.365**): low, meaning the model is uncertain HOW the two proteins dock — not that they don't interact. This is expected for a novel phage-host pair with no PDB template.
- Key question for PI: does the predicted interface region match the known TonB-ExbB binding face from *E. coli* literature?

### 2. Structural Confidence Numbers
**File:** `05_structure_prediction/outputs/boltz2/.../affinity.json`
```
confidence_score: 0.683   (overall complex quality)
interface_ipTM:  0.365    (binding interface confidence — low but expected)
chain_A_ptm:     0.808    (rbp_01 monomer prediction — high)
chain_B_ptm:     0.494    (TonB monomer prediction — medium)
```
**What to discuss:** The low ipTM is not a failure — it tells us we need wet lab data. The high chain A ptm means rbp_01 is structurally well-constrained, making it a reliable starting point for variant design.

### 3. RBP Candidate List with Confidence Scores
**File:** `03_rbp_identification/outputs/EU717894.1_rbp_candidates.csv`

| Rank | Protein | Length | Evidence | Score |
|------|---------|--------|----------|-------|
| 1 | EU717894.1_rbp_01 | **712 aa** | HMM: Tail_spike_N + unknown_C54 | 1.0 |
| 2 | EU717894.1_rbp_02 | 918 aa | HMM: collagen-like repeat | 1.0 |
| 3 | EU717894.1_rbp_03 | 224 aa | HMM: unknown_C294 | 1.0 |

**What to discuss:** rbp_01 is the primary candidate (Tail_spike_N domain = canonical tail spike anchor). rbp_02 and rbp_03 are backup candidates. Do we want to also express rbp_02 in Cycle 0 as a chimera source?

### 4. Interaction Matrix (Ground Truth)
**File:** `01_data_ground_truth/outputs/interaction_matrix.csv`
- **2,236 phage-host pairs** total (315 positive, 1,920 negative)
- **1 ground-truth row** for our reference system: phiL7 × Xcc (confirmed by Hung 2003, PMID 12646254)
- This is the training data baseline before any wet lab data comes in
**What to discuss:** The 1 confirmed pair is thin. The Boltz-2 ipTM and ESM-2 embeddings provide the "synthetic prior" to compensate before Cycle 0 data arrives.

### 5. Protein Embeddings (Module 04)
**File:** `04_protein_embedding/outputs/embeddings_esm2_t6_8M_phiL7_rbps.npz`
- 3 RBPs × 320-dim vectors (ESM-2 8M, local version)
- Production version (1280-dim, ESM-2 650M) to be run on Laguna
**What to discuss:** These vectors are the input to the deep ensemble model. The 8M version is a proof-of-concept; before Cycle 0, Laguna should run the 650M version.

### 6. Validation Strategy Recommendation
**Reference:** "Validation Strategy Analysis" section above + `docs/planning/progress_report_2026-05-11.md`
**What to discuss with PI:**
- Recommend **Tier 3**: ELISA + Plaque + ΔtonB/ΔexbB/ΔexbD1 knockouts
- Start pK18mobsacB knockouts on **May 17** alongside strain isolation
- ΔexbD2 is NOT required for phiL7 infection (Hung 2003) — a free built-in negative control
- If knockouts fail, fall back to Tier 2 (ELISA + plaque, no knockouts)

### 7. Literature Audit Corrections (PI Should Know)
**File:** `docs/reference/paper_reading_notes.md`
Five errors corrected in project documents after reading 19 papers:
1. ExbD2 NOT essential — only TonB/ExbB/ExbD1 (Hung 2003)
2. Boltz-2 affinity head is small molecule-only (protein-protein → NaN)
3. Greenman 2025 is PLoS Comput Biol, not NAR; "no single best UQ method"
4. Hie 2024 uses ESM-1b/1v, not ESM-2
5. Lee 2009 never names a tail spike — rbp_01 identified computationally

### 8. Open Questions Requiring PI Input
| Question | Urgency | Why it matters |
|----------|---------|----------------|
| pK18mobsacB vs CRISPRi for receptor knockouts? | 🔴 Before May 17 | Affects which reagents to order this week |
| AF3 model weights application — PI or Alex? | 🔴 This week | 1–7 day review; needed for variant design |
| Sewage/runoff source for phage enrichment? | 🟡 Before June 1 | Broader phage diversity → better chance of lytic phage isolation |
| Parallel manuscript submission target? | 🟡 Before Cycle 0 | Affects data documentation granularity from the start |
| Is *exbD2* data sufficient to claim "negative control"? | 🟡 PI judgment | Hung 2003 is a single study; PI may want confirmatory experiment |

---
---

# 简体中文版

## 这几天做了什么

### 5月7–8日：项目方向转型 + 隔夜平行构建

项目从「6因子手工特征加权 + ML 分类」转型为**闭环主动学习系统**。7个 AI agent 在隔夜（5/7→5/8）并行构建，每个负责一个模块。早晨起来，Module 00–06 全部有工作代码、测试套件和 README。

5月10日修复：
- 重新下载 phiL7 基因组（headers 污染问题修复）
- 替换两个无效细菌 accession（KY000037、PY746849）
- 630MB 基因组 binary gitignore（可通过 fetch 脚本重新下载）
- 设置 CARC Laguna HPC（boltz2 conda 环境，历经 8 次调试解决 torch/CUDA/trifast 版本冲突）

### 5月10日：第一次 Boltz-2 GPU 跑

在 CARC Laguna（NVIDIA L40S，job 59949，47秒）成功跑出。**但**：用了错误的蛋白——phiL7 P25（85 aa）而非 rbp_01（712 aa）。ipTM = 0.345 结果作废。

### 5月11日：文献审查 + 正确的 Boltz-2 结果

读完全部 19 篇核心论文，逐条核查 project documents 中的每个定量声明。发现并修正 5 个事实性错误：

1. **ExbD2 不是必需的**（Hung 2003）：只有 TonB、ExbB、ExbD1 是 essential，已修正所有规划文档
2. **Boltz-2 的 affinity head 只支持小分子-蛋白**（Passaro 2025）：蛋白-蛋白对输出 NaN，我们用 ipTM 作 proxy
3. **Greenman 2025 期刊是 PLoS Comput Biol**（不是 NAR Genomics）；结论是「没有单一最佳 UQ 方法」
4. **Hie 2024 用的是 ESM-1b/1v**（不是 ESM-2）；每个抗体 ~20 个 variant（不是 ~50）
5. **Lee 2009 从未指定哪个蛋白是 tail spike**——rbp_01 的 RBP 身份来自 PhageRBPdetect HMM（Tail_spike_N 结构域），不是 Lee 的注释

清除 Boltz-2 的 processed/ cache 后重新跑（job 59986）。结果：**ipTM = 0.365，chain A ptm = 0.808**。Module 04 用真实 3 个 RBP 序列重新跑（不是 5 个含 mock 的）。所有论文的阅读笔记已写入 `docs/reference/paper_reading_notes.md`。

---

## 每个 Module 的作用

| 模块 | 在 Pipeline 中的作用 | 关键输出 |
|------|---------------------|---------|
| **00 原始数据** | 为训练数据提供基因组库 | 777 phage + 34 bacteria 基因组 |
| **01 Ground Truth** | 已知 phage-host 相互作用对（系统唯一的「真实标签」） | interaction_matrix.csv：2,236 对 |
| **02 注释** | 把原始 DNA 翻译成蛋白质序列 | phiL7：80 个 ORF；Xcc：4,344 个 ORF |
| **03 RBP 识别** | 找到每个噬菌体上决定宿主范围的「钥匙蛋白」（RBP） | rbp_01（712 aa）：主要尾刺蛋白候选 |
| **04 Embedding** | 把蛋白质序列转成神经网络能处理的向量 | 1280 维 ESM-2 embedding（Laguna 上） |
| **05 结构预测** | 预测 RBP 与受体结合的 3D 结构；给模型提供「结构先验」 | PDB 文件 + ipTM 分数 |
| **06 深度集成** | 训练 5 个独立神经网络；它们的分歧 = 不确定性 | 每个 variant 的（预测分数，不确定性）tuple |
| **07 BALD** | 用不确定性决定下一个实验：选模型最不确定的 variant | 给 wet lab 的 ranked variant list |
| **08 循环数据** | 摄入 wet lab ELISA 数据，触发重训练，闭合循环 | 每轮 ELISA 数据 + 模型 checkpoint |

---

## 验证策略分析

实验验证分为三个实际可行的层次：

### 第一层：纯 ELISA

**具体操作：** 表达 His6-RBP variant，在 heat-inactivated Xanthomonas（10⁸ CFU/孔）上测结合，4PL 拟合 EC50。

**对应的 Dry Lab：** Module 06（预测结合分数）+ Module 07（BALD 选哪个 variant 测）。

**能回答的生物问题：** 「哪些 rbp_01 variant 比 wild-type 与 Xcc 细胞结合更强或更弱？」

**回答不了的：** 结合能否转化为真实感染。ELISA 检测的是结合，不是感染。T4 噬菌体 RBP 与 85% 的 *E. coli* 菌株结合，但只在 11% 上形成噬斑（Farquharson 2021）——ELISA 单独无法区分。

**时间/成本：** 每轮 1–2 周，试剂约 $400/批。**这是 Cycle 0–2 实际跑的内容。**

**故事完整性：** ⭐⭐⭐ — 「我们找到了结合更好的 variant。」定量干净，但机制层面留有空白。

---

### 第二层：ELISA + 噬斑测定（WT 菌株）

**具体操作：** 同第一层 ELISA，加上用 phiL7 在 Xanthomonas 菌坪上做 spot assay 验证裂解性感染。

**对应的 Dry Lab：** 同第一层；增加一个二元「感染/不感染」数据层。

**能回答的生物问题：** 「更高的 ELISA 结合是否与真实感染相关？」开始解决 binding ≠ infection 这个混杂因素。

**时间/成本：** 每批增加 1–2 天（噬斑测定简单，Benchling 上已有 protocol）。边际成本几乎为零。

**故事完整性：** ⭐⭐⭐⭐ — 「我们找到了结合更好的 variant，并且确认了感染。」达到 reviewer 认可的验证水准。**这是可信 iGEM 故事的最低要求。**

---

### 第三层：ELISA + 噬斑 + WT / ΔtonB / ΔexbB / ΔexbD1

**具体操作：** 用 pK18mobsacB markerless deletion 生成受体敲除株（4–6 周）。在 WT 和三种敲除株上同时做 ELISA + 噬斑测定。

**对应的 Dry Lab：** 第一层 + Layer 2 因果分解。比较 WT 与 ΔReceptor 的结合率和感染率差值。

**能回答的生物问题：**
- X − Y = 受体特异性贡献（WT 结合 - ΔtonB 结合）
- 模型的结合分数是否真的捕捉到了**受体特异性**亲和力，而非非特异性表面粘附
- 定量分解：「X% 的结合信号来自 TonB 特异性」

**ExbD2 的价值：** Hung 2003 已确认 ΔexbD2 不影响 phiL7 感染——这天然成为一个内置阴性对照，用来验证整个实验体系。

**时间/成本：** 生成敲除株 4–6 周（电穿孔 + 双交叉选择 + PCR 验证）。pK18mobsacB 质粒 ~$95。中等技术难度——Xanthomonas 电穿孔 protocol 已在 Benchling 上。

**故事完整性：** ⭐⭐⭐⭐⭐ — 「我们不只预测和测量了结合，还因果地证明了该结合信号是受体特异性的。」这是**Nature 级别的机制性 claim**，绝大多数 iGEM 团队做不到。

---

### 第四层（可选）：加入防御系统敲除

**具体操作：** 额外敲除 Xcc 的抗噬菌体防御系统（CRISPR、限制性修饰系统）。

**能回答的生物问题：** 完整分解：受体贡献 vs 防御系统贡献。可以声称「X% 来自受体识别，Y% 来自防御逃避」。

**现实评估：** Xcc 的 CRISPR 和 RM 基因座没有 tonB/exbB 那么有据可查。额外工作量大，回报不确定。**不推荐纳入 iGEM 2026 时间线**——留给后续论文。

---

### 推荐方案

| 层次 | iGEM 2026 可行性 | 故事完整性 | 建议 |
|------|----------------|-----------|------|
| 纯 ELISA | ✅ 容易 | ⭐⭐⭐ | 最低保底 |
| + 噬斑测定 | ✅ 容易（+2 天） | ⭐⭐⭐⭐ | **必做** |
| + ΔtonB/ΔexbB/ΔexbD1 | 🟡 中等（4–6 周） | ⭐⭐⭐⭐⭐ | **如果 5/17 启动就能做** |
| + 防御系统敲除 | ❌ 难，不确定 | — | 仅限后续论文 |

**结论：** 5/17 wet lab 启动时就开始做受体敲除。如果 7 月初前成功，你们有五星级的故事。失败了就退回到第二层，依然是强有力的提交。ELISA + 噬斑 + 受体敲除的组合正是 Hung 2003 用的方法——那篇文章发表在 *Biochemical and Biophysical Research Communications*，证明这个实验方案在 Xanthomonas 上完全可行。

---

*Report generated: 2026-05-11 | Branch: active-learning-pipeline*

---

## Dry Lab ↔ Wet Lab Interface Protocol

### The Two-Way Handoff

Each cycle has exactly two handoff events:

```
DRY LAB OUTPUT                      WET LAB OUTPUT
────────────────                     ──────────────
07_acquisition_function/             08_cycle_data/
outputs/cycle_N/                     outputs/cycle_N/
  recommendations.csv    ──────▶       elisa_processed.csv
  primer_sequences.txt                 plaque_results.csv
  uncertainty_bands.png                qc_report.md
  random_control.csv
  (hidden from wet lab)
```

### Dry Lab → Wet Lab Deliverables (48h SLA)

After receiving wet lab ELISA data, dry lab must within 48 hours produce:

| File | Content | Who uses it |
|------|---------|------------|
| `recommendations.csv` | Ranked list: variant ID, mutation, BALD score, predicted EC50 ± uncertainty | Wet lab cloning lead |
| `primer_sequences.txt` | NEB Q5-compatible primers for each SDM variant (auto-generated via NEBaseChanger API) | Wet lab orders from IDT |
| `uncertainty_bands.png` | Calibration plot: predicted vs measured from previous cycle | PI review |
| `safe_pick_list.csv` | Pre-selected fallback variants (expert intuition, no model needed) | Used ONLY if dry lab misses 48h SLA |

**Cycle 0 exception:** No ELISA data yet → variant design is structure-based (Boltz-2 ipTM + expert picks) + gene synthesis order (not SDM).

### Wet Lab → Dry Lab Deliverables (end of each cycle)

| File | Format | Required columns |
|------|--------|-----------------|
| `elisa_processed.csv` | CSV | variant_id, receptor_id, ec50_nM, hill_slope, r2, plate_id, date |
| `plaque_results.csv` | CSV | variant_id, strain_id, pfu_per_ml, plaque_morphology, date |
| `qc_report.md` | Markdown | SDS-PAGE image path, Bradford concentration, notes on expression issues |

**Minimum for model retraining:** At least 3 valid EC50 measurements (R² > 0.9) per new variant. If a variant failed expression or ELISA, wet lab marks `ec50_nM = NaN` with a `failed_reason` column — the model handles missing data gracefully.

### SDM as the Physical Execution of Dry Lab Recommendations

Cycles 1+ use site-directed mutagenesis instead of gene synthesis:

```
Dry lab recommends:           Wet lab executes:
"Mutate K450A in rbp_01"  →  NEB Q5 SDM (4 days)
+ primer sequences               ↓
                              Expression in BL21 (overnight, 18°C)
                                 ↓
                              Ni-NTA purification (1 day)
                                 ↓
                              ELISA + plaque assay (2 days)
                                 ↓
                              Data back to dry lab
```

**Total wet lab cycle time:** ~10–14 days from receiving recommendations to delivering ELISA data.

**Total dry lab response time:** 48 hours (retraining 5 ensemble members on Laguna + acquisition function + generating recommendations).

### Control Arm: The Hidden Random Pick

Every batch includes 1 randomly-selected variant that is:
- **NOT** the highest BALD score
- Selected by uniform random draw from the unmeasured pool
- Added to `recommendations.csv` without marking it as random (wet lab doesn't know which one is random)
- Recorded separately in `random_control.csv` (dry lab eyes only)

This design allows retrospective comparison: "if we had always picked randomly, would we have learned as fast?"

### Quality Gates and Failure Modes

| Scenario | Action |
|----------|--------|
| Dry lab misses 48h SLA | Wet lab uses `safe_pick_list.csv` (pre-approved by PI, no model needed) |
| Variant fails expression (insoluble) | Mark NaN, wet lab attempts backup truncation variant |
| ELISA R² < 0.9 (noisy data) | Mark as low-confidence, model down-weights this point |
| All 5 BALD recommendations fail QC | Wet lab runs 1 expert pick + 2 random; dry lab retrains on partial data |
| Model calibration ECE > 0.1 | Apply temperature scaling before next BALD run |
