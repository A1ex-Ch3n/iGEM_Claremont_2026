# iGEM Claremont 2026 — PI Briefing & Project Status
**Date:** 2026-05-12 (created 2026-05-11, updated 2026-05-12) | **By:** Alex Chen | **For:** Prof. J. Cesar Ignacio-Espinoza

---

# ENGLISH VERSION

## TL;DR

The **entire computational pipeline (Modules 00–07) is fully built and tested.** We have the first Boltz-2 3D structure of phiL7 rbp_01 × Xcc TonB (ipTM = 0.365), a calibrated 5-member deep ensemble, and a working BALD acquisition function that produces ranked variant recommendations in under 1 second. The system will begin active learning as soon as wet lab delivers the first ELISA measurements (~June 1). **No critical dry-lab work remains before May 17 wet lab launch.**

Five factual corrections were made after a full literature audit (19 papers). All documented in `docs/reference/paper_reading_notes.md`.

---

## Pipeline Status

| Module | Status | Key fact |
|--------|--------|---------|
| 00 Raw Data | ✅ | 777 phage + 34 bacteria genomes |
| 01 Interaction Matrix | ✅ | 2,236 pairs; 1 confirmed (phiL7 × Xcc) |
| 02 Annotation | ✅ | phiL7: 80 ORFs; Xcc: 4,344 ORFs |
| 03 RBP Identification | ✅ | 3 candidates; rbp_01 (712 aa) primary |
| 04 Embeddings | ✅ | ESM-2 vectors ready (650M version pending Laguna) |
| 05 Structure | ✅ | rbp_01 × TonB PDB + ipTM = 0.365 |
| 06 Deep Ensemble | ✅ | 5-member MLP, calibrated, outputs epistemic_std |
| 07 BALD | ✅ | Scorer + CLI orchestrator; 18 tests pass; first cycle run |
| 08 Cycle Data | ⏳ Cycle 0 starts ~6/1 | Waiting for ELISA measurements |

---

## Computational Results

### 1. phiL7 rbp_01 × Xcc TonB — First 3D Complex Prediction (Boltz-2, job 59986)

| Metric | Value | Interpretation |
|--------|-------|----------------|
| `interface_ipTM` | **0.365** | Low — model uncertain about HOW they dock. Expected for a novel system with no PDB template. ELISA will resolve this. |
| `chain A ptm` | **0.808** | High — rbp_01 monomer structure is well-predicted. Reliable basis for variant design. |
| `confidence_score` | **0.683** | Overall complex quality — moderate. |

The low ipTM is not a failure — it defines the experiment. That uncertainty IS the question the ELISA + active learning loop is designed to answer. The high chain A ptm (0.808) means rbp_01 is structurally well-constrained — good news for expression and stability.

**File:** `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/` (view in PyMOL or ChimeraX)

### 2. RBP Candidates from phiL7

| Candidate | Length | Domain hit | Priority |
|-----------|--------|-----------|---------|
| rbp_01 | **712 aa** | Tail_spike_N + C-terminal binding | **Primary — Cycle 0 target** |
| rbp_02 | 918 aa | Collagen-like repeat | Backup / chimera source |
| rbp_03 | 224 aa | Short C-terminal domain | Low priority |

Note: rbp_01 is computationally identified. Lee et al. 2009 (phiL7 genome) suggests p20 (1105 aa) for host range but does not name a tail spike. Our identification is based on the Tail_spike_N HMM hit from PhageRBPdetect.

### 3. BALD Acquisition — First Cycle Recommendations (synthetic prior)

Before any ELISA data arrives, BALD ranks by epistemic uncertainty (Var_k[μ_k], variance of ensemble member means). On the synthetic Cycle 0 predictions:

| Priority | Variant | Receptor | BALD score | Rationale |
|----------|---------|----------|------------|-----------|
| 1 (BALD top) | rbp_07 | rec_02 | 0.218 | Highest ensemble disagreement |
| 2 | rbp_03 | rec_01 | 0.197 | — |
| 3 | rbp_05 | rec_02 | 0.197 | — |
| 4 | rbp_01 | rec_02 | 0.190 | — |
| 5 (random) | rbp_03 | rec_03 | 0.127 | Control arm for retrospective comparison |

These are synthetic placeholders — actual Cycle 1 recommendations will be based on real rbp_01 variants × TonB after Cycle 0 ELISA data is in.

---

## What Each Module Does

| Module | Role | Key output |
|--------|------|-----------|
| 00 Raw Data | Genome library for training data | 777 phage + 34 bacteria genomes |
| 01 Ground Truth | Labeled phage-host pairs (the only real labels in the system) | interaction_matrix.csv: 2,236 pairs |
| 02 Annotation | Translates raw DNA → protein sequences | phiL7: 80 ORFs; Xcc: 4,344 ORFs |
| 03 RBP ID | Finds the "key" protein determining host range | rbp_01 (712 aa): primary tail spike candidate |
| 04 Embedding | Converts protein sequences → neural network input vectors | 1280-dim ESM-2 embeddings |
| 05 Structure | 3D structure of RBP × receptor complex; structural confidence prior | PDB + ipTM score |
| 06 Ensemble | 5 independent neural networks; their disagreement = what the model doesn't know | (predicted score, epistemic_std) per variant |
| 07 BALD | Uses uncertainty to pick the next experiment: select the variant the model is most uncertain about | Ranked variant CSV for wet lab |
| 08 Cycle Data | Ingests ELISA results, triggers retraining, closes the loop | Per-cycle data + model checkpoints |

---

## Validation Strategy

Three practical tiers, all using the same workflow (ELISA-first, then add layers):

| Tier | What you do | Scientific story | Feasibility |
|------|------------|-----------------|-------------|
| 1 — ELISA only | Binding curves on WT Xcc, 4–6 variants/cycle | "We found variants that bind better." | ✅ Guaranteed — this is Cycle 0–2 |
| 2 — + Plaque assay | +2 days: spot assay on WT confirms lytic infection | "Binding → infection confirmed." | ✅ Near-zero incremental cost |
| 3 — + ΔtonB/ΔexbB/ΔexbD1 | Markerless deletions (pK18mobsacB, 4–6 weeks) | "Binding is causally receptor-specific." This is a **mechanistic paper-quality claim**. | 🟡 If knockouts start May 17 |

**Built-in negative control:** Hung 2003 confirms ΔexbD2 does NOT affect phiL7 infection — generating this strain validates the entire knockout system for free. If ΔexbD2 still allows infection and ΔtonB does not, the receptor specificity claim is experimentally supported.

**Recommendation:** Commit to Tier 3 now. Start pK18mobsacB knockouts May 17. If knockouts succeed by early July, the result is a five-star iGEM story and paper-grade mechanistic claim.

---

## Corrected Biological Facts (Literature Audit)

After reading Hung et al. 2003 (BBRC 302:878–884, PMID 12646254) directly:

**ExbD2 is NOT required for phiL7 infection.** Only TonB, ExbB, ExbD1 are essential. The exbD2 knockout (CH620) retained full sensitivity. All planning documents have been corrected.

Four additional corrections from the 19-paper audit:
1. Boltz-2 affinity head = small molecule–protein only. Protein-protein pairs output NaN. We use ipTM as structural confidence proxy — NOT a binding affinity.
2. Greenman 2025 journal = *PLoS Comput Biol* 21(1):e1012639 (not NAR Genomics); conclusion: "no single best UQ method."
3. Hie 2024 uses ESM-1b/1v (not ESM-2); ~20 variants per antibody (not ~50).
4. Lee 2009 never names a tail spike — rbp_01 identified computationally by PhageRBPdetect (Tail_spike_N HMM domain).

---

## Key Decisions Needed from PI

### 1. Receptor knockout system — pK18mobsacB or CRISPRi? 🔴 Before May 17
- **pK18mobsacB (our plan):** Markerless deletion, proven in *Xanthomonas* (Hung 2003). 4–6 weeks, permanent clean knockouts. Plasmid on Addgene (#87097, ~$95).
- **CRISPRi:** Knockdown, reversible, faster (~2 weeks). Less established in Xcc; guides need design.
- **Recommendation:** pK18mobsacB for tonB/exbB/exbD1 in parallel. Start May 17.

### 2. AlphaFold 3 model weights application 🔴 This week
- AF3 gives higher-quality structures than Boltz-2, including trimer predictions.
- Requires Google form (academic use). 1–7 day approval.
- Should Alex or PI submit? (Institutional email preferred.)

### 3. Validation tier commitment 🟡 Before Cycle 0 (June 1)
See Validation Strategy section above. Recommend Tier 3 if knockouts start May 17.

### 4. Phage enrichment source 🟡 Before June 1
Co-isolation requires enrichment substrate broader than crop tissue alone. What agricultural runoff or sewage sources does the lab have access to in LA County?

### 5. Manuscript ambition 🟡 Before Cycle 0
Aim for concurrent submission (*Bioinformatics* / *NAR Genomics & Bioinformatics*) alongside iGEM wiki? This sets data documentation granularity from the first ELISA measurement.

---

## Discussion-Ready Outputs (for PI meeting)

| # | Item | File | Key talking point |
|---|------|------|-------------------|
| 1 | rbp_01 × TonB 3D structure | `05_structure_prediction/outputs/boltz2/.../predictions/*.pdb` | View in PyMOL/ChimeraX. Does predicted interface match known TonB-ExbB face from *E. coli* literature? |
| 2 | Structural confidence numbers | `...affinity.json`: ipTM=0.365, chain_A_ptm=0.808 | Low ipTM = experiment needed, not model failure. High chain A = reliable for variant design. |
| 3 | RBP candidate list | `03_rbp_identification/outputs/EU717894.1_rbp_candidates.csv` | rbp_01 primary. Should we express rbp_02 in Cycle 0 as chimera source? |
| 4 | Interaction matrix | `01_data_ground_truth/outputs/interaction_matrix.csv` | 2,236 pairs; 1 confirmed. Thin baseline compensated by Boltz-2 structural prior. |
| 5 | ESM-2 embeddings | `04_protein_embedding/outputs/*.npz` | 8M proof-of-concept; 650M needed on Laguna before Cycle 0. |
| 6 | BALD recommendations | `07_acquisition_function/outputs/cycle_1/recommendations.csv` | 4 BALD picks + 1 random control. Synthetic now; real picks after Cycle 0 ELISA. |
| 7 | Validation strategy comparison | This document, "Validation Strategy" section | Recommend Tier 3. Hung 2003 used the same approach in Xcc — proven feasible. |
| 8 | Literature audit corrections | `docs/reference/paper_reading_notes.md` | 5 errors corrected. Most important: ExbD2 is NOT required → free negative control. |

---

## Dry-Lab → Wet-Lab Handoff Protocol

Each cycle has two handoffs:

**Dry lab → Wet lab (48-hour SLA after ELISA data arrives):**
| File | Content |
|------|---------|
| `recommendations.csv` | Ranked variants: ID, mutation, BALD score, predicted EC50 ± uncertainty |
| `primer_sequences.txt` | NEB Q5-compatible primers for SDM variants |
| `uncertainty_bands.png` | Calibration plot: predicted vs measured from previous cycle |
| `safe_pick_backup.csv` | Pre-selected fallback list for PI/team if 48-h SLA is missed |

**Wet lab → Dry lab (end of each cycle):**
| File | Required columns |
|------|-----------------|
| `elisa_processed.csv` | variant_id, receptor_id, ec50_nM, hill_slope, r2, plate_id, date |
| `plaque_results.csv` | variant_id, strain_id, pfu_per_ml, date |
| `qc_report.md` | SDS-PAGE image path, concentration, expression notes |

Minimum for retraining: ≥3 valid EC50 measurements (R² > 0.9) per new variant. Failed variants marked `ec50_nM = NaN` with `failed_reason` — the model handles missing data.

**Cycle 0 exception:** No ELISA yet — variant design is structure-based (Boltz-2 interface + expert picks) + gene synthesis order (not SDM).

---

## Timeline

```
2026-05-11   Module 07 BALD implemented ✅ (completed same day as briefing)
2026-05-12   Full pipeline (00–07) complete and committed

2026-05-17   Wet lab launches:
             • Brassica sampling (LA County)
             • Receptor knockout plasmid construction (pK18mobsacB)
             • Cycle 0 variants ordered (gene synthesis, IDT/Twist)

2026-05-17 → 2026-06-01
             Strain + phage isolation; Cycle 0 variants expressed in BL21

2026-06-01 → 2026-06-14
             Cycle 0: ELISA optimization + first binding measurements

2026-06-14 → 2026-06-28
             Cycle 1: Ensemble retrained; BALD recommends next variants (SDM)

2026-06-28 → 2026-07-12
             Cycle 2: Round 2 + receptor knockouts complete + causal validation
```

---

*Full pipeline: `docs/planning/iGEM_2026_Project_Plan.md` (EN) | `docs/planning/iGEM_2026_项目大纲_中文版.md` (ZH)*
*Boltz-2 structure: `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/`*
*Literature audit: `docs/reference/paper_reading_notes.md`*

---
---

# 中文版（PI 简报）

**日期：** 2026-05-12（2026-05-11 撰写，2026-05-12 更新）| **撰写：** Alex Chen | **呈送：** Prof. J. Cesar Ignacio-Espinoza

---

## 三句话总结

**整个计算 pipeline（Module 00–07）已全部构建完成并通过测试。** 我们已得到 phiL7 rbp_01 × Xcc TonB 复合体的首个 Boltz-2 3D 结构（ipTM = 0.365），一个经过校准的 5 成员深度集成模型，以及能在 1 秒内输出 variant 排名的 BALD 采集函数。系统已准备好在 wet lab 交出第一批 ELISA 数据（~6/1）后即刻启动主动学习循环。**5/17 wet lab 启动前，干实验室已无关键待办事项。**

完整阅读 19 篇核心论文后，对规划文档进行了 5 处事实性修正，全部记录于 `docs/reference/paper_reading_notes.md`。

---

## Pipeline 状态

| 模块 | 状态 | 关键事实 |
|------|------|---------|
| 00 原始数据 | ✅ | 777 phage + 34 bacteria 基因组 |
| 01 相互作用矩阵 | ✅ | 2,236 对；1 个已确认（phiL7 × Xcc）|
| 02 注释 | ✅ | phiL7: 80 ORF；Xcc: 4,344 ORF |
| 03 RBP 识别 | ✅ | 3 个候选；rbp_01（712 aa）为主 |
| 04 Embedding | ✅ | ESM-2 向量就绪（650M 版本待 Laguna）|
| 05 结构预测 | ✅ | rbp_01 × TonB PDB + ipTM = 0.365 |
| 06 深度集成 | ✅ | 5 成员 MLP，已校准，输出 epistemic_std |
| 07 BALD | ✅ | 采集函数 + CLI 流程脚本；18 个测试通过；首轮跑通 |
| 08 循环数据 | ⏳ Cycle 0 ~6/1 | 等待 ELISA 测量结果 |

---

## 计算结果

### 1. phiL7 rbp_01 × Xcc TonB — 首个 3D 复合体预测（Boltz-2，job 59986）

| 指标 | 数值 | 解读 |
|------|------|------|
| `interface_ipTM` | **0.365** | 低——模型对两个蛋白如何对接不确定。对于无 PDB 模板的全新系统，这是预期结果。ELISA 数据将填补这个空白。 |
| `chain A ptm` | **0.808** | 高——rbp_01 单体结构预测可信。是 variant design 的可靠基础。 |
| `confidence_score` | **0.683** | 整体复合体质量——中等。 |

低 ipTM 不是失败——它定义了实验目标。这种不确定性正是 ELISA + 主动学习循环要解答的问题。高 chain A ptm（0.808）说明 rbp_01 是结构约束良好的蛋白——对表达和稳定性是好消息。

**文件：** `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/`（用 PyMOL 或 ChimeraX 打开）

### 2. phiL7 的 RBP 候选

| 候选 | 长度 | 结构域 | 优先级 |
|------|------|--------|--------|
| rbp_01 | **712 aa** | Tail_spike_N + C 端结合域 | **主要候选——Cycle 0 目标** |
| rbp_02 | 918 aa | Collagen-like repeat | 备选 / Chimera 来源 |
| rbp_03 | 224 aa | 短 C 端结构域 | 低优先级 |

注意：rbp_01 是计算识别的。Lee et al. 2009 建议 p20（1105 aa）可能参与宿主范围决定，但未明确命名任何蛋白为 tail spike。我们的识别依据是 Tail_spike_N HMM。

### 3. BALD 采集函数——第一轮推荐（合成先验）

当前基于合成 Cycle 0 预测（真实 ELISA 数据到位后将替换）：

| 优先级 | Variant | Receptor | BALD 分数 | 说明 |
|-------|---------|----------|----------|-----|
| 1（BALD 最优） | rbp_07 | rec_02 | 0.218 | 集成成员分歧最大 |
| 2 | rbp_03 | rec_01 | 0.197 | — |
| 3 | rbp_05 | rec_02 | 0.197 | — |
| 4 | rbp_01 | rec_02 | 0.190 | — |
| 5（随机对照） | rbp_03 | rec_03 | 0.127 | 用于事后学习曲线对比 |

---

## 每个 Module 的作用

| 模块 | 作用 | 关键输出 |
|------|------|---------|
| 00 原始数据 | 基因组库——训练数据来源 | 777 phage + 34 bacteria 基因组 |
| 01 Ground Truth | 已知 phage-host 相互作用标签（系统中唯一的真实标签） | interaction_matrix.csv：2,236 对 |
| 02 注释 | 原始 DNA → 蛋白质序列 | phiL7：80 ORF；Xcc：4,344 ORF |
| 03 RBP 识别 | 找到决定宿主范围的「钥匙蛋白」 | rbp_01（712 aa）：主要尾刺蛋白候选 |
| 04 Embedding | 蛋白质序列 → 神经网络输入向量 | 1280 维 ESM-2 embedding |
| 05 结构预测 | RBP × 受体复合体 3D 结构；提供结构先验 | PDB + ipTM 分数 |
| 06 深度集成 | 5 个独立神经网络；分歧 = 不确定性 = 模型不知道的地方 | 每个 variant 的（预测分数，epistemic_std）|
| 07 BALD | 用不确定性选下一个实验：选模型最不确定的 variant | 给 wet lab 的 ranked variant CSV |
| 08 循环数据 | 摄入 ELISA 数据，触发重训练，闭合循环 | 每轮 ELISA 数据 + 模型 checkpoint |

---

## 验证策略

三个实际可行的层次，用同一个工作流（ELISA 为基础，逐层叠加）：

| 层次 | 具体操作 | 科学故事 | 可行性 |
|------|---------|---------|--------|
| 第一层：纯 ELISA | WT Xcc 上测结合曲线，每轮 4–6 个 variant | 「我们找到了结合更好的 variant。」 | ✅ 有保障——这是 Cycle 0–2 实际内容 |
| 第二层：+ 噬斑测定 | +2 天：WT 上做 spot assay 确认裂解性感染 | 「结合 → 感染已确认。」 | ✅ 几乎零边际成本 |
| 第三层：+ ΔtonB/ΔexbB/ΔexbD1 | Markerless deletion（pK18mobsacB，4–6 周） | 「结合信号是受体特异性的。」这是**论文级别的机制性 claim**。 | 🟡 若 5/17 启动敲除 |

**内置阴性对照：** Hung 2003 已确认 ΔexbD2 不影响 phiL7 感染——生成此菌株即可免费验证整个敲除实验体系。若 ΔexbD2 仍允许感染而 ΔtonB 不允许，则受体特异性 claim 得到实验支持。

**建议：** 现在就承诺第三层。5/17 启动 pK18mobsacB 敲除。若 7 月初前成功，你们有五星级 iGEM 故事和论文级机制性 claim。

---

## 文献核查修正的事实

直接阅读 Hung et al. 2003（BBRC 302:878–884，PMID 12646254）后：

**ExbD2 不是 phiL7 感染的必需基因。** 只有 TonB、ExbB、ExbD1 是 essential。exbD2 敲除株（CH620）保留对 phiL7 的完全敏感性。所有规划文档已修正。

另外 4 处修正：
1. Boltz-2 affinity head 只支持小分子-蛋白。蛋白-蛋白对输出 NaN。用 ipTM 作结构信心 proxy——**不是**结合亲和力。
2. Greenman 2025 期刊是 *PLoS Comput Biol* 21(1):e1012639（不是 NAR Genomics）；结论：「没有单一最佳 UQ 方法」。
3. Hie 2024 用的是 ESM-1b/1v（不是 ESM-2）；每个抗体 ~20 个 variant（不是 ~50）。
4. Lee 2009 从未指定哪个蛋白是 tail spike——rbp_01 由 PhageRBPdetect HMM（Tail_spike_N 结构域）计算识别。

---

## 需要 PI 做决定的事项

### 1. 受体敲除系统——pK18mobsacB 还是 CRISPRi？🔴 5/17 前
- **pK18mobsacB（当前计划）：** Markerless deletion，已在 *Xanthomonas* 上验证（Hung 2003）。4–6 周，永久干净的敲除株。质粒在 Addgene（#87097，~$95）。
- **CRISPRi：** 基因沉默，可逆，更快（~2 周）。在 Xcc 上不如 pK18mobsacB 成熟。
- **建议：** 对 tonB/exbB/exbD1 并行使用 pK18mobsacB。5/17 启动。

### 2. AlphaFold 3 模型权重申请 🔴 本周
- AF3 提供比 Boltz-2 更高质量的结构，含 trimer 预测。
- 需通过 Google 表单申请（学术用途）。1–7 天审批。
- 由 Alex 还是 PI 提交？（学术机构 email 优先。）

### 3. 验证层次承诺 🟡 6/1 前
见上方「验证策略」。若 5/17 启动敲除，建议承诺第三层。

### 4. 噬菌体富集来源 🟡 6/1 前
共分离需要多样性更丰富的富集底物。实验室在 LA 县有哪些农业污水或灌溉水来源？

### 5. 论文投稿意向 🟡 Cycle 0 前
是否计划与 iGEM wiki 并行投稿（*Bioinformatics* / *NAR Genomics & Bioinformatics*）？这会影响从第一个 ELISA 测量开始的数据文档化粒度。

---

## 时间线

```
2026-05-11   Module 07 BALD 实现完成 ✅
2026-05-12   全 pipeline（00–07）完成并提交

2026-05-17   Wet lab 启动：
             • 十字花科蔬菜采样（LA 县）
             • 受体敲除质粒构建（pK18mobsacB）
             • Cycle 0 variants 下单（基因合成，IDT/Twist）

2026-05-17 → 2026-06-01
             菌株 + 噬菌体分离；Cycle 0 variants 在 BL21 中表达

2026-06-01 → 2026-06-14
             Cycle 0：ELISA 优化 + 第一批结合测量

2026-06-14 → 2026-06-28
             Cycle 1：集成重训，BALD 推荐下一批 variants（SDM）

2026-06-28 → 2026-07-12
             Cycle 2：第二轮 + 受体敲除完成 + 因果验证
```

---

*完整 pipeline：`docs/planning/iGEM_2026_Project_Plan.md`（EN）| `docs/planning/iGEM_2026_项目大纲_中文版.md`（ZH）*
*Boltz-2 结构：`05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/`*
*文献审查笔记：`docs/reference/paper_reading_notes.md`*
