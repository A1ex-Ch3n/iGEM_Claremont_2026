# iGEM Claremont 2026 — PI Briefing
**Date:** May 11, 2026 | **Prepared by:** Alex Chen | **For:** Prof. J. Cesar Ignacio-Espinoza

---

# ENGLISH VERSION

## TL;DR

The computational pipeline (Modules 00–06) is **fully built and tested**. We have our first Boltz-2 3D structure of the phiL7 tail spike × Xcc TonB receptor complex. The system is ready to begin active learning as soon as wet lab delivers the first ELISA measurements (~June 1). **Module 07 (BALD acquisition function) is the only critical remaining piece and is being built this week.**

Five factual corrections were made to planning documents after a full literature audit (19 papers). All are documented in `docs/reference/paper_reading_notes.md`.

---

## What We Have Now

### Computational Results

**1. phiL7 rbp_01 × Xcc TonB — First 3D Complex Prediction (Boltz-2, job 59986)**

The 3D structure of phiL7's putative tail spike (rbp_01, 712 aa) bound to the TonB receptor (604 aa) has been predicted. Key numbers:

| Metric | Value | Interpretation |
|--------|-------|----------------|
| `interface_ipTM` | **0.365** | Low — model uncertain about HOW they dock. Expected for a novel system with no PDB template. Wet lab ELISA will resolve this. |
| `chain A ptm` | **0.808** | High — rbp_01 monomer structure is well-predicted. Reliable basis for variant design. |
| `confidence_score` | **0.683** | Overall complex quality — moderate. |

**What this means:** The model cannot yet tell us where exactly rbp_01 touches TonB. That is the core scientific unknown that the ELISA + active learning loop will map. The high chain A ptm (0.808) tells us rbp_01 is a structurally well-constrained protein — good news for expression and stability.

**2. RBP Candidates from phiL7**

PhageRBPdetect (HMM + XGBoost) identified 3 receptor-binding protein candidates from phiL7's 80 ORFs:

| Candidate | Length | Domain hit | Priority |
|-----------|--------|-----------|---------|
| rbp_01 | **712 aa** | Tail_spike_N + C-terminal binding | **Primary — Cycle 0 target** |
| rbp_02 | 918 aa | Collagen-like repeat | Backup / chimera source |
| rbp_03 | 224 aa | Short C-terminal domain | Low priority |

Note: rbp_01 is identified computationally. Lee et al. 2009 (phiL7 genome paper) suggests p20 (tail protein III, 1105 aa) as a possible host range determinant, but does not explicitly name a tail spike. Our identification is based on the Tail_spike_N HMM domain hit.

**3. Pipeline Status**

| Module | Status | Key fact |
|--------|--------|---------|
| 00 Raw Data | ✅ | 777 phage + 34 bacteria genomes |
| 01 Interaction Matrix | ✅ | 2,236 pairs; 1 confirmed (phiL7 × Xcc) |
| 02 Annotation | ✅ | phiL7: 80 ORFs; Xcc: 4,344 ORFs |
| 03 RBP Identification | ✅ | 3 candidates; rbp_01 primary |
| 04 Embeddings | ✅ | ESM-2 vectors ready (650M version pending Laguna) |
| 05 Structure | ✅ | rbp_01 × TonB PDB + ipTM scores |
| 06 Ensemble | ✅ (architecture) | Waiting for real ELISA data (~June 1) |
| 07 BALD | ❌ Building now | Must be done before May 17 |

---

## Corrected Biological Facts (Literature Audit)

After reading Hung et al. 2003 (BBRC 302:878–884, PMID 12646254) directly:

**ExbD2 is NOT required for phiL7 infection.**
Only TonB, ExbB, ExbD1 are essential. The exbD2 knockout (CH620) retained full sensitivity to phiL7. All planning documents have been corrected.

**Practical implication:** We should still generate a ΔexbD2 strain — it serves as a **built-in negative control** to validate the entire knockout experiment. If ΔexbD2 still allows infection and ΔtonB does not, that validates our receptor specificity claim.

---

## Key Decisions Needed from PI

### 1. Receptor knockout system — pK18mobsacB or CRISPRi? 🔴 Before May 17
- **pK18mobsacB (our current plan):** Markerless deletion, proven in Xanthomonas (Hung 2003 used it). Requires 4–6 weeks but gives permanent, clean knockouts. Plasmid available from Addgene (#87097, ~$95).
- **CRISPRi:** Knockdown, reversible, faster (~2 weeks). Less established in Xcc; guides need to be designed.
- **Recommendation:** pK18mobsacB for tonB/exbB/exbD1 in parallel. Start May 17.

### 2. AlphaFold 3 model weights application 🔴 This week
- AF3 provides higher-quality structures than Boltz-2, including trimer prediction.
- Requires application via Google form (academic use). 1–7 day approval.
- Should Alex or PI submit? (Academic institutional email preferred.)

### 3. Validation tier — how far do we commit? 🟡 Before Cycle 0 (June 1)

| Tier | What we do | Story | Feasibility |
|------|-----------|-------|-------------|
| ELISA only | Binding curves, 4-6 variants/cycle | ⭐⭐⭐ | ✅ Guaranteed |
| + Plaque assay | Add infection confirmation (+2 days) | ⭐⭐⭐⭐ | ✅ Easy |
| + ΔtonB/ΔexbB/ΔexbD1 | Receptor-specific causal proof | ⭐⭐⭐⭐⭐ | 🟡 If knockouts start May 17 |

**Recommendation:** Commit to Tier 3 now. Start knockouts May 17. If they succeed by July, we have a mechanistic paper-quality story.

### 4. Phage enrichment source 🟡 Before June 1
For co-isolation of lytic phages, we need an enrichment substrate with broader phage diversity than the crop tissue alone. What agricultural runoff or sewage sources does the lab have access to in LA County?

### 5. Manuscript ambition 🟡 Before Cycle 0
Does the team aim for concurrent manuscript submission (*Bioinformatics* / *NAR Genomics & Bioinformatics*) alongside the iGEM wiki? This affects data documentation granularity from the first ELISA measurement.

---

## Timeline (Next 6 Weeks)

```
May 11 (today)    Module 07 BALD implementation (Alex, this week)
May 17            Wet lab launches:
                  • Brassica sampling (LA county)
                  • Receptor knockout plasmid construction (pK18mobsacB)
                  • Cycle 0 variants ordered (gene synthesis, IDT/Twist)
May 17–June 1     Strain + phage isolation; Cycle 0 variants expressed
June 1–14         Cycle 0: ELISA optimization + first binding measurements
June 14–28        Cycle 1: Model retrained, BALD recommends next variants (SDM)
June 28–July 12   Cycle 2: Round-2 + receptor knockouts complete + causal validation
```

---
---

# 中文版（PI 简报）

**日期：** 2026-05-11 | **撰写：** Alex Chen | **呈送：** Prof. J. Cesar Ignacio-Espinoza

---

## 三句话总结

计算 pipeline（Module 00–06）**全部构建完成并通过测试**。我们已经得到了 phiL7 尾刺蛋白 × Xcc TonB 受体复合物的第一个 Boltz-2 3D 结构预测。系统已准备好在 wet lab 交出第一批 ELISA 数据（~6/1）后即刻启动主动学习循环。**Module 07（BALD 采集函数）是唯一剩余的关键模块，正在本周构建中。**

完整阅读 19 篇核心论文后，对规划文档进行了 5 处事实性修正，全部记录于 `docs/reference/paper_reading_notes.md`。

---

## 当前已有的成果

### 计算结果

**1. phiL7 rbp_01 × Xcc TonB — 首个 3D 复合体预测（Boltz-2，job 59986）**

phiL7 的尾刺蛋白候选（rbp_01，712 aa）与 TonB 受体（604 aa）的复合体结构已经预测完成。关键数字：

| 指标 | 数值 | 解读 |
|------|------|------|
| `interface_ipTM` | **0.365** | 低——模型对两个蛋白如何对接不确定。对于没有 PDB 模板的全新系统，这是预期结果。ELISA 数据将填补这个空白。 |
| `chain A ptm` | **0.808** | 高——rbp_01 单体结构预测可信。是 variant design 的可靠基础。 |
| `confidence_score` | **0.683** | 整体复合体质量——中等。 |

**含义：** 模型目前还无法告诉我们 rbp_01 在哪个位置接触 TonB。这正是 ELISA + 主动学习循环要解答的核心科学问题。高 chain A ptm（0.808）说明 rbp_01 是一个结构约束良好的蛋白——对表达和稳定性是好消息。

**2. phiL7 的 RBP 候选**

PhageRBPdetect（HMM + XGBoost）从 phiL7 的 80 个 ORF 中识别出 3 个受体结合蛋白候选：

| 候选 | 长度 | 结构域 | 优先级 |
|------|------|--------|--------|
| rbp_01 | **712 aa** | Tail_spike_N + C 端结合域 | **主要候选——Cycle 0 目标** |
| rbp_02 | 918 aa | Collagen-like repeat | 备选 / Chimera 来源 |
| rbp_03 | 224 aa | 短 C 端结构域 | 低优先级 |

注意：rbp_01 是计算识别的。Lee et al. 2009（phiL7 基因组文章）建议 p20（tail protein III，1105 aa）可能参与宿主范围决定，但未明确命名任何蛋白为 tail spike。我们的识别依据是 Tail_spike_N HMM domain hit。

**3. Pipeline 状态**

| 模块 | 状态 | 关键事实 |
|------|------|---------|
| 00 原始数据 | ✅ | 777 phage + 34 bacteria 基因组 |
| 01 相互作用矩阵 | ✅ | 2,236 对；1 个已确认（phiL7 × Xcc）|
| 02 注释 | ✅ | phiL7: 80 ORF；Xcc: 4,344 ORF |
| 03 RBP 识别 | ✅ | 3 个候选；rbp_01 为主 |
| 04 Embedding | ✅ | ESM-2 向量就绪（650M 版本待 Laguna）|
| 05 结构预测 | ✅ | rbp_01 × TonB PDB + ipTM 分数 |
| 06 深度集成 | ✅（架构） | 等待真实 ELISA 数据（~6/1）|
| 07 BALD | ❌ 构建中 | 5/17 前必须完成 |

---

## 文献核查修正的生物学事实

直接阅读 Hung et al. 2003（BBRC 302:878–884，PMID 12646254）后：

**ExbD2 不是 phiL7 感染的必需基因。**
只有 TonB、ExbB、ExbD1 是 essential。exbD2 敲除株（CH620）保留了对 phiL7 的完全敏感性。所有规划文档已修正。

**实践含义：** 我们仍应生成 ΔexbD2 菌株——它可以作为整个敲除实验的**内置阴性对照**。如果 ΔexbD2 仍然允许感染而 ΔtonB 不允许，这将验证我们的受体特异性 claim。

---

## 需要 PI 做决定的事项

### 1. 受体敲除系统——pK18mobsacB 还是 CRISPRi？🔴 5/17 前
- **pK18mobsacB（当前计划）：** Markerless deletion，已在 Xanthomonas 上验证（Hung 2003 用的就是这个）。需要 4–6 周，但产生永久、干净的敲除株。质粒可从 Addgene 购买（#87097，~$95）。
- **CRISPRi：** 基因沉默，可逆，更快（~2 周）。在 Xcc 上不如 pK18mobsacB 成熟；需要设计引导 RNA。
- **建议：** 对 tonB/exbB/exbD1 并行使用 pK18mobsacB。5/17 启动。

### 2. AlphaFold 3 模型权重申请 🔴 本周
- AF3 提供比 Boltz-2 更高质量的结构，包括 trimer 预测。
- 需要通过 Google 表单申请（学术用途）。1–7 天审批。
- 由 Alex 还是 PI 提交？（学术机构 email 优先。）

### 3. 验证层次——我们承诺走到哪一步？🟡 Cycle 0 前（6/1）

| 层次 | 具体操作 | 故事完整性 | 可行性 |
|------|---------|-----------|--------|
| 纯 ELISA | 结合曲线，4-6 个 variant/cycle | ⭐⭐⭐ | ✅ 有保障 |
| + 噬斑测定 | 加感染确认（+2 天） | ⭐⭐⭐⭐ | ✅ 容易 |
| + ΔtonB/ΔexbB/ΔexbD1 | 受体特异性因果证明 | ⭐⭐⭐⭐⭐ | 🟡 若 5/17 启动敲除 |

**建议：** 现在就承诺第三层。5/17 启动敲除。如果 7 月前成功，我们有论文级别的机制性故事。

### 4. 噬菌体富集来源 🟡 6/1 前
为了共分离裂解性噬菌体，需要一个多样性更丰富的富集底物。实验室在 LA 县有哪些农业污水或灌溉水来源可以使用？

### 5. 论文投稿意向 🟡 Cycle 0 前
团队是否计划与 iGEM wiki 并行投稿（目标期刊：*Bioinformatics* / *NAR Genomics & Bioinformatics*）？这会影响从第一个 ELISA 测量开始的数据文档化粒度。

---

## 时间线（未来 6 周）

```
5/11（今天）      Module 07 BALD 实现（Alex，本周）
5/17              Wet lab 启动：
                  • 卷心菜采样（LA 县）
                  • 受体敲除质粒构建（pK18mobsacB）
                  • Cycle 0 variants 下单（基因合成，IDT/Twist）
5/17–6/1         菌株 + 噬菌体分离；Cycle 0 variants 表达
6/1–6/14         Cycle 0：ELISA 优化 + 第一批结合测量
6/14–6/28        Cycle 1：模型重训，BALD 推荐下一批 variants（SDM）
6/28–7/12        Cycle 2：第二轮 + 受体敲除完成 + 因果验证
```

---

*Full pipeline details: `docs/planning/iGEM_2026_Project_Plan.md` (EN) | `docs/planning/iGEM_2026_项目大纲_中文版.md` (ZH)*
*Boltz-2 structure: `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/`*
*Literature audit notes: `docs/reference/paper_reading_notes.md`*
