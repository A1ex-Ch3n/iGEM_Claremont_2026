# CHEAT_SHEET_zh.md — 一页安全网

**打印它。放讲台上。不要照念——但知道它在那里会让你冷静。**

---

## 关键数字 — 冷记

| 项目 | 值 |
|------|-----|
| Phage 基因组 (Module 00) | **777** |
| 细菌基因组 (Module 00) | **34** |
| Interaction matrix 对数 | **2,236** (315 + / 1,921 −) |
| rbp_01 长度 | **712 aa** |
| rbp_01 HMM bit score (Tail_spike_N) | **342** |
| phiL7 RBP candidate 数 | **3** |
| ESM-2 8M / 650M / 3B 维度 | **320 / 1280 / 2560** |
| Boltz-2 ipTM (rbp_01 × TonB, job 59986) | **0.365** |
| Boltz-2 chain_A_ptm | **0.808** |
| Boltz-2 confidence_score | **0.683** |
| Ensemble 成员 | **5** |
| 测试通过:M01/M02/M03/M04/M05/M06/M07 | 22/26/25+/17/28/9/**18** |
| 每 cycle picks | **4 BALD + 1 random** |
| 干实验 SLA | **48 小时** |
| 湿实验 cycle | **10–14 天** |
| SDM 成本 / 时间 | **$50 / 4 天** |
| 基因合成成本 / 时间 | **$150 / 2 周** |

---

## 引用 — 简短形

| 标签 | Claim |
|------|-------|
| **Lee 2009** *AEM* 75:7828 | phiL7 基因组;主动搜 OP1 ORF25 同源物,找不到 |
| **Hung 2003** *BBRC* 302:878 PMID 12646254 | TonB/ExbB/ExbD1 essential;**ExbD2 NOT** 必需 |
| **Boeckaerts 2022** *Viruses* 14:1329 | PhageRBPdetect HMM + XGBoost;PR AUC 93.8% |
| **Boeckaerts 2024** *Nat Commun* 15:4768 | PhageHostLearn AUC 至 0.82(属内) |
| **Mutalik 2025** *bioRxiv* | PAML benchmark;跨属 AUC 0.67–0.70 |
| **Lin 2023** *Science* 379:1123 | ESM-2;6500 万蛋白,masked LM |
| **Liu 2025** *Nat Commun* 16:64512 | PLM-interact;+16–28% AUPR 跨物种 PPI |
| **Abramson 2024** *Nature* 630:493 | AlphaFold 3 |
| **Passaro 2025** *bioRxiv* | Boltz-2;affinity head **只小分子** |
| **Lakshminarayanan 2017** *NeurIPS* arXiv:1612.01474 | Deep ensemble + Gaussian NLL |
| **Ovadia 2019** *NeurIPS* | UQ under shift;ensemble > MC Dropout 校准 |
| **Greenman 2025** *PLoS Comp Biol* 21:e1012639 | "无单一最佳 UQ 方法";**不是** *NAR Genomics* |
| **Houlsby 2011** arXiv:1112.5745 | BALD;最初用于 GPC 分类 |
| **Hie 2024** *Nat Biotechnol* 42:275 | 抗体亲和成熟;用 **ESM-1b/1v**(不是 ESM-2);~20 个 variant |
| **Yang 2025** *Nat Commun* 16:55987 | ALDE;**Thompson sampling** + one-hot(**不是** BALD/ESM-2) |
| **McNair 2019** *Bioinformatics* 35:4537 | PHANOTATE (phage ORF, DP) |
| **Hyatt 2010** *BMC Bioinf* 11:119 | Prodigal(细菌 ORF)——我们用 pyrodigal 绑定 |
| **Schäfer 1994** *Gene* 145:69 | pK18mobsacB(suicide vector + sacB 反向选择) |
| **Farquharson 2021** | T4 × E. coli:结合 ≠ 感染(85% 结合,11% 噬菌斑) |

---

## ⚠️ 三句必须原话讲

> **1.** "Boltz-2 的 affinity head 只针对小分子——蛋白-蛋白对输出 NaN。我们用 ipTM 作为结构置信度 proxy,不是结合亲和力。"

> **2.** "我们 deep-ensemble 回归的 BALD = `epistemic_std = Std_k[μ_k]` ——ensemble 成员均值的标准差。Houlsby 2011 在 GPC 分类上推导,deep-ensemble 回归形式是我们的扩展。"

> **3.** "Lee 2009 用 BLAST 明确写明找不到 OP1 ORF25 同源物。我们用 Hidden Markov Model ——更敏感的工具——找到了 rbp_01。两个结果是补充,不是矛盾。"

---

## ⚠️ 三句绝不能说

> ❌ "Boltz-2 给我们 zero-shot 结合亲和力。"
> ❌ "ALDE (Yang 2025) 验证了 BALD。"(他们用 Thompson sampling。)
> ❌ "Lee 2009 漏掉了" 或 "Hie 2024 用 ESM-2"(ESM-1b/1v;Lee 主动搜索得到负面结果)。

---

## 卡住了 — 三句救场

> 🎯 **"让我把代码调出来…"** ——打开 `bald.py` 或 `ensemble.py`。代码本身回答大多数 ML 问题。从自己 repo 念看起来像有能力,不像慌乱。

> 🎯 **"这是开放问题——我记下来。"** ——拿出笔,可见地写在纸上。不要瞎吹。诚实 + 记录 比错答好。

> 🎯 **"短答案是 X。长答案在 `docs/onboarding/guide_zh.md` 第 Y 节——会后我可以详走。"** ——保持节奏,延后但不躲避。

---

## 听众资料 — 分享路径

- Slides (EN): `docs/onboarding/slides_en.pdf`
- Slides (ZH): `docs/onboarding/slides_zh.pdf`
- 配套指南 (EN): `docs/onboarding/guide_en.md`
- 配套指南 (ZH): `docs/onboarding/guide_zh.md`
- 演示 runbook: `docs/onboarding/DEMO.md`
- 论文审计: `docs/reference/paper_reading_notes.md`
- 项目计划 (EN): `docs/planning/iGEM_2026_Project_Plan.md`
- PI briefing 快照: `docs/planning/PI_briefing_2026-05-11.md`

---

## 演示应急 — 最小可行演示

如果全坏了:
1. `cat 07_acquisition_function/outputs/cycle_1/recommendations.csv` ——展示实际 recommendation 文件。
2. `cat 06_uncertainty_model/outputs/cycle_0/predictions.csv | head` ——展示带 epistemic_std 的预测。
3. 打开 `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/affinity.json` ——展示 ipTM 0.365。

这三个文件证明管道可工作。不必现场跑。

---

## 双语切换 — 短提示

- "下面用中文复述一下…" / "Let me say that in Chinese…"
- "回到英文 / Back to English."
- "用 Sarah 听得最舒服的语言…" / "In the language Sarah hears most easily…"

句中切换被允许——听众是双语的。
