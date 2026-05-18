# REHEARSAL_zh.md — 自测题(不要偷看答案)

**用法:** 今晚睡前(冷测)、明天咖啡时(暖测)。每题标 ✅ / 🤔 / ❌。❌ → 重读相关 slide / concept。答案在 `REHEARSAL_ANSWERS_zh.md`。

---

## Tier 1 — 10 秒内冷答(20 题)

### 数字事实
1. Module 00 多少个 phage 基因组?
2. Module 00 多少个细菌基因组?
3. Module 01 interaction matrix 多少对?(总 / 阳性 / 阴性)
4. rbp_01 多少个氨基酸长?
5. rbp_01 对 Tail_spike_N 的 HMM bit score?
6. Module 03 给 phiL7 找到多少 RBP candidate?
7. ESM-2 8M embedding 维度?
8. ESM-2 650M embedding 维度?
9. ESM-2 3B embedding 维度?
10. rbp_01 × TonB Boltz-2 跑 (job 59986) 的 ipTM?
11. 同一跑的 chain_A_ptm?
12. Module 06 ensemble 几个成员?
13. Module 07 多少测试通过?
14. 每 cycle 默认 n_bald 和 n_random 各多少?
15. 干实验 cycle 周转 SLA(小时)?

### 工具-模块映射
16. 哪个工具注释 phage ORF?哪个工具注释细菌 ORF?
17. 哪个工具找 RBP candidate?两个子轨道?
18. Module 05 用的两个结构预测工具?

### 论文-claim 映射
19. 哪篇论文确立了 phiL7 受体系统?它发现什么 essential / 不 essential?
20. 哪篇论文引入 BALD?用在什么类型模型上?

---

## Tier 2 — 1 分钟解释(10 题)

对应 CONCEPT_DEEPDIVE_zh.md 的 7 个概念加 3 道延伸。

21. 为什么主动学习解决数据稀缺(对比随机 / 监督)。
22. ESM-2 embedding 实际是什么——向湿实验队员解释。
23. ipTM 衡量什么,为什么**不是**结合亲和力。
24. 为什么 deep ensemble 给 epistemic 不确定性。
25. 推导(或速写)为什么 BALD ≈ Var_k[μ_k] 在我们的设置下。
26. 为什么 48 小时 SLA 重要,拖了什么会坏。
27. 为什么 Lee 2009 + 我们的 HMM rbp_01 是补充而非矛盾。
28. ΔexbD2 在敲除组里干什么。
29. 用一口气解释 epistemic 和 aleatoric 不确定性的区别。
30. 随机对照臂的用途,一句话。

---

## Tier 3 — 必须辩护设计选择(5 题)

31. 为什么 BALD 不用 Thompson sampling(像 ALDE)?
32. 为什么 ESM-2 不用 ProtBERT(或 ProtTrans,或手工特征集)?
33. 为什么用 Boltz-2 而不只是 AlphaFold 3?
34. 为什么用 ipTM 而不是 pLDDT 作为结构置信 proxy?
35. 为什么 deep ensemble 而不是 MC Dropout 或 Bayesian 神经网络?

---

## 延伸(只有 Tier 1–3 全绿时)

36. ECE 是什么,温度缩放怎么修。
37. `safe_pick_backup.csv` 是什么,湿实验什么时候用。
38. 为什么湿实验不知道 recommendations.csv 哪行是随机对照。
39. "model_version = aa99d51_cycle_0" 告诉你什么,为什么对复现重要。
40. 有人让你把管道扩展到预测 T7 × E. coli K-12 RBP × receptor 结合,最小改动是什么?

---

## 打分(诚实)

| Tier | 答对 | 评语 |
|------|------|------|
| Tier 1 (20) | __/20 | <16 → 重浏览 slide;<12 → 取消讲演 |
| Tier 2 (10) | __/10 | <7 → 重读 CONCEPT_DEEPDIVE |
| Tier 3 (5)  | __/5  | <3 → Q&A 时预期可信度损失 |

**硬性规则:** 如果 Lee 2009 / BALD 数学 / Boltz-2 ipTM / ExbD2 四题答不对,**不要讲** ——这四个是可信度杀手。
