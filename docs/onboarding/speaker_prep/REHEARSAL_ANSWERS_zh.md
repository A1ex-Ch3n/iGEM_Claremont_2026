# REHEARSAL_ANSWERS_zh.md — 答案

尝试 REHEARSAL_zh.md 后再打开。自打分。

---

## Tier 1 — 数字 / 映射

1. **777** 个 phage 基因组。
2. **34** 个细菌基因组。
3. **2,236** 对总 = **315 阳性 + 1,921 阴性 + 1 ground-truth**。
4. **712** 氨基酸 (rbp_01 = orf_00001)。
5. HMM bit score **342** 对 `Tail_spike_N` (`unknown_C54`)。
6. **3** 个 RBP candidate (rbp_01: 712 aa @ 342; rbp_02: 918 aa @ 235.1; rbp_03: 224 aa @ 56.7)。
7. ESM-2 8M → **320-D**。
8. ESM-2 650M → **1280-D**。
9. ESM-2 3B → **2560-D**。
10. ipTM = **0.365**。
11. chain_A_ptm = **0.808**。(Confidence score = 0.683。)
12. **5** 个成员 (`n_members = 5` 默认,Lakshminarayanan 2017 配方)。
13. **18/18** 测试通过 (Module 07)。(Module 06: 9/9。)
14. **n_bald=4, n_random=1** → 每 cycle 5 个 pick。
15. **48 小时** SLA (干实验 cycle)。
16. **PHANOTATE** phage ORF (McNair 2019);**pyrodigal**(Prodigal 的 Python 绑定,Hyatt 2010)细菌。绝不互换。
17. **PhageRBPdetect** (Boeckaerts 2022)。HMM 轨道 + XGBoost ML 轨道。我们用 HMM。
18. **AlphaFold 3** (Abramson 2024) + **Boltz-2** (Passaro 2025)。
19. **Hung 2003 *BBRC* 302:878–884 PMID 12646254**。Essential: **TonB, ExbB, ExbD1**。**ExbD2 NOT 必需** (ΔexbD2 保留感染)。
20. **Houlsby 2011 arXiv:1112.5745**。最初应用于**高斯过程分类 (GPC)**。我们延伸到 deep-ensemble 回归。

---

## Tier 2 — 1 分钟解释

(完整答案在 `CONCEPT_DEEPDIVE_zh.md`;快速核对如下。)

21. **AL 解决数据稀缺** ——被动 ML 在模型已自信的地方浪费测量。AL 挑模型最不确定的地方 → 每次 ELISA 最大化减少不确定。邻近领域 2–5× 加速 (Hie 2024 抗体,Yang 2025 酶)。

22. **ESM-2 embedding** ——一个 transformer 在 6500 万蛋白上学填空,内部表征(650M 模型 1280 维)捕捉结构 + 功能上下文。对湿实验说:"这是模型对每个蛋白长相的压缩理解。"

23. **ipTM = interface predicted TM-score**,0–1,对接几何**置信度**。**不是亲和力**。Boltz-2 affinity head 只小分子 → 蛋白-蛋白 NaN。ipTM 0.365 = "模型不确定怎么对接",不是"它们不结合"。

24. **Deep ensemble** ——5 个网络,同架构,不同种子 → 训练好的在数据告诉的地方一致,其他地方分歧。预测散度 = epistemic 不确定性(模型不知道什么,随数据缩小)。

25. **BALD ≈ Var_k[μ_k]** ——Houlsby 互信息目标 `H[y|x,D] - E_θ[H[y|x,θ]]`,在高斯 ensemble 同方差下化为 `½ log(1 + Var_k[μ_k]/σ²)`,关于 `Var_k[μ_k]` 单调。用 std 是单位一致。

26. **48-h SLA** ——湿实验 cycle 10–14 天,干实验 48h。任一滑期 → (a) 湿实验空转(浪费一周)或 (b) 过期推荐(浪费信息)。Safe-pick backup 是降落伞。三次滑期 = wiki freeze 风险。

27. **Lee 2009 + HMM** ——Lee 2009 BLAST 搜 OP1 ORF25 在 phiL7 同源物,找不到。我们 HMM 搜 Tail_spike_N 家族 profile,找到 rbp_01。HMM 在低序列同一性下比 BLAST 敏感。不矛盾:不同工具,不同敏感度。

28. **ΔexbD2** ——Hung 2003 显示 ExbD2 不必需;我们在做 ΔtonB/ΔexbB/ΔexbD1 同时做 ΔexbD2,它应该**保留**感染。这免费验证我们的敲除系统——若 ΔexbD2 不感染,上游某处坏了。

29. **Epistemic** = 可减少的模型不确定性(随更多数据缩小);**aleatoric** = 不可减少的测量噪声(ELISA 移液)。总方差 = epistemic + aleatoric。BALD 只针对 epistemic。

30. **随机对照臂** ——5 个 picks 里 1 个随机采样,湿实验盲法;项目结束时让我们证明 BALD 确实胜过随机选择(Hie 2024 标准)。

---

## Tier 3 — 辩护设计选择

31. **BALD 而非 Thompson sampling**:BALD 是信息论的(最大化不确定性减少);Thompson 偏 exploitation(从后验采样,行动)。早期 cycle 受益于 exploration;BALD 是规范的 exploration acquisition。ALDE (Yang 2025) 用 Thompson 因为他们在定向进化中更靠后,exploitation 更重要。不同问题,不同选择。

32. **ESM-2 而非 ProtBERT 等**:ESM-2 在 6500 万蛋白上训练 (UniRef50/90, masked LM),33 层 transformer 650M;在大多数 benchmark 上胜过 ProtBERT (Lin 2023 *Science*)。加上 PLM-interact (Liu 2025) 提供清晰的 PPI fine-tune 路径。ESM 家族有最强社区 + 工具。手工特征基本上在 2020 以来所有蛋白 ML benchmark 上输给学习的 embedding。

33. **Boltz-2 而非只用 AF3**:AF3 也预测复合物但权重获取门槛 + 跑得慢。Boltz-2 开源权重 (Github jwohlwend/boltz),给同样 ipTM + PAE 家族输出,L40S ~3 分钟跑完。我们**两者都用** ——AF3 做高质量静态结构,Boltz-2 在迭代 RBP 变体时速度 + 开源权重重要。

34. **ipTM 而非 pLDDT**:pLDDT 是**每残基局部**置信(backbone 位置)。ipTM 是**界面**置信(链怎么对接)。对 RBP × receptor 对,我们关心界面,不是每条单体的残基局部几何。单体 pLDDT(均 0.76)告诉我们每条链预测好;ipTM 告诉我们真正关心的对接。

35. **Deep ensemble 而非 MC Dropout / BNN**:校准。Ovadia 2019 显示 MC Dropout 在 shift 下校准差;完整 BNN 在我们规模算力难。Deep ensemble 简单,能扩展到 ESM-2 1280 维输入,蛋白 ML 有发表先例 (ALDE)。Greenman 2025 警告没单一方法普遍最好,但 ensemble 是可辩护的默认,加温度缩放后备。

---

## 延伸答案

36. **ECE (Expected Calibration Error)** ——把预测概率分箱,每箱看真实频率,平均差距。ECE > 0.1 = 过度自信。**温度缩放** = logits 除以学习的标量 `T`;保留排名,重校准置信度。单参数后处理修复。

37. **`safe_pick_backup.csv`** = 按 BALD 分数排序的 top-10 BALD picks,**只在**干实验错过 48h SLA 时用。预先经 PI / 干实验审核,湿实验能立即推进不等重跑。

38. **盲法** ——保留项目结束 AL-vs-random 对比的统计有效性。如果湿实验知道哪行是随机,他们可能下意识在 BALD picks 上做更好的工作(霍桑效应),或在随机上更差(跳对照)。

39. **model_version = git_sha_cycle_N** ——每个预测从这个字符串可复现:`git checkout <sha>` + 载入 `ensemble_member_*.pt` + 同 `seed=42` → 一样预测。"Cycle 3 为什么挑这个变体?" debug 关键。

40. **扩展到 T7 × E. coli K-12**:最小改动是 Module 01 输入(加 T7 × *E. coli* K-12 ground-truth 对)和 Module 04(嵌入 T7 RBP + *E. coli* 受体)。Module 02 (PHANOTATE)、03 (PhageRBPdetect)、06 (ensemble)、07 (BALD) 不变——这就是数据契约设计的整个意义。Cycle 0 重训,挑新变体。无代码更改。
