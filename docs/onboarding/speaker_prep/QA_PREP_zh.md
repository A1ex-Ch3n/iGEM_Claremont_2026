# QA_PREP_zh.md — 预期问题与示范回答 + 兜底

三个听众分组。每个问题:**自信 3–5 句回答** + Alex 不确定时的**兜底**。

---

## [湿实验 Q&A] — Sarah、Olivia、Weitao、Carol

### W1. "`recommendations.csv` 到了我具体做什么?"
**自信回答:** 从 `07_acquisition_function/outputs/cycle_<N>/` 拉 CSV。里面 5 行——4 个 BALD picks + 1 个随机对照,**打乱过的,你不知道哪个是哪个**。每行 `rbp_id`、`receptor_id`、预测分。配套 `primer_sequences.txt` (NEB Q5 兼容),按 Benchling protocol 跑 SDM。ELISA 完成后交回 `elisa_processed.csv`。提前做完了就做 `safe_pick_backup.csv` 里的额外 picks——已经审核过的备选。
**兜底:** "好问题——我把 protocol 调出来 / 我们查 `08_cycle_data/README.md`。"

### W2. "变体在 BL21 表达不出来怎么办?"
**自信回答:** 标 `ec50_nM = NaN` 加 `failed_reason = "insoluble"`(或 "low_expression")。模型能处理缺失数据,不会卡 cycle。备份方案:超过 ~30% 不溶时改用 GCN4 三聚体融合标签——tail spike 天然要三聚化。PI briefing 把这个标为中等概率/中等影响,有明确升级路径。
**兜底:** "PI briefing 风险表里有备份计划——我查具体升级路径。"

### W3. "BALD 是聪明 pick,为什么还要随机对照?"
**自信回答:** 三个原因。**第一**,项目结束时回答"BALD 比随机好吗"的唯一办法——Hie 2024 把这立为标准。**第二**,校准合理性检查——随机持续胜过 BALD = 模型坏了。**第三**,你不知道哪个是对照,保留对比有效性。不是浪费——是我们能宣称 AL 有效(如果有效)的整个依据。
**兜底:** "标准对照臂模式——给你看 Hie 2024 的设置。"

### W4. "整批 ELISA 失败怎么办?"
**自信回答:** 触发质量门第 5 条——模型用部分数据(前一 cycle 累积集)重训,湿实验用 `safe_pick_backup.csv`,失败 cycle 视为时间线上的"损失一周"。Cycle 0 专门预留两周用于 ELISA 优化,所以单批失败不致命。升级路径在 `docs/onboarding/guide_en.md` §quality-gates。
**兜底:** "让我查——每个失败模式我们都写了升级路径。"

### W5. "Cycle 0 基因合成要 2 周,Cycle 1+ SDM 只要 4 天,为什么?"
**自信回答:** Cycle 0 从零开始——IDT/Twist 全长 DNA。150 美元 + 2 周。Cycle 1+ 拿 Cycle 0 的质粒做单点突变(NEB Q5 SDM)。50 美元 + 4 天。**3 倍便宜 3.5 倍快**。小数据 BALD 倾向挑点突变,正好匹配 SDM 经济性。BALD 哪天推荐巨大多位点重写,那个变体回到合成。
**兜底:** "成本和速度——克隆 slide 上有数字。"

### W6. "为什么 ELISA 用热灭活 *Xanthomonas* 而不是活细胞?"
**自信回答:** 两个原因。**安全/操作** ——超过 BSL-1 菌株时,活细胞在 ELISA 孔里加防控负担。**重复性** ——热灭活给稳定可存储的板,表面展示一致。代价是测的是变性受体表面的结合,筛选 OK,但 Tier 3 plaque assay 用活细胞确认受体特异感染。这是 Boeckaerts 2024 + Latka 2021 protocol;我们按已发表方法。
**兜底:** "Protocol 来自 Boeckaerts 2024——查论文具体理由。"

### W7. "ΔexbD2 是负对照——如果它也不感染怎么办?"
**自信回答:** 那就知道敲除系统本身有问题——可能删除的极性效应破坏邻近基因,或抗生素选择引入脱靶突变。ΔexbD2 的整个意义是**系统验证**:Hung 2003 显示它不重要,所以如果在我们手里重要了,其他敲除也可疑。要 PCR + 互补救援验证才信任 ΔtonB / ΔexbB / ΔexbD1 数据。
**兜底:** "好——这种情况我会想一起 debug。让我想想极性效应的可能。"

### W8. "怎么知道我的 Cycle 0 ELISA 板数据好?"
**自信回答:** Protocol 内置三个检查。**R² > 0.9** 4PL 拟合(低则 EC50 标低置信度)。**WT-RBP 阳性对照** 在历史 EC50 基线 2 倍以内(板间归一化)。**BSA 单空白背景** 低于阈值(~0.1 OD~450~)。三个都在 `qc_report.md` 模板里。任何一个失败,模型重训时降权那块板的数据。
**兜底:** "QC 模板在 `08_cycle_data/templates/` ——调出来。"

### W9. "我能加自选变体吗?"
**自信回答:** 可以——5 个 pick 之外的都标 `selection_reason = "expert_pick"`。模型当额外训练数据。请在 `qc_report.md` 注理由,这样我们能分析专家挑选是否系统性同意或偏离 BALD ——这也是可发表的。
**兜底:** "可以——具体标签我们需要规范,我跟进。"

### W10. "一个 cycle 你们绝对最低需要什么实验输出?"
**自信回答:** **3 个有效 EC50,R² > 0.9**。低于这个,模型拿不到足够新信号有意义更新。如果只 2 个有效,我们还重训但 PI briefing 标低置信 cycle,跳过 BALD 选择——下一 cycle 用 safe-pick backup。
**兜底:** "让我查——最低数据阈值在 cycle infra doc 里。"

---

## [干实验 Q&A] — Ryan、Leah、ESM 经验队员

### D1. "为什么 deep ensemble 不用 MC Dropout?"
**自信回答:** 校准。Ovadia 2019 (NeurIPS) 显示 MC Dropout 在 dataset shift 下校准差,deep ensemble 保持更好校准的预测区间。AL 来说,校准比纯精度更重要——校准差的 UQ 让 BALD 挑错变体。代价 5× 训练算力,我们这规模微不足道(Cycle 0 端到端 3 秒)。MC Dropout 唯一真正优势——单模型放进内存——在我们规模无关紧要。
**兜底:** "校准是头条差异——要细节我调 Ovadia 2019 的校准图。"

### D2. "ESM-2 8M 本地、650M Laguna ——规模买了什么?"
**自信回答:** Lin 2023 显示 embedding 结构质量随参数大致 sublinearly 提升,650M 之后回报递减。8M = **概念验证** ——无 GPU 跑通 embedding 管道。650M (1280 维 vs 320 维) 给真正更好的结构特征,是生产目标。3B 是 Cycle 4+ 最终 benchmark。PLM-interact 论文 (Liu 2025) 更有趣——在 PPI 数据上 fine-tune 650M 胜过原生 3B 的交互预测,提示 scale-vs-fine-tune 是我们想探索的轴。
**兜底:** "经验上的——让我调 Lin 2023 的 scaling 图。"

### D3. "没有真 ELISA 标签的冷启动怎么办?"
**自信回答:** 两个回答。**今天** ——Cycle 0 在合成数据上训练,结构有文档:低秩目标加噪声,只用来验证管道跑通且 BALD 评分排序符合预期。合成数据图在 calibration.png。**6 月 1 日起** ——首批真 ELISA 测量(4 个基因合成的结构 picks)成为种子训练集。Cycle 1 在那上面重训。我们事实上从一个 Cycle 0 batch bootstrap ——小但充分,因为每个成员被 Gaussian NLL regularize。
**兜底:** "从 Cycle 0 bootstrap ——让我走具体数据流。"

### D4. "模块之间的数据契约?"
**自信回答:** 每个模块文件夹有 `inputs/`、`processes/`、`outputs/`。`inputs/` 只读指向上游 `outputs/`。代码**只**在 `processes/` 里。Outputs 是规范化的,配 `MANIFEST.csv`(SHA-256 + 大小 + 记录数)。Schema 在 `INTERFACE.md`。最重要契约:predictions.csv 必须有 `rbp_id, receptor_id, predicted_score, std, epistemic_std` ——BALD 直接读 `epistemic_std`。
**兜底:** "`INTERFACE.md` 里——给你看。"

### D5. "为什么 BALD 用 `epistemic_std` 不用 total `std`?"
**自信回答:** Total `std` 包括 **aleatoric** 不确定性——不可减少的测量噪声(ELISA 移液方差)。按 total `std` 排,BALD 会偏向测量噪的对。这正好与我们要的相反——我们要**模型不确定**的对,不是**实验不确定**的对。`epistemic_std` 是 variance-of-member-means 分量;这是真正随数据收缩的。这在 `bald.py:38–52` 和 Lakshminarayanan 2017 eq. 3 分解里。
**兜底:** "Aleatoric vs epistemic ——分解在 `ensemble.py:286–294`。让我调出来。"

### D6. "测试覆盖怎么样?怎么知道管道没默默坏?"
**自信回答:** 每个模块有 pytest 套——模块 00–07 总共 140+ 测试。Module 07 (BALD) 18 测试,覆盖空池、全零分、已测对排除、种子可重现等边界。Module 06 9 测试,含 ensemble 多样性断言 (`frac_diverse > 0.5` 在 `run_cycle0.py:194`)。Module 05 28 测试 + 1 expected skip 等完成 GPU 运行。通过率在 "Tests — current pass rates" slide。
**兜底:** "每模块测试数在测试 slide 上;愿意走特定套。"

### D7. "为什么 Module 07 是 `.py` 其他是 notebook?"
**自信回答:** 生产编排,不是探索。Module 07 跑在 48 小时 SLA 上——需要是确定性、幂等的 CLI,能从 cron 或 shell 脚本调。Notebook 在 `输入 → 调试 → 输出` 探索上很好,在 `cron → CI/CD → 部署` 上糟糕。CLAUDE.md 明确标这例外。其他模块稳定后都会从 notebook 毕业到 `.py`;Module 07 一出生就是生产。
**兜底:** "CLAUDE.md 里的故意例外——生产代码路径需要 CI 友好。"

### D8. "重训时旧模型 checkpoint 怎么办?"
**自信回答:** 每 cycle 有自己子文件夹:`06_uncertainty_model/outputs/cycle_<N>/`。predictions.csv 里 model_version = `<git_sha>_cycle_<N>`,所以每个预测溯源到特定 git commit + cycle。旧 `ensemble_member_*.pt` 无限期保留(小文件;每 cycle <10 MB)。每 cycle 校准图让 cycle 间 diff 容易。MLflow tracking 将在 Cycle 1 接入。
**兜底:** "溯源通过 git SHA + cycle ——给你看 model_version 字段。"

### D9. "怎么扩展到通用噬菌体-宿主模型?"
**自信回答:** Module 04 embedding 步本来就是通用的——任何噬菌体 RBP 和任何细菌受体都同样 ESM-2 处理。要扩展,把 Module 06 ensemble 在完整 2236 对 Module 01 interaction matrix 上训(目前 Cycle 0 用合成)。Mutalik 2025 PAML benchmark 给我们 0.67–0.70 cross-genus AUC 做基线——超过就是贡献。Post-Jamboree 研究 roadmap 明确包括"在土壤来源噬菌体基因组上训练"为 Phase 2。
**兜底:** "项目自然延伸——给你看 post-Jamboree 计划。"

### D10. "你真的端到端跑过吗?"
**自信回答:** 跑过——模块 00–07 在笔记本上 5 分钟以内端到端,加上 Laguna 3 分钟 Boltz-2 跑。`DEMO.md` 里的 runbook 就是精确序列。唯一需要外部数据的步骤是 Module 08 (ELISA),6 月 1 日开始产出。所有输出 CSV 都已 commit(Module 06 起合成),所以契约锁定。
**兜底:** "演示马上来——任何模块愿意走。"

---

## [PI Q&A] — Prof. Ignacio-Espinoza、Prof. Libeskind-Hadas

### P1. "iGEM 奖项相关的差异化在哪?"
**自信回答:** 三个。**第一**,没有已发表系统把 ESM-2 + deep ensemble UQ + BALD 闭环 ELISA 配对到噬菌体 RBP × 细菌受体——Hie 2024 做抗体,Yang 2025 做酶,Boeckaerts 2024 做宿主范围预测不带 AL。我们第一个在这个领域整合。**第二**,Tier 3 受体敲除验证给的是因果而非相关——这是让它超出 iGEM 达到 paper 级别的部分。**第三**,in-house *Xanthomonas* 和噬菌体分离是给 iGEM Registry + NCBI 的社区贡献。
**兜底:** "组合新——让我列最接近的前作。"

### P2. "你怎么知道模型真的在学,不是随机?"
**自信回答:** 随机对照臂。每 cycle 5 个 picks 里 1 个是随机从未测池采样,对湿实验盲法。项目结束时,比较 BALD 的 test R² 和信息增益轨迹 vs 随机臂轨迹。如果无法区分,BALD 没帮。我们也每 cycle 跟踪 ECE 和后验 KL 散度——这些是模型内部合理性检查。整个比较预注册(在 PLAN.md + 风险 slide)。
**兜底:** "对照臂 + 预注册指标——愿意走实验设计。"

### P3. "什么失败模式会推翻整个方法?"
**自信回答:** 两个场景。**一** ——ELISA 定量信噪比太低,无法把有意义的结合增益从移液方差中区分(aleatoric 主导 epistemic)。Cycle 0 专门留两周 ELISA 优化对冲这个。**二** ——ESM-2 embedding 缺受体相关信号(比如特定表面 loop 重要,mean-pool 冲掉)。备份是 per-residue 或 attention-pooled embedding,cycle 间换。任一失败如实报告本身就是可发表负面结果。
**兜底:** "两大风险——ELISA 噪声 + ESM-2 信号。让我详化对策。"

### P4. "为什么这些论文不是别的?"
**自信回答:** 每个引用映射特定设计决定。Lakshminarayanan 2017 → ensemble;Houlsby 2011 → BALD 目标;Lin 2023 → ESM-2;Boeckaerts 2022 → PhageRBPdetect HMM;Hung 2003 → 定义 Tier 3 对照的受体系统。我们**不**靠 ALDE (Yang 2025) 作为 BALD 验证——他们用 Thompson sampling,不同 acquisition。slide 和 `paper_reading_notes.md` 都明确区分。5 月 2026 论文审计修正了早期草稿里五处错引——已修。
**兜底:** "每篇论文支撑一个具体决定——愿意逐引用走。"

### P5. "你的 iGEM 金牌关键路径?"
**自信回答:** 三个交付驱动奖。**Best Composite Part** ——wiki freeze 之前注册至少 4–6 个 RBP-His6 构建。**Best Model** ——闭环 AL 管道 + 至少 2 个湿实验 cycle 显示 BALD 行为(正面或如实负面)。**Best Agriculture Project** ——来自加州十字花科的 in-house *Xanthomonas* + 裂解噬菌体,测序沉积。Tier 3 受体敲除验证是 Model 赛道金牌级差异化。时间线紧但 5 月 17 开始克隆能完成。
**兜底:** "三个交付——让我走时间线。"

### P6. "iGEM 之后这代码怎么办?"
**自信回答:** 开源在现有 GitHub repo (active-learning-pipeline 分支是规范的——main 故意空)。5 月审计 + 双语注释 + per-module README 让它真的可 onboard。已讨论并发投稿 *Bioinformatics* 或 *NAR Genomics*(Cycle 0 前 PI 决定)。Repo、in-house 分离物序列、训练模型 checkpoint 都将是公共交付。
**兜底:** "稿件时间是待 PI 决定——我标记。"

### P7. "你跟团队外谁讨论过方法?"
**自信回答:** 还没正式。方法学综述在 `docs/planning/iGEM_2026_Project_Plan.md`,5 月 7 日 PI 咨询时 review。非正式接触到 ALDE 组(Caltech Frances Arnold 实验室)的已发表代码,需要时可借力。提交稿件之前,我想要至少一次外部方法学 review ——很乐意请你建议 reviewer。
**兜底:** "还没——很想要你建议的外部 reviewer。"

### P8. "让你失眠最大的风险?"
**自信回答:** 老实说:**5–6 月菌株分离失败**。下游所有——Tier 3 验证、composite part、农业项目——都链到有 in-house *Xanthomonas* 分离物加裂解噬菌体。我们有双源计划(brassica + citrus)和 Phage Directory 备份噬菌体源,但 5 月采样轮失败了,我们就得重推湿实验故事。这是桌面上唯一"没好备份"的风险。干实验能自跑;奖牌故事需要湿实验。
**兜底:** "菌株分离是关键路径风险——让我详化双源计划。"

---

## 加值:敌意听众问题(约 10 分钟焦虑模拟)

### H1. "ipTM 0.365 你说中等——这难道不意味着 Boltz-2 基本上在界面预测上失败了吗?"
**回答:** 是,对——Boltz-2 在这个界面上置信低。这是没被晶体学解析的新颖噬菌体-受体对的**预期**结果。模型训练数据里没有同源复合物。这次跑的价值不是自信结构——是**单体折叠**(chain pTM 0.808)给我们可信的 SDM 设计骨架,以及**界面不确定的承认** ——正好是主动学习要解决的。

### H2. "你的合成 Cycle 0 校准看起来完美——这不就是自洽剧场吗?"
**回答:** 正确——合成数据校准是合理性检查管道不破,不是模型质量验证。它告诉我们 `epistemic_std` 导出正确、`select_batch` 排除已测对、BALD 排名符合规范。真正校准测试是带真实 ELISA 数据的 Cycle 1+。我们在校准 slide 和 predictions.csv `data_source = synthetic_fallback_random` metadata 字段明确标记。

### H3. "为什么 PI 该信你对 Lee 2009 的解读高于已发表共识?"
**回答:** phiL7 tail spike 的"已发表共识"恰恰是 Lee 2009 说的——他们用对 OP1 ORF25 的 BLAST 同源找不到。我们不矛盾;我们说 HMM 搜索(更敏感方法,2022 Boeckaerts 发表)找到候选,他们用 2009 工具找不到。"rbp_01 = phiL7 tail spike" 是待测试假设——ELISA + ΔtonB 上的 plaque assay 直接测试。我们不宣称已确认。

### H4. "主动学习 1956 年起就有,本科 ML 课就教。这里到底什么是新的?"
**回答:** 组合 + 应用。BALD on deep ensemble regression with ESM-2 embeddings,闭环带定量 ELISA,在噬菌体 RBP × 细菌受体上,带受体敲除因果验证——这个精确栈我们查不到已发表先例。每个组件已知;集成和噬菌体系统的湿实验落地是贡献。reviewer 指出我们漏的前作,会立即更新文献综述。

### H5. "ELISA EC50 本身来自带置信区间的噪 4PL 拟合,你怎么处理 epistemic-aleatoric 分解?"
**回答:** 目前把 EC50 点估当 label,R² > 0.9 当质量过滤。4PL 拟合的参数 SE 可作为 Gaussian NLL 的观测噪声接入——这是 Cycle 2+ 精化。当前我们对噪板降权(R² 0.7–0.9 半权进入;<0.7 排除)。Per-measurement aleatoric 在 post-Cycle-1 roadmap。

### H6. "Hie 2024 每个抗体只用 20 个变体——你为什么承诺每 cycle 4–5 个 × 多个 cycle?"
**回答:** 不同问题。Hie 2024 是语言模型 likelihood 过滤——一批候选变体排序一次,不是闭环。我们的 ALDE 类比 (Yang 2025) 用更大 batch(~90 / 轮,2–3 轮)。我们居中——每 cycle 小 batch 因为 Cycle 0 ELISA 优化限制吞吐,多 cycle 后总量大。湿实验吞吐增长,batch 也增。

### H7. "你的 model_version 是 git SHA + cycle。半年后想复现 Cycle 1 怎么办?"
**回答:** `git checkout <sha>` 重建代码。`06_uncertainty_model/outputs/cycle_1/ensemble_member_*.pt` 有训练权重。合成数据用种子 RNG seed=42 在 model_meta.json 记录。真实数据复现需要 `08_cycle_data/outputs/cycle_1/elisa_processed.csv` 的 ELISA CSV,gitignore 但 MANIFEST.csv 有 checksum。半年后从这三个 artifact 复现。

### H8. "你说湿实验不知道哪个 pick 是随机。实际怎么强制?"
**回答:** `recommendations.csv` 行打乱过且 `selection_reason` 列**不**发给湿实验——只发 `rbp_id, receptor_id, predicted_score, primer_sequence`。带 `selection_reason` 的全 CSV 归档在干实验 outputs 给回顾对比。是荣誉系统——Sarah、Olivia、Weitao、Carol 不能偷看——但文件分离让操作上容易。
