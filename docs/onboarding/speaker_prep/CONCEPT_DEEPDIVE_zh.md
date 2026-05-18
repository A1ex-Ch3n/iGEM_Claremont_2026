# CONCEPT_DEEPDIVE_zh.md — Alex 必须吃透的七个概念

每个概念:30 秒电梯版(PI/湿实验)· 2 分钟版(干实验)· 完整机制 · 类比 · 一句必须背下来 · 已知失败模式。

---

## 1. 为什么主动学习能解决数据稀缺(对比随机 / 监督学习)

### 30 秒电梯版
每次 ELISA 花 50 美元 + 4 天。随机挑变体时,大多数测量落在模型已经懂的区域。主动学习相反——挑模型**最不确定**的变体,这些测量信息含量最高。在邻近领域(Hie 2024 抗体、Yang 2025 酶)已经证明:每次实验进步 2–5 倍。

### 2 分钟版
固定数据集上训练的监督模型,被该数据集的偏差所困。主动学习把环路关上:模型用**自己的不确定性**作为下一次测量的信号。形式上,有了数据 `D` 后得到参数后验 `p(θ | D)`。在输入 `x` 上的下一次测量,如果 `y(x)` 能降低后验熵就有信息价值。Lindley 1956 的最优实验设计目标就是这个。BALD 把它操作化:挑 `x` 最大化互信息 `I(y; θ | x, D)`。实践中就是"挑 ensemble 成员分歧最大的点"。对比随机:变体空间巨大(712 aa 蛋白 = 20^712)时,随机击中信息丰富的变体概率极低。主动学习按设计击中——Yang 2025 两个 cycle 把 yield 从 12% 提到 93%。

### 完整机制
1. 在 `D` 上训练 → 后验 `p(θ | D)`。
2. 池中每个候选 `x`,计算 `acq(x) = I(y; θ | x, D)`,高斯下有闭式。
3. 按 `acq(x)` 取 top-K。
4. 湿实验测量 → `y_obs`。
5. 更新 `D ← D ∪ {(x, y_obs)}`,重训练,迭代。
"闭环" = 迭代;"主动" = 模型(不是人)挑 `x`。

### 类比
- **医生诊断**:不要把所有检查都开,只开那个结果最能缩小诊断范围的。
- **二十问游戏**:每个问题应当把剩余可能性砍半,不要测试随意假设。
- **雾中声纳**:打你最不确定的方向,不是你已经看清的方向。

### 必背一句
> "主动学习挑下一个最大化模型不确定性减少的实验——所以每次 ELISA 都买到最多的信息。"

### 失败模式
- **校准错误** ——模型在不该自信的地方自信,BALD 挑错地方。对策:温度缩放,每个 cycle 监控 ECE。
- **模型偏差是根本性的** ——特征空间错(比如 ESM-2 漏掉受体特异信号)。BALD 那就在错的问题上高效问。对策:随机对照臂 + Tier 3 因果验证。
- **Aleatoric 噪声主导** ——ELISA 噪声太大,epistemic 信号被淹没。对策:3 个技术重复 + ECE 守门。

---

## 2. ESM-2 ——蛋白质语言模型 embedding 实际是什么

### 30 秒电梯版
ESM-2 是个神经网络,在约 6500 万蛋白序列上做"完形填空"训练——类似 GPT,但针对蛋白质。训完后,任何新蛋白扔进去,提取模型内部的数字——一个 1280 维的 **embedding**,捕捉了"在见过 6500 万其他蛋白的模型眼中"这个蛋白长什么样。功能相似的蛋白在这 1280 维空间里彼此靠近,即便它们的氨基酸字母不匹配。

### 2 分钟版
ESM-2 (Lin 2023 *Science*) 是个 transformer,在 UniRef50/90 上做 **masked language modelling**:随机遮 15% 氨基酸,让模型根据上下文预测。要做好这个,模型必须在内部学到生化——疏水簇、二级结构 motif、进化保守性。训完后,我们不用预测的氨基酸,用的是 **per-residue hidden states** ——1280 维向量(650M 模型),编码了位置感知的结构和功能上下文。**Mean-pool** 跨残基得到每个蛋白一个向量。在我们的管道里,这个向量替代手工特征(电荷、疏水性等),作为 deep ensemble 的输入。

### 完整机制
架构:33 层 transformer (650M 参数), 1280 维 hidden state。训练目标:masked AA 预测 cross-entropy。训练数据:UniRef50/90,~6500 万独立序列。推理:encoder forward pass,提取 `hidden_state[layer=33]` 形状 `(L, 1280)`。Mean-pool over L:`embedding = hidden_state.mean(dim=0) → (1280,)`。然后拼接 `embedding_RBP ⊕ embedding_receptor → (2560,)` 作为 MLP 输入。

### 类比
- **Spotify 推荐 embedding**:你从没听过的歌,基于模型从几百万首学到的模式,落在相似歌附近。
- **词向量 (word2vec)**:"king − man + woman ≈ queen" ——空间几何编码了语义。
- **声音频谱图**:原始波形看不出东西;频率-时间表示让模式可见。

### 必背一句
> "ESM-2 是个 transformer,训练目标是在 6500 万蛋白上填空——我们用的是内部表征,不是预测,因为表征捕捉了结构和功能上下文。"

### 失败模式
- **Domain shift**:ESM-2 训练数据主要是细胞蛋白。噬菌体 tail spike 是奇特结构(β-helix、三聚体交锁)——embedding 可能低估关键部分。
- **Mean-pooling 丢失位置信息**:距离结合位点 200 aa 的点突变和结合位点上的点突变全局均值相同。可以用 per-residue 或 attention-weighted 替代方案。
- **规模重要但 sublinear**:8M → 650M 跨越大;650M → 3B 小。PLM-interact (Liu 2025) 暗示在 PPI 数据上 fine-tune 比纯放大更有效。

---

## 3. Boltz-2 ipTM ——最大的误解风险

### 30 秒电梯版
Boltz-2 预测蛋白复合物的 3D 结构,给一个叫 **ipTM** (interface predicted TM-score) 的分数,0–1 之间。**ipTM 是置信度衡量** ——模型对两条链怎么对接有多确定。它**不是**结合亲和力。Boltz-2 确实有 affinity head,但训练数据是**小分子**(PubChem、ChEMBL、BindingDB)——蛋白-蛋白对 affinity 输出 `NaN`。我们 rbp_01 × TonB 的 ipTM = 0.365:界面几何置信度低。这不是说"它们不结合",而是模型不知道怎么结合。

### 2 分钟版
TM-score (template modeling score) 衡量结构相似性,0–1,>0.5 通常意味着同 fold。**pTM** 是结构模型内部置信度预测的 TM-score——模型自评"假设我有真实结构,我的预测有多像?"。**ipTM** 把这个限定在两条链的**界面区域** ——对接几何被预测得多好。Boltz-2 (Passaro 2025) 在 `affinity.json` 报告 ipTM。模型主要在 PDB + distillation 上训练;蛋白-蛋白产生带 ipTM 的结构预测。另外,Boltz-2 有一个 **affinity head** 预测结合自由能,但这个 head 在小分子结合数据库上训练。喂蛋白-蛋白时,affinity head 输出 `NaN`,因为没有学过这种映射。我们的 `predicted_dG = null`。Boltz-2 论文和我们的 build report 都确认。

### 完整机制
- 输入:两条氨基酸序列(或一条 + ligand)。
- 模型:diffusion-based 原子结构预测 + 置信头。
- 输出:
  - `*.pdb` ——复合物原子坐标
  - `chain_X_ptm` ——单体 pTM (rbp_01: 0.808 = 折叠可信)
  - `interface_ipTM` ——界面置信度 (我们 0.365 = 低)
  - `confidence_score` ——加权总体 (我们 0.683)
  - `pae_*.npz` ——pairwise alignment error 矩阵
  - `plddt_*.npz` ——per-residue 局部置信度 (0.30–0.98, 均值 0.76)
  - `predicted_dG = null` 蛋白-蛋白

对于我们的 rbp_01 × TonB:ipTM 0.365 = 模型确实不确定 rbp_01 在 TonB 哪里对接。单体预测好(chain pTM 高)。这和生物学一致——Hung 2003 已知 phiL7 通过 TonB 进入,但具体结合界面从未被晶体学解析过。

### 类比
- **GPS 置信圆圈**:ipTM 是"95% 把握你在这里面"圆圈的半径——小圆圈 = 高 ipTM。圆圈不告诉你能不能找到停车位。
- **天气预报概率**:"降雨概率 70%" 是置信度,不是"会下多少雨"。ipTM 是预测置信度,不是结合有多强。
- **认人队列**:ipTM 0.365 = "我觉得是这三张脸之一" ——模型缩小但没定论。

### 必背一句
> "Boltz-2 的 affinity head 只针对小分子——蛋白-蛋白对输出 NaN。我们用 ipTM 作为结构置信度 proxy,不是结合亲和力。"

### 失败模式
- **把 ipTM 当亲和力说** ——最大的可信度风险;反复背这句。
- **低 ipTM ≠ 不结合** ——Hung 2003 已经确认 phiL7 通过 TonB 进入。低 ipTM 是预测不出几何,不是没有相互作用。
- **高 ipTM ≠ 强结合** ——即便 0.9 ipTM 也只说对接置信高,不说亲和力。置信度和强度是正交的。

---

## 4. Deep ensemble ——为什么 5 个网络给 epistemic 不确定性

### 30 秒电梯版
我们用不同随机种子把同样的神经网络架构训 5 遍。5 个模型在数据告诉它们怎么做的地方一致,在数据不告诉它们的地方**分歧**。预测的散度 = epistemic 不确定性 = "模型不知道什么"。这是神经网络量化模型无知最简单最可靠的方法。

### 2 分钟版
单个神经网络给点预测——没有诚实的不确定性信号。**贝叶斯神经网络**给权重分布但计算难。**MC Dropout** 近似但校准差(Ovadia 2019)。**Deep ensembles** (Lakshminarayanan 2017) 极简:用 K 个随机种子训 K 个模型,预测之间的方差近似后验预测方差。我们的设置中,每个成员输出**高斯** ——`μ_k(x)` 和 `σ_k(x)`。**Epistemic** 方差 = `Var_k[μ_k(x)]` ——成员均值的散度。**Aleatoric** 方差 = `E_k[σ_k²(x)]` ——成员内方差的平均。Total = epistemic + aleatoric (Lakshminarayanan eq. 3,在 `ensemble.py:286–294` 实现)。K=5 选择:原论文 5 个之后回报递减;ALDE 也用 5。Greenman 2025 提醒:ensemble 通常 accuracy 最好但 calibration 最差——ECE > 0.1 时用温度缩放补救。

### 完整机制
对每个成员 k ∈ {0, …, 4}:
1. 随机种子 = k(控制权重 init 和 DataLoader shuffle)。
2. 初始化 3 层 MLP (256/256/128, ReLU, dropout 0.1)。
3. 两输出头:`head_mean → μ_k`、`head_log_sigma → log σ_k`,log_sigma 截到 [-7, 7]。
4. Gaussian NLL loss 训练,max 200 epochs,val NLL early stop (patience 10)。
5. 保存 state_dict。

预测 (`DeepEnsemble.predict()`):
- 5 个成员 forward → `means[5, N]`, `sigmas[5, N]`。
- `epistemic_var = means.var(axis=0)`
- `aleatoric_var = (sigmas**2).mean(axis=0)`
- `total_std = sqrt(epistemic_var + aleatoric_var)`
- `epistemic_std = sqrt(epistemic_var)` ← BALD 用这个。

### 类比
- **五位医生第二意见**:全一致 = 诊断定;分歧 = 分歧本身就是信号说要更多检查。
- **五个 GPS 来自不同卫星**:一致 = 位置低不确定;分歧 = 你在覆盖差的区域。
- **Bootstrap 重采样**:经典统计从单一数据集得置信区间的方法。

### 必背一句
> "5 个 MLP,同样架构,不同随机种子——它们在数据告诉它们怎么做的地方一致,不告诉的地方分歧。这个分歧就是 epistemic 不确定性。"

### 失败模式
- **成员塌缩到同一解** ——在 `run_cycle0.py:194–199` 通过 `frac_diverse > 0.5` 检查。塌缩则 BALD 全 0。
- **校准漂移** ——ensemble 在 out-of-distribution 输入上倾向过度自信(Greenman 2025)。
- **Aleatoric vs epistemic 混淆** ——如果你不小心用 `total_std` 当 BALD,会挑"测量噪声大"的变体,不是"模型不确定"的。

---

## 5. BALD 数学 ——回归下的 `Std_k[μ_k]`(推导)

### 30 秒电梯版
BALD = "Bayesian Active Learning by Disagreement"。挑给你最多模型参数信息的实验。对高斯输出的 deep ensemble 回归,这化简为:挑 ensemble 成员均值分歧最大的输入。数学上:`BALD(x) ≈ Var_k[μ_k(x)] → score = epistemic_std`。

### 2 分钟版
Houlsby et al. 2011 在高斯过程**分类**上引入 BALD。目标:`argmax_x I(y ; θ | x, D)`,`I` 是 label 和模型参数的互信息。论文 equation 2 改写为 `H[y | x, D] − E_{θ ~ p(θ|D)}[H[y | x, θ]]`。**第一项** = 预测熵(总不确定性)。**第二项** = 期望的 per-parameter-sample 熵(只有 aleatoric——固定 θ 移除参数不确定性)。差 = epistemic 不确定性。

对我们的设置——高斯 deep ensemble 回归——代入高斯熵。Per-member 高斯熵 `½ log(2π e σ_k²)`。混合高斯(K 个成员平均)预测方差 `Var_k[μ_k] + E_k[σ_k²]`。代入:
- `H[y | x, D] ≈ ½ log(2π e (Var_k[μ_k] + E_k[σ_k²]))`
- `E_θ[H[y | x, θ]] = ½ log(2π e) + ½ E_k[log σ_k²]`
- 差(忽略 `log` Jensen 不等式松弛)≈ `½ log(1 + Var_k[μ_k] / E_k[σ_k²])`,当 aleatoric 在 `x` 间大致恒定时,关于 `Var_k[μ_k]` 单调。

我们的体制下,`E_k[σ_k²]` 在池中大致稳定(同方差 per-pair 噪声),所以按 `Var_k[μ_k]` 排名等价于按完整 BALD 分数排名。我们取**标准差** (`epistemic_std`) ——单调,与预测同单位。

**从 GPC 分类延伸到 deep ensemble 回归是我们的工作** ——引用 Houlsby 作为原始目标,化简视作我们的贡献。

### 完整推导(显式高斯)
Per-member predictive:`y | x, θ_k ~ N(μ_k(x), σ_k²(x))`。
混合 predictive:`y | x, D ~ (1/K) Σ_k N(μ_k(x), σ_k²(x))`。
- `Var[y | x, D] = E_k[σ_k²] + Var_k[μ_k]`(总方差定律)。
- 混合的 predictive entropy **上界**为同方差高斯熵:`H[y | x, D] ≤ ½ log(2π e (E_k[σ_k²] + Var_k[μ_k]))`。
- `E_θ[H[y | x, θ]] = E_k[½ log(2π e σ_k²)]`。
- `BALD(x) ≈ ½ log( (E_k[σ_k²] + Var_k[μ_k]) / exp(E_k[log σ_k²]) )`。

当 `σ_k(x) ≈ σ`(per-x 同方差):
- 分子 → `σ² + Var_k[μ_k]`。
- 分母 → `σ²`。
- BALD → `½ log(1 + Var_k[μ_k] / σ²)`,关于 `Var_k[μ_k]` 单调。所以按 `Var_k[μ_k]` 排 = 按 `epistemic_std` 排。QED。

### 类比
- **众包**:问群体最分裂的问题——在那上面测量真相对所有人信念偏移最大。
- **投资**:把仓位再平衡到你最不确定价值的那笔——边际价格发现在那里最重要。
- **三角测量**:测量员不在两条视线已经交叉的地方测;在交叉最模糊的地方测。

### 必背一句
> "我们 deep-ensemble 回归的 BALD = `epistemic_std = Std_k[μ_k]` ——ensemble 成员均值的标准差。Houlsby 2011 在 GPC 分类上推导,deep-ensemble 回归形式是我们的扩展。"

### 失败模式
- **异方差噪声** ——`σ_k(x)` 变化大,同方差假设破,按 `Var_k[μ_k]` 排不再等价于完整 BALD 排名。
- **后验近似** ——deep ensemble 是真贝叶斯后验的粗糙近似。校准审计(Ovadia 2019, Greenman 2025)指出。
- **批内相关** ——top-4 BALD picks 在相似输入上 epistemic 不确定性相关,可能冗余。Batch-BALD (Kirsch 2019) 解决——我们暂未用。

---

## 6. 48 小时闭环 ——为什么 SLA 重要

### 30 秒电梯版
湿实验 cycle 端到端 10–14 天(克隆→表达→ELISA)。干实验 cycle 48 小时(在新 ELISA 上重训 → BALD → 下一轮推荐)。干实验拖,湿实验要么空转(浪费时间)要么用旧推荐(浪费信息)。48 h SLA 让闭环真"闭" ——不只是顺序的。

### 2 分钟版
天真的想法是"模型想花多久就花多久"。但湿实验按固定 cycle 运行——SDM 引物订购、质粒转化、诱导、纯化、ELISA、回归。整个流水线 ~10–14 天。新模型推荐**在**湿实验计划下一轮**之前**到达就进 Cycle N+1。**之后**到达就两个坏选项:湿实验空转(损失一周吞吐)或用过期推荐(损失刚结束 cycle 的信息)。48 小时是让湿实验不停又不用过期数据的缓冲。操作上 SLA 由这些强制:预计算 `safe_pick_backup.csv`(上一 cycle top-10 BALD picks——只在 SLA miss 时用)、MLflow 跟踪训练时间、文档化的升级路径。Cycle 1 拖 → 下游拖。三次 cycle 拖 → 错过 iGEM wiki freeze。

### 完整机制
Cycle N → Cycle N+1 时间线:
- **Day 0**:ELISA 数据 finalized,CSV 上传到 `08_cycle_data/outputs/cycle_<N>/elisa_processed.csv`。
- **Day 0+0h**:干实验拉数据,若新变体则重生 ESM-2 embedding,重训 ensemble (Module 06)。
- **Day 0+24h**:重训完,校准检查 (ECE ≤ 0.1 或应用温度缩放)。全池预测导出。
- **Day 0+30h**:BALD 打分 + 选择 (`run_bald.py`)。
- **Day 0+36h**:PI 交叉检查(top picks 合理性 review)。
- **Day 0+48h** (SLA):`recommendations.csv` + `primer_sequences.txt` 交付湿实验。
- **Day 0+48h+**:湿实验订 SDM 引物 (NEB Q5 设计工具,接进 `primer_sequences.txt`),开始 Cycle N+1。

质量门:ECE > 0.1 → 温度缩放;上一 cycle SDM 失败率 > 50% → 备份 picks;BALD 分数方差缩小快于预期 → flag 给 PI(模型可能过拟合)。

### 类比
- **准时制制造**:Toyota kanban——下游需要时才拉,不积压陈旧信息。
- **电视烹饪节目两版同时做**:第一版反馈来之前先做第二版;反馈晚到第二版就是盲做的。
- **CI/CD 流水线**:代码 commit → 测试 → 部署。测试比下次 commit 久就要排队或跳。

### 必背一句
> "48 小时 SLA 让闭环关上——湿实验不等我们,我们不推过期推荐。Safe-pick backup 是降落伞,不是默认。"

### 失败模式
- **重训不收敛** ——载入旧模型,新超参重试,推到 PI 之前升级。
- **ECE 爆炸** ——温度缩放;还差就退到 safe-pick。
- **ELISA 数据格式错** ——schema validation 直接在管道里拒掉,清楚的报错。
- **累积拖延** ——三次 cycle 拖 = wiki freeze 风险。我们视拖延为关键路径 PI 升级。

---

## 7. Lee 2009 HMM vs BLAST ——为什么 rbp_01 是补充,不是矛盾

### 30 秒电梯版
Lee 2009 ——phiL7 基因组论文——明确搜索 OP1 ORF25 的 tail fiber 同源物,**在文中写明找不到**用序列同源。我们的管道用 Hidden Markov Model 发现 rbp_01。HMM 抓 BLAST 漏的——分歧到字母不匹配但结构模式还在的蛋白。Lee 2009 用了 2009 年合适的工具;我们有更敏感的工具。结果是补充,不是矛盾。

### 2 分钟版
BLAST(和 Lee 2009 用的搜索工具)通过**序列相似性**找蛋白同源物——氨基酸字母对齐,用 BLOSUM62 这类替换矩阵打分。序列 >25–30% 一致时敏感;以下急剧失敏。HMM (Hidden Markov Models, profile HMM 通过 HMMER 工具) 反过来,从多序列比对构建一个**蛋白家族**的概率模型——捕捉哪些位置保守、哪些容忍替换、哪些插入/删除 gap。序列分歧但结构保守的蛋白,HMM 在 BLAST 失明之后仍然能检测到。PhageRBPdetect (Boeckaerts 2022) 提供 tail-spike domain HMM profile——包括 `Tail_spike_N`。rbp_01 hit 这个 profile 分数 342(很高)。它可能就是 Lee 在找的那一类 tail spike,只是序列分歧到 2009 BLAST 看不到。**说"Lee 2009 漏掉了"是错的** ——他们用了那时合适的工具。说"rbp_01 推翻 Lee 2009"也是错的——Lee 找的是 OP1-ORF25 *同源物* 找不到;我们没说 rbp_01 是那个同源物,只是说它是用更敏感方法找到的 tail spike 候选。

### 完整机制
BLAST:在对齐残基上 `S = Σ_i s(a_i, b_i)` 打分,然后对随机数据库算 E-value。阈值典型 E < 0.01。序列同一性低于 ~30% 时敏感度急降。

Profile HMM (HMMER `hmmsearch`):从 N 个同源物构建 profile,每个位置有 emission 概率(允许哪些氨基酸)和 transition 概率(插入/删除)。分数 = log-likelihood 对 profile,校准到 bit score。Domain-level hit:profile 捕捉到的短保守 motif 也能注册。

关键:profile HMM 用 **N 个同源物的信息**,不只查询对单一目标。Boeckaerts 2022 从几百种噬菌体的人工策选 tail-spike domain 对齐建 `Tail_spike_N`。我们 rbp_01 hit 分 342——比可信阈值高数量级。

Lee 2009 把 OP1 ORF25 当模板跑 BLAST——单一模板。rbp_01 vs OP1 ORF25:不 BLAST hit。rbp_01 vs Tail_spike_N HMM(多同源物建):hit。两者不冲突。

### 类比
- **警察素描师 vs 人脸识别**:素描是一张图(BLAST),人脸识别用上千张训练人脸(HMM)——能认出素描认不出的伪装的人。
- **品尝从没喝过的葡萄酒**:知道一款特定 Bordeaux (BLAST 模板) vs 喝过 1000 款酒 (HMM profile) ——后者能抓家族相似性即便葡萄比例新颖。
- **翻译**:知道一句双语句子 vs 知道语法——语法能推广到没见过的句子。

### 必背一句
> "Lee 2009 用 BLAST 明确写明找不到 OP1 ORF25 同源物。我们用 Hidden Markov Model ——更敏感的工具——找到了 rbp_01。两个结果是补充,不是矛盾。"

### 失败模式
- **HMM false positive**:rbp_01 可能偶然 hit profile 而不是真 tail spike。对策:Boltz-2 单体 pTM 0.808 支持合理 fold;ELISA 直接测对 TonB 的结合。
- **rbp_01 不是结合用的 RBP**:phiL7 可能有多个 host-range 蛋白;rbp_01 可能是不结合 TonB 的结构 tail 蛋白。在 ΔtonB 上的 plaque assay 测试这个。
- **序列分歧真实但无关**:rbp_01 是对的家族但不特异结合 *Xanthomonas* TonB。整个管道就是设计来测试这点。
