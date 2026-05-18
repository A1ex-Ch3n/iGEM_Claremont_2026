---
title: "针对黄单胞菌生物防治的噬菌体主动学习工程"
subtitle: "iGEM Claremont 2026 — 团队入职培训"
author: "核心工程师：Alex Chen · 导师：J. Cesar Ignacio-Espinoza · 学术顾问：Ran Libeskind-Hadas"
date: "2026 年 5 月"
mainfont: "Hiragino Sans GB"
CJKmainfont: "Hiragino Sans GB"
monofont: "Menlo"
fontsize: 10pt
aspectratio: 169
classoption:
  - "aspectratio=169"
theme: "default"
colortheme: "seahorse"
innertheme: "rectangles"
outertheme: "infolines"
slide_level: 2
header-includes:
  - \usepackage{xeCJK}
  - \setCJKmainfont{Hiragino Sans GB}
  - \usepackage{fancyvrb}
  - \usepackage{booktabs}
  - \usepackage{xcolor}
  - \usepackage{mdframed}
  - \definecolor{cmcblue}{HTML}{1A4480}
  - \definecolor{wetlab}{HTML}{A01A4D}
  - \setbeamercolor{frametitle}{bg=cmcblue,fg=white}
  - \setbeamercolor{title}{fg=cmcblue}
  - \setbeamertemplate{navigation symbols}{}
  - \newenvironment{noteblock}{\begin{mdframed}[backgroundcolor=cmcblue!8,linecolor=cmcblue,linewidth=0.5pt,roundcorner=2pt,innertopmargin=4pt,innerbottommargin=4pt]}{\end{mdframed}}
  - \newenvironment{warnblock}{\begin{mdframed}[backgroundcolor=red!5,linecolor=red!60,linewidth=0.5pt,roundcorner=2pt,innertopmargin=4pt,innerbottommargin=4pt]}{\end{mdframed}}
---

# 第 0 部分 —— 路线图

## 三句话总结

- 针对 *Xanthomonas campestris* pv. *campestris* 的噬菌体 RBP 工程的**闭环主动学习**系统。
- **Module 00–07 已全部完成**（计算部分）；Module 08（湿实验室数据）大约 6 月 1 日启动。
- 首个 Boltz-2 结果 `rbp_01 (712 aa) × TonB`：**ipTM = 0.365，chain_A_ptm = 0.808**。
- BALD 在**不到 1 秒**内推荐下一轮 4 个 variant + 1 个随机对照。

\vspace{0.5em}
\small *目标：Best Agriculture Project · Best Model · Best Composite Part。*

## 路线图

1. **科学背景** —— 黄单胞菌、噬菌体宿主范围、数据稀缺问题。  `[湿实验室] [PI]`
2. **流程架构** —— Module 00–08，数据约定。                     `[干实验室] [PI]`
3. **ML 核心** —— 深度集成、BALD、代码与输出样本。             `[干实验室] [PI]`
4. **当前 Boltz-2 结果** —— 0.365 ipTM 意味着什么。              全员
5. **48 小时循环** —— 干湿对接、ELISA、敲除。                   全员
6. **复现与演示** —— 快速开始、测试、Laguna。                   `[干实验室] [湿实验室]`
7. **风险、决策与请求**                                          `[PI]`
8. **参考文献与附录**                                            全员

## 如何阅读本演示

每个章节标题后的标签告诉你这部分给谁看：

- **`[湿实验室]`** —— 生物学动机、实验流程、需要执行的对接。
- **`[干实验室]`** —— 代码、算法、文件约定、可重现性。
- **`[PI]`** —— 整体框架、需要决策的事项、风险。

未标记的部分面向全员。

## 推荐分段（presenter cuts）

- **30 分钟湿实验室入门：** slide 1–4，7–12，17–22（模块速览），37–44（48 小时循环与实验流程）。
- **45 分钟 PI / 顾问简报：** slide 1–4，5–12（科学背景），13–14，23，28–32（ML 核心），33–36（Boltz-2），51–55（风险与决策），56（交付物）。
- **90 分钟全员深入：** 全部 slide。

演示（约 15 分钟）独立于此 deck，可在第 48 张后或最后运行（见 `DEMO.md`）。

# 第 1 部分 —— 科学背景 · `[湿实验室] [PI]`

## *Xanthomonas* —— 是什么、为何重要

- 植物病原性 γ-变形菌属；侵染 >400 种宿主植物（Ryan 2011，*Nat Rev Microbiol*）。
- 加州相关致病变种：
  - **Xcc** —— 十字花科黑腐病（甘蓝、西兰花、羽衣甘蓝）。
  - *X. citri* subsp. *citri* —— 柑橘溃疡病。
  - *X. perforans* / *X. euvesicatoria* —— 番茄细菌斑病。
- 当前防治依赖**铜制剂**；耐铜性已广泛出现（Aiello 2019，*Plant Dis*）。
- 田间试验显示噬菌体防治可达铜制剂效果（Iriarte 2018；Holtappels 2022）。

## 噬菌体生物防治 —— 宿主范围问题

- 噬菌体**宿主特异性强** —— 既是优点（无脱靶微生物）也是缺点（每个菌株需要对应的噬菌体）。
- 预测哪种噬菌体能感染哪种细菌仍不可靠：
  - **PhageHostLearn**（Boeckaerts 2024，*Nat Commun*）**菌种内 AUC 0.82**。
  - **跨菌属 AUC 降至 0.67–0.70**（PAML benchmark，Mutalik 2025）。
- 瓶颈：有标签的（噬菌体，宿主）相互作用数据稀缺、昂贵、物种受限。

## 我们的参考骨架

\begin{center}
\begin{tabular}{ll}
\toprule
宿主 & \textbf{Xcc ATCC 33913} \quad NCBI：GCF\_000007145.1 \\
噬菌体 & \textbf{phiL7} \quad NCBI：EU717894 \\
来源 & da Silva 2002（\emph{Nature}）；Lee 2009（\emph{AEM}） \\
\bottomrule
\end{tabular}
\end{center}

\vspace{0.4em}

- 干实验室使用上述公共参考。
- **湿实验室从加州十字花科作物自分离** —— 绕过 4 个月的 USDA APHIS PPQ-526 许可证（PI 协商，2026-05-07）。

## phiL7 受体系统 —— Hung 2003

**Hung, C.-H. 等 (2003) *BBRC* 302:878–884，PMID 12646254。**

- phiL7 通过 **TonB-ExbB-ExbD1** 复合体进入 Xcc（能量耦合的外膜转运）。
- Tn5 突变 + 互补实验：**TonB、ExbB、ExbD1 必需**。
- ΔexbD2（菌株 CH620）保留完全敏感性 → **ExbD2 不是必需**。

\begin{noteblock}
\textbf{免费的阴性对照。} 2026 年 5 月文献核查修正了早先错误将 ExbD2 列为必需的引用。构建 $\Delta$exbD2 与 $\Delta$tonB / $\Delta$exbB / $\Delta$exbD1 平行 —— 可验证整个敲除系统：$\Delta$exbD2 应当\emph{仍允许}感染。
\end{noteblock}

## phiL7 的 RBP —— Lee 2009 与我们的 HMM 重新发现

**Lee, C.N. 等 (2009) *AEM* 75:7828** —— phiL7 基因组论文。

- Lee 2009 **主动搜索**了 OP1 tail fiber（ORF25）的同源物，**通过序列相似性找不到**：

  > *"We were unable to identify a homolog of the OP1 tail fiber protein (ORF25)…"*

- p20（1105 aa）被建议与宿主范围相关；论文中没有任何蛋白被命名为"tail spike"。
- **我们的 rbp_01（712 aa）**来自 PhageRBPdetect 的 Tail_spike_N HMM —— 一种结构 profile 方法，能找到序列已分歧到 BLAST 找不到的蛋白。
- **与 Lee 2009 不矛盾 —— 这是用更敏感工具的互补发现。**

## 数据稀缺瓶颈

- 我们的系统**只有 1 个**实验确认的（噬菌体，宿主）相互作用：phiL7 × Xcc。
- 相互作用矩阵（Module 01）有 2,236 个文献整理对 —— 但没有定量结合亲和数据。
- 对于 712 aa 蛋白，变种空间为 $20^{712}$ —— 远超任何湿实验预算。
- **需要一种让每次昂贵 ELISA 测量都尽可能多产生信息的学习策略。**

## 为什么主动学习是正确的框架

- **主动学习（AL）：** 模型通过最大化预期信息增益来选择下一个训练点（Lindley 1956；Settles 2009）。
- 在相邻领域已证明 2–5× 优于随机选择：
  - **Hie 2024**（*Nat Biotechnol*）—— ESM-1b/1v + AL 进行抗体亲和成熟（每个抗体约 20 个 variant）。
  - **Yang 2025 ALDE**（*Nat Commun*）—— DNN 集成 + **Thompson 采样**（**不是** BALD）做酶定向进化，产率 12% → 93%，2 轮 wet lab。
- **没有人将其应用于噬菌体 RBP × 细菌受体。** 这正是我们的方法学贡献。

## 对应到 iGEM 的赛道

\begin{tabular}{ll}
\toprule
\textbf{赛道} & \textbf{我们的交付物} \\
\midrule
Best Agriculture Project & 针对 Xcc 的噬菌体生物防治流程 \\
Best Model               & 闭环 AL + UQ + 因果验证 \\
Best Composite Part      & 注册的 RBP-His6 表达库 \\
\bottomrule
\end{tabular}

\vspace{0.6em}

外加：自分离的 *Xanthomonas* + 噬菌体菌株，测序并存入 NCBI。

# 第 2 部分 —— 流程架构 · `[干实验室] [PI]`

## 三层架构

\begin{center}
\begin{tabular}{ll}
\toprule
\textbf{第 0 层 —— 先验}      & Boltz-2 ipTM、ESM-2、PLM-interact（可选） \\
\textbf{第 1 层 —— AL 循环}   & 集成 $\to$ BALD $\to$ ELISA $\to$ 重训 \\
\textbf{第 2 层 —— 因果性}    & $\Delta$tonB / $\Delta$exbB / $\Delta$exbD1 敲除 \\
\bottomrule
\end{tabular}
\end{center}

\vspace{0.6em}

结合（第 1 层）是感染的必要但不充分条件（Farquharson 2021）。第 2 层将受体特异性结合与防御系统贡献解耦。

## Pipeline 一览 —— Module 00 → 08

\begin{center}
\includegraphics[width=\textwidth]{docs/onboarding/figures/pipeline_flow.png}
\end{center}

\small 颜色编码：数据（蓝）· 特征（黄）· ML（绿）· 湿实验室（粉）。

## 数据约定

\begin{center}
\includegraphics[width=0.85\textwidth]{docs/onboarding/figures/data_contract.png}
\end{center}

\small
- `inputs/` = 指向上游 outputs 的指针（只读）。
- `processes/` = **唯一**写代码的地方。
- `outputs/` = 经典产物 + `MANIFEST.csv`（大文件 gitignore）。

完整规范：`INTERFACE.md`（路径、标识符、FASTA 头、MANIFEST schema）。

## Notebook 优先 + 双语注释（CLAUDE.md）

- 新代码先以 **Jupyter notebook** 形式编写（`<NN>_<short_name>.ipynb`）。
- 当 notebook (a) 能端到端运行 (b) 通过验证 cell → **冻结为 `.py`**。
- **双语注释**为强制规范 —— 简短注释 `# English / 中文` 同行。
- Module 07 例外 —— 生产流程，从第一天就写 `.py`。
- `nbstripout` 全仓库启用，notebook 输出不会污染 git diff。

## Module 00 —— 原始数据

- **777 个噬菌体 + 34 个细菌**完整基因组（NCBI）。
- 基因组二进制 **gitignore**（约 630 MB）；可通过下列命令重新下载：
  - `python 00_raw_data/processes/fetch_phages.py`
  - `python 00_raw_data/processes/fetch_bacteria.py`
- `MANIFEST.csv` 提供每个基因组的 SHA-256 + 大小，可重现。
- 标准参考对：`EU717894.1.fna`（phiL7）+ `GCF_000007145.1.fna`（Xcc）。

## Module 01 —— Ground Truth 相互作用矩阵

`01_data_ground_truth/outputs/interaction_matrix.csv` —— 前 4 行：

\scriptsize
\begin{tabular}{lllll}
\toprule
phage\_acc & host\_acc & label & source & notes \\
\midrule
EU717894.1 & GCF\_000007145.1 & 1 & literature\_curated & Hung 2003（受体）；Lee 2009（基因组） \\
NC\_054459.1 & NZ\_CP150073    & 1 & literature\_curated & X. oryzae pv. oryzae \\
ON758385.1 & — & 1 & literature\_curated & \emph{Xanthomonas} sp. \\
ON711490.1 & — & 1 & literature\_curated & Xcc XC114 \\
\bottomrule
\end{tabular}

\normalsize
- 2,236 对（315 阳性 + 1,921 阴性 + 1 ground truth）。
- 22/22 测试通过。

## Module 02 —— 注释

- **噬菌体：** PHANOTATE（McNair 2019）—— 动态规划处理重叠 ORF。
- **细菌：** pyrodigal 绑定 Prodigal（Hyatt 2010）。
- **绝不互换。** Prodigal 假设 ORF 不重叠 → 丢失约 10–15 % 噬菌体基因。
- 可选第二轮：pharokka（PHROG / VFDB / CARD）做功能分类。

\begin{tabular}{lll}
\toprule
基因组 & 工具 & 检出 ORF 数 \\
\midrule
phiL7（EU717894.1） & PHANOTATE & 80 \\
Xcc（GCF\_000007145.1） & pyrodigal & 4,344 \\
\bottomrule
\end{tabular}

## Module 03 —— RBP 识别（PhageRBPdetect，Boeckaerts 2022）

`03_rbp_identification/outputs/EU717894.1_rbp_candidates.csv` —— 80 个 ORF 中的 top 3：

\scriptsize
\begin{tabular}{llllll}
\toprule
orf\_id & length\_aa & hmm\_score & hmm\_match & combined\_score & rank \\
\midrule
EU717894.1\_orf\_00001 & \textbf{712} & 342.0 & unknown\_C54  & 1.000 & 1 \\
EU717894.1\_orf\_00021 & 918 & 235.1 & unknown\_C112 & 1.000 & 2 \\
EU717894.1\_orf\_00003 & 224 &  56.7 & unknown\_C294 & 1.000 & 3 \\
\bottomrule
\end{tabular}

\normalsize
- **rbp_01**（= orf_00001，712 aa，HMM score 342）—— Cycle 0 主要目标。
- 命中 Tail_spike_N HMM —— 正是 Lee 2009 用 BLAST 找不到的那种蛋白。
- 25+ 测试通过（本地需 `hmmpress`，否则 2 个会失败）。

## Module 04 —— 蛋白嵌入（ESM-2，Lin 2023）

- ESM-2 = 在约 65 M 蛋白序列上预训练的 masked language model（UniRef50/90）。
- 残基级嵌入 mean-pooling → 序列级向量。

\begin{tabular}{lll}
\toprule
变体 & 维度 & 运行位置 \\
\midrule
\texttt{esm2\_t6\_8M\_UR50D}    & 320  & 本地 CPU（概念验证，当前 outputs） \\
\texttt{esm2\_t33\_650M\_UR50D} & 1280 & Laguna A100 / L40S（生产目标） \\
\texttt{esm2\_t36\_3B\_UR50D}   & 2560 & 仅 Laguna（最终基准） \\
\bottomrule
\end{tabular}

\vspace{0.5em}
\small 可选层：PLM-interact（Liu 2025）—— 在人类 PPI 上 fine-tune 的 ESM-2。迁移到 mouse / fly / worm / yeast / *E. coli* PPI 后 AUPR +16–28 %。**尚未在噬菌体-细菌上测试** —— 我们可能是第一个。

## Module 05 —— 结构预测（Boltz-2 + AF3）

- **Boltz-2**（Passaro 2025）—— 预测复合体 3D 结构 + ipTM 界面信心分数。
- **AF3**（Abramson 2024）—— 更高质量静态结构；支持 trimer。

\begin{warnblock}
\textbf{关键警语 —— Boltz-2 的 affinity head 只支持小分子。} 训练数据来自 PubChem / ChEMBL / BindingDB。对\textbf{蛋白-蛋白对}（如 RBP × TonB）affinity head 输出 \texttt{NaN}。使用 \textbf{ipTM} 作为\emph{结构信心 proxy}，\textbf{不是}定量亲和力。
\end{warnblock}

# 第 3 部分 —— ML 核心 · `[干实验室] [PI]`

## Module 06 —— 深度集成用于预测不确定性

- **Lakshminarayanan 2017**（*NeurIPS*）—— 5 个 MLP，独立训练，输出预测均值 + 不确定性。
- 每个成员通过 Gaussian NLL loss 输出高斯分布 $(\mu_k, \sigma_k^2)$。

**为什么用 deep ensemble？**
- 校准性优于 MC Dropout（Ovadia 2019）。
- 可扩展到 ESM-2 的 1280 维输入（GP 做不到）。

\begin{noteblock}
\textbf{Greenman 2025 核查：} \emph{"没有单一最佳 UQ 方法"}（蛋白工程基准下）。集成通常 accuracy 最高但 calibration 最差。我们选集成的理由：(a) 可扩展性，(b) ALDE 的先例，(c) ECE 可以用 temperature scaling 修正。
\end{noteblock}

## ensemble.py —— 单成员预测（代码片段）

\scriptsize
```python
def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Returns (mean, sigma) — both shape (batch,).
    返回 (均值, 标准差)，各形状为 (batch,)。
    """
    h = self.net(x)
    mean = self.head_mean(h).squeeze(-1)
    # Clamp log_sigma to [-7, 7] for numerical stability
    # 将 log_sigma 截断到 [-7, 7] 防止数值爆炸
    log_sigma = self.head_log_sigma(h).squeeze(-1).clamp(-7.0, 7.0)
    sigma = torch.exp(log_sigma)
    return mean, sigma
```
\normalsize
\small 来源：`06_uncertainty_model/processes/ensemble.py:57–68`。

## Epistemic vs Aleatoric —— 关键分解

对 $K$ 成员的高斯集成，总预测方差分解为：

$$\sigma_{\text{total}}^2(x) \;=\; \underbrace{\mathrm{Var}_k[\mu_k(x)]}_{\text{epistemic（可约）}} \;+\; \underbrace{\mathbb{E}_k[\sigma_k^2(x)]}_{\text{aleatoric（噪声）}}$$

- **Epistemic** = \emph{模型}不知道的部分 → 数据越多越小 → **BALD 的目标**。
- **Aleatoric** = 测量 / 实验噪声 → 本质上不可约。

\vspace{0.4em}
`predictions.csv` 同时导出：`std`（总）和 `epistemic_std`（BALD 输入）。

## Module 06 输出 —— Cycle 0（合成数据）

`06_uncertainty_model/outputs/cycle_0/predictions.csv` —— 部分行：

\scriptsize
\begin{tabular}{llrrr}
\toprule
rbp\_id & receptor\_id & predicted\_score & std & epistemic\_std \\
\midrule
rbp\_00 & rec\_00 & 4.872 & 0.736 & 0.125 \\
rbp\_00 & rec\_01 & 4.959 & 0.755 & 0.114 \\
rbp\_01 & rec\_02 & \textbf{5.128} & 0.726 & \textbf{0.190} \\
rbp\_01 & rec\_03 & 5.094 & 0.715 & 0.049 \\
rbp\_07 & rec\_02 & 5.177 & 0.721 & \textbf{0.218} \\
\bottomrule
\end{tabular}

\normalsize
- 80 个 (RBP × receptor) 对。
- `epistemic_std` 在候选池中范围 0.04–0.22。
- `model_version = aa99d51_cycle_0`（git SHA + 循环）。

## Calibration 校准图 —— Cycle 0（合成数据）

\begin{center}
\includegraphics[width=0.62\textwidth]{06_uncertainty_model/outputs/cycle_0/calibration.png}
\end{center}

\small 可信区间覆盖率 vs 观察值。良好校准 = 点落在对角线上。合成数据下接近对角线 —— Cycle 0 ELISA 后将产生真正可信的版本。

## Module 07 —— BALD 直观示意

\begin{center}
\includegraphics[width=0.85\textwidth]{docs/onboarding/figures/bald_intuition.png}
\end{center}

\small 选集成成员分歧最大的点 —— 这些测量最快降低模型不确定性。

## BALD 数学（回归延伸）

原始版本（Houlsby 2011，GP 分类）：
$$\mathrm{BALD}(x) \;=\; I(y;\theta \mid x, D) \;=\; H[y \mid x, D] - \mathbb{E}_{\theta \sim p(\theta \mid D)}\!\bigl[H[y \mid x, \theta]\bigr]$$

对高斯深度集成 + 回归目标：
$$\mathrm{BALD}(x) \;\approx\; \mathrm{Var}_k[\mu_k(x)] \;\;\Longrightarrow\;\; \text{score} = \text{epistemic\_std}$$

\begin{noteblock}
\textbf{延伸说明。} Houlsby 2011 原始论文应用于高斯过程\emph{分类器}。把信息论目标延伸到深度集成\emph{回归}是合理的，但属于我们的延伸 —— \textbf{不是}原文直接引用。学术场合呈现时应注明这一点。
\end{noteblock}

## bald.py —— `bald_score()` + `select_batch()`

\scriptsize
```python
def bald_score(epistemic_std: np.ndarray) -> np.ndarray:
    """BALD score = epistemic_std（与 epistemic 方差单调等价）。"""
    if np.any(epistemic_std < 0):
        raise ValueError("epistemic_std must be non-negative; got negatives.")
    return epistemic_std.copy()

# select_batch() 内部：按 BALD 分数降序排名，取 top-K
ranked = pool.sort_values(bald_col, ascending=False).reset_index(drop=True)
ranked["bald_rank"] = ranked.index + 1
bald_picks = ranked.iloc[:n_bald].copy()
bald_picks["selection_reason"] = [f"bald_top_{i+1}" for i in range(len(bald_picks))]

# 随机对照：从未选中的 BALD top-K 之外抽样
remaining = ranked[~ranked.index.isin(set(bald_picks.index))]
random_picks = ranked.loc[_random.Random(seed).sample(list(remaining.index), k=n_random)]
random_picks["selection_reason"] = "random_control"
```
\normalsize
\small 来源：`07_acquisition_function/processes/bald.py:38–135`。

## Module 07 输出 —— Cycle 1 推荐

`07_acquisition_function/outputs/cycle_1/recommendations.csv`：

\scriptsize
\begin{tabular}{llrrrl}
\toprule
rbp\_id & receptor\_id & predicted\_score & epistemic\_std & bald\_rank & selection\_reason \\
\midrule
rbp\_07 & rec\_02 & 5.177 & \textbf{0.218} & 1 & bald\_top\_1 \\
rbp\_03 & rec\_01 & 4.810 & 0.197 & 2 & bald\_top\_2 \\
rbp\_05 & rec\_02 & 4.906 & 0.196 & 3 & bald\_top\_3 \\
rbp\_01 & rec\_02 & 5.128 & 0.190 & 4 & bald\_top\_4 \\
rbp\_03 & rec\_03 & 5.128 & 0.127 & 19 & random\_control \\
\bottomrule
\end{tabular}

\normalsize
- 4 个 BALD pick + 1 个随机对照（Hie 2022 对照臂设计）。
- 这是**合成占位数据**。真正的 Cycle 1 = rbp_01 variants × TonB（Cycle 0 ELISA 后）。

## ALDE 警语 —— Yang 2025 不验证 BALD

- **Yang 2025 ALDE**（*Nat Commun*）是最接近的公开工作 —— AL + DNN 集成 + 蛋白工程。
- 但是：ALDE 用的是 **Thompson 采样**作为采集函数，**one-hot encoding** 作为特征。不是 BALD，不是 ESM-2。
- 因此 ALDE 验证的是：*"AL + UQ + DNN 集成在蛋白工程中可行"*。
- BALD 本身仍需另立引用链（Houlsby 2011 + 我们的延伸）。

# 第 4 部分 —— 当前 Boltz-2 结果

## 我们跑了什么（job 59986）

\begin{tabular}{ll}
\toprule
Chain A & rbp\_01（712 aa 尾刺候选，EU717894.1） \\
Chain B & TonB（604 aa，GCF\_000007145.1） \\
GPU     & NVIDIA L40S（CARC Laguna） \\
运行时间 & 约 3 分钟 \\
\bottomrule
\end{tabular}

\vspace{0.6em}

历史注：早先一次（job 59949）用错了蛋白 —— 85 aa P25 ORF 被误标为 `rbp_01`。**job 59986 用的是正确的 712 aa rbp_01**（来自 `03_rbp_identification/outputs/EU717894.1_rbp_sequences.faa`）。

## 三个关键数字

\begin{tabular}{lrl}
\toprule
\textbf{指标} & \textbf{数值} & \textbf{解读} \\
\midrule
\texttt{interface\_ipTM} & \textbf{0.365} & 低 —— 模型对如何对接不确定。 \\
\texttt{chain\_A\_ptm}   & \textbf{0.808} & 高 —— rbp\_01 单体结构良好预测。 \\
\texttt{confidence\_score} & 0.683 & 中等的总体复合体质量。 \\
\bottomrule
\end{tabular}

\vspace{0.6em}

低 ipTM **不是失败** —— 它定义了实验。这种不确定性正是 ELISA + 主动学习循环要回答的问题。高 chain A ptm 说明 rbp_01 结构约束良好 → 是 variant design 的可靠基础。

## PAE 热图 —— 界面块

\begin{center}
\includegraphics[width=0.55\textwidth]{pae_heatmap.png}
\end{center}

\scriptsize 1316×1316 PAE 矩阵（rbp_01 残基 0–711，TonB 712–1315）。低 PAE（深色）= 相对置信度高。**非对角块**（行 712–1315 × 列 0–711）= rbp_01 × TonB 界面 —— 浅色 = 低信心。

通过 `np.load(.../pae_*.npz)['pae']` 加载（完整路径见下一张）。

## Boltz-2 输出位置

\scriptsize
\texttt{05\_structure\_prediction/outputs/boltz2/\\
EU717894.1\_rbp\_01\_\_GCF\_000007145.1\_tonB/}

目录内：

- `affinity.json` —— interface_ipTM、chain_A_ptm、confidence_score。
- `boltz_results_.../predictions/.../*.pdb` —— 完整原子坐标（PyMOL / ChimeraX 打开）。
- `pae_*.npz` —— PAE 矩阵。
- `plddt_*.npz` —— 逐残基 pLDDT（范围 0.30–0.98，均值 0.76）。

完整路径见 `docs/planning/PI_briefing_2026-05-11.md`（无省略号）。

\normalsize
**Boltz-2 提醒：** `predicted_dG = null`，因 affinity head 仅支持小分子。

# 第 5 部分 —— 48 小时循环 · `[湿实验室] [干实验室]`

## 循环结构

\begin{center}
\includegraphics[width=0.95\textwidth]{docs/onboarding/figures/cycle_48h.png}
\end{center}

\small 干实验室在 **48 小时**内把 ELISA → 推荐；湿实验室循环（SDM → 表达 → ELISA）需 **10–14 天**。

## 湿实验室 → 干实验室对接

每轮结束，湿实验室提交到 `08_cycle_data/outputs/cycle_<N>/`：

\begin{tabular}{lll}
\toprule
文件 & 必需列 & 说明 \\
\midrule
\texttt{elisa\_processed.csv} & variant\_id, receptor\_id, ec50\_nM, hill\_slope, r2, plate\_id, date & \\
\texttt{plaque\_results.csv}  & variant\_id, strain\_id, pfu\_per\_ml, date & WT 和 $\Delta$Receptor \\
\texttt{qc\_report.md}        & SDS-PAGE 图像路径、浓度、表达备注 & 每轮一份 \\
\bottomrule
\end{tabular}

\vspace{0.4em}
**重训最低门槛：** 每个 variant ≥3 个有效 EC50 且 R² > 0.9。失败 variant 标记 `ec50_nM = NaN` + `failed_reason` —— 模型可处理缺失数据。

## 干实验室 → 湿实验室对接（48 小时 SLA）

收到 ELISA 数据后：

\begin{tabular}{ll}
\toprule
文件 & 用途 \\
\midrule
\texttt{recommendations.csv}    & 4 BALD + 1 随机；主要任务列表 \\
\texttt{primer\_sequences.txt}  & NEB Q5 兼容的 SDM 引物 \\
\texttt{uncertainty\_bands.png} & 校准图：上一轮 预测 vs 实测 \\
\texttt{safe\_pick\_backup.csv} & Top-10 BALD，48 小时 SLA 未达时手动选择 \\
\texttt{run\_meta.json}         & git SHA + 时间戳 + 池统计（溯源） \\
\bottomrule
\end{tabular}

\vspace{0.4em}
湿实验室**不知道**哪一行是随机对照 —— 保留盲法基准，用于项目结束时的 AL vs random 事后比较。

## Cloning 执行

- **Cycle 0** —— 基因合成（IDT/Twist），4–6 个 variant，基于结构设计。约 2 周交付，~$150/variant。
- **Cycle 1+** —— 定点突变（NEB Q5），在已有构建体上做点突变。约 4 天，~$50/variant（3× 更便宜，3.5× 更快）。
- BALD 在小数据场景下倾向选点突变 → SDM 是自然的执行方法。

载体：pET-28a（Addgene 69864-3）。宿主：BL21(DE3)。诱导：0.5 mM IPTG，18 °C 过夜（利于可溶 trimer 装配）。

## ELISA 协议 —— 全菌结合

来自 Boeckaerts 2024 + Latka 2021：

1. 96 孔板包被热失活 *Xanthomonas*（10⁸ CFU/孔）。
2. 3 % BSA 封闭 1 h。
3. 系列稀释 His6-RBP（1 nM – 1 µM）。
4. HRP-anti-His6 检测，TMB 显色，读 OD~450~。
5. **4PL 拟合 → EC50** 是主动学习的目标变量。

每板对照：BSA-only（背景）；WT-RBP（板间归一化）；热变性 RBP（折叠特异性）；每个浓度 3 个技术重复。

## 受体敲除 —— pK18mobsacB（Schäfer 1994）

无标记敲除（自杀载体 + 蔗糖反选）：

1. 构建敲除质粒：靶基因上下游各约 500 bp + pK18mobsacB。
2. 电转入 Xcc 分离株；卡那霉素正选（单交换）。
3. 蔗糖反选（sacB 致死非解离体）；PCR 确认。

**靶基因（基于 Hung 2003）：**

- ΔtonB、ΔexbB、ΔexbD1 → 预期**阻断** phiL7 感染。
- ΔexbD2 → 预期**保留**感染（内置阴性对照）。

时间线：每个基因 4–6 周（可并行）。

## 验证层次

\begin{tabular}{lll}
\toprule
\textbf{层次} & \textbf{测什么} & \textbf{故事} \\
\midrule
1 & 仅 ELISA（WT 宿主） & "我们找到了结合更好的 variant。" \\
2 & + Plaque 测定（WT） & "结合 → 感染已确认。" \\
3 & + $\Delta$tonB / $\Delta$exbB / $\Delta$exbD1 & "受体特异性因果"（论文级）。 \\
\bottomrule
\end{tabular}

\vspace{0.6em}

**建议（PI 简报）：** 若 5/17 启动敲除，承诺第三层。ΔexbD2 阴性对照免费验证整个敲除系统。

## 质量门与失败模式

\scriptsize
\begin{tabular}{ll}
\toprule
\textbf{失败} & \textbf{处理} \\
\midrule
干实验室未达 48-h SLA   & 湿实验室使用 \texttt{safe\_pick\_backup.csv} \\
Variant 不溶解         & 标 NaN；尝试备用 truncation \\
ELISA R² $< 0.9$       & 重训时下调权重 \\
5 个 pick 全部 QC 失败  & 1 个专家选 + 2 个随机；用部分数据重训 \\
Calibration ECE $> 0.1$ & 下次 BALD 前做 temperature scaling \\
\bottomrule
\end{tabular}

# 第 6 部分 —— 复现与演示 · `[干实验室] [湿实验室]`

## 快速开始

\scriptsize
```bash
# 克隆 + 切分支
git clone https://github.com/A1ex-Ch3n/iGEM_Claremont_2026.git
cd iGEM_Claremont_2026
git checkout active-learning-pipeline

# 环境（一次性）
conda env create -f shared/env/environment.yml
conda activate igem2026
pip install nbstripout && nbstripout --install

# 本地开发最小基因组（共约 5 MB）
python 00_raw_data/processes/fetch_phages.py   --accession EU717894.1
python 00_raw_data/processes/fetch_phages.py   --accession NC_001604.1
python 00_raw_data/processes/fetch_bacteria.py --accession GCF_000007145.1

# 启动 JupyterLab
jupyter lab
```

\normalsize
完整逐模块清单：`GETTING_STARTED.md`。

## 每个模块的入口

\scriptsize
\begin{tabular}{ll}
\toprule
Module & 从这里打开 \\
\midrule
00 & \texttt{00\_raw\_data/processes/01\_verify\_dataset.ipynb} \\
01 & \texttt{01\_data\_ground\_truth/processes/01\_fetch\_reference\_genomes.ipynb} \\
02 & \texttt{02\_annotation/processes/01\_run\_phanotate.ipynb} \\
03 & \texttt{03\_rbp\_identification/processes/01\_run\_phagerbpdetect.ipynb} \\
04 & \texttt{04\_protein\_embedding/processes/01\_embed\_esm2.ipynb} \\
05 & \texttt{05\_structure\_prediction/processes/01\_run\_boltz2.ipynb}（生产用 Laguna） \\
06 & \texttt{06\_uncertainty\_model/processes/run\_cycle0.py} \\
07 & \texttt{07\_acquisition\_function/processes/run\_bald.py} \\
\bottomrule
\end{tabular}

## 输出文件位置

\scriptsize
\begin{tabular}{ll}
\toprule
文件 & 路径（相对仓库根） \\
\midrule
相互作用矩阵          & \texttt{01\_data\_ground\_truth/outputs/interaction\_matrix.csv} \\
RBP 候选              & \texttt{03\_rbp\_identification/outputs/EU717894.1\_rbp\_candidates.csv} \\
ESM-2 embedding       & \texttt{04\_protein\_embedding/outputs/embeddings\_esm2\_t6\_8M\_phiL7\_rbps.npz} \\
Boltz-2 PDB / ipTM    & \texttt{05\_structure\_prediction/outputs/boltz2/EU717894.1\_rbp\_01\_\_GCF\_000007145.1\_tonB/} \\
集成预测              & \texttt{06\_uncertainty\_model/outputs/cycle\_0/predictions.csv} \\
校准图                & \texttt{06\_uncertainty\_model/outputs/cycle\_0/calibration.png} \\
BALD 推荐             & \texttt{07\_acquisition\_function/outputs/cycle\_1/recommendations.csv} \\
PAE 热图              & \texttt{pae\_heatmap.png}（仓库根） \\
\bottomrule
\end{tabular}

## 现场演示计划（完整 runbook → `DEMO.md`）

1. **Module 03** —— `pytest 03_rbp_identification/processes/tests/ -v`（~15 s）。
2. **Module 04** —— 检查已缓存的 `.npz`（shape `(3, 320)`，seq_ids 对应 rbp_01/02/03）（<1 s）。
3. **Module 06** —— `python 06_uncertainty_model/processes/run_cycle0.py`（~3 s，已验证）。
4. **Module 07** —— `python 07_acquisition_function/processes/run_bald.py --cycle 1 --n_bald 4 --n_random 1`（<1 s）。
5. **Module 05（只读）** —— 打开 `pae_heatmap.png` 与 PDB（PyMOL）。

\vspace{0.4em}
\small 每步都有"说什么"和"失败怎么办"，见 `DEMO.md`。

## Laguna HPC —— 何时切换到 GPU

\scriptsize
\begin{tabular}{lll}
\toprule
任务 & 本地？ & 为什么 Laguna？ \\
\midrule
ESM-2 8M（3 个 RBP） & 是 & CPU 秒级 \\
ESM-2 650M / 3B（777 phages）   & 否 & GPU 显存 + 时间 \\
Boltz-2 蛋白-蛋白       & 否 & 约 15 分钟/对，仅 A100/L40S \\
AF3 batch                       & 否 & GPU + 权重审批 \\
集成训练 + BALD                 & 是 & CPU 秒级 \\
\bottomrule
\end{tabular}

\normalsize
环境配方 + SLURM 模板：`LAGUNA.md`（CUDA 12.4，torch 2.5.1+cu121）。

## 测试 —— 当前通过率

\begin{tabular}{ll}
\toprule
Module & 测试 \\
\midrule
00 原始数据             & 15+ 通过（3 个预期失败 —— GCF/T7 不在 seed list） \\
01 Ground Truth         & 22/22 \\
02 注释                 & 26/26 \\
03 RBP ID               & 25+ 通过（2 个预期失败 —— HMM 需要本地 \texttt{hmmpress}） \\
04 Embedding            & 17/17 \\
05 结构                 & 28/28（1 个预期跳过 —— 需要完成的 GPU 跑） \\
06 不确定性             & 9/9 \\
07 BALD                 & 18/18 \\
\bottomrule
\end{tabular}

# 第 7 部分 —— 风险、决策与请求 · `[PI]`

## 待 PI 决策事项

\scriptsize
\begin{tabular}{lll}
\toprule
\textbf{决策} & \textbf{截止} & \textbf{默认} \\
\midrule
敲除系统 pK18mobsacB vs CRISPRi   & 5/17 & \textbf{pK18mobsacB}（Hung 2003 在 Xcc 上已用过） \\
AlphaFold 3 权重 —— 谁去申请？     & 本周 & Alex 或 PI；学术机构 email 优先 \\
验证层次承诺                       & 6/1 前 & \textbf{第三层}（若 5/17 启动敲除） \\
LA 县噬菌体富集来源                & 6/1 前 & 待定 —— 需 PI 实验室渠道 \\
论文投稿意向                       & Cycle 0 前 & 与 iGEM wiki 并行投 \emph{Bioinformatics} / \emph{NAR Genomics} \\
\bottomrule
\end{tabular}

## 关键风险与缓解

\scriptsize
\begin{tabular}{lll}
\toprule
\textbf{风险} & \textbf{可能性 / 影响} & \textbf{缓解} \\
\midrule
菌株分离失败                  & 低 / 高 & 双源（十字花科 + 柑橘）；2 轮采样 \\
噬菌体分离无裂解性            & 中 / 高 & 污水 / 农业径流富集；Phage Directory 备用 \\
RBP 表达不溶解                & 中 / 中 & 先用 T7 gp17 阳性对照；GCN4 trimer tag 备选 \\
ELISA 动态范围不足            & 中 / 高 & Cycle 0 前 2 周专门优化 ELISA \\
AL 不如随机                   & 低 / 高（可报告） & 标准对照臂；诚实报告 \\
干实验室未达 48-h SLA          & 中 / 中 & PI 预批的 \texttt{safe\_pick\_backup.csv} \\
受体敲除失败                  & 中 / 中 & 多靶基因并行；CRISPRi 备用 \\
\bottomrule
\end{tabular}

## "AL 不如随机"的样子

- 在 Cycle 2 后：BALD 轨迹的测试集 R² 与随机选择**统计等价**。
- 我们承诺诚实报告：
  - 用与正面结果同等严谨度记录该负面结果。
  - 对领域是有用贡献 —— *Hie 2024 行，ALDE 行，phage-RBP（暂时）不行*。
- 每轮归档的具体指标：
  - 测试集 R²（在 held-out variant 上）。
  - Calibration ECE。
  - 每次实验的信息增益（后验 vs 先验 KL 散度）。

## 时间线

\scriptsize
\begin{tabular}{ll}
\toprule
日期 & 里程碑 \\
\midrule
2026-05-17 & 湿实验室启动：十字花科采样、pK18mobsacB 构建、Cycle 0 合成订单 \\
2026-06-01 & Cycle 0 启动：ELISA 优化 + 首批结合测量 \\
2026-06-14 & Cycle 1：集成重训，首批 BALD picks，SDM \\
2026-06-28 & Cycle 2：第二轮，受体敲除完成 \\
2026-07-12 & 分析与拓展：训练模型应用于土壤源噬菌体基因组 \\
2026-09-12 & Wiki 冻结准备 \\
2026-10-21 & Wiki 冻结 \\
2026-11-13 & Grand Jamboree \\
\bottomrule
\end{tabular}

## iGEM 交付物清单

- **Wiki：** Engineering DBTL · Modeling · Hardware/Software · Human Practices。
- **Composite Part：** ≥4–6 个 RBP-His6 构建体（来自自分离噬菌体）。
- **软件仓库：** 开源闭环 AL pipeline（即本仓库）。
- **宣传视频：** 一个闭环周期的演示。
- **社区贡献：** 自分离 *Xanthomonas* + 噬菌体 → iGEM Registry；测序基因组 → NCBI。

## 给团队的开放性问题

1. 湿实验室 —— LA 县的污水 / 农业径流渠道？
2. 湿实验室 —— 在投 variant 库前是否先做 T7 gp17 阳性表达对照？
3. 干实验室 —— 谁负责 ESM-2 3B Laguna 任务（1280 → 2560 维，最终基准）？
4. 全员 —— 论文投稿与 iGEM wiki 冻结的时间关系？
5. PI —— 若分离株鉴定结果是检疫级 pathovar，APHIS 许可的姿态？

# 第 8 部分 —— 参考文献与附录

## 各模块关键论文 —— 紧凑表

\scriptsize
\begin{tabular}{lll}
\toprule
Module & 论文 & 角色 \\
\midrule
00 / 01 & da Silva 2002 \emph{Nature}；Lee 2009 \emph{AEM}；Hung 2003 \emph{BBRC} & 参考基因组 + 受体 \\
02 & McNair 2019 \emph{Bioinformatics}（PHANOTATE）；Hyatt 2010（Prodigal）；Bouras 2023（pharokka） & 基因预测 \\
03 & Boeckaerts 2022 \emph{Viruses}（PhageRBPdetect） & RBP HMM + XGBoost \\
04 & Lin 2023 \emph{Science}（ESM-2）；Liu 2025 \emph{Nat Commun}（PLM-interact） & embedding + PPI 先验 \\
05 & Abramson 2024 \emph{Nature}（AF3）；Passaro 2025（Boltz-2） & 结构 + ipTM \\
06 & Lakshminarayanan 2017 \emph{NeurIPS}；Greenman 2025 \emph{PLoS Comput Biol} & 深度集成 + UQ 核查 \\
07 & Houlsby 2011（BALD）；Yang 2025 \emph{Nat Commun}（ALDE）；Hie 2024 \emph{Nat Biotechnol} & 采集 \\
08 & Schäfer 1994（pK18mobsacB）；Gibson 2009；Latka 2021 \emph{mBio} & 湿实验室方法 \\
\bottomrule
\end{tabular}

\vspace{0.4em}
\small 完整带注释阅读指南：`docs/reference/papers.md`（19 篇）。

## 五项文献核查修正（2026 年 5 月）

\scriptsize
\begin{tabular}{ll}
\toprule
\textbf{原说法} & \textbf{修正后} \\
\midrule
"TonB-ExbB-ExbD1-ExbD2 全必需"（Wang 2003）            & TonB-ExbB-ExbD1 必需；ExbD2 不是（Hung 2003 BBRC） \\
"Boltz-2 affinity head 用于 RBP × 受体"                & affinity head 仅小分子 $\to$ \texttt{NaN}；用 ipTM \\
"Greenman 2025 \emph{NAR Genomics}；深度集成最佳 UQ"   & Greenman 2025 \emph{PLoS Comput Biol}；"无单一最佳 UQ" \\
"Hie 2024 用 ESM-2，约 50 个抗体 variant"             & Hie 2024 用 ESM-1b/1v，每个抗体 ~20 个 variant \\
"Lee 2009 命名了 tail spike"                           & Lee 2009 \emph{主动搜索并找不到} OP1 ORF25 同源物 \\
\bottomrule
\end{tabular}

\vspace{0.5em}
\small 完整核查：`docs/reference/paper_reading_notes.md`。

## 导航参考

\small
\begin{tabular}{ll}
\toprule
想知道 …… & 看这里 \\
\midrule
每个模块的产出   & \texttt{<module>/README.md} \\
端到端状态快照   & \texttt{docs/planning/PI\_briefing\_2026-05-11.md} \\
项目计划（EN）   & \texttt{docs/planning/iGEM\_2026\_Project\_Plan.md} \\
项目计划（ZH）   & \texttt{docs/planning/iGEM\_2026\_项目大纲\_中文版.md} \\
论文列表（注释） & \texttt{docs/reference/papers.md} \\
文献修正记录     & \texttt{docs/reference/paper\_reading\_notes.md} \\
数据约定         & \texttt{INTERFACE.md} \\
HPC 设置         & \texttt{LAGUNA.md} \\
湿实验室 SOP     & \texttt{docs/protocols/}（5 份 Benchling PDF） \\
\bottomrule
\end{tabular}

## 术语表 —— 快速查询

\scriptsize
\begin{tabular}{ll}
\toprule
术语 & 含义 \\
\midrule
RBP        & Receptor-binding protein —— 噬菌体的"钥匙" \\
HMM        & Hidden Markov Model —— 序列 profile 方法 \\
ESM-2      & 进化尺度蛋白语言模型（Lin 2023） \\
ipTM       & Interface predicted TM-score（0–1，结构信心） \\
ELISA      & 酶联免疫吸附测定（结合读数） \\
EC50       & 半数有效浓度（4PL 拟合） \\
SDM        & 定点突变（NEB Q5 试剂盒） \\
BALD       & Bayesian Active Learning by Disagreement（Houlsby 2011） \\
Epistemic  & 可约的模型不确定性（= BALD 目标） \\
Aleatoric  & 不可约的测量噪声 \\
\bottomrule
\end{tabular}

\vspace{0.4em}
\small 完整术语表：`docs/reference/glossary.md`。

## 谢谢 / Thank you

\vspace{0.6em}
\Large
**欢迎提问、评论、质疑。**

\normalsize
\vspace{0.5em}

- 干实验室：Alex Chen — `CChen29@cmc.edu`
- PI：Prof. J. Cesar Ignacio-Espinoza
- 学术顾问：Prof. Ran Libeskind-Hadas
- 湿实验室负责人：Sarah、Olivia、Weitao、Carol
- 贡献者：Ryan、Leah

\vspace{0.6em}
\small 接下来现场演示（15 分钟），见 `docs/onboarding/DEMO.md`。
