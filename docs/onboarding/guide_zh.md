# 入职指南 —— iGEM Claremont 2026

**针对黄单胞菌生物防治的噬菌体主动学习工程**

> 配套 `slides_zh.pdf`。如果错过了演示，或要找到某张 slide 背后的具体文件路径、命令、引用，请读本文。

| | |
|---|---|
| **核心工程师** | Alex Chen |
| **湿实验室负责人** | Sarah、Olivia、Weitao、Carol |
| **贡献者** | Ryan、Leah |
| **PI** | Prof. J. Cesar Ignacio-Espinoza |
| **学术顾问** | Prof. Ran Libeskind-Hadas |
| **分支** | `active-learning-pipeline`（main 是空的 —— 不要在 main 上工作） |
| **当前状态** | Module 00–07 已全部完成并通过测试。Module 08 从首次 ELISA 交付（约 6 月 1 日）开始。 |

---

## 1. 项目一览

我们构建了一个**闭环主动学习 pipeline**，将噬菌体 RBP（受体结合蛋白）× 细菌受体相互作用的机器学习模型与迭代式湿实验室验证结合起来。每轮 ELISA 之后，模型重训，BALD 采集函数按 epistemic 不确定性对所有未测试 variant 打分排序，湿实验室测下一批 4–5 个。整个循环的设计目标是：每一次昂贵 ELISA 测量都产生尽可能多的信息 —— 直击噬菌体-宿主预测的核心痛点：**数据稀缺**。

参考干实验室骨架：
- 宿主：*Xanthomonas campestris* pv. *campestris*（Xcc）ATCC 33913 —— NCBI `GCF_000007145.1`。
- 噬菌体：phiL7 —— NCBI `EU717894`。
- 受体系统：TonB-ExbB-ExbD1 必需，ExbD2 **不**必需（Hung 2003，*BBRC* 302:878–884，PMID 12646254）。

湿实验室从加州十字花科作物自分离 *Xanthomonas* + 噬菌体（PI 协商 2026-05-07），绕过数月的 USDA APHIS PPQ-526 许可证。

iGEM 目标赛道：**Best Agriculture Project · Best Model · Best Composite Part**。

---

## 2. 科学背景 —— 为什么做这个

### *Xanthomonas* 与宿主范围问题

*Xanthomonas* 属是 >400 种植物的病原（Ryan 2011，*Nat Rev Microbiol*）。加州相关致病变种包括 Xcc（十字花科黑腐病，商业蔬菜生产中普遍存在）、*X. citri* subsp. *citri*（柑橘溃疡病）、*X. perforans* / *X. euvesicatoria*（番茄细菌斑病）。当前防治依赖铜制剂，但耐铜性已广泛出现（Aiello 2019，*Plant Dis*）。

噬菌体生物防治是有吸引力的替代方案 —— 宿主特异性、自扩增、环境可降解。田间试验在番茄细菌斑病上达到铜制剂效果（Iriarte 2018），在甘蓝 Xcc 上成功抑制（Holtappels 2022）。但同样的特异性让部署变难：每个菌株都需要对应的噬菌体。

预测哪种噬菌体感染哪种细菌仍不可靠。当前最强模型 PhageHostLearn（Boeckaerts 2024，*Nat Commun*）**菌种内 AUC 约 0.82**，**跨菌属降至 0.67–0.70**（PAML benchmark，Mutalik 2025）。根本瓶颈是数据稀缺 —— 定量（噬菌体，宿主）结合标签生成慢且贵。

### phiL7 + Xcc —— 已知什么

phiL7 是一个 44,080 bp 的裂解性 Siphoviridae，Lee 等 2009 年（*AEM* 75:7828）做了基因组特征描述。它通过 **TonB-ExbB-ExbD1** 外膜转运复合体感染 *X. campestris*。Hung 等 2003 年通过 Tn5 突变实验确认：TonB、ExbB、ExbD1 对噬菌体侵入必需；尽管 ExbD2 与它们共转录，但**不是必需的**（CH620 株，ΔexbD2，保留完全敏感性）。

最后这一点比表面看到的更有用。构建 ΔexbD2 与 ΔtonB / ΔexbB / ΔexbD1 平行，我们就有了**免费的阴性对照**：如果敲除系统正常工作，ΔexbD2 应当仍允许感染。这是对整个敲除流程的内置验证。

### rbp_01 —— Lee 2009 与我们的 HMM 重新发现

Lee 2009 主动搜索 phiL7 中 OP1 ORF25 的 tail-fiber 同源物，并**通过序列相似性找不到**：

> "We were unable to identify a homolog of the OP1 tail fiber protein (ORF25) thought to be involved in host range determination…"

p20（1105 aa，tail protein III）被建议与宿主范围相关，但没有任何蛋白被命名为"tail spike"。我们的 **rbp_01（712 aa）**来自 PhageRBPdetect 的 Tail_spike_N HMM（Boeckaerts 2022，*Viruses*）—— 一种结构 profile 方法，能识别序列已分歧到 BLAST 无法找到的蛋白。这与 Lee 2009 不矛盾；HMM 与 BLAST 灵敏度不同，rbp_01 是用更敏感工具的互补发现。

### 为什么主动学习

主动学习（AL）是对数据稀缺的数学回应（Lindley 1956；Settles 2009）。AL 系统不是被动地在已有数据上训练，而是用当前不确定性来选择哪个实验若执行后能最大程度地减少这种不确定性。BALD（Houlsby 2011）等采集函数将其形式化为信息增益最大化。

近期示范：
- **Hie 2024**（*Nat Biotechnol* 42:275）—— 用 ESM-1b/1v 语言模型 likelihood + AL 做抗体亲和成熟，每个抗体约 20 个 variant，最高 160 倍亲和力改善。
- **Yang 2025 ALDE**（*Nat Commun*）—— DNN 集成 + Thompson 采样（注意：**不是** BALD）做酶定向进化；产率 12% → 93%，2 轮 wet lab。

**没有人将其应用到噬菌体 RBP × 细菌受体。** 这正是我们的方法学贡献。

---

## 3. Pipeline 架构

### 三层

```
第 0 层 —— 先验             Boltz-2 ipTM、ESM-2 嵌入、PLM-interact（可选）
                            （任何 wet lab 数据之前的信息先验）
                                       │
                                       ▼
第 1 层 —— AL 循环          集成 → BALD → recommendations.csv
                            → 湿实验室 ELISA → 重训 → 下一批推荐
                                       │
                                       ▼
第 2 层 —— 因果性           ΔtonB / ΔexbB / ΔexbD1 敲除 + ΔexbD2 阴性对照
                            （解耦受体特异性结合与防御系统贡献）
```

结合（第 1 层）是有效感染的必要但**不充分**条件（Farquharson 2021，T4 × *E. coli*：RBP 结合 85 % 的菌株，但只在 11 % 上形成噬斑）。第 2 层量化模型结合信号有多少能转化为感染。

### 模块图（00 → 08）

![Pipeline 流程](figures/pipeline_flow.png)

| Module | 工具 | 状态 | 关键输出 |
|--------|------|------|---------|
| 00 原始数据 | NCBI fetch | ✅ | 777 噬菌体 + 34 细菌；`MANIFEST.csv`（二进制 gitignore） |
| 01 Ground Truth | curated | ✅ | `interaction_matrix.csv`：2,236 phage–host 对 |
| 02 注释 | PHANOTATE / pyrodigal | ✅ | phiL7：80 ORF；Xcc：4,344 ORF |
| 03 RBP ID | PhageRBPdetect HMM + XGBoost | ✅ | `EU717894.1_rbp_candidates.csv`；rbp_01 = 712 aa |
| 04 Embedding | ESM-2（本地 8M；Laguna 上 650M / 3B） | ✅ | `embeddings_esm2_*.npz` |
| 05 结构 | Boltz-2（Laguna）；AF3 待批 | ✅ | `affinity.json`（`ipTM = 0.365`）、PDB、PAE 矩阵 |
| 06 不确定性 | 5-MLP 深度集成 | ✅ | `predictions.csv` 含 `std` 和 **`epistemic_std`** |
| 07 BALD | Var_k[μ_k] 采集 | ✅ | `recommendations.csv` + `safe_pick_backup.csv` |
| 08 循环数据 | 湿实验室 ELISA | ⏳ ~6/1 | `elisa_processed.csv` |

### 数据约定

![数据约定](figures/data_contract.png)

每个模块都有同样的三个子目录：

- **`inputs/`** —— 指向上游 `outputs/` 的只读指针（或外部 seed，如 NCBI accession 列表）。**绝不在此写生成数据。**
- **`processes/`** —— 唯一写代码的地方。脚本与 notebook 从 `inputs/` 读取，写入自己的 `outputs/`。
- **`outputs/`** —— 下游模块消费的经典产物。大文件树 gitignore —— 一个含 SHA-256 + size 的 `MANIFEST.csv` 让数据集可重现。

完整规范 —— 包括标识符格式、FASTA 头约定、MANIFEST schema —— 见仓库根的 `INTERFACE.md`。

### Notebook 优先工作流 + 双语注释

来自 `CLAUDE.md`：

- 新代码先以 **Jupyter notebook** 形式编写（`<NN>_<short_name>.ipynb`）用于探索。
- 当 notebook (a) 端到端运行 (b) 通过验证 cell，**冻结为 `.py`** 在同一文件夹；把 notebook 改名为 `<NN>_<short_name>__frozen.ipynb`，加上指向 `.py` 的指针。
- **双语注释**强制 —— 简短内联：`# English / 中文` 同行；长注释：在 markdown cell 用单独段落。
- Module 07 例外 —— 生产编排代码，从第一天就写 `.py`。
- `nbstripout` 在仓库范围启用，notebook 输出不会污染 git diff。每次新 clone 后跑一次 `pip install nbstripout && nbstripout --install`。

---

## 4. ML 核心深度解析

### Module 06 —— 深度集成（Lakshminarayanan 2017，*NeurIPS*）

5 个独立初始化的 MLP（3 层隐藏，ReLU + dropout，两个输出头分别预测 `mean` 和 `log_sigma`）在 ESM-2 嵌入上训练，预测 ELISA 结合分数。每个成员输出 Gaussian $(\mu_k, \sigma_k^2)$；训练损失是 Gaussian 负对数似然。

```python
# 来自 06_uncertainty_model/processes/ensemble.py:57–68
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

总预测方差分解：

$$\sigma^2_{\text{total}}(x) \;=\; \underbrace{\mathrm{Var}_k[\mu_k(x)]}_{\text{epistemic（可约）}} + \underbrace{\mathbb{E}_k[\sigma_k^2(x)]}_{\text{aleatoric（噪声）}}$$

`predictions.csv` 同时导出两者：`std`（总）和 `epistemic_std`（BALD 输入）。

**为什么用深度集成？** 比 MC Dropout 校准更好（Ovadia 2019）；可扩展到 ESM-2 的 1280 维输入（GP 做不到）。Greenman 2025 核查（*PLoS Comput Biol* 21:e1012639）很重要：

> *"没有单一最佳 UQ 方法"*（蛋白工程基准下）。集成通常 accuracy 最高但 calibration 最差。

我们选集成是因为可扩展性 + ALDE 的先例，如果 ECE 漂到 0.1 以上就用 temperature scaling 修补。

#### Cycle 0 输出预览

`06_uncertainty_model/outputs/cycle_0/predictions.csv`（部分行）：

| rbp_id | receptor_id | predicted_score | std | epistemic_std |
|--------|-------------|-----------------|-----|---------------|
| rbp_00 | rec_00      | 4.872           | 0.736 | 0.125 |
| rbp_01 | rec_02      | **5.128**       | 0.726 | **0.190** |
| rbp_07 | rec_02      | 5.177           | 0.721 | **0.218** |

80 个 pair 评分；`epistemic_std` 在池中范围 0.04–0.22。合成数据；Cycle 0 后将被真实 ELISA 数据替换。

校准图（`outputs/cycle_0/calibration.png`）在合成数据上接近对角线 —— Cycle 0 ELISA 后才会出现有意义的版本。

### Module 07 —— BALD 采集函数（Houlsby 2011 延伸）

原始 Houlsby 2011 公式（高斯过程分类器）：

$$\mathrm{BALD}(x) \;=\; I(y;\theta \mid x, D) \;=\; H[y \mid x, D] - \mathbb{E}_{\theta \sim p(\theta \mid D)}\![H[y \mid x, \theta]]$$

对高斯深度集成回归目标：

$$\mathrm{BALD}(x) \;\approx\; \mathrm{Var}_k[\mu_k(x)] \;\;\Longrightarrow\;\; \text{score} = \text{epistemic\_std}$$

**延伸说明。** Houlsby 2011 原本是分类 + GP。把信息论目标延伸到深度集成回归很自然，但这是我们的延伸 —— *不是*原文直接引用。学术评审时应明确说明。

```python
# 来自 07_acquisition_function/processes/bald.py:38–135（精简）
def bald_score(epistemic_std: np.ndarray) -> np.ndarray:
    """BALD score = epistemic_std（与 epistemic 方差单调等价）。"""
    if np.any(epistemic_std < 0):
        raise ValueError("epistemic_std must be non-negative; got negatives.")
    return epistemic_std.copy()

# select_batch() 内部：
ranked = pool.sort_values(bald_col, ascending=False).reset_index(drop=True)
ranked["bald_rank"] = ranked.index + 1
bald_picks = ranked.iloc[:n_bald].copy()
bald_picks["selection_reason"] = [f"bald_top_{i+1}" for i in range(len(bald_picks))]

# 随机对照：从未选中的 BALD top-K 之外抽样
remaining = ranked[~ranked.index.isin(set(bald_picks.index))]
random_picks = ranked.loc[_random.Random(seed).sample(list(remaining.index), k=n_random)]
random_picks["selection_reason"] = "random_control"
```

#### Cycle 1 输出预览

`07_acquisition_function/outputs/cycle_1/recommendations.csv`：

| rbp_id | receptor_id | predicted_score | epistemic_std | bald_rank | selection_reason |
|--------|-------------|-----------------|---------------|-----------|------------------|
| rbp_07 | rec_02      | 5.177           | **0.218**     | 1         | bald_top_1       |
| rbp_03 | rec_01      | 4.810           | 0.197         | 2         | bald_top_2       |
| rbp_05 | rec_02      | 4.906           | 0.196         | 3         | bald_top_3       |
| rbp_01 | rec_02      | 5.128           | 0.190         | 4         | bald_top_4       |
| rbp_03 | rec_03      | 5.128           | 0.127         | 19        | random_control   |

4 个 BALD pick + 1 个随机对照（Hie 2022 模式）。这是合成占位数据 —— 真正的 Cycle 1 = rbp_01 variants × TonB（Cycle 0 ELISA 后）。

### ALDE 警语 —— Yang 2025 不验证 BALD

Yang 2025 ALDE 是离我们项目最近的已发表工作但**不验证 BALD**。它使用 Thompson 采样和 one-hot encoding（不是 ESM-2）。因此 ALDE 支持更宽泛的 claim *"AL + UQ + DNN 集成在蛋白工程中可行"*；对 BALD 本身，引用链是 Houlsby 2011 → 我们对集成的延伸。

### 一个数据约定的代码示例

`07_acquisition_function/processes/run_bald.py` 从 `inputs/` 指针读取并写入 `outputs/`：

```python
# 输入（只读指针）
predictions_path = (
    REPO_ROOT
    / "06_uncertainty_model" / "outputs" / f"cycle_{cycle-1}"
    / "predictions.csv"
)
measured_pairs = load_measured_pairs(args.measured_csv) if args.measured_csv else set()

df = pd.read_csv(predictions_path)
df["bald_score"] = bald_score(df["epistemic_std"].values)

recommendations, random_replay, safe_pick_backup = select_batch(
    df, n_bald=args.n_bald, n_random=args.n_random,
    measured_pairs=measured_pairs, seed=args.seed,
)

# 输出（经典产物 + 溯源）
out_dir = REPO_ROOT / "07_acquisition_function" / "outputs" / f"cycle_{cycle}"
out_dir.mkdir(parents=True, exist_ok=True)
recommendations.to_csv(out_dir / "recommendations.csv", index=False)
safe_pick_backup.to_csv(out_dir / "safe_pick_backup.csv", index=False)
random_replay.to_csv(out_dir / "random_replay.csv", index=False)
(out_dir / "run_meta.json").write_text(json.dumps(run_meta, indent=2))
```

无硬编码绝对路径。`REPO_ROOT` 通过 `Path(__file__).resolve().parents[2]` 锚定。

---

## 5. 当前 Boltz-2 结果

我们在 Laguna 上跑了 Boltz-2（Passaro 2025，job 59986，NVIDIA L40S）：

- **Chain A：** rbp_01，712 aa，来自 `03_rbp_identification/outputs/EU717894.1_rbp_sequences.faa`
- **Chain B：** TonB，604 aa

`affinity.json` 三个关键数字：

| 指标 | 数值 | 解读 |
|------|------|------|
| `interface_ipTM` | **0.365** | 低。模型对 rbp_01 与 TonB 如何对接不确定。对没有 PDB 模板的全新系统来说是预期结果；ELISA 来解决。 |
| `chain_A_ptm` | **0.808** | 高。rbp_01 单体结构约束良好 —— 是 variant design 的可靠基础。 |
| `confidence_score` | 0.683 | 整体复合体质量，中等。 |
| `predicted_dG` | `null` | **affinity head 仅小分子。** |

低 ipTM 不是模型失败 —— 它定义了实验。这种结构不确定性正是 ELISA + 主动学习循环要解答的问题。

**Boltz-2 affinity head 警语（Passaro 2025）。** Affinity head 在小分子 × 蛋白结合数据（PubChem、ChEMBL、BindingDB）上训练。对蛋白-蛋白对（RBP × TonB），输出 `NaN`。始终用 **ipTM** 作为结构信心 proxy —— 绝不要将 Boltz-2 描述为对蛋白-蛋白提供"zero-shot affinity prior"。

### PAE 热图

![rbp_01 × TonB PAE](../../pae_heatmap.png)

由 `pae_*.npz`（1316×1316 float32）生成。残基 0–711 = Chain A（rbp_01）；712–1315 = Chain B（TonB）。**非对角块**（行 712–1315 × 列 0–711）是界面区域 —— 低值（深蓝）表示模型对该残基对的相对位置有信心。块内的浅带对应界面不确定性 —— ELISA 最能信息性地约束的部分。

如何从原始 `.npz` 重新生成热图：见 `docs/planning/PI_briefing_2026-05-11.md` 底部的代码片段。

### 文件位置

所有 Boltz-2 输出位于：

```
05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/
├── affinity.json                                          ← 汇总分数
└── boltz_results_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/
    └── predictions/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/
        ├── boltz_input_..._model_0.pdb                    ← 3D 原子结构
        ├── pae_..._model_0.npz                            ← PAE 矩阵
        ├── plddt_..._model_0.npz                          ← 逐残基 pLDDT
        └── confidence_..._model_0.json                    ← 逐链置信度
```

用 PyMOL（`pymol.org`）或 ChimeraX（`rbvi.ucsf.edu/chimerax/`）打开 PDB 检查预测复合体。PI_briefing.md 含完整未截断路径。

---

## 6. 48 小时循环 —— 干湿对接

![48 小时循环](figures/cycle_48h.png)

### 湿实验室 → 干实验室（每轮结束）

湿实验室提交至 `08_cycle_data/outputs/cycle_<N>/`：

| 文件 | 必需列 |
|------|--------|
| `elisa_processed.csv` | `variant_id, receptor_id, ec50_nM, hill_slope, r2, plate_id, date` |
| `plaque_results.csv`  | `variant_id, strain_id, pfu_per_ml, plaque_morphology, date`（WT 与 ΔReceptor） |
| `qc_report.md`        | SDS-PAGE 图像路径、Bradford 浓度、表达备注 |

**重训最低门槛：** 每个 variant ≥3 个有效 EC50 且 R² > 0.9。失败 variant 标 `ec50_nM = NaN` + `failed_reason` —— 模型可处理缺失数据。

### 干实验室 → 湿实验室（48 小时 SLA）

收到 ELISA 数据后 48 小时内，干实验室在 `07_acquisition_function/outputs/cycle_<N+1>/` 生成：

| 文件 | 用途 |
|------|------|
| `recommendations.csv` | 4 BALD + 1 随机；主要任务列表 |
| `primer_sequences.txt` | NEB Q5 兼容的 SDM 引物（自动生成） |
| `uncertainty_bands.png` | 校准图：上一轮 预测 vs 实测 |
| `safe_pick_backup.csv` | Top-10 BALD；仅当 48 小时 SLA 未达时使用 |
| `run_meta.json` | 溯源：git SHA、时间戳、池大小、分数统计 |

**盲法对照：** 湿实验室**不知道**哪一行是随机对照 —— recommendations CSV 在交给湿实验室时不会单独标 `random_control` 列。这保留了 AL vs random 的事后比较。

### Cloning 执行

- **Cycle 0** —— 基因合成（IDT/Twist），4–6 个 variant 由结构设计 + 专家选定。约 2 周交付，~$150/variant。
- **Cycle 1+** —— 定点突变（NEB Q5）在已有构建体上。约 4 天，~$50/variant —— 比基因合成便宜 3 倍、快 3.5 倍。

载体：pET-28a（Addgene 69864-3）。宿主：BL21(DE3)。诱导：0.5 mM IPTG，18 °C 过夜（利于可溶 trimer 装配，Studier & Moffatt 1986）。

### ELISA 协议（全菌结合，Boeckaerts 2024 + Latka 2021）

1. 96 孔板包被热失活 *Xanthomonas*（10⁸ CFU/孔）。
2. 3% BSA 封闭 1 h。
3. 系列稀释 His6-RBP（1 nM – 1 µM）。
4. HRP-anti-His6 检测，TMB 显色，读 OD₄₅₀。
5. **4PL 拟合 → EC50** 是主动学习的目标变量。

每板对照：BSA-only（背景）、WT-RBP 固定浓度（板间归一化）、热变性 RBP（折叠特异性结合）、每个浓度 3 个技术重复。

### 受体敲除 —— pK18mobsacB（Schäfer 1994）

无标记敲除（自杀载体 + 蔗糖反选）：

1. 构建敲除质粒：靶基因上下游各约 500 bp + pK18mobsacB（Addgene 87097）。
2. 电转入 Xcc 分离株；卡那霉素正选（单交换）。
3. 蔗糖反选；sacB 致死非解离体。PCR + 测序确认。

靶基因与预期表型（Hung 2003）：

| 菌株 | 预期 phiL7 结果 |
|------|----------------|
| WT | 敏感 |
| ΔtonB | 抗性 |
| ΔexbB | 抗性 |
| ΔexbD1 | 抗性 |
| **ΔexbD2** | **仍敏感**（阴性对照） |

时间线：每个基因 4–6 周，可并行。

### 验证层次

| 层次 | 测什么 | 故事 |
|------|--------|------|
| 1 | 仅 ELISA（WT 宿主） | "我们找到了结合更好的 variant。" |
| 2 | + Plaque 测定（WT） | "结合 → 感染已确认。" |
| 3 | + Δreceptor 组合 | **"受体特异性因果。"**（论文级） |

建议（PI 简报 2026-05-11）：若 5/17 启动敲除，承诺第三层。

### 失败模式与质量门

| 失败 | 处理 |
|------|------|
| 干实验室未达 48-h SLA | 湿实验室使用 `safe_pick_backup.csv`（PI 预批） |
| Variant 不溶解 | 标 NaN；尝试备用 truncation |
| ELISA R² < 0.9 | 重训时下调权重 |
| 5 个 BALD pick 全部 QC 失败 | 1 个专家选 + 2 个随机；用部分数据重训 |
| Calibration ECE > 0.1 | 下次 BALD 前做 temperature scaling |

---

## 7. 如何复现 —— 环境、命令、测试

### 快速开始（共约 10 分钟）

```bash
# 1. 克隆 + 切分支（active-learning-pipeline，不是 main）
git clone https://github.com/A1ex-Ch3n/iGEM_Claremont_2026.git
cd iGEM_Claremont_2026
git checkout active-learning-pipeline

# 2. Conda 环境（一次性，约 5 分钟）
conda env create -f shared/env/environment.yml
conda activate igem2026
pip install nbstripout && nbstripout --install  # 每次新 clone 一次

# 3. 本地开发最小基因组（约 5 MB）
python 00_raw_data/processes/fetch_phages.py   --accession EU717894.1
python 00_raw_data/processes/fetch_phages.py   --accession NC_001604.1     # T7，被测试使用
python 00_raw_data/processes/fetch_bacteria.py --accession GCF_000007145.1

# 4. 一次性 PhageRBPdetect 数据 + HMM press（Module 03）
bash 03_rbp_identification/inputs/setup_inputs.sh

# 5. 启动 JupyterLab
jupyter lab
```

完整 777 噬菌体 + 34 细菌数据集（约 630 MB）已 gitignore。只在 Laguna 上批处理任务前获取。

### 每模块入口

| Module | 从这里打开 |
|--------|-----------|
| 00 | `00_raw_data/processes/01_verify_dataset.ipynb` |
| 01 | `01_data_ground_truth/processes/01_fetch_reference_genomes.ipynb` |
| 02 | `02_annotation/processes/01_run_phanotate.ipynb`（噬菌体）+ `02_run_prodigal.ipynb`（宿主） |
| 03 | `03_rbp_identification/processes/01_run_phagerbpdetect.ipynb` |
| 04 | `04_protein_embedding/processes/01_embed_esm2.ipynb` |
| 05 | `05_structure_prediction/processes/01_run_boltz2.ipynb`（生产用 Laguna） |
| 06 | `06_uncertainty_model/processes/run_cycle0.py`（CPU 上 3 秒；2026-05-17 验证） |
| 07 | `07_acquisition_function/processes/run_bald.py`（<1 秒） |

### 运行测试

```bash
pytest 00_raw_data/processes/tests/             -v    # 15+ 通过（3 个预期失败）
pytest 01_data_ground_truth/processes/tests/    -v    # 22/22
pytest 02_annotation/processes/tests/           -v    # 26/26
pytest 03_rbp_identification/processes/tests/   -v    # 25+ 通过（2 个预期失败 —— hmmpress）
pytest 04_protein_embedding/processes/tests/    -v    # 17/17
pytest 05_structure_prediction/processes/tests/ -v    # 28/28（1 个预期跳过 —— GPU run）
pytest 06_uncertainty_model/processes/tests/    -v    # 9/9
pytest 07_acquisition_function/processes/tests/ -v    # 18/18
```

### Laguna HPC —— 何时

| 任务 | 本地？ | 为什么 Laguna？ |
|------|--------|----------------|
| ESM-2 8M（3 个 RBP） | 是 | CPU 秒级 |
| ESM-2 650M / 3B（777 noisemic phages） | 否 | GPU 显存 + 时间 |
| Boltz-2 蛋白-蛋白 | 否 | ~15 分钟/对，仅 A100/L40S |
| AF3 批处理 | 否 | GPU + 权重审批（Google 表单） |
| 集成重训 + BALD | 是 | CPU 秒级 |

完整 Laguna 配方（CUDA 12.4，torch 2.5.1+cu121，conda `boltz2` env，NVIDIA L40S，OnDemand 门户）见 `LAGUNA.md`。boltz 安装会升级 torch 至 cu130 —— `LAGUNA.md` 中有一行修复将其固定回 cu121。

### 重新构建本入职 deck

```bash
# 图表
python docs/onboarding/figures/make_figures.py

# 英文 slides
pandoc docs/onboarding/slides_en.md \
  -t beamer --pdf-engine=xelatex --slide-level=2 \
  -o docs/onboarding/slides_en.pdf

# 中文 slides
pandoc docs/onboarding/slides_zh.md \
  -t beamer --pdf-engine=xelatex --slide-level=2 \
  -o docs/onboarding/slides_zh.pdf
```

一次性 TeX 依赖：`tlmgr install beamer booktabs mdframed xecjk ctex pgfplots fontaxes needspace zref`。Mac CJK 字体：`Hiragino Sans GB`（系统自带）。

---

## 8. 约定参考

### 标识符格式（来自 `INTERFACE.md`）

| 概念 | 格式 | 示例 |
|------|------|------|
| 噬菌体 NCBI accession | `<base>.<version>` | `EU717894.1` |
| 细菌组装 | `GCF_*` / `GCA_*` | `GCF_000007145.1` |
| ORF id（每个基因组） | `<acc>_orf_<5-digit>` | `EU717894.1_orf_00031` |
| RBP 候选 id | `<acc>_rbp_<2-digit>` | `EU717894.1_rbp_01` |
| 受体 id | `<host_acc>_<gene_name>` | `GCF_000007145.1_tonB` |
| Variant id | `<parent_rbp>_<change>_<idx>` | `EU717894.1_rbp_01_trunc_03` |

### 路径锚定

- `.ipynb` 内：`REPO_ROOT = Path.cwd().resolve().parents[1]`（notebook 在 `<module>/processes/`）。
- `.py` 内：`REPO_ROOT = Path(__file__).resolve().parents[2]`。
- **绝不**硬编码绝对路径。**绝不**用 `~`。**绝不**假设 CWD。

### 循环版本

每个模型 checkpoint、predictions CSV、ELISA 数据集都打 `cycle_<N>` 标签并放在对应目录。溯源（`run_meta.json`）含跑时 git SHA。MLflow 跟踪 run（计划 Cycle 0 起启用）。

### 工具分工（Module 02）

PHANOTATE 用于**噬菌体** ORF 调用；Prodigal / pyrodigal 用于**细菌**宿主。**绝不互换。** Prodigal 假设 ORF 不重叠，会丢失约 10–15 % 噬菌体基因；PHANOTATE 用动态规划处理重叠 ORF。

---

## 9. 参考文献 —— 各模块紧凑

| Module | 论文 | 角色 |
|--------|------|------|
| 00 / 01 | da Silva 2002 *Nature*；Lee 2009 *AEM*；Hung 2003 *BBRC* | 参考基因组 + 受体 |
| 02 | McNair 2019 *Bioinformatics*（PHANOTATE）；Hyatt 2010（Prodigal）；Bouras 2023（pharokka） | 基因预测 |
| 03 | Boeckaerts 2022 *Viruses*（PhageRBPdetect） | RBP HMM + XGBoost |
| 04 | Lin 2023 *Science*（ESM-2）；Liu 2025 *Nat Commun*（PLM-interact） | embedding + PPI 先验 |
| 05 | Abramson 2024 *Nature*（AF3）；Passaro 2025（Boltz-2） | 结构 + ipTM |
| 06 | Lakshminarayanan 2017 *NeurIPS*；Greenman 2025 *PLoS Comput Biol* | 深度集成 + UQ 核查 |
| 07 | Houlsby 2011（BALD）；Yang 2025 *Nat Commun*（ALDE）；Hie 2024 *Nat Biotechnol* | 采集 |
| 08 | Schäfer 1994（pK18mobsacB）；Gibson 2009；Latka 2021 *mBio* | 湿实验室方法 |

带 🔴/🟡/⚪ 优先级标签的完整阅读指南：`docs/reference/papers.md`（19 篇）。

### 五项文献核查修正（2026 年 5 月）

这些是 Alex 在 2026-05-11 完整阅读 19 篇核心论文期间发现并修正的。不要再引入这些错误。

1. **ExbD2 不是必需** —— Hung 2003（真实来源）显示 ΔexbD2 保留 phiL7 敏感性。先前 "Wang 2003 *Mol Microbiol*" 引用是幻觉。
2. **Boltz-2 affinity head 仅小分子** —— 对蛋白-蛋白对输出 NaN。用 ipTM 作为结构信心 proxy，不作为亲和先验。
3. **Greenman 2025 期刊是 *PLoS Comput Biol***（不是 *NAR Genomics*）；结论是"无单一最佳 UQ 方法"，不是"深度集成胜出"。
4. **Hie 2024 用 ESM-1b/1v**（不是 ESM-2）；每个抗体看 **~20 个 variant**（不是 ~50）；机制是语言模型 likelihood 过滤，不是 BALD 闭环。
5. **Lee 2009 主动搜索并明确找不到** phiL7 中 OP1 ORF25 的同源物。我们的 rbp_01 是基于 HMM 的互补发现 —— *不*与 Lee 2009 矛盾。HMM 能识别结构相似但序列分歧到 BLAST 找不到的蛋白。

完整核查：`docs/reference/paper_reading_notes.md`。

---

## 10. 术语表 + 进一步阅读

### 快速查询

| 术语 | 含义 |
|------|------|
| RBP | Receptor-binding protein —— 噬菌体进入宿主细胞的"钥匙" |
| HMM | Hidden Markov Model —— 序列 profile 方法（Boeckaerts 2022） |
| ESM-2 | Evolutionary-Scale Modeling v2 —— 蛋白语言模型（Lin 2023） |
| ipTM | Interface predicted TM-score（0–1，结构信心） |
| ELISA | 酶联免疫吸附测定（结合读数） |
| EC50 | 半数有效浓度（4PL 拟合） |
| SDM | 定点突变（NEB Q5） |
| BALD | Bayesian Active Learning by Disagreement（Houlsby 2011） |
| epistemic | 可约的模型不确定性（= BALD 目标） |
| aleatoric | 不可约的测量噪声 |
| pK18mobsacB | 无标记敲除的自杀载体（Schäfer 1994） |

### 进一步阅读

- `docs/reference/glossary.md` —— 完整术语表。
- `docs/planning/iGEM_2026_Project_Plan.md` —— 给 PI 看的英文项目计划。
- `docs/planning/iGEM_2026_项目大纲_中文版.md` —— 中文版，含更深的干实验室模块机制 + 6 层数据稀缺策略。
- `docs/planning/PI_briefing_2026-05-11.md` —— 双语状态简报，含每个输出路径与 5 月 7–12 日工作记录。
- `docs/protocols/` —— 湿实验室 Benchling SOP（培养、转化、噬斑、感染曲线、裂解液扩增）。
- `INTERFACE.md` —— 模块间锁定的数据约定。
- `LAGUNA.md` —— HPC 设置、SLURM 模板、CUDA 注意事项。

---

## 附录 —— 复现一行命令

```bash
# 从一次新 clone —— 端到端跑 Module 06 和 07（约 5 秒）
conda activate igem2026
python 06_uncertainty_model/processes/run_cycle0.py
python 07_acquisition_function/processes/run_bald.py --cycle 1 --n_bald 4 --n_random 1
# 查看：
head 06_uncertainty_model/outputs/cycle_0/predictions.csv
head 07_acquisition_function/outputs/cycle_1/recommendations.csv
```

这就是主动学习闭环的两行命令。难的部分从 6 月 1 日开始。
