# 项目方向调整说明（团队同步文档）

## 📌 项目一句话梳理

**iGEM Claremont 2026** 做的是用噬菌体生物防治农业植物病原 *Xanthomonas*——目标拿 **Best Agriculture Project** + **Best Model**。原 proposal 的 dry lab 是「用 6 个手工生化特征 + ML 二元分类预测哪个 phage 能感染 Xanthomonas」；现在 pivot 成「**主动学习 (Active Learning / Bayesian Optimal Experimental Design) 闭环系统**——dry lab 模型告诉 wet lab 下一个该做哪个 RBP variant 实验，每轮迭代收敛」。Wet lab 端从 proposal 的「噬菌体 vs 整菌」拉到 motif-level「**RBP-His6 蛋白纯化 + ELISA 定量结合 + 受体敲除因果验证**」。研究系统选定 **Xanthomonas campestris pv. campestris (Xcc) ATCC 33913 + 噬菌体 phiL7 (NCBI EU717894)**——后者的受体已知是 **TonB-dependent system**，让 Step 8 受体敲除有明确目标。时间线极紧：今天 5/7，wet lab 启动 5/17，cycle 跑到 7/12，wiki freeze 10/21，jamboree 11/13。**最大 blocker 是 USDA APHIS plant pest permit（平均 127 天）**——必须 48 小时内跟 PI 确认实验室既有库存或 permit 状态，决定走 Path A（PI stock）/ B（本地分离）/ C（向作者要 phage）哪条路。

---

> **目的**：与全体队友同步 dry lab 方向的调整、与原 proposal 的差异、以及未来两个月的关键时间节点。
>
> **状态**：方向初步定案，等待团队 + PI 最终确认
>
> **写作日期**：2026-05-07
>
> **作者**：Alex（dry lab core）

---

## TL;DR（30 秒看懂）

1. 我们把 dry lab 的核心思路从「**用 6 个手工特征加权 → 训一个 ML 模型预测感染**」升级成「**让 ML 模型告诉 wet lab 下一个该做哪个实验，每轮迭代逐步收敛**」。
2. **wet lab 计划基本不动**（约 80% 沿用 proposal）。变的是 dry lab 的数学基础和迭代方式。
3. 这个改动让我们从「我们做了一个 ML pipeline」升级成「我们做了一个数学最优实验设计系统」，对 **Best Model 奖项**是质的提升。
4. 时间线很紧但可行：今天 5/7 → 5/17 wet lab 启动，dry lab 必须在这 10 天内把基础架构跑起来。

---

## 一、新方向到底是什么？

### 1.1 核心想法

我们要建立一个 **闭环主动学习系统（Active Learning / Bayesian Optimal Experimental Design, BOED）**：

```
模型当前的"信念" → 数学上算出"下一个最有信息量的实验是什么"
        ↓
   wet lab 去做这个实验
        ↓
   实验结果回灌进模型
        ↓
   模型更新 → 推荐下一个实验
        ↓
        （循环 2-3 轮）
```

### 1.2 学术基础（这不是我们拍脑袋想的）

这个领域有约 70 年的数学积累可以站：

- **Lindley 1956** — Bayesian experimental design 的开山之作
- **Settles 2009** — Active learning 的经典 survey
- **Frazier 2018** — Bayesian optimization 教学论文
- **蛋白质工程脉络**：Romero & Arnold 2009、Wittmann et al. Cell Systems 2021、**Hie et al. Cell 2022**（用蛋白质语言模型做高效抗体进化 — 跟我们要做的最接近）

### 1.3 为什么这个方向真的"新"

**关键问题**：phage-host interaction 领域有人做过 closed-loop active learning 吗？

**据我所知没有**：
- Yehl 2019 是 one-shot DMS（不是迭代）
- Mutalik 2020 是 CRISPRi screen（不是模型驱动）
- GenoPHI 2026 是 retrospective genome screens（回头看，不是前瞻设计）

**这个 niche 是空的**。而且 iGEM 评审会很买帐这个 framing —— 因为「数据稀缺」是 phage ML 公认的 #1 痛点，而我们的答案不是「收更多数据」，是「**让每个数据点的信息量最大化**」。

### 1.4 这个方向同时优雅地解决一个老问题

之前一直困扰我们的一个事实：**蛋白质结合 ≠ 真实感染**（Farquharson 2021 显示 T4 RBP 跟 85% 菌株结合，但只在 11% 上形成噬斑）。

在主动学习框架下，这变成 **feature 而不是 bug**：

> 模型对"结合机制 vs 防御机制谁主导"有不确定性 → acquisition function 自动推荐能区分这两者的实验（例如同一菌株的 WT + ΔReceptor + ΔDefense 三角对照）

最终 paper 的 claim 变成：「我们的系统在 N 轮 cycle 之后，把结合贡献跟防御贡献的相对权重学出来了」。这个故事比单纯的 binding atlas 强一个档次。

---

## 二、跟原 proposal 比，改动多大？

### ✅ Proposal 里已经写进去、跟新方向完全一致的部分

| Proposal 步骤 | 跟新方向的关系 |
|---|---|
| **Step 5**：His-tagged RBP + 免疫印迹对 *Xanthomonas* 细胞 | 这就是我们要的定量结合实验（换个 readout 而已） |
| **Step 7**：RBP 结合 loop 的定向突变 | 这就是 motif-level variants |
| **Step 8**：宿主受体敲除 + Western blot | 这就是因果验证（KD validation）！ |
| **Step 10**：「模型通过实验反馈精修」 | 这就是闭环的雏形 |
| **Highlights → Model**：「闭环 DBTL，digital twin」 | **PI 已经 buy in 这个大方向了！** |
| **Highlights → Composite Part**：RBP-His6 表达载体库 | 完全对得上 |

**结论**：proposal 的 wet lab 计划 ~80% 不用动。变的主要是 dry lab 的数学基础。

### ❌ 真正要改的 3 件事

#### 改动 1：把 Step 3（六因子特征加权）整个 retire

我们原本用 6 个手工特征：GC content、pI/acidity、protein length、GRAVY、sequence similarity、CAI。

**坦白说，这 6 个特征学术 support 都很薄弱：**

| Feature | 文献 support 强度 | 备注 |
|---|---|---|
| GC content | 弱 | Kirchberger 2020 显示 phage GC ≈ host GC，但不是 strain-level predictor |
| pI/acidity | 很弱 | 影响蛋白质表达/溶解度，不是 specificity feature |
| Protein length | 很弱 | Tail fiber 长度跨家族差异大，但同家族内没 predictive power |
| GRAVY | 很弱 | 没有 phage-host literature 支持 |
| Sequence similarity | 中等 | BLAST-based 方法有，但被 PLM 全面压制 |
| CAI | 弱 | Lucks et al. 显示 phage codon ≈ host，但 strain-level 没 resolution |

**对照 PhageHostLearn (Boeckaerts 2024 *Nature Communications*)**：完全不用这 6 个 features，纯 PLM embedding，AUC 0.82。

→ **这一步换成 ESM-2 embedding + 不确定性感知回归器（uncertainty-aware regressor），文献 support 直接从「弱」跳到「SOTA-aligned」**。

#### 改动 2：Step 10 从 passive retrain → active selection

Proposal 写的是「实验数据回来 → recalibrate 模型 → 改 weights」。这是 **被动**学习。

新方向是：「模型 → **告诉 wet lab 下一个该做哪个实验** → 拿到数据 → retrain」。这是 **主动**学习。

**这个改动对 "Best Model" 奖是质变**：从「我们有个 ML pipeline」升级成「我们有个数学最优的实验设计系统」。

#### 改动 3：每一步补 literature anchor

Proposal 现在每步都写得像 statement of intent，没 cite paper。后面我们每个 sub-step 都该有 ≥1 篇 anchor paper 撑着。这是 wiki / presentation video 阶段的硬功夫。

---

## 三、时间线（必须严格执行）

```
今天 = 2026-05-07
   ↓
   ↓ 10 天 dry lab 冲刺窗口（baseline + acquisition function + cycle infrastructure）
   ↓ 同步：gene synthesis 必须在 5/17 前下单（synthesis 本身要 1-2 周）
   ↓
2026-05-17  Wet Lab 启动（Phase 3）
   ↓
   ↓ ~4 周 cloning + expression（cycle 0 = seed 4-6 个 variants）
   ↓
2026-06-14  Dose-Response Loops 开始（cycle 1：模型推荐 next 4-6 个）
   ↓
   ↓ ~4 周 ELISA + cycle iteration（cycle 2：模型推荐 next 4-6 个）
   ↓
2026-07-12  Capsule/Gel UV Tests（wet lab 核心冲刺终点，~1.5 个月）
   ↓
2026-08-09  Phase 4 Metagenomics + writeup
   ↓
2026-10-21  Wiki Freeze
   ↓
2026-11-13  Grand Jamboree
```

**核心节奏**：8 周 wet lab 内必须跑完至少 **2 轮 active learning cycle**（cycle 1 + cycle 2），才能在 paper / wiki 上 claim「我们的系统 work」。

---

## 四、wet lab 工作量怎么排

3 位生物背景队友 + PI 指导。基于 cloning / expression / ELISA 三条 parallel pipeline 估算：

| 阶段 | 时间 |
|------|------|
| Cloning + small-scale expression test | 1-2 周/批（gene synthesis 平行化） |
| Large-scale purification + ELISA | 1 周/批 |
| 3 人 parallel + 经验 OK | **12-18 个 RBP variants** 是 8 周内的 realistic 上限 |

→ 不是原 brief 里写的 25-50 个。我们要 **scope 收紧**，把质量做扎实。

---

## 五、必须诚实点出的 3 个 critical risks

### Risk 1：必须有 control arm 证明「主动学习真的比随机选好」

iGEM 评审第一个会问：「你怎么知道你的 acquisition function 真的比随便挑实验好？」

**对策**：
- Retrospectively simulate（用最终所有数据回头测 random vs 我们的 acquisition function）
- 或每轮 cycle 推荐里塞 1 个 random variant 当对照

**Alex 的看法**：
这是**最致命**的 risk，因为它是回头补救不了的方法学问题。但好消息是修复 wet lab 成本几乎为零。我的具体提议：

1. **Cycle 0-N 全程保存三件东西**：(a) 每轮模型完整 state，(b) acquisition function 推荐的 top-K，(c) 平行计算「随机挑会选到什么」的 random-K 列表
2. **每轮 wet lab batch 偷塞 1 个 random pick**：4-6 个 actual ELISA targets 里 1 个是随机选的，给我们 in-distribution 直接对照
3. **最终 retrospective replay**：用所有 final data 重跑「如果一开始就用 random 选，learning curve 会长怎样」
4. 这是 **Hie et al. 2022 的标准做法**，不是我们发明的，reviewer 看了会服

→ 必须从 cycle 0 就 commit 这件事，不能等到 wiki freeze 前才补。

### Risk 2：seed data 太少 → 不确定性估计全错 → 推荐变 garbage

4-6 个 seed variants 对高斯过程 (GP) 来说是极小的数据量。

**对策**：Phase 0 用 public phage RBP-receptor data 做 informative prior，模型不是从零学起。

**Alex 的看法**：
其实可以叠 **3 层 prior** 来撑过冷启动：

| Prior 层 | 来源                                                                         | 成本            |
| ------- | -------------------------------------------------------------------------- | ------------- |
| Layer 1 | PLM-interact 预训练（人类 PPI 数据迁移）                                              | 0 — 直接下载权重    |
| Layer 2 | 我们 00_raw_data 已经有的 777 个 phage 基因组的 RBP，做 unsupervised contrastive        | 1-2 天 dry lab |
| Layer 3 | Boltz-2 对所有候选 variants 跑 zero-shot affinity prediction → 当 synthetic prior | 一晚 GPU        |

**真正的 danger 不是 prior 不够，是过度自信**。Cycle 0-1 model 实质上还在猜。我的承诺：每轮 cycle 在团队 sync meeting 里都会画 uncertainty bands 给大家看，明确讲「这一轮模型还在乱猜，wet lab 的人凭直觉觉得不对就 push back」。Trust 会随 cycle 累积，但不能假装第一轮就准。

### Risk 3：wet/dry sync 卡住 critical path

如果 cycle 1 蛋白质表达失败 1-2 个 variants（很常见），整个时间表就塌。

**对策**：每轮推荐 6-8 个 variants 但只 commit 4-6 个进 ELISA，留 redundancy buffer。

**Alex 的看法**：
这本质上是 project management risk 不是技术 risk。修复方法是 process：

1. **硬性 SLA**：Alex 承诺收到 cycle N data 后 **48 小时内**出 cycle N+1 推荐。错过 SLA → 自动 trigger 备案
2. **备案：Safe pick backup list**：开始就由 PI + 生物 teammates 凭直觉手挑 ~20 个候选 variants 排成优先级队列。如果某轮 dry lab 来不及，wet lab 直接从 backup list 顺位往下做，**不等**
3. **Buffer day**：每轮 wet lab schedule 内建 1 天 buffer，可以拿来跑 replicate ELISA、additional controls、prep 下一批 cell pellets — 永远不会让 wet lab idle

核心原则：**wet lab 永远不会因为 dry lab 而停下来**。这是说服 PI 接受闭环模式的关键。

---

## 六、接下来 10 天 dry lab 必须交付的 4 件事

1. **Phase 0 baseline 模型** — ESM-2 embedding + uncertainty-aware regressor（先用 deep ensemble，1-2 天能动；GP 等 cycle 2 再升级）
2. **Acquisition function** — 不重新发明，直接用 BALD 或简化版 EIG
3. **Public data 上的 transfer learning prior** — PhageRBPdetect 5k 序列 + PDB phage complex 结构
4. **Cycle infrastructure** — 一个 pipeline 接收新 ELISA 数据 → 自动 retrain → 输出下一轮推荐

---

## 七、团队确认结果与下一步

下面三件事原本是需要团队确认的 open questions，目前已有初步答案。把当前状态记录下来，team 可以看到后再 pushback。

### 议题 1：*Xanthomonas* 株 + phage scaffold 选哪个？

**当前状态**：完全没选，Benchling 上现有提到的 Xcc / X. citri 都只是 placeholder。

> ⚠️ **本文件 v1.0 / v1.1 都推荐错了**：之前 Xp10 (X. oryzae host 错配)、CP2 (filamentous Inoviridae) 都不对。这一版基于 May 2026 的文献调研重写。

#### **Alex 的推荐（v1.2）**：**Xcc (ATCC 33913) + phage phiL7 (NCBI EU717894)**

| 项目 | 选择 | 理由 |
|---|---|---|
| **Host strain** | *Xanthomonas campestris* pv. *campestris*  ATCC 33913 (= DSM 3586, NCPPB 528, LMG 568, type strain) | BSL-1，模式菌株，全基因组已测 (NCBI: AE008922)，2002 年起广泛用于学术 plant pathology research，Brassica 黑腐病的 causal agent — 跟「Best Agriculture」叙事完美契合（加州 brassica 产业大） |
| **Phage** | phiL7 (NCBI: **EU717894**, 44 kb) | T7-like Siphoviridae **lytic phage**（对的，不是 podovirus；icosahedral capsid + long noncontractile tail），感染 Xcc 的 reference phage，2009 年由 Academia Sinica (Taiwan) 完整测序与功能注释 |
| **Receptor** | **TonB + ExbB + ExbD1**（ExbD2 co-transcribed 但非必需）operon | 这是 **本研究的金子** — phiL7 的 receptor 已经由 Tn5 mutagenesis 鉴定（Hung et al. 2003, *BBRC* 302:878–884, PMID 12646254），TonB 是 inner membrane protein，TonB-dependent receptor system 是 well-defined target，**Step 8 受体敲除直接知道要 KO 哪个基因**（Hung 2003 確認 TonB、ExbB、ExbD1 三個 essential；exbD2 突變保留 phiL7 敏感性——是負對照） |

**为什么这个组合是 winning combo**：

1. **Receptor 已知** → 整个 pipeline 不用瞎猜。Layer 1 ELISA target 用 TonB-dependent OMP。Layer 2 KO 直接 targeted to tonB / exbB / exbD1（exbD2 作為負對照）
2. **Strain + phage 都是 reference standards** → reviewer 不会 question host range 之类基础问题
3. **Phage 是 lytic, T7-like** → 跟蛋白质工程 literature（T7 gp17 family）有 cross-reference 可用
4. **ATCC 33913 + EU717894 都在 NCBI/PubMed 全公开** → dry lab Phase 0 可以立即开始，不用等 wet lab

#### ⚠️ **CRITICAL BLOCKER：USDA APHIS PPQ-526 permit**

**ATCC 明确说明 33913 需要 USDA-APHIS Permit to Move Live Plant Pests 才能 ship**。

USDA APHIS 数据：
- 平均处理时间 **127 天**（约 4.2 个月）
- APHIS 建议 **40 周（280 天）前**就要送件，特别是需要 facility inspection 的话
- 我们今天送件 = 最快 **2026 年 9 月** 才能拿到 = **远超过 wet lab 核心期 (5/17 - 7/12)**

**这意味着 ATCC 直接订购路线不可行**。必须走下面三条 alternative path 之一。

#### 三条 alternative path

**Path A（最理想，必须先确认）：PI 实验室现有库存**

PI Cesar Ignacio-Espinoza 的研究就是环境微生物 + phage，他实验室 **极有可能已经有 Xanthomonas + 相关 phage 的 glycerol stocks**。如果是，permit 不是问题（既有材料不需要新 permit）。

**MUST ASK PI 的 3 个问题**：
1. 实验室有没有任何 Xanthomonas 株（任何 pathovar）的 glycerol stock？
2. 实验室有没有 standing PPQ-526 permit（如果之前做过 plant pathogen 工作很可能已有）？
3. 有没有合作 lab 可以 transfer Xanthomonas 株（interstate movement 通常比 ATCC import 快）？

**Path B（fallback 1）：Local isolation from California 农产品/土壤**

环境采样 + 自行分离通常 **不需要 PPQ-526 permit**（permit 主要管跨州 import / commercial transfer，不管 in-state 环境采样）。

具体做法（proposal Step 6 已经规划了）：
- 从 farmer's market / 加州农场买带黑腐病症状的 brassica（卷心菜、花椰菜、芥菜）
- Standard enrichment protocol：磨碎组织 + Xcc-selective medium (NYG + cycloheximide) → 30°C → 单菌落分离
- 用 16S rRNA + MALDI-TOF 确认 ID
- 时间：~1-2 周
- **同时 isolate phage**：把环境样品过 0.22 µm 滤膜 + Xcc lawn → plaque

这条 path 的 bonus：**直接产出新颖 phage isolates** 给 iGEM Part Registration（"Best Composite Part" 加分）

**Path C（fallback 2）：直接联系 phiL7 原作者**

phiL7 paper 的通讯作者是 NCHU（台灣國立中興大學）Tseng lab 的研究员（Hung et al. 2003, BBRC + Lee et al. 2009, AEM）。学术 phage 通常作者会乐意 share。

问题：
- 国际 transfer 仍需 USDA permit（material 从台湾入境美国）
- Email + 等回复可能 1-3 周
- 时间风险高

**推荐的 path**：**Path A 优先确认 → 失败立刻 Path B（local isolation）**。Path C 当 long-shot backup。

#### dry lab 端可以立即开始（不依赖 wet lab）

**关键 insight**：dry lab 的 motif-level RBP 工作 **完全不需要实体材料**。我们需要的是：

1. ✅ **phiL7 完整基因组** — 已在 NCBI EU717894，立即可下载
2. ✅ **Xcc ATCC 33913 完整基因组** — 已在 NCBI AE008922，立即可下载（含 tonB, exbB, exbD1, exbD2 序列）
3. ✅ **phiL7 RBP 候选基因** — 用 PHANOTATE 在 EU717894 上识别 tail spike / tail fiber gene
4. ✅ **AlphaFold 预测 RBP 3D 结构** — 没有 PDB structure for phiL7 RBP，必须靠 AF
5. ✅ **TonB-dependent receptor 3D 结构** — Xcc 自身 TonB / ExbB / ExbD1 也可以 AF 预测（ExbD2 作為負對照一並預測）

→ Alex 1 周 dry lab sprint 可以全部基于上述公开数据完成，**5/17 wet lab 启动延迟不影响 dry lab 进度**。

### 议题 2：Benchling 集成方案

**当前状态**：Sarah 已经把现有 Benchling protocols 直接下载成 PDF 放进 `docs/`。这一步**优先级降级了**——dry lab 不需要 SDK 拉数据，可以直接读 `docs/` 里的 PDF / 我等下整理出来的 markdown。

**已下载的 5 份 protocols**（详见下面第八节的 gap 分析）：

```
docs/
├── Xanthomonas Cultivation Protocol · Benchling.pdf
├── Xanthomonas Transformation Protocol · Benchling.pdf
├── Phage Propagation Plaque Assay (Small Scale) Protocol · Benchling.pdf
├── Phage Infection Curves · Benchling.pdf
└── Phage Propagation Liquid Amplification (Large-Scale Lysate) Protocol · Benchling.pdf
```

**Benchling SDK 是否还需要？**

短期内（5/17 wet lab 启动前）**不需要**。中期（cycle 之间）如果团队要快速 sync 新 protocol 版本到 dry lab，可以再考虑接 SDK。

如果之后要接，最简版本如下（已经测试过这种 pattern 在 benchling-integration skill 里 work）：

```python
from benchling_sdk.benchling import Benchling
from benchling_sdk.auth.api_key_auth import ApiKeyAuth

benchling = Benchling(
    url="https://claremont.benchling.com",
    auth_method=ApiKeyAuth(os.environ["BENCHLING_API_KEY"])
)
```

API key 申请方式：Benchling Profile Settings → API Tokens。**绝对不要 commit 进 repo**，存到 `.env`。

**下一步**：Alex 这周不碰 Benchling SDK，先用 `docs/` 里的 PDF 推进。

### 议题 3：dry lab 分工

**当前状态**：明确了，结构如下。

| 角色 | 人 | 职责 |
|---|---|---|
| **Dry lab core (coding)** | Alex | 全部 model / pipeline / acquisition function 实作。**主负责，避免沟通成本** |
| **Dry lab review** | 全员 | 每周一次 code/methodology review |
| **ESM 顾问** | 队内有 ESM 经验的同学 | 验证 embedding pipeline、PLM 选型咨询 |
| **HPC 操作** | 队内会用 Laguna 超算的同学 | ESM-2 pretrain / large-scale embedding / Boltz-2 batch inference |

**这个分工的关键 implication**：

1. **Laguna access 是大杀器**：原本 12GB GPU 限制 ESM-2 650M，有 Laguna 后 ESM-2 3B 甚至 15B 都能用，embedding 质量直接升一个档次。Boltz-2 batch inference 也不用排队
2. **ESM 顾问的存在** = 我们不用从头摸 PLM 工程化的坑
3. **Alex 单点的 risk**：如果 Alex 病倒/被功课压垮，整个 dry lab 停摆。**Mitigation**：每周一次的 review session 必须包含 walkthrough，让至少一位生物背景队友能看懂 codebase 大致结构（不需要会改，但要能 debug 输入输出）

---

## 八、现有 Protocols 盘点 + Gap 分析

Sarah 已经把 Benchling 上的 5 份 protocols 下到 `docs/`。我把每一份读完了，下面是诚实的 coverage 评估。

### 8.1 已有 protocols（5 份）的内容跟覆盖范围

| Protocol | 关键参数 | 在新方向中的角色 | 评估 |
|---|---|---|---|
| **Xanthomonas Cultivation** | NYG medium, 30°C, 200 rpm, OD600, kanamycin 20 µg/mL | Wet lab 全程基础——所有实验都要先把 Xcc 养起来 | ✅ **可直接使用** |
| **Xanthomonas Transformation** | Xcc 电转 (14 kV/cm, 10 µF), 28°C 恢复, kanamycin 选择 | 受体敲除（Step 8）的核心工具 | ✅ **可直接使用**——但需要补 knockout-specific protocol（pK18mobsacB 或类似系统） |
| **Plaque Assay (Small Scale)** | 100 µL log-phase Xcc + 系列稀释 phage, NYG soft agar 0.6-0.7%, 30°C 隔夜 | Spot test / 感染验证 | ✅ **可直接使用** |
| **Phage Infection Curves** | OD600 每 6 小时一次，36 小时，CFU 计数 | 名义上是 phage infection，**实际 protocol 写的是 bacterial density 对生长的影响**——没真的加 phage | ⚠️ **协议命名跟内容不符**，需要重写一版真的「phage challenge curve」protocol（详见下方 Gap 1） |
| **Liquid Amplification (Large Scale)** | MOI 0.01-0.1, 30°C 200 rpm, 4-8 小时 lysis, 0.22 µm 过滤 | Phage stock 制备 | ✅ **可直接使用** |

### 8.2 缺失的 protocols（按优先级排序）

新方向（active learning + motif-level RBP binding）必须新增的 protocols：

| Priority | 缺失 protocol | 用途 | 估计撰写时间 | 谁写 |
|---|---|---|---|---|
| 🔴 **P0** | **RBP-His6 表达 in BL21** | Cycle 0-2 全部 RBP variants 的核心表达系统 | 0.5 天（adapt 标准 pET 协议） | Alex + PI |
| 🔴 **P0** | **Ni-NTA 纯化 + SDS-PAGE QC** | 拿到纯蛋白才能做 ELISA | 0.5 天 | Sarah + Alex |
| 🔴 **P0** | **RBP-Xanthomonas 全细胞结合 ELISA / dot blot** | **整个新方向最核心的 readout**—— ground truth label 的来源 | 1-2 天，需要 PI 协助 optimization | 全员 |
| 🟡 P1 | **真正的 Phage Infection Curve（with phage challenge）** | Step 6 / Layer 2 calibration 的定量 readout | 0.5 天（基于现有 protocol 改） | Sarah |
| 🟡 P1 | **Spot test on receptor-KO panel** | Step 8 因果验证 | 0.5 天 | Sarah |
| 🟡 P1 | **Receptor knockout protocol（具体哪个系统）** | Step 8 的实操方案——markerless deletion / Tn-mutant / CRISPRi | 1 天，PI 选定系统后写 | PI 主导 |
| 🟢 P2 | **Western blot for RBP QC** | 表达验证 | 0.5 天（standard） | Sarah |
| 🟢 P2 | **Cloning workflow（Gibson 或 Golden Gate）** | RBP variants 入 pET vector | 0.5 天 | Sarah |

**估计总工作量**：~6-8 个 protocol working days，分散到 wet lab team 三个人 + Alex 协助起草，**5/17 前能写完 P0 的 3 份 + P1 的 2 份**。

### 8.3 关键注意点（对 PI / wet lab 队友的几个 flag）

**Flag 1：Phage Infection Curves protocol 名实不符**

现有 protocol 是测 *bacterial density* 对生长的影响，**没加 phage**。新方向需要的是「固定 MOI 下的 phage challenge over time」，跟现有 protocol 起点不同。需要 Sarah 重写一版。

**Flag 2：RBP 表达 host 是 BL21 不是 Xcc**

Transformation protocol 是 Xcc 电转——这个在做受体敲除（Step 8）时用得上，**但 RBP-His6 表达不在 Xcc 做，要在 *E. coli* BL21(DE3) 做**（标准 pET 系统）。需要确认 PI 实验室有 BL21 跟 IPTG 跟 Ni-NTA。

**Flag 3：ELISA optimization 是 wet lab 最大不确定性**

「RBP 蛋白 vs Xanthomonas 全细胞 / 纯化 LPS / 纯化 OMP」的 binding assay format 在文献上有几种做法（plate-coating with cells、dot blot、SPR），**每种都要 1-2 周 optimization 才稳定**。建议 cycle 0 的前两周先用 P2 gpH 或 T7 gp17 这类有 published binding data 的 RBP 当 positive control，确认 assay 能 work，再放 Xanthomonas RBP variants 进去。

**Flag 4：Step 7 / Step 8 是 critical path**

Proposal Step 7（directed mutagenesis）跟 Step 8（receptor KO）是 Active Learning 真正能发挥的地方—— acquisition function 推荐的 variants 进 Step 7，KO panel 当 Layer 2 ground truth。这两步的 protocol **必须在 cycle 1 启动前（~6/14）就 work**。

---

## 九、接下来 48-72 小时具体动作（Alex's checklist）

> **本节根据 v1.2 strain/phage 选择更新**。USDA permit 风险让「跟 PI 确认材料」变成 **#1 最高优先级**——其他事都次于此。

### 🔴 P0 - 必须 48 小时内完成（材料获取 critical path）

1. ☐ **跟 PI 紧急 sync 3 个问题**（决定 wet lab 能不能按时启动）：
   - 实验室有没有 Xanthomonas glycerol stock（任何 pathovar 都问）？
   - 实验室有没有 standing USDA APHIS PPQ-526 permit？
   - 有没有合作 lab（Pomona / HMC / KGI / 别的 PI）可以 transfer Xanthomonas？

2. ☐ **如果 PI 答案全是 No**：立即启动 Path B，去 farmer's market / 加州本地农场买带黑腐病症状的卷心菜 / 花椰菜 / 芥菜，回来分离 Xcc + 同时 plaque 出 phage（Sarah / 生物组主导，1-2 周）

3. ☐ **同时 fire-and-forget**：发 email 给 phiL7 paper 通讯作者（Academia Sinica）询问 phage 是否可 share，作为 long-shot Path C

### 🟡 P1 - 1 周内（dry lab 可独立推进，不依赖材料）

**Dry lab 端**：

4. ☐ 下载 phiL7 完整基因组（NCBI EU717894）跟 Xcc ATCC 33913 基因组（AE008922）
5. ☐ 用 PHANOTATE 注释 phiL7，识别 tail spike / tail fiber RBP 候选基因
6. ☐ AlphaFold 跑 phiL7 RBP 3D 结构预测（Laguna 上跑，跟 ESM 顾问 sync）
7. ☐ AlphaFold 跑 Xcc TonB / ExbB / ExbD1 结构预测（ExbD2 保留作為負對照）
8. ☐ 写 baseline ML pipeline：ESM-2 embedding + deep ensemble + BALD acquisition function（先简化版）
9. ☐ Design seed variants（4-6 个 truncations + 2 chimeras + 2 point mutants）—— 这次基于 phiL7 RBP 真实序列设计

**Wet lab 端（请 Sarah / 队友主导）**：

10. ☐ 把缺失的 P0 三份 protocol（RBP-His6 expression、Ni-NTA 纯化、binding ELISA）的初稿在 Benchling 起头
11. ☐ 重写 Phage Infection Curve protocol（加入 phage challenge 步骤）
12. ☐ 跟 PI 确认 BL21(DE3) + IPTG + Ni-NTA + IDA / NTA resin 系统齐备性
13. ☐ 跟 PI 确认 receptor knockout 用什么系统（pK18mobsacB markerless / Tn-mutant 库 / CRISPRi）

### 🟢 P2 - 1-2 周内

14. ☐ 全员 30 分钟 sync meeting，对齐方向 + 时间表 + permit 现况
15. ☐ Alex 写 1-page summary 给 PI 拿到 sign-off
16. ☐ 如果走 Path B，开始 16S rRNA / MALDI-TOF identification of isolated Xcc strains
17. ☐ 准备 5/17 前 gene synthesis 下单（针对最终选定的 phiL7 RBP variants）

### 关键时间节点

```
今天 (5/7) ──────► 48h (5/9) PI 确认材料状态
                      ↓
        Path A 成功 ─────────► 5/17 wet lab 按计划启动
                      ↓ Path A 失败
        启动 Path B 本地分离 ─► 5/24 拿到 in-house Xcc isolate
                      ↓
                              5/30 ELISA 准备 → 6/1 cycle 0 启动
```

**最坏情况评估**：如果 PI 没库存、没 permit、没合作 lab，全部走 Path B local isolation，wet lab 启动会 delay 1-2 周，cycle 数从 2-3 轮压缩到 1-2 轮，但 **不会 kill 项目**。Dry lab 端可以平行推进不浪费时间。

---

## 八、关于「文献支持不足」这件事的承诺

Alex 自己的判断：原 proposal 中很多 dry lab assumption 都是用比较直觉的数学逻辑去推的，并没有学术理论支持。这是真实的 weakness。

**新方向往前走的核心原因之一就是要把这块补扎实**。后续每一步——从 PLM 选型、uncertainty 估计方法、acquisition function 设计、到 cycle 设计——我都会要求每个 sub-step 至少有 1 篇 anchor paper。会建立一个 living `literature_anchors.md` 文档跟 wiki 同步。

---

## 附录：关键参考文献（先看这几篇）

### 主动学习 / 实验设计
- Settles (2009) Active learning literature survey
- Frazier (2018) A tutorial on Bayesian optimization
- Hie et al. (2022) *Cell* — Efficient evolution of human antibodies from general protein language models

### Phage-host prediction（同领域 SOTA）
- Boeckaerts et al. (2024) *Nature Communications* — PhageHostLearn
- Liu et al. (2025) *Nature Communications* — PLM-interact
- Coelho et al. (2025) *Scientific Reports* — PPI-based phage-host prediction

### Phage RBP 生物学
- Yehl et al. (2019) *Cell* — T7 deep mutational scanning
- Latka et al. (2021) *mBio* — Klebsiella RBP modularity
- Garcia-Doval & van Raaij (2012) *PNAS* — T7 gp17 tip 结构

### 不确定性估计 (uncertainty quantification)
- Lakshminarayanan et al. (2017) NeurIPS — Deep ensembles
- Gal & Ghahramani (2016) ICML — MC Dropout

### **本项目核心 references（v1.2 strain/phage 选择依据）**
- **Lee et al. (2009)** Genomic Characterization of the Intron-Containing T7-Like Phage phiL7 of Xanthomonas campestris. *Applied and Environmental Microbiology* 75(24):7828. — phiL7 完整测序，NCBI EU717894
- **Hung, C.-H. et al. (2003)** Involvement of tonB-exbBD1D2 operon in infection of Xanthomonas campestris phage phiL7. *Biochem Biophys Res Commun* 302(4):878–884. PMID 12646254 — **phiL7 receptor 鉴定为 TonB-dependent system**
- **da Silva et al. (2002)** Comparison of the genomes of two Xanthomonas pathogens with differing host specificities. *Nature* 417:459. — Xcc ATCC 33913 全基因组（NCBI AE008922）
- **Holtappels et al. (2022)** The potential of bacteriophages to control Xanthomonas campestris pv. campestris at different stages of disease development. *Microbial Biotechnology*. — Xcc + phage biocontrol 综述
- **Nakayinga et al. (2021)** Xanthomonas bacteriophages: a review of their biology and biocontrol applications in agriculture. *BMC Microbiology* 21:291. — Xanthomonas phage 全面综述

### 关键 NCBI accession（dry lab 立即可下载）
- Xcc ATCC 33913 全基因组：[AE008922](https://www.ncbi.nlm.nih.gov/nuccore/AE008922)
- phiL7 全基因组：[EU717894](https://www.ncbi.nlm.nih.gov/nuccore/EU717894)
- Xcc 在 BacDive 的全 cross-reference：[BacDive 17584](https://bacdive.dsmz.de/strain/17584)

### USDA / 监管 references
- USDA APHIS [PPQ Form 526 申请指南](https://www.aphis.usda.gov/sites/default/files/ppq-526-guidance.pdf) — 平均 127 天处理
- ATCC [33913 产品页](https://www.atcc.org/products/33913) — $486, 需 PPQ-526 permit

---

**文档版本**：v1.2
**最后更新**：2026-05-07
**变更日志**：
- v1.0 → v1.1：修正 Xp10 (X. oryzae host 错配) → 改推 CP2 / Xp15
- v1.1 → v1.2：基于 web 调研重写 — CP2 实为 filamentous Inoviridae 不适用；最终推荐 **Xcc ATCC 33913 + phiL7 (EU717894)**，识别 USDA permit 为关键 blocker，提出 Path A/B/C 应对策略

**下一次更新**：48 小时后 PI 回复 material 状态后
