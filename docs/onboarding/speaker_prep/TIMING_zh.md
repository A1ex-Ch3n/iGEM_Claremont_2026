# TIMING_zh.md — 分钟级时间表

**目标总长:** 75 分钟讲演 + 15 分钟 Q&A = 90 分钟从头到尾。
**硬底线:** 60 分钟(赶时跳过 OPTIONAL)。
**演示 (15 分):** 独立 — 在 Part 6 之后或最后 Q&A 缓冲跑。

---

## 各 Part 分配

| Part | 分 | 累计 | 标签 | 听众优先 |
|------|-----|------|------|---------|
| 0. Roadmap | 4 | 4 | [CORE] | 全员 |
| 1. 科学 | 10 | 14 | [CORE] | 湿实验 + PI ★★★ |
| 2. 管道架构 | 12 | 26 | [CORE 除 Module 04 备路径] | 干实验 + PI |
| 3. ML 核心 | 15 | 41 | [CORE] | 干实验 + PI ★★★ |
| 4. Boltz-2 结果 | 8 | 49 | [CORE] | 全员 |
| 5. 48-h 闭环 | 12 | 61 | [CORE] | 湿实验 + 干实验 ★★★ |
| 6. 复现 | 6 | 67 | [OPTIONAL — 拖时跳] | 干实验 + 湿实验 |
| 7. 风险 & 决策 | 6 | 73 | [CORE for PI] | PI ★★★ |
| 8. 引用 & 致谢 | 2 | 75 | [CORE 只致谢] | 全员 |
| Q&A 缓冲 | 15 | 90 | — | 全员 |
| 演示(可选) | 15 | +15 | [OPTIONAL] | 干 + 湿 |

---

## 检查点 — 落后则调整

### 第 15 分钟 — 应在 Part 1 末
- **准时:** Part 1 结束,过渡到架构。
- **落后 5 分:** 跳 "Suggested cuts" + "How to read";压缩 vocab block;Part 1 续接。
- **落后 >10 分:** 砍"数据稀缺瓶颈" slide(与开场重复);跳到"为什么主动学习"。

### 第 30 分钟 — 应在 Part 2 中(Module 03 / 04)
- **准时:** Module 03 RBP slide 结束,启动 embedding。
- **落后 5 分:** 跳 "Module 02 Annotation" 和 "Module 04 PLM-interact 段";跳到 Module 05 Boltz-2。
- **落后 >10 分:** Module 00–04 压成一分钟,跳到 Module 05。

### 第 45 分钟 — 应在 Part 3 中(BALD 数学)
- **准时:** BALD 数学 slide。
- **落后 5 分:** 跳 "ensemble.py snippet" + "calibration plot";跳到 BALD 数学。
- **落后 >10 分:** 完全靠"BALD 直觉"图;跳过推导;跳到 bald.py snippet。

### 第 60 分钟 — 应在 Part 5 中(handoff)
- **准时:** 在"湿→干 handoff"或"干→湿 handoff"。
- **落后 5 分:** 压缩受体敲除 slide;跳到验证 tier。
- **落后 >10 分:** 跳过 cloning 和 ELISA protocol slide(guide 里有);跳到验证 tier + 质量门。

### 第 75 分钟 — 应在致谢
- **准时:** Part 7/8 收尾,开 Q&A。
- **超时:** 立即开 Q&A。要求时才跑演示。

---

## 按听众优先压缩

**湿实验偏多** :
- 压缩:Part 3 BALD 数学 + ensemble 代码 snippet(砍到"5 网络投票,挑分歧大的")。
- 展开:Part 5 handoff 细节,ELISA protocol,敲除策略。
- 时间:Part 3 → Part 5 转 5 分。

**干实验/PI 偏多**:
- 压缩:Part 5 cloning + ELISA protocol(一张 slide 概览)。
- 展开:Part 3 BALD 数学 + ALDE caveat + Greenman 2025 caveat。
- 时间:Part 5 → Part 3 转 5 分。

**只 PI(45 分钟)**:
- 砍 Part 6。
- 压缩 Part 5 为"48h SLA 下两个 handoff"。
- 以 Part 7 风险 + 决策开头。
- 顺序:1, 5 (Hung 2003 受体), 3 (简), 4 (Boltz-2 ipTM), 5 (48h cycle + Tier 3), 7 (决策)。

---

## 可砍而不毁讲演的 slide

任何优先模式下都可跳:
- "How to read this deck" (Part 0)
- "Suggested cuts" (Part 0)
- "Notebook-first workflow + bilingual comments" (Part 2)
- "Calibration plot — Cycle 0" (Part 3) ——口头提一下
- "ensemble.py snippet" (Part 3) ——纯视觉
- "Where the Boltz-2 outputs live" (Part 4) ——指 PI briefing
- "Quality gates and failure modes" (Part 5) ——挥一下
- "Per-module entry points" (Part 6) ——`GETTING_STARTED.md` 有
- 四张引用附录 slide (Part 8) ——指 `papers.md`

全砍能省 ~10 分钟。

---

## 即便迟到也**不要**跳

- Module 03 (rbp_01 发现) ——你的关键干实验结果
- Module 05 Boltz-2 (ipTM caveat) ——最大可信度风险 slide
- BALD 数学 + ALDE caveat ——你的方法学贡献 + 准确性守门
- Boltz-2 三个数字 slide (ipTM 0.365) ——你的实时结果
- 受体敲除 (Tier 3 + ΔexbD2 负对照) ——PI 相关
- 风险 + 决策 slide ——PI 来这就为决定这些

---

## Q&A 拖长怎么办

- 演示可推到讲演后(小组)。
- 一直答到 95 分钟,然后宣布:"我会留到所有人都问完——但我们正式结束,让需要走的人走。"
- 不要在问题中途突然结束。在自然节点结束。

---

## EN 镜像 — TIMING.md

英文版同样计划在 `TIMING.md`。讲演时主用一种语言;Q&A 时遇到 ZH 听众可即时切换。
