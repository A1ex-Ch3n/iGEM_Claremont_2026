"""
make_figures.py — Generate onboarding diagrams.
生成入职演示所需的流程图与示意图。

Run from repo root:
    python docs/onboarding/figures/make_figures.py
"""
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

OUT = Path(__file__).resolve().parent

plt.rcParams.update({
    "font.family": "Hiragino Sans GB",
    "axes.spines.top": False,
    "axes.spines.right": False,
})


# ---------------------------------------------------------------------------
# 1. pipeline_flow.png — Modules 00..08 with data-flow arrows
# 1. pipeline_flow.png — 模块 00..08 数据流总图
# ---------------------------------------------------------------------------

def pipeline_flow():
    modules = [
        ("00", "Raw Data\n原始基因组",       "777 phage + 34 bacteria",        "#E8F1FF"),
        ("01", "Ground Truth\n相互作用矩阵", "2,236 phage–host pairs",         "#E8F1FF"),
        ("02", "Annotation\n基因注释",       "PHANOTATE / pyrodigal\n.faa proteins", "#E8F1FF"),
        ("03", "RBP ID\n受体结合蛋白识别",   "PhageRBPdetect HMM\nrbp_01 = 712 aa", "#FFF4E0"),
        ("04", "Embedding\n蛋白嵌入",       "ESM-2 (Lin 2023)\n1280-dim vector", "#FFF4E0"),
        ("05", "Structure\n结构预测",       "Boltz-2 ipTM\n(not affinity!)",  "#FFF4E0"),
        ("06", "Deep Ensemble\n深度集成",   "5 MLPs · epistemic_std\nLakshmin. 2017", "#E8F8E8"),
        ("07", "BALD\n采集函数",            "Var_k[μ_k]\nrecommendations.csv", "#E8F8E8"),
        ("08", "Cycle Data\n循环数据",      "ELISA → retrain loop\n(~June 1)", "#FFE0E8"),
    ]

    fig, ax = plt.subplots(figsize=(16, 4.2))
    ax.set_xlim(0, len(modules))
    ax.set_ylim(0, 4)
    ax.axis("off")

    for i, (num, name, descr, color) in enumerate(modules):
        x = i + 0.5
        box = FancyBboxPatch(
            (i + 0.05, 1.2), 0.9, 2.0,
            boxstyle="round,pad=0.02,rounding_size=0.08",
            linewidth=1.4, edgecolor="#2c3e50", facecolor=color,
        )
        ax.add_patch(box)
        ax.text(x, 2.85, num, ha="center", va="center",
                fontsize=15, fontweight="bold", color="#1a4480")
        ax.text(x, 2.25, name, ha="center", va="center",
                fontsize=9.5, fontweight="bold")
        ax.text(x, 1.55, descr, ha="center", va="center", fontsize=7.5,
                color="#34495e")

        if i < len(modules) - 1:
            arr = FancyArrowPatch(
                (i + 0.95, 2.2), (i + 1.05, 2.2),
                arrowstyle="->", mutation_scale=14, color="#2c3e50", linewidth=1.4,
            )
            ax.add_patch(arr)

    # Stripe labels
    stripe_y = 0.55
    ax.add_patch(mpatches.Rectangle((0.05, stripe_y), 2.9, 0.4, facecolor="#E8F1FF", edgecolor="none"))
    ax.text(1.5, stripe_y + 0.2, "DATA  数据", ha="center", va="center", fontsize=9, fontweight="bold")
    ax.add_patch(mpatches.Rectangle((3.05, stripe_y), 2.9, 0.4, facecolor="#FFF4E0", edgecolor="none"))
    ax.text(4.5, stripe_y + 0.2, "FEATURES  特征", ha="center", va="center", fontsize=9, fontweight="bold")
    ax.add_patch(mpatches.Rectangle((6.05, stripe_y), 1.9, 0.4, facecolor="#E8F8E8", edgecolor="none"))
    ax.text(7.0, stripe_y + 0.2, "MODEL  建模", ha="center", va="center", fontsize=9, fontweight="bold")
    ax.add_patch(mpatches.Rectangle((8.05, stripe_y), 0.9, 0.4, facecolor="#FFE0E8", edgecolor="none"))
    ax.text(8.5, stripe_y + 0.2, "WET LAB", ha="center", va="center", fontsize=9, fontweight="bold")

    ax.text(len(modules) / 2, 3.7,
            "iGEM Claremont 2026 — Active-Learning Pipeline / 主动学习流程",
            ha="center", va="center", fontsize=13, fontweight="bold", color="#1a4480")

    fig.savefig(OUT / "pipeline_flow.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# 2. data_contract.png — inputs / processes / outputs convention
# 2. data_contract.png — 模块内部数据流约定
# ---------------------------------------------------------------------------

def data_contract():
    fig, ax = plt.subplots(figsize=(12, 5.2))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 6)
    ax.axis("off")

    def folder(x, y, w, h, label_en, label_zh, descr, color):
        box = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.04,rounding_size=0.1",
            linewidth=1.6, edgecolor="#2c3e50", facecolor=color,
        )
        ax.add_patch(box)
        ax.text(x + w / 2, y + h - 0.45, label_en, ha="center", va="center",
                fontsize=12, fontweight="bold", color="#1a4480")
        ax.text(x + w / 2, y + h - 0.85, label_zh, ha="center", va="center",
                fontsize=9.5, color="#1a4480")
        ax.text(x + w / 2, y + h / 2 - 0.5, descr, ha="center", va="center",
                fontsize=8.5, color="#34495e")

    # Module N
    ax.text(6, 5.6, "Module N — data contract / 模块内数据流约定",
            ha="center", fontsize=13, fontweight="bold")

    folder(0.3, 1.0, 3.4, 3.5, "inputs/", "输入（只读）",
           "• Pointers to upstream\n  outputs (or seeds)\n"
           "• READ-ONLY\n• 指向上游 outputs\n  或外部 seed", "#E8F1FF")
    folder(4.3, 1.0, 3.4, 3.5, "processes/", "处理脚本与笔记本",
           "• .ipynb (exploration)\n• .py (frozen)\n"
           "• bilingual comments\n• .ipynb 笔记本 +\n  冻结的 .py 模块",  "#FFF4E0")
    folder(8.3, 1.0, 3.4, 3.5, "outputs/", "输出（可消费）",
           "• Canonical artifacts\n• Large files gitignored\n  → MANIFEST.csv\n"
           "• 下游模块消费\n• 大文件 gitignore\n  → MANIFEST.csv", "#E8F8E8")

    for x0, x1 in [(3.7, 4.3), (7.7, 8.3)]:
        arr = FancyArrowPatch((x0, 2.75), (x1, 2.75),
                              arrowstyle="->", mutation_scale=18,
                              color="#2c3e50", linewidth=1.6)
        ax.add_patch(arr)

    # Bottom — flow into next module
    ax.text(6, 0.45,
            "outputs/  →  Module (N+1) inputs/   (READ-ONLY POINTER, no copy)",
            ha="center", fontsize=10, fontweight="bold", color="#1a4480")

    fig.savefig(OUT / "data_contract.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# 3. cycle_48h.png — Wet ↔ dry 48-hour handoff
# 3. cycle_48h.png — 干湿实验室 48 小时往返循环
# ---------------------------------------------------------------------------

def cycle_48h():
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 7)
    ax.axis("off")

    ax.text(5.5, 6.55, "The 48-hour Active-Learning Cycle / 48 小时主动学习循环",
            ha="center", fontsize=13.5, fontweight="bold", color="#1a4480")

    def step(x, y, w, h, title, body, color):
        box = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.04,rounding_size=0.1",
            linewidth=1.4, edgecolor="#2c3e50", facecolor=color,
        )
        ax.add_patch(box)
        ax.text(x + w / 2, y + h - 0.35, title, ha="center", va="center",
                fontsize=10.5, fontweight="bold", color="#1a4480")
        ax.text(x + w / 2, y + h / 2 - 0.15, body, ha="center", va="center",
                fontsize=8.5, color="#34495e")

    # Dry lab top, wet lab bottom
    step(0.5, 4.4, 2.3, 1.5, "① Ensemble retrain\n① 集成重训",
         "5 MLPs on cumulative\nELISA data (~25 s)\n5 个 MLP 在累计数据\n上独立训练", "#E8F8E8")
    step(3.4, 4.4, 2.3, 1.5, "② BALD score\n② BALD 打分",
         "Var_k[μ_k] per\nunmeasured pair\n每个未测对的\nepistemic_std", "#E8F8E8")
    step(6.3, 4.4, 2.3, 1.5, "③ Pick batch\n③ 选下批",
         "Top-4 BALD\n+ 1 random control\n4 个 BALD 最优\n+ 1 个随机对照", "#E8F8E8")
    step(8.7, 4.4, 2.0, 1.5, "④ Recommendations\n④ 推荐文件",
         "recommendations.csv\nprimer_sequences.txt\nsafe_pick_backup.csv", "#E8F8E8")

    step(0.5, 0.8, 2.3, 1.5, "⑧ Process ELISA\n⑧ 整理数据",
         "Fit 4PL → EC50\n4 参数曲线拟合\nelisa_processed.csv", "#FFE0E8")
    step(3.4, 0.8, 2.3, 1.5, "⑦ ELISA + plaque\n⑦ ELISA + 噬斑",
         "Cell-based ELISA\n(WT and ΔReceptor)\n基础 ELISA + 受体\n敲除株对照", "#FFE0E8")
    step(6.3, 0.8, 2.3, 1.5, "⑥ Express + purify\n⑥ 表达纯化",
         "BL21(DE3) + Ni-NTA\nHis6-RBP variants\n～10–14 d (SDM)", "#FFE0E8")
    step(8.7, 0.8, 2.0, 1.5, "⑤ SDM / synthesis\n⑤ SDM 或合成",
         "NEB Q5 SDM\n($50/variant)\nor gene synthesis\n(Cycle 0 only)", "#FFE0E8")

    # Lane labels
    ax.text(0.3, 5.15, "DRY LAB\n干实验室", ha="right", va="center",
            fontsize=10, fontweight="bold", color="#1a4480", rotation=0)
    ax.text(0.3, 1.55, "WET LAB\n湿实验室", ha="right", va="center",
            fontsize=10, fontweight="bold", color="#a01a4d", rotation=0)

    # Flow arrows (dry → wet → dry)
    for (x0, y0, x1, y1) in [
        (2.8, 5.15, 3.4, 5.15), (5.7, 5.15, 6.3, 5.15), (8.6, 5.15, 8.7, 5.15),
        (2.8, 1.55, 3.4, 1.55), (5.7, 1.55, 6.3, 1.55), (8.6, 1.55, 8.7, 1.55),
    ]:
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="->",
                                     mutation_scale=14, color="#2c3e50", linewidth=1.4))
    # Down/up bridges
    ax.add_patch(FancyArrowPatch((9.7, 4.3), (9.7, 2.4), arrowstyle="->",
                                 mutation_scale=18, color="#1a4480", linewidth=2.0))
    ax.text(10.15, 3.35, "48-h SLA\n48 小时", fontsize=9, color="#1a4480", va="center")

    ax.add_patch(FancyArrowPatch((1.0, 2.4), (1.0, 4.3), arrowstyle="->",
                                 mutation_scale=18, color="#a01a4d", linewidth=2.0))
    ax.text(0.55, 3.35, "10–14 d\n湿实验室周期", fontsize=9,
            color="#a01a4d", va="center", ha="right")

    fig.savefig(OUT / "cycle_48h.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# 4. bald_intuition.png — toy plot of BALD picks
# 4. bald_intuition.png — BALD 直观示意（玩具数据）
# ---------------------------------------------------------------------------

def bald_intuition():
    rng = np.random.default_rng(42)
    n = 60
    preds = rng.normal(5.0, 0.6, size=n)
    epis = np.abs(rng.normal(0.0, 0.08, size=n)) + 0.02
    # Inject a few high-uncertainty points
    high_idx = rng.choice(n, size=4, replace=False)
    epis[high_idx] = rng.uniform(0.18, 0.24, size=4)

    top4 = np.argsort(epis)[-4:]
    random_pick = rng.choice([i for i in range(n) if i not in top4])

    fig, ax = plt.subplots(figsize=(10, 5.4))
    colors = ["#7f8c8d"] * n
    for idx in top4:
        colors[idx] = "#1a4480"
    colors[random_pick] = "#a01a4d"

    ax.scatter(preds, epis, c=colors, s=80, edgecolor="white", linewidth=1.0, zorder=3)
    ax.set_xlabel("Predicted binding score / 预测结合得分", fontsize=11)
    ax.set_ylabel("Epistemic std (BALD score) / epistemic 标准差", fontsize=11)
    ax.set_title("BALD selection intuition — pick where the model disagrees most\n"
                 "BALD 选取直观：选模型成员分歧最大的点",
                 fontsize=12, color="#1a4480", pad=12)
    ax.axhline(np.sort(epis)[-5], color="#1a4480", linestyle="--", alpha=0.4, lw=1)
    ax.text(ax.get_xlim()[1] - 0.05, np.sort(epis)[-5] + 0.005,
            "BALD top-4 cutoff / BALD top-4 阈值",
            color="#1a4480", ha="right", fontsize=9)

    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#7f8c8d",
                   markersize=10, label="Candidate / 候选"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#1a4480",
                   markersize=10, label="BALD top-4 (recommend)"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#a01a4d",
                   markersize=10, label="Random control (mandatory)"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=10, frameon=True)
    ax.grid(True, alpha=0.25, linestyle=":")

    fig.tight_layout()
    fig.savefig(OUT / "bald_intuition.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    pipeline_flow()
    data_contract()
    cycle_48h()
    bald_intuition()
    print("All 4 figures written to:", OUT)
