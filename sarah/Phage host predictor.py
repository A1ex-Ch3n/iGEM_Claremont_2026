"""
Phage-Host Interaction Predictor
==================================
支持中等规模数据（20-100 phages × 50-200 bacteria）

功能：
  1. 数据输入：手动 dict 或 CSV/Excel 文件
  2. 进化距离：通过 NCBI Taxonomy 层级自动估算
  3. NaN 预测：距离加权计算 predicted infection probability
  4. Confidence：每个 phage 的已知数据覆盖率
  5. 输出：CSV、Excel（多 sheet）、Heatmap（PNG）

依赖：
  pip install pandas numpy openpyxl matplotlib seaborn requests
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import requests
import time
import os
from typing import Optional


# ══════════════════════════════════════════════════════════
#  1. 示例数据（手动输入模式）
# ══════════════════════════════════════════════════════════

def create_example_data() -> pd.DataFrame:
    """
    手动定义 phage-host 交互矩阵。
    - 行   = bacteria（物种全名，用于 NCBI 查询）
    - 列   = phage
    - 值   = 1（感染）/ 0（不感染）/ np.nan（无记录）

    ← 替换成你自己的数据即可 →
    """
    bacteria = [
        "Escherichia coli",
        "Salmonella enterica",
        "Klebsiella pneumoniae",
        "Pseudomonas aeruginosa",
        "Staphylococcus aureus",
        "Bacillus subtilis",
        "Listeria monocytogenes",
        "Vibrio cholerae",
    ]
    data = {
        "PhageA": [1,       1,       0,       np.nan,  0,       np.nan,  np.nan,  1      ],
        "PhageB": [0,       np.nan,  1,       1,       np.nan,  0,       np.nan,  np.nan ],
        "PhageC": [np.nan,  0,       np.nan,  0,       1,       1,       np.nan,  0      ],
        "PhageD": [1,       np.nan,  np.nan,  1,       0,       np.nan,  1,       np.nan ],
        "PhageE": [np.nan,  1,       0,       np.nan,  np.nan,  0,       1,       1      ],
    }
    return pd.DataFrame(data, index=bacteria)


# ══════════════════════════════════════════════════════════
#  2. 从文件读取数据
# ══════════════════════════════════════════════════════════

def load_from_file(filepath: str) -> pd.DataFrame:
    """
    从 CSV 或 Excel 读取 phage-host 矩阵。

    文件格式：
      - 第一列 = bacteria names（行索引）
      - 其余列 = phage names
      - 值     = 1 / 0 / 空白（空白 → NaN）

    示例 CSV：
        bacteria,PhageA,PhageB
        Escherichia coli,1,0
        Salmonella enterica,1,
        Klebsiella pneumoniae,,1
    """
    ext = os.path.splitext(filepath)[-1].lower()
    if ext == ".csv":
        df = pd.read_csv(filepath, index_col=0)
    elif ext in (".xlsx", ".xls"):
        df = pd.read_excel(filepath, index_col=0)
    else:
        raise ValueError(f"不支持的文件格式：{ext}，请使用 .csv 或 .xlsx")

    df = df.apply(pd.to_numeric, errors="coerce")
    print(f"✅ 已读取文件：{filepath}  →  {df.shape[0]} bacteria × {df.shape[1]} phages")
    return df


# ══════════════════════════════════════════════════════════
#  3. NCBI Taxonomy 层级距离估算
# ══════════════════════════════════════════════════════════

# 每个 taxonomy 层级的距离权重（层级越高距离越大）
TAXONOMY_RANKS = ["species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"]
RANK_DISTANCE  = {r: w for r, w in zip(TAXONOMY_RANKS, [0.0, 0.1, 0.3, 0.5, 0.7, 0.85, 0.95, 1.0])}


def fetch_taxonomy(species_name: str, retries: int = 3) -> dict:
    """
    通过 NCBI E-utilities 获取物种的 taxonomy 层级信息。
    返回 dict: {rank: name}，例如 {"genus": "Escherichia", "family": "Enterobacteriaceae", ...}
    """
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    # Step 1: 搜索 TaxID
    for attempt in range(retries):
        try:
            r = requests.get(
                base + "esearch.fcgi",
                params={"db": "taxonomy", "term": species_name, "retmode": "json"},
                timeout=10,
            )
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            if ids:
                break
        except Exception:
            time.sleep(1)
    else:
        print(f"  ⚠️  NCBI 查询失败（无结果）：{species_name}")
        return {}

    if not ids:
        print(f"  ⚠️  未找到 TaxID：{species_name}")
        return {}

    tax_id = ids[0]
    time.sleep(0.35)  # NCBI 限速：~3 req/s

    # Step 2: 获取 lineage 信息
    for attempt in range(retries):
        try:
            r = requests.get(
                base + "efetch.fcgi",
                params={"db": "taxonomy", "id": tax_id, "retmode": "json"},
                timeout=10,
            )
            result = r.json()
            break
        except Exception:
            time.sleep(1)
    else:
        print(f"  ⚠️  NCBI efetch 失败：{species_name}")
        return {}

    result_section = result.get("result", {})
    tax_data = result_section.get(tax_id, {})
    if not isinstance(tax_data, dict):
        print(f"  ⚠️  NCBI 返回格式异常：{species_name}")
        return {}
    lineage = tax_data.get("lineage", "")
    lineageex = tax_data.get("lineageex", [])

    taxonomy = {}
    for item in lineageex:
        rank = item.get("rank", "").lower()
        name = item.get("scientificname", "")
        if rank in RANK_DISTANCE:
            taxonomy[rank] = name

    # 加入 species 本身
    taxonomy["species"] = species_name
    return taxonomy


def build_taxonomy_cache(bacteria_list: list, cache_file: str = "taxonomy_cache.csv") -> dict:
    """
    批量获取所有细菌的 taxonomy，带本地缓存（避免重复请求 NCBI）。
    返回 dict: {bacteria_name: {rank: name}}
    """
    # 读取已有缓存
    cache = {}
    if os.path.exists(cache_file):
        df_cache = pd.read_csv(cache_file, index_col=0)
        for bact, row in df_cache.iterrows():
            cache[bact] = row.dropna().to_dict()
        print(f"📂 读取缓存：{cache_file}（{len(cache)} 条记录）")

    # 查询未缓存的物种
    to_query = [b for b in bacteria_list if b not in cache]
    if to_query:
        print(f"🔍 向 NCBI 查询 {len(to_query)} 个物种的 taxonomy...")
        for i, bact in enumerate(to_query, 1):
            print(f"  [{i}/{len(to_query)}] {bact}")
            cache[bact] = fetch_taxonomy(bact)

        # 保存缓存
        rows = []
        for bact, tax in cache.items():
            row = {"bacteria": bact}
            row.update(tax)
            rows.append(row)
        pd.DataFrame(rows).set_index("bacteria").to_csv(cache_file)
        print(f"💾 缓存已保存：{cache_file}")

    return cache


def taxonomy_distance(tax_a: dict, tax_b: dict) -> float:
    """
    根据两个物种的 taxonomy 层级，计算进化距离。

    算法：找到最近公共祖先（LCA）所在的层级，
    距离 = 该层级对应的权重值。
    层级越高（越远离 species），距离越大。
    """
    if not tax_a or not tax_b:
        return 1.0  # 无法确定时，视为最远

    # 从 species 往上找 LCA
    for rank in TAXONOMY_RANKS:
        name_a = tax_a.get(rank)
        name_b = tax_b.get(rank)
        if name_a and name_b and name_a == name_b:
            return RANK_DISTANCE[rank]

    return 1.0  # 完全不同


def build_distance_matrix(bacteria_list: list, taxonomy_cache: dict) -> pd.DataFrame:
    """
    利用 taxonomy 信息构建细菌间的进化距离矩阵（对称，对角线=0）。
    """
    n = len(bacteria_list)
    dist = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            d = taxonomy_distance(
                taxonomy_cache.get(bacteria_list[i], {}),
                taxonomy_cache.get(bacteria_list[j], {}),
            )
            dist[i, j] = d
            dist[j, i] = d

    return pd.DataFrame(dist, index=bacteria_list, columns=bacteria_list)


# ══════════════════════════════════════════════════════════
#  4. 核心预测算法
# ══════════════════════════════════════════════════════════

def compute_confidence(interaction_matrix: pd.DataFrame) -> pd.Series:
    """
    Confidence = 非 NaN 条目数 / 总细菌数（按 phage 列计算）
    值域 [0, 1]，越高代表该 phage 的已知数据越充分。
    """
    n_total = len(interaction_matrix)
    n_known = interaction_matrix.notna().sum(axis=0)
    conf = (n_known / n_total).rename("confidence")
    return conf


def predict_infection_probability(
    interaction_matrix: pd.DataFrame,
    distance_matrix: pd.DataFrame,
    distance_scale: float = 0.5,
) -> pd.DataFrame:
    """
    对所有 NaN 位置计算 predicted infection probability。

    算法（距离加权 KNN 思想）：
      对于 phage p、bacteria b（值为 NaN）：
        1. 取 p 的所有已知宿主 b_i（值为 0 或 1）
        2. 计算 b 与每个 b_i 的进化距离 d_i
        3. 高斯权重：w_i = exp(−d_i² / distance_scale²)
        4. 归一化：w_i_norm = w_i / Σw
        5. predicted prob = Σ(w_i_norm × y_i)

    distance_scale：越小 → 只信任近邻；越大 → 远近均等参考
    """
    bacteria = interaction_matrix.index.tolist()
    phages   = interaction_matrix.columns.tolist()

    predicted = pd.DataFrame(np.nan, index=bacteria, columns=phages, dtype=float)

    # 检测整个矩阵是否只有 1（没有任何 0 记录）
    all_values = interaction_matrix.stack()
    has_zeros  = (all_values == 0).any()
    has_ones   = (all_values == 1).any()

    if has_zeros and has_ones:
        mode = "weighted_average"   # 标准模式：0/1 加权平均
    elif has_ones and not has_zeros:
        mode = "similarity_only"    # 只有 1：用进化相似度估算感染概率
    else:
        mode = "weighted_average"   # fallback

    print(f"  📊 预测算法：{'加权平均 (0+1)' if mode == 'weighted_average' else '进化相似度 (仅1，无0)'}")

    for phage in phages:
        col = interaction_matrix[phage]
        known_mask   = col.notna()
        known_bact   = col[known_mask].index.tolist()
        known_values = col[known_mask].values.astype(float)

        if len(known_bact) == 0:
            continue  # 该 phage 完全无记录

        for bact in bacteria:
            if not pd.isna(col[bact]):
                continue  # 已有记录，跳过

            dists = np.array(
                [distance_matrix.loc[bact, kb] for kb in known_bact], dtype=float
            )

            # 高斯核权重（距离越近权重越大）
            weights = np.exp(-(dists ** 2) / (distance_scale ** 2))
            w_sum   = weights.sum()

            if w_sum == 0:
                predicted.loc[bact, phage] = np.nan
                continue

            if mode == "weighted_average":
                # 标准模式：加权平均 0/1 值
                predicted.loc[bact, phage] = float(np.dot(weights / w_sum, known_values))

            else:
                # 相似度模式：只看已知感染宿主（值=1）的距离
                # 只用 1 的宿主计算权重，距离越近概率越高
                # 再用全局平均距离做归一化，压缩到 [0, 1]
                infect_bact  = [b for b, v in zip(known_bact, known_values) if v == 1]
                if not infect_bact:
                    continue

                infect_dists = np.array(
                    [distance_matrix.loc[bact, ib] for ib in infect_bact], dtype=float
                )
                # 相似度 = exp(-d²/scale²)，取所有已知感染宿主的加权平均相似度
                similarities = np.exp(-(infect_dists ** 2) / (distance_scale ** 2))

                # 归一化：除以自身距离为 0 时的最大可能相似度（即 1.0）
                # 取平均相似度作为概率估计
                predicted.loc[bact, phage] = float(similarities.mean())

    return predicted


# ══════════════════════════════════════════════════════════
#  5. 输出：CSV + Excel + Heatmap
# ══════════════════════════════════════════════════════════

def save_csv(
    interaction_matrix: pd.DataFrame,
    predicted_prob: pd.DataFrame,
    confidence: pd.Series,
    out_dir: str = "output",
) -> None:
    os.makedirs(out_dir, exist_ok=True)

    interaction_matrix.to_csv(f"{out_dir}/original_matrix.csv")
    predicted_prob.to_csv(f"{out_dir}/predicted_probability.csv", float_format="%.4f")
    confidence.to_frame().to_csv(f"{out_dir}/confidence.csv", float_format="%.4f")

    # 综合矩阵（原始值 + 预测值）
    combined = _build_combined(interaction_matrix, predicted_prob)
    combined.to_csv(f"{out_dir}/combined_matrix.csv")

    print(f"📄 CSV 已保存至：{out_dir}/")


def save_excel(
    interaction_matrix: pd.DataFrame,
    predicted_prob: pd.DataFrame,
    confidence: pd.Series,
    distance_matrix: pd.DataFrame,
    out_dir: str = "output",
    filename: str = "phage_host_results.xlsx",
) -> None:
    os.makedirs(out_dir, exist_ok=True)
    filepath = f"{out_dir}/{filename}"

    combined = _build_combined(interaction_matrix, predicted_prob)

    conf_df = confidence.to_frame()
    conf_df["known_count"] = interaction_matrix.notna().sum(axis=0)
    conf_df["total_bacteria"] = len(interaction_matrix)
    conf_df["unknown_count"] = conf_df["total_bacteria"] - conf_df["known_count"]

    with pd.ExcelWriter(filepath, engine="openpyxl") as writer:
        interaction_matrix.to_excel(writer, sheet_name="Original_Matrix")
        predicted_prob.round(4).to_excel(writer, sheet_name="Predicted_Probability")
        combined.to_excel(writer, sheet_name="Combined_Matrix")
        conf_df.to_excel(writer, sheet_name="Confidence")
        distance_matrix.round(4).to_excel(writer, sheet_name="Distance_Matrix")

    print(f"📊 Excel 已保存：{filepath}")


def save_heatmaps(
    interaction_matrix: pd.DataFrame,
    predicted_prob: pd.DataFrame,
    confidence: pd.Series,
    out_dir: str = "output",
) -> None:
    os.makedirs(out_dir, exist_ok=True)

    # ── Heatmap 1：原始矩阵 ────────────────────────────────
    _plot_original_heatmap(interaction_matrix, f"{out_dir}/heatmap_original.png")

    # ── Heatmap 2：预测概率矩阵（仅 NaN 位置）──────────────
    _plot_predicted_heatmap(predicted_prob, f"{out_dir}/heatmap_predicted.png")

    # ── Heatmap 3：综合矩阵（原始 + 预测合并）─────────────
    _plot_combined_heatmap(interaction_matrix, predicted_prob, confidence,
                           f"{out_dir}/heatmap_combined.png")

    print(f"🎨 Heatmap 已保存至：{out_dir}/")


def _build_combined(
    interaction_matrix: pd.DataFrame,
    predicted_prob: pd.DataFrame,
) -> pd.DataFrame:
    """构建展示用综合矩阵（字符串格式，预测值加括号）"""
    combined = interaction_matrix.copy().astype(object)
    for phage in interaction_matrix.columns:
        for bact in interaction_matrix.index:
            if pd.isna(interaction_matrix.loc[bact, phage]):
                prob = predicted_prob.loc[bact, phage]
                combined.loc[bact, phage] = (
                    f"[{prob:.3f}]" if not pd.isna(prob) else "N/A"
                )
    return combined


def _plot_original_heatmap(matrix: pd.DataFrame, savepath: str) -> None:
    n_bact  = len(matrix.index)
    n_phage = len(matrix.columns)
    fig_w   = max(8, n_phage * 0.9 + 2)
    fig_h   = max(6, n_bact  * 0.5 + 2)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # 数值矩阵用于颜色
    num = matrix.copy().astype(float)

    # 自定义颜色：0=蓝, 1=红, NaN=灰
    cmap = matplotlib.colors.ListedColormap(["#4575b4", "#d73027"])
    bounds = [-0.5, 0.5, 1.5]
    norm   = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    im = ax.imshow(num.values, cmap=cmap, norm=norm, aspect="auto")

    # NaN 位置覆灰
    nan_mask = num.isna().values
    nan_overlay = np.where(nan_mask, 0.5, np.nan)
    gray_cmap = matplotlib.colors.ListedColormap(["#cccccc"])
    ax.imshow(nan_overlay, cmap=gray_cmap, aspect="auto", vmin=0, vmax=1)

    # 标注数值
    for i in range(n_bact):
        for j in range(n_phage):
            val = matrix.iloc[i, j]
            if pd.isna(val):
                ax.text(j, i, "?", ha="center", va="center",
                        fontsize=9, color="#888888")
            else:
                ax.text(j, i, str(int(val)), ha="center", va="center",
                        fontsize=10, color="white", fontweight="bold")

    ax.set_xticks(range(n_phage));  ax.set_xticklabels(matrix.columns, rotation=45, ha="right")
    ax.set_yticks(range(n_bact));   ax.set_yticklabels(matrix.index, fontsize=8)
    ax.set_title("Original Phage-Host Interaction Matrix", fontsize=13, fontweight="bold", pad=12)

    patches = [
        mpatches.Patch(color="#d73027", label="1 = Infects"),
        mpatches.Patch(color="#4575b4", label="0 = Does not infect"),
        mpatches.Patch(color="#cccccc", label="? = No data (NaN)"),
    ]
    ax.legend(handles=patches, loc="upper right", bbox_to_anchor=(1.28, 1), fontsize=9)

    plt.tight_layout()
    plt.savefig(savepath, dpi=150, bbox_inches="tight")
    plt.close()


def _plot_predicted_heatmap(predicted_prob: pd.DataFrame, savepath: str) -> None:
    n_bact  = len(predicted_prob.index)
    n_phage = len(predicted_prob.columns)
    fig_w   = max(8, n_phage * 0.9 + 2)
    fig_h   = max(6, n_bact  * 0.5 + 2)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    cmap = sns.diverging_palette(220, 20, as_cmap=True)  # blue→white→red

    # 仅绘制有预测值的位置
    data = predicted_prob.values.astype(float)
    masked = np.ma.masked_invalid(data)

    im = ax.imshow(masked, cmap=cmap, vmin=0, vmax=1, aspect="auto")

    # NaN（既无原始也无预测）覆灰
    nan_mask    = predicted_prob.isna().values
    nan_overlay = np.where(nan_mask, 0.5, np.nan)
    ax.imshow(nan_overlay, cmap=matplotlib.colors.ListedColormap(["#eeeeee"]),
              aspect="auto", vmin=0, vmax=1)

    # 标注概率值
    for i in range(n_bact):
        for j in range(n_phage):
            val = predicted_prob.iloc[i, j]
            if not pd.isna(val):
                color = "white" if (val > 0.65 or val < 0.35) else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=8, color=color)
            else:
                ax.text(j, i, "—", ha="center", va="center",
                        fontsize=9, color="#aaaaaa")

    plt.colorbar(im, ax=ax, label="Predicted Infection Probability", shrink=0.8)
    ax.set_xticks(range(n_phage));  ax.set_xticklabels(predicted_prob.columns, rotation=45, ha="right")
    ax.set_yticks(range(n_bact));   ax.set_yticklabels(predicted_prob.index, fontsize=8)
    ax.set_title("Predicted Infection Probability (NaN positions only)",
                 fontsize=13, fontweight="bold", pad=12)

    plt.tight_layout()
    plt.savefig(savepath, dpi=150, bbox_inches="tight")
    plt.close()


def _plot_combined_heatmap(
    interaction_matrix: pd.DataFrame,
    predicted_prob: pd.DataFrame,
    confidence: pd.Series,
    savepath: str,
) -> None:
    """
    综合热图：
      - 已知 1   → 深红
      - 已知 0   → 深蓝
      - 预测高概率 → 浅红
      - 预测低概率 → 浅蓝
      - 无数据     → 灰
    底部附 confidence bar。
    """
    n_bact  = len(interaction_matrix.index)
    n_phage = len(interaction_matrix.columns)
    fig_w   = max(10, n_phage * 0.9 + 3)
    fig_h   = max(7,  n_bact  * 0.5 + 3)

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs  = fig.add_gridspec(2, 1, height_ratios=[n_bact, 2], hspace=0.35)
    ax_main = fig.add_subplot(gs[0])
    ax_conf = fig.add_subplot(gs[1])

    # 构建颜色矩阵（RGBA）
    rgba = np.ones((n_bact, n_phage, 4))
    for i, bact in enumerate(interaction_matrix.index):
        for j, phage in enumerate(interaction_matrix.columns):
            known = interaction_matrix.iloc[i, j]
            pred  = predicted_prob.iloc[i, j]
            if not pd.isna(known):
                if known == 1:
                    rgba[i, j] = [0.84, 0.19, 0.15, 1.0]   # 深红
                else:
                    rgba[i, j] = [0.27, 0.46, 0.71, 1.0]   # 深蓝
            elif not pd.isna(pred):
                r = pred
                # 概率 → 红蓝渐变（浅色表示预测）
                rgba[i, j] = [0.9 + 0.05 * r, 0.9 - 0.5 * r, 0.9 - 0.6 * r, 0.85]
            else:
                rgba[i, j] = [0.85, 0.85, 0.85, 1.0]       # 灰

    ax_main.imshow(rgba, aspect="auto")

    # 标注
    for i in range(n_bact):
        for j in range(n_phage):
            known = interaction_matrix.iloc[i, j]
            pred  = predicted_prob.iloc[i, j]
            if not pd.isna(known):
                ax_main.text(j, i, str(int(known)), ha="center", va="center",
                             fontsize=10, color="white", fontweight="bold")
            elif not pd.isna(pred):
                ax_main.text(j, i, f"{pred:.2f}", ha="center", va="center",
                             fontsize=7.5, color="#333333")
            else:
                ax_main.text(j, i, "?", ha="center", va="center",
                             fontsize=9, color="#aaaaaa")

    ax_main.set_xticks(range(n_phage))
    ax_main.set_xticklabels(interaction_matrix.columns, rotation=45, ha="right")
    ax_main.set_yticks(range(n_bact))
    ax_main.set_yticklabels(interaction_matrix.index, fontsize=8)
    ax_main.set_title("Combined Matrix: Known (bold) + Predicted (decimal)",
                      fontsize=12, fontweight="bold", pad=10)

    # ── Confidence bar ────────────────────────
    conf_vals = [confidence[p] for p in interaction_matrix.columns]
    colors    = ["#2ecc71" if c >= 0.7 else "#f39c12" if c >= 0.4 else "#e74c3c"
                 for c in conf_vals]
    bars = ax_conf.bar(range(n_phage), conf_vals, color=colors, edgecolor="white", linewidth=0.8)
    ax_conf.set_xticks(range(n_phage))
    ax_conf.set_xticklabels(interaction_matrix.columns, rotation=45, ha="right", fontsize=8)
    ax_conf.set_ylim(0, 1.15)
    ax_conf.set_ylabel("Confidence", fontsize=9)
    ax_conf.set_title("Phage Confidence Score  (green ≥0.7 | orange ≥0.4 | red <0.4)",
                      fontsize=9)
    ax_conf.axhline(0.7, color="#2ecc71", linestyle="--", linewidth=0.8, alpha=0.6)
    ax_conf.axhline(0.4, color="#e74c3c", linestyle="--", linewidth=0.8, alpha=0.6)
    for bar, val in zip(bars, conf_vals):
        ax_conf.text(bar.get_x() + bar.get_width() / 2, val + 0.02,
                     f"{val:.0%}", ha="center", va="bottom", fontsize=7.5)

    # 图例
    patches = [
        mpatches.Patch(color="#d73027", label="Known: 1 (infects)"),
        mpatches.Patch(color="#4575b4", label="Known: 0 (no infect)"),
        mpatches.Patch(color="#f5c0b0", label="Predicted: high prob"),
        mpatches.Patch(color="#b0c8e8", label="Predicted: low prob"),
        mpatches.Patch(color="#cccccc", label="No data"),
    ]
    ax_main.legend(handles=patches, loc="upper right",
                   bbox_to_anchor=(1.32, 1.0), fontsize=8, title="Legend")

    plt.savefig(savepath, dpi=150, bbox_inches="tight")
    plt.close()


# ══════════════════════════════════════════════════════════
#  6. 终端报告
# ══════════════════════════════════════════════════════════

def print_report(
    interaction_matrix: pd.DataFrame,
    predicted_prob: pd.DataFrame,
    confidence: pd.Series,
    distance_matrix: pd.DataFrame,
) -> None:
    sep = "─" * 72

    print(f"\n{'═'*72}")
    print("   PHAGE-HOST INTERACTION ANALYSIS REPORT")
    print(f"{'═'*72}\n")

    print("【1】Original Interaction Matrix  (1=infects | 0=no | NaN=unknown)")
    print(sep)
    print(interaction_matrix.to_string())

    print(f"\n【2】Evolutionary Distance Matrix  (NCBI Taxonomy-based)")
    print(sep)
    print(distance_matrix.round(3).to_string())

    print(f"\n【3】Confidence Score per Phage")
    print(sep)
    conf_df = confidence.to_frame()
    conf_df["known"]   = interaction_matrix.notna().sum(axis=0)
    conf_df["total"]   = len(interaction_matrix)
    conf_df["unknown"] = conf_df["total"] - conf_df["known"]
    print(conf_df.to_string())

    print(f"\n【4】Predicted Infection Probability  (NaN positions only)")
    print(sep)
    print(predicted_prob.to_string(float_format="{:.4f}".format))

    print(f"\n【5】Combined Matrix  (known values | predicted in [brackets])")
    print(sep)
    combined = _build_combined(interaction_matrix, predicted_prob)
    print(combined.to_string())

    print(f"\n{'═'*72}\n")


# ══════════════════════════════════════════════════════════
#  7. 主流程
# ══════════════════════════════════════════════════════════

def run_analysis(
    interaction_matrix: pd.DataFrame,
    distance_scale: float = 0.5,
    out_dir: str = "output",
    taxonomy_cache_file: str = "taxonomy_cache.csv",
    use_ncbi: bool = True,
) -> dict:
    """
    完整分析流程。

    参数
    ────
    interaction_matrix   : phage-host 矩阵（行=bacteria, 列=phage）
    distance_scale       : 高斯核宽度，控制距离权重衰减（默认 0.5）
    out_dir              : 输出目录
    taxonomy_cache_file  : NCBI taxonomy 本地缓存路径
    use_ncbi             : True=联网查询 NCBI；False=使用模拟距离（离线测试用）

    返回
    ────
    dict 含 interaction_matrix, distance_matrix, predicted_prob, confidence
    """
    bacteria_list = interaction_matrix.index.tolist()

    # ── Step 1: 进化距离矩阵 ──────────────────
    if use_ncbi:
        taxonomy_cache = build_taxonomy_cache(bacteria_list, taxonomy_cache_file)
        distance_matrix = build_distance_matrix(bacteria_list, taxonomy_cache)
    else:
        print("⚠️  离线模式：使用随机模拟距离矩阵（仅用于测试）")
        np.random.seed(42)
        n = len(bacteria_list)
        rand = np.random.uniform(0.1, 0.9, (n, n))
        rand = (rand + rand.T) / 2
        np.fill_diagonal(rand, 0)
        distance_matrix = pd.DataFrame(rand, index=bacteria_list, columns=bacteria_list)

    # ── Step 2: Confidence ────────────────────
    confidence = compute_confidence(interaction_matrix)

    # ── Step 3: 预测 ──────────────────────────
    print("🧮 计算 predicted infection probability...")
    predicted_prob = predict_infection_probability(
        interaction_matrix, distance_matrix, distance_scale
    )

    # ── Step 4: 报告 ──────────────────────────
    print_report(interaction_matrix, predicted_prob, confidence, distance_matrix)

    # ── Step 5: 保存 ──────────────────────────
    save_csv(interaction_matrix, predicted_prob, confidence, out_dir)
    save_excel(interaction_matrix, predicted_prob, confidence, distance_matrix, out_dir)
    save_heatmaps(interaction_matrix, predicted_prob, confidence, out_dir)

    return {
        "interaction_matrix": interaction_matrix,
        "distance_matrix":    distance_matrix,
        "predicted_prob":     predicted_prob,
        "confidence":         confidence,
    }


# ══════════════════════════════════════════════════════════
#  入口
# ══════════════════════════════════════════════════════════

if __name__ == "__main__":

    # ── 选项 A：手动输入示例数据 ─────────────────────────
    print("=" * 55)
    print("  模式：手动输入示例数据 + 离线模拟距离（测试用）")
    print("=" * 55)
    # matrix = create_example_data()

    # ── 选项 B：从文件读取（取消注释使用）────────────────
    # matrix = load_from_file("/Users/mac/Documents/Mudd/iGEM/dataCleaningCode/phage_host_matrix_with_ids.csv")
    # matrix = matrix[~matrix.index.duplicated(keep='first')]  # 去重
    df_raw = pd.read_csv("/Users/mac/Documents/Mudd/iGEM/dataCleaningCode/phage_host_matrix_with_ids.csv")
    df_raw = df_raw[df_raw['Host_Name'].notna() & ~df_raw['Host_Name'].isin(['Unknown', 'ERROR', 'No Complete Genome Found', 'N/A'])]
    df_raw['Affinity'] = pd.to_numeric(df_raw['Affinity'], errors='coerce')
    df_raw = df_raw.dropna(subset=['Affinity'])
    matrix = df_raw.pivot_table(index='Host_Name', columns='Phage', values='Affinity', aggfunc='max')
    matrix = matrix[~matrix.index.duplicated(keep='first')]
    print(f"✅ 转换完成：{matrix.shape[0]} bacteria × {matrix.shape[1]} phages")

    results = run_analysis(
        interaction_matrix=matrix,
        distance_scale=0.5,       # 可调节：越小越依赖近邻
        out_dir="/Users/mac/Documents/Mudd/iGEM/dataCleaningCode/results_26-4-12",
        use_ncbi=False,           # 改为 True 则真实查询 NCBI
    )

    print("\n✅ 分析完成！输出文件：")
    print("   output/original_matrix.csv")
    print("   output/predicted_probability.csv")
    print("   output/combined_matrix.csv")
    print("   output/confidence.csv")
    print("   output/phage_host_results.xlsx  （多 sheet）")
    print("   output/heatmap_original.png")
    print("   output/heatmap_predicted.png")
    print("   output/heatmap_combined.png")