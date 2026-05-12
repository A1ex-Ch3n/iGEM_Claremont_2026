"""
bald.py — BALD acquisition function for deep ensemble regression.
BALD 采集函数：基于深度集成的贝叶斯主动学习（回归版本）。

BALD(x) = I(y; θ | x, D)
         = H[y | x, D] − E_{p(θ|D)}[ H[y | x, θ] ]
         ≈ Var_k[ μ_k(x) ]  (epistemic variance of ensemble member means)

For a Gaussian deep ensemble with K members:
  - Each member predicts (μ_k, σ_k).
  - Epistemic variance = Var_k[μ_k] — the part BALD targets.
  - Aleatoric variance = E_k[σ_k²] — irreducible noise, not targeted.

We use epistemic_std (square root of epistemic variance) as the BALD score
because it has the same units as the prediction and is monotone with rank.

Reference / 参考文献:
  Houlsby et al. (2011) "Bayesian Active Learning for Classification and
  Preference Learning." arXiv:1112.5745. (Original GPC formulation.)
  Lakshminarayanan et al. (2017) "Simple and Scalable Predictive Uncertainty
  Estimation Using Deep Ensembles." NeurIPS 2017. (Ensemble architecture.)
"""

from __future__ import annotations

import random as _random
from typing import Optional, Set, Tuple

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# BALD scoring
# BALD 打分
# ---------------------------------------------------------------------------

def bald_score(epistemic_std: np.ndarray) -> np.ndarray:
    """
    BALD score = epistemic_std (monotone rank-equivalent to epistemic variance).
    BALD 分数 = epistemic_std（与 epistemic 方差单调等价，可直接用于排名）。

    Args:
        epistemic_std: (N,) array from DeepEnsemble.predict(return_epistemic=True).

    Returns:
        (N,) BALD scores — higher is more informative.
        (N,) BALD 分数——越高越有信息价值。
    """
    if np.any(epistemic_std < 0):
        raise ValueError("epistemic_std must be non-negative; got negatives.")
    return epistemic_std.copy()


# ---------------------------------------------------------------------------
# Batch selection
# 批次选择
# ---------------------------------------------------------------------------

def select_batch(
    df: pd.DataFrame,
    bald_col: str = "bald_score",
    n_bald: int = 4,
    n_random: int = 1,
    measured_pairs: Optional[Set[Tuple[str, str]]] = None,
    rbp_col: str = "rbp_id",
    rec_col: str = "receptor_id",
    seed: int = 42,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Select n_bald top-BALD + n_random random from unmeasured candidates.
    从未测量的候选对中选出 n_bald 个 BALD 最优 + n_random 个随机对照。

    Args:
        df:             Candidate pool DataFrame — must include bald_col, rbp_col, rec_col.
                        候选池 DataFrame，必须包含 bald_col、rbp_col、rec_col 列。
        bald_col:       Column name with BALD scores (descending = more informative).
        n_bald:         Number of top-BALD picks. / BALD 最优选取数量。
        n_random:       Number of random-control picks. / 随机对照选取数量。
        measured_pairs: Set of (rbp_id, receptor_id) already measured — excluded.
                        已测量的 (rbp_id, receptor_id) 集合，排除在外。
        rbp_col:        Column name for RBP identifier.
        rec_col:        Column name for receptor identifier.
        seed:           Random seed for reproducible random pick.

    Returns:
        recommendations: Combined picks with 'selection_reason' column.
                         合并后的推荐，含 'selection_reason' 列。
        random_replay:   What pure-random selection would have picked (same batch size).
                         纯随机选择的结果（供回顾性对比）。
        safe_pick_backup: Extended top-BALD pool for manual PI/team selection.
                          扩展的 BALD 最优列表，供 PI/团队手动选择（SLA 未达时使用）。
    """
    if measured_pairs is None:
        measured_pairs = set()

    # -----------------------------------------------------------------------
    # Filter out measured pairs / 过滤掉已测量的对
    # -----------------------------------------------------------------------
    pair_key = list(zip(df[rbp_col], df[rec_col]))
    mask_unmeasured = pd.Series(
        [pair not in measured_pairs for pair in pair_key], index=df.index
    )
    pool = df[mask_unmeasured].copy()

    n_needed = n_bald + n_random
    if len(pool) < n_needed:
        raise ValueError(
            f"Candidate pool has only {len(pool)} unmeasured pairs, "
            f"need at least {n_needed} (n_bald={n_bald} + n_random={n_random})."
        )

    # -----------------------------------------------------------------------
    # Rank by BALD score (descending) / 按 BALD 分数降序排名
    # -----------------------------------------------------------------------
    ranked = pool.sort_values(bald_col, ascending=False).reset_index(drop=True)
    ranked["bald_rank"] = ranked.index + 1

    # -----------------------------------------------------------------------
    # Top-BALD picks / 顶部 BALD 选取
    # -----------------------------------------------------------------------
    bald_picks = ranked.iloc[:n_bald].copy()
    bald_picks["selection_reason"] = [
        f"bald_top_{i+1}" for i in range(len(bald_picks))
    ]

    # -----------------------------------------------------------------------
    # Random control: sample from remaining unmeasured (not in BALD top-K)
    # 随机对照：从 BALD 未选中的未测量候选中抽取
    # -----------------------------------------------------------------------
    bald_indices = set(bald_picks.index)
    remaining = ranked[~ranked.index.isin(bald_indices)]

    rng = _random.Random(seed)
    random_idx = rng.sample(list(remaining.index), k=n_random)
    random_picks = ranked.loc[random_idx].copy()
    random_picks["selection_reason"] = "random_control"

    # -----------------------------------------------------------------------
    # Combine into recommendations / 合并成推荐列表
    # -----------------------------------------------------------------------
    recommendations = pd.concat([bald_picks, random_picks], ignore_index=True)

    # -----------------------------------------------------------------------
    # Random replay: what would pure-random have picked? (for retrospective)
    # 随机重放：纯随机策略会选什么？（供 Cycle 复盘用）
    # -----------------------------------------------------------------------
    rng_replay = _random.Random(seed + 1)  # different seed for independence
    replay_idx = rng_replay.sample(list(ranked.index), k=n_needed)
    random_replay = ranked.loc[replay_idx].copy()
    random_replay["selection_reason"] = "random_replay"

    # -----------------------------------------------------------------------
    # Safe-pick backup: top-2*(n_bald+n_random) BALD picks for PI manual cut
    # 安全备选：BALD 最优前 2*(n_bald+n_random) 个，供 PI 手动挑选（SLA 未达时）
    # -----------------------------------------------------------------------
    n_backup = min(2 * n_needed, len(ranked))
    safe_pick_backup = ranked.iloc[:n_backup].copy()
    safe_pick_backup["selection_reason"] = [
        f"safe_pick_bald_top_{i+1}" for i in range(len(safe_pick_backup))
    ]

    return recommendations, random_replay, safe_pick_backup
