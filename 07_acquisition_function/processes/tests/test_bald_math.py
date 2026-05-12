"""
test_bald_math.py — Unit tests for BALD scoring and batch selection math.
BALD 打分和批次选择数学的单元测试。

Tests:
  1. BALD score is monotone with epistemic_std.
  2. High-epistemic candidate ranked first AND lowest-epistemic NOT in top-K.
  3. Random-control pick is not in BALD top-K.
  4. Measured pairs are excluded from selection.
  5. Cycle 0 empty-measured path works.
  6. Exact batch size (n_bald + n_random).
  7. Non-negative epistemic_std enforced.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PROCESSES_DIR = Path(__file__).resolve().parents[1]
if str(PROCESSES_DIR) not in sys.path:
    sys.path.insert(0, str(PROCESSES_DIR))

from bald import bald_score, select_batch


# ---------------------------------------------------------------------------
# Fixtures / 测试夹具
# ---------------------------------------------------------------------------

def make_pool(n: int = 10, seed: int = 0) -> pd.DataFrame:
    """
    Synthetic candidate pool DataFrame with known epistemic_std values.
    生成已知 epistemic_std 值的合成候选池 DataFrame。
    """
    rng = np.random.default_rng(seed)
    epistemic = rng.uniform(0.01, 1.0, n).astype(np.float32)
    df = pd.DataFrame({
        "rbp_id":          [f"rbp_{i:02d}" for i in range(n)],
        "receptor_id":     ["rec_00"] * n,
        "predicted_score": rng.uniform(3.0, 7.0, n).tolist(),
        "std":             (epistemic + rng.uniform(0.0, 0.3, n)).tolist(),
        "epistemic_std":   epistemic.tolist(),
        "bald_score":      bald_score(epistemic).tolist(),
    })
    return df


# ---------------------------------------------------------------------------
# 1. bald_score is monotone with epistemic_std
# bald_score 与 epistemic_std 单调一致
# ---------------------------------------------------------------------------

def test_bald_score_monotone():
    """bald_score preserves rank ordering from epistemic_std."""
    epist = np.array([0.1, 0.5, 0.3, 0.9, 0.05], dtype=np.float32)
    scores = bald_score(epist)
    assert np.allclose(scores, epist), "bald_score must equal epistemic_std"
    assert np.all(np.argsort(scores) == np.argsort(epist)), "Rank order must be preserved"


def test_bald_score_rejects_negative():
    """bald_score raises ValueError on negative epistemic_std."""
    with pytest.raises(ValueError, match="non-negative"):
        bald_score(np.array([-0.1, 0.5], dtype=np.float32))


# ---------------------------------------------------------------------------
# 2. Highest-epistemic is ranked first; lowest is NOT in top-K
# 最高 epistemic 排第一；最低的 NOT 在 top-K 内
# ---------------------------------------------------------------------------

def test_top_pick_has_highest_bald():
    """
    The single highest-epistemic_std candidate must be the first BALD pick.
    Also verify the lowest-epistemic candidate is not selected as a BALD pick
    (guards against trivially-passing one-sided assertions).
    """
    # Construct pool where one candidate has clearly the highest epistemic_std
    # 构造一个候选池，其中一个候选有明显最高的 epistemic_std
    n = 8
    rng = np.random.default_rng(1)
    epistemic = rng.uniform(0.1, 0.5, n).astype(np.float32)
    epistemic[3] = 2.0   # clear winner / 明显的最高值
    epistemic[7] = 0.001  # clear loser / 明显的最低值

    df = pd.DataFrame({
        "rbp_id":      [f"rbp_{i:02d}" for i in range(n)],
        "receptor_id": ["rec_00"] * n,
        "predicted_score": rng.uniform(3.0, 7.0, n).tolist(),
        "std":         (epistemic + 0.1).tolist(),
        "epistemic_std": epistemic.tolist(),
        "bald_score":  bald_score(epistemic).tolist(),
    })

    recs, _, _ = select_batch(df, n_bald=3, n_random=1, seed=99)

    # Top BALD pick must be rbp_03 (index 3, highest epistemic)
    bald_picks = recs[recs["selection_reason"].str.startswith("bald_top")]
    assert bald_picks.iloc[0]["rbp_id"] == "rbp_03", (
        f"Expected rbp_03 as top BALD pick, got {bald_picks.iloc[0]['rbp_id']}"
    )

    # rbp_07 (lowest epistemic) must NOT appear as a BALD pick
    bald_rbp_ids = set(bald_picks["rbp_id"])
    assert "rbp_07" not in bald_rbp_ids, (
        "Lowest-epistemic candidate rbp_07 should not be in BALD top-K picks."
    )


# ---------------------------------------------------------------------------
# 3. Random-control pick is not in BALD top-K
# 随机对照不在 BALD top-K 内
# ---------------------------------------------------------------------------

def test_random_control_not_in_bald_top():
    """Random-control pick must not duplicate a BALD top pick."""
    df = make_pool(n=10)
    recs, _, _ = select_batch(df, n_bald=4, n_random=1, seed=42)

    bald_picks = recs[recs["selection_reason"].str.startswith("bald_top")]
    random_pick = recs[recs["selection_reason"] == "random_control"]

    bald_pairs = set(zip(bald_picks["rbp_id"], bald_picks["receptor_id"]))
    for _, row in random_pick.iterrows():
        pair = (row["rbp_id"], row["receptor_id"])
        assert pair not in bald_pairs, (
            f"Random control {pair} duplicates a BALD top pick."
        )


# ---------------------------------------------------------------------------
# 4. Measured pairs are excluded
# 已测量的对被排除
# ---------------------------------------------------------------------------

def test_measured_pairs_excluded():
    """Pairs in measured_pairs must not appear in any output."""
    df = make_pool(n=10)

    # Mark the top-2 bald candidates as already measured
    # 将前 2 名 BALD 候选标记为已测量
    sorted_df = df.sort_values("bald_score", ascending=False)
    measured = set(zip(sorted_df["rbp_id"].iloc[:2], sorted_df["receptor_id"].iloc[:2]))

    recs, replay, backup = select_batch(df, n_bald=3, n_random=1, measured_pairs=measured, seed=42)

    all_outputs = pd.concat([recs, replay, backup])
    for _, row in all_outputs.iterrows():
        pair = (row["rbp_id"], row["receptor_id"])
        assert pair not in measured, (
            f"Measured pair {pair} appeared in BALD output — should have been excluded."
        )


# ---------------------------------------------------------------------------
# 5. Cycle 0 empty-measured path
# Cycle 0 空测量集路径
# ---------------------------------------------------------------------------

def test_cycle0_empty_measured():
    """
    With measured_pairs=None (Cycle 0), the full pool is eligible.
    All n_bald + n_random candidates come from the full pool.
    """
    df = make_pool(n=12)
    recs, replay, backup = select_batch(df, n_bald=4, n_random=1, measured_pairs=None, seed=0)

    assert len(recs) == 5, f"Expected 5 recommendations, got {len(recs)}"
    assert len(replay) == 5, f"Expected 5 replay picks, got {len(replay)}"

    # All recommended pairs must be in the original pool
    pool_pairs = set(zip(df["rbp_id"], df["receptor_id"]))
    for _, row in recs.iterrows():
        assert (row["rbp_id"], row["receptor_id"]) in pool_pairs


# ---------------------------------------------------------------------------
# 6. Exact batch size
# 精确批次大小
# ---------------------------------------------------------------------------

def test_batch_size():
    """Total recommendations = n_bald + n_random exactly."""
    df = make_pool(n=20)
    for n_bald, n_random in [(4, 1), (3, 2), (5, 0), (1, 1)]:
        recs, replay, _ = select_batch(df, n_bald=n_bald, n_random=n_random, seed=7)
        assert len(recs) == n_bald + n_random, (
            f"n_bald={n_bald}, n_random={n_random}: expected {n_bald+n_random} recs, got {len(recs)}"
        )
        assert len(replay) == n_bald + n_random, (
            f"random_replay size mismatch: expected {n_bald+n_random}, got {len(replay)}"
        )


# ---------------------------------------------------------------------------
# 7. Selection reason labels are correct
# selection_reason 标签正确
# ---------------------------------------------------------------------------

def test_selection_reason_labels():
    """recommendations must have correct selection_reason values."""
    df = make_pool(n=10)
    recs, replay, backup = select_batch(df, n_bald=3, n_random=1, seed=0)

    bald_reasons = [f"bald_top_{i+1}" for i in range(3)]
    rec_reasons = list(recs["selection_reason"])
    for expected in bald_reasons:
        assert expected in rec_reasons, f"Missing reason '{expected}' in recommendations."

    assert "random_control" in rec_reasons, "Missing 'random_control' in recommendations."

    replay_reasons = set(replay["selection_reason"])
    assert replay_reasons == {"random_replay"}, (
        f"random_replay has unexpected reasons: {replay_reasons}"
    )

    backup_reasons = set(backup["selection_reason"])
    assert all(r.startswith("safe_pick_bald_top_") for r in backup_reasons), (
        f"safe_pick_backup has unexpected reasons: {backup_reasons}"
    )
