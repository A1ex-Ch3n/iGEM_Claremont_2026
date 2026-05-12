"""
test_schema.py — Schema validation tests for Module 07 output CSVs.
Module 07 输出 CSV 的列结构验证测试。

Tests:
  1. recommendations.csv has all required columns.
  2. random_replay.csv has all required columns.
  3. safe_pick_backup.csv has all required columns.
  4. bald_score and epistemic_std are non-negative in all outputs.
  5. selection_reason values are from the expected vocabulary.
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


def make_pool(n: int = 15, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    epistemic = rng.uniform(0.01, 1.0, n).astype(np.float32)
    df = pd.DataFrame({
        "rbp_id":          [f"rbp_{i:02d}" for i in range(n)],
        "receptor_id":     [f"rec_{i % 4:02d}" for i in range(n)],
        "predicted_score": rng.uniform(3.0, 7.0, n).tolist(),
        "std":             (epistemic + rng.uniform(0.0, 0.3, n)).tolist(),
        "epistemic_std":   epistemic.tolist(),
        "bald_score":      bald_score(epistemic).tolist(),
    })
    return df


REQUIRED_COLS_RECS = {
    "rbp_id", "receptor_id", "predicted_score", "std", "epistemic_std",
    "bald_score", "bald_rank", "selection_reason",
}

REQUIRED_COLS_REPLAY = {
    "rbp_id", "receptor_id", "bald_score", "bald_rank", "selection_reason",
}


def test_recommendations_schema():
    """recommendations DataFrame must contain all required columns."""
    df = make_pool()
    recs, _, _ = select_batch(df, n_bald=4, n_random=1)
    missing = REQUIRED_COLS_RECS - set(recs.columns)
    assert not missing, f"recommendations.csv missing columns: {missing}"


def test_random_replay_schema():
    """random_replay DataFrame must contain all required columns."""
    df = make_pool()
    _, replay, _ = select_batch(df, n_bald=4, n_random=1)
    missing = REQUIRED_COLS_REPLAY - set(replay.columns)
    assert not missing, f"random_replay.csv missing columns: {missing}"


def test_safe_pick_backup_schema():
    """safe_pick_backup DataFrame must contain all required columns."""
    df = make_pool()
    _, _, backup = select_batch(df, n_bald=4, n_random=1)
    missing = REQUIRED_COLS_REPLAY - set(backup.columns)
    assert not missing, f"safe_pick_backup.csv missing columns: {missing}"


def test_bald_score_non_negative_in_outputs():
    """bald_score in all outputs must be >= 0."""
    df = make_pool()
    recs, replay, backup = select_batch(df, n_bald=3, n_random=1)
    for name, frame in [("recommendations", recs), ("random_replay", replay), ("safe_pick_backup", backup)]:
        assert (frame["bald_score"] >= 0).all(), (
            f"Negative bald_score found in {name}"
        )


def test_selection_reason_vocabulary():
    """selection_reason values must be from the expected vocabulary."""
    df = make_pool(n=15)
    recs, replay, backup = select_batch(df, n_bald=4, n_random=1)

    # Recommendations: bald_top_<N> or random_control
    # 推荐列表：bald_top_<N> 或 random_control
    for reason in recs["selection_reason"]:
        assert reason.startswith("bald_top_") or reason == "random_control", (
            f"Unexpected selection_reason in recommendations: '{reason}'"
        )

    # random_replay: all must be random_replay
    assert all(r == "random_replay" for r in replay["selection_reason"])

    # safe_pick_backup: all must start with safe_pick_bald_top_
    assert all(r.startswith("safe_pick_bald_top_") for r in backup["selection_reason"])


def test_no_duplicate_recommendations():
    """No (rbp_id, receptor_id) pair should appear twice in recommendations."""
    df = make_pool(n=15)
    recs, _, _ = select_batch(df, n_bald=4, n_random=1)
    pairs = list(zip(recs["rbp_id"], recs["receptor_id"]))
    assert len(pairs) == len(set(pairs)), "Duplicate pairs in recommendations."
