"""
test_integration.py — End-to-end integration tests for run_bald.py.
run_bald.py 端到端集成测试。

Tests:
  1. Cycle 0 (no measured_csv): runs end-to-end, all 3 output CSVs written.
  2. Cycle 1 (with measured_csv): measured pairs excluded, outputs correct.
  3. Auto-detection of predictions_csv from Module 06 outputs.
  4. Fallback: missing epistemic_std falls back to total std.
  5. run_meta.json written with correct fields.
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PROCESSES_DIR = Path(__file__).resolve().parents[1]
if str(PROCESSES_DIR) not in sys.path:
    sys.path.insert(0, str(PROCESSES_DIR))

from run_bald import main as run_bald_main


# ---------------------------------------------------------------------------
# Helpers / 工具函数
# ---------------------------------------------------------------------------

def make_predictions_csv(path: Path, n: int = 20, seed: int = 0) -> Path:
    """Write a synthetic predictions.csv to path. 写入合成 predictions.csv。"""
    rng = np.random.default_rng(seed)
    n_rbp = 5
    n_rec = 4
    rbp_ids = [f"rbp_{i:02d}" for i in range(n_rbp) for _ in range(n_rec)]
    rec_ids = [f"rec_{j:02d}" for _ in range(n_rbp) for j in range(n_rec)]
    n_pairs = n_rbp * n_rec

    epistemic = rng.uniform(0.05, 0.5, n_pairs).astype(np.float32)
    aleatoric = rng.uniform(0.0, 0.3, n_pairs).astype(np.float32)
    total_std = np.sqrt(epistemic**2 + aleatoric**2)
    means = rng.uniform(3.0, 7.0, n_pairs).astype(np.float32)

    df = pd.DataFrame({
        "rbp_id":          rbp_ids,
        "receptor_id":     rec_ids,
        "predicted_score": means.tolist(),
        "std":             total_std.tolist(),
        "epistemic_std":   epistemic.tolist(),
        "lower_95":        (means - 1.96 * total_std).tolist(),
        "upper_95":        (means + 1.96 * total_std).tolist(),
        "model_version":   "test_cycle_0",
    })
    df.to_csv(path, index=False)
    return path


def make_measured_csv(path: Path, pairs: list) -> Path:
    """Write a measured-pairs CSV. 写入已测量 pair 的 CSV。"""
    df = pd.DataFrame(pairs, columns=["rbp_id", "receptor_id"])
    df.to_csv(path, index=False)
    return path


# ---------------------------------------------------------------------------
# 1. Cycle 0 end-to-end (no measured_csv)
# Cycle 0 端到端（无 measured_csv）
# ---------------------------------------------------------------------------

def test_cycle0_end_to_end(tmp_path):
    """
    run_bald.py with no measured_csv produces 3 output CSVs + run_meta.json.
    无 measured_csv 时，run_bald 生成 3 个 CSV + run_meta.json。
    """
    pred_csv = make_predictions_csv(tmp_path / "predictions.csv")
    out_base = tmp_path / "07_outputs"

    # Patch MODULE_DIR to write under tmp_path
    # 将 MODULE_DIR 临时重定向到 tmp_path
    import run_bald as rb
    orig_module_dir = rb.MODULE_DIR
    rb.MODULE_DIR = out_base

    try:
        run_bald_main(
            cycle=1,
            predictions_csv=str(pred_csv),
            measured_csv=None,
            n_bald=4,
            n_random=1,
            seed=42,
        )
    finally:
        rb.MODULE_DIR = orig_module_dir

    cycle_dir = out_base / "outputs" / "cycle_1"
    assert (cycle_dir / "recommendations.csv").exists()
    assert (cycle_dir / "random_replay.csv").exists()
    assert (cycle_dir / "safe_pick_backup.csv").exists()
    assert (cycle_dir / "run_meta.json").exists()

    recs = pd.read_csv(cycle_dir / "recommendations.csv")
    assert len(recs) == 5  # 4 BALD + 1 random

    replay = pd.read_csv(cycle_dir / "random_replay.csv")
    assert len(replay) == 5


# ---------------------------------------------------------------------------
# 2. Cycle 1 with measured_csv
# Cycle 1 含 measured_csv
# ---------------------------------------------------------------------------

def test_cycle1_measured_excluded(tmp_path):
    """Measured pairs do not appear in any output when measured_csv is provided."""
    pred_csv = make_predictions_csv(tmp_path / "predictions.csv")

    # Mark the first 2 pairs as already measured
    # 标记前 2 个 pair 已测量
    measured_pairs = [("rbp_00", "rec_00"), ("rbp_00", "rec_01")]
    meas_csv = make_measured_csv(tmp_path / "measured.csv", measured_pairs)

    import run_bald as rb
    orig_module_dir = rb.MODULE_DIR
    out_base = tmp_path / "07_outputs"
    rb.MODULE_DIR = out_base

    try:
        run_bald_main(
            cycle=1,
            predictions_csv=str(pred_csv),
            measured_csv=str(meas_csv),
            n_bald=3,
            n_random=1,
            seed=0,
        )
    finally:
        rb.MODULE_DIR = orig_module_dir

    cycle_dir = out_base / "outputs" / "cycle_1"
    recs = pd.read_csv(cycle_dir / "recommendations.csv")
    replay = pd.read_csv(cycle_dir / "random_replay.csv")
    backup = pd.read_csv(cycle_dir / "safe_pick_backup.csv")

    all_out = pd.concat([recs, replay, backup])
    measured_set = set(measured_pairs)
    for _, row in all_out.iterrows():
        pair = (row["rbp_id"], row["receptor_id"])
        assert pair not in measured_set, (
            f"Measured pair {pair} appeared in output — should be excluded."
        )


# ---------------------------------------------------------------------------
# 3. Fallback: missing epistemic_std uses total std
# 回退：缺少 epistemic_std 时使用 total std
# ---------------------------------------------------------------------------

def test_fallback_to_total_std(tmp_path, caplog):
    """If predictions.csv lacks epistemic_std, falls back to total std with a warning."""
    pred_csv = make_predictions_csv(tmp_path / "predictions.csv")

    # Remove epistemic_std column / 删除 epistemic_std 列
    df = pd.read_csv(pred_csv)
    df = df.drop(columns=["epistemic_std"])
    df.to_csv(pred_csv, index=False)

    import run_bald as rb
    import logging
    orig_module_dir = rb.MODULE_DIR
    out_base = tmp_path / "07_outputs"
    rb.MODULE_DIR = out_base

    try:
        with caplog.at_level(logging.WARNING):
            run_bald_main(
                cycle=1,
                predictions_csv=str(pred_csv),
                measured_csv=None,
                n_bald=4,
                n_random=1,
                seed=42,
            )
    finally:
        rb.MODULE_DIR = orig_module_dir

    assert "falling back to total std" in caplog.text.lower(), (
        "Expected a warning about falling back to total std."
    )

    recs = pd.read_csv(out_base / "outputs" / "cycle_1" / "recommendations.csv")
    assert len(recs) == 5


# ---------------------------------------------------------------------------
# 4. run_meta.json has required fields
# run_meta.json 包含必要字段
# ---------------------------------------------------------------------------

def test_run_meta_fields(tmp_path):
    """run_meta.json must contain all required provenance fields."""
    pred_csv = make_predictions_csv(tmp_path / "predictions.csv")

    import run_bald as rb
    orig_module_dir = rb.MODULE_DIR
    out_base = tmp_path / "07_outputs"
    rb.MODULE_DIR = out_base

    try:
        run_bald_main(
            cycle=2,
            predictions_csv=str(pred_csv),
            measured_csv=None,
            n_bald=3,
            n_random=1,
            seed=7,
        )
    finally:
        rb.MODULE_DIR = orig_module_dir

    meta_path = out_base / "outputs" / "cycle_2" / "run_meta.json"
    with open(meta_path) as f:
        meta = json.load(f)

    required_fields = {
        "cycle", "n_bald", "n_random", "seed",
        "n_candidates", "n_unmeasured", "bald_score_col",
        "top_bald_picks", "repo_commit_sha", "timestamp_utc",
    }
    missing = required_fields - set(meta.keys())
    assert not missing, f"run_meta.json missing fields: {missing}"

    assert meta["cycle"] == 2
    assert meta["n_bald"] == 3
    assert meta["n_random"] == 1
    assert len(meta["top_bald_picks"]) == 3
