"""
test_schema.py — Validate output file schemas against INTERFACE.md §Module 06.
验证输出文件 schema 是否符合 INTERFACE.md §Module 06 规范。
"""

import json
import os
from pathlib import Path

import pandas as pd
import pytest

# Locate outputs relative to this test file
# 相对于测试文件定位 outputs 目录
TESTS_DIR = Path(__file__).resolve().parent
MODULE_DIR = TESTS_DIR.parent.parent          # 06_uncertainty_model/
OUTPUTS_DIR = MODULE_DIR / "outputs" / "cycle_0"

REQUIRED_PREDICTION_COLS = {
    "rbp_id", "receptor_id", "predicted_score", "std",
    "lower_95", "upper_95", "model_version",
}

REQUIRED_META_KEYS = {
    "arch", "hidden_dim", "train_size", "val_size",
    "repo_commit_sha", "timestamp_utc",
}


@pytest.mark.skipif(
    not (OUTPUTS_DIR / "predictions.csv").exists(),
    reason="predictions.csv not yet generated (run notebook first)"
)
def test_predictions_csv_schema():
    """
    predictions.csv must have exactly the columns in INTERFACE.md §Module 06.
    predictions.csv 必须包含 INTERFACE.md §Module 06 规定的全部列。
    """
    df = pd.read_csv(OUTPUTS_DIR / "predictions.csv")
    missing = REQUIRED_PREDICTION_COLS - set(df.columns)
    assert not missing, f"predictions.csv missing columns: {missing}"

    # No NaN in core numeric columns
    # 核心数值列不能有 NaN
    for col in ("predicted_score", "std", "lower_95", "upper_95"):
        assert df[col].notna().all(), f"NaN found in column '{col}'"

    # std must be non-negative / std 必须非负
    assert (df["std"] >= 0).all(), "std column contains negative values"

    # lower_95 <= predicted_score <= upper_95
    # 置信区间方向正确
    assert (df["lower_95"] <= df["predicted_score"]).all(), \
        "lower_95 > predicted_score for some rows"
    assert (df["predicted_score"] <= df["upper_95"]).all(), \
        "predicted_score > upper_95 for some rows"


@pytest.mark.skipif(
    not (OUTPUTS_DIR / "model_meta.json").exists(),
    reason="model_meta.json not yet generated"
)
def test_model_meta_json_schema():
    """
    model_meta.json must contain all required keys.
    model_meta.json 必须包含所有必要字段。
    """
    with open(OUTPUTS_DIR / "model_meta.json") as f:
        meta = json.load(f)
    missing = REQUIRED_META_KEYS - set(meta.keys())
    assert not missing, f"model_meta.json missing keys: {missing}"


@pytest.mark.skipif(
    not (OUTPUTS_DIR.parent / "MANIFEST.csv").exists(),
    reason="MANIFEST.csv not yet generated"
)
def test_manifest_schema():
    """
    MANIFEST.csv must have columns required by universal INTERFACE.md convention.
    MANIFEST.csv 必须包含 INTERFACE.md 通用规范要求的列。
    """
    REQUIRED = {
        "filename", "sha256", "bytes", "n_records",
        "created_utc", "source_acc", "source_module", "notes",
    }
    df = pd.read_csv(OUTPUTS_DIR.parent / "MANIFEST.csv")
    missing = REQUIRED - set(df.columns)
    assert not missing, f"MANIFEST.csv missing columns: {missing}"
