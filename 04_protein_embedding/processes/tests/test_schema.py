"""
test_schema.py — verify INTERFACE §Module 04 .npz keys and CSV schemas.
验证 INTERFACE §Module 04 要求的 .npz 键名和 CSV 列名。
"""
import json
import pytest
import numpy as np
import pandas as pd
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[2]
OUTPUTS = MODULE_DIR / "outputs"


# ── .npz schema ─────────────────────────────────────────────────────────────

REQUIRED_NPZ_KEYS = {"seq_ids", "array", "lengths", "meta"}
REQUIRED_META_KEYS = {"model", "pooling", "date_utc", "repo_commit_sha"}


def _check_npz(path: Path) -> None:
    assert path.exists(), f"NPZ not found: {path}"
    data = np.load(path, allow_pickle=True)

    # Required keys / 必需的键名
    missing = REQUIRED_NPZ_KEYS - set(data.files)
    assert not missing, f"Missing keys in {path.name}: {missing}"

    # meta must parse as JSON / meta 必须能解析为 JSON
    meta = json.loads(str(data["meta"]))
    missing_meta = REQUIRED_META_KEYS - set(meta.keys())
    assert not missing_meta, f"Missing meta fields: {missing_meta}"

    # array dtype and shape / array 的数据类型和形状
    arr = data["array"]
    assert arr.dtype == np.float32, f"array dtype should be float32, got {arr.dtype}"
    assert arr.ndim == 2, f"array should be 2-D (mean-pooled), got {arr.ndim}-D"
    N = arr.shape[0]

    # consistent N across keys / 各键的 N 须一致
    assert data["seq_ids"].shape[0] == N
    assert data["lengths"].shape[0] == N

    # lengths dtype / lengths 数据类型
    assert data["lengths"].dtype in (np.int32, np.int64), \
        f"lengths dtype should be int32/int64, got {data['lengths'].dtype}"

    # no NaN / 不含 NaN
    assert not np.any(np.isnan(arr)), "NaN found in embedding array"


def test_phiL7_rbp_npz_schema():
    """phiL7 RBP embedding .npz has all required keys and valid types."""
    _check_npz(OUTPUTS / "embeddings_esm2_t6_8M_phiL7_rbps.npz")


def test_xcc_receptors_npz_schema():
    """Xcc receptor embedding .npz has all required keys and valid types."""
    _check_npz(OUTPUTS / "embeddings_esm2_t6_8M_xcc_receptors.npz")


# ── embedding_index.csv schema ───────────────────────────────────────────────

REQUIRED_INDEX_COLS = {"seq_id", "source_file", "model", "n_residues", "embed_dim", "npz_file"}


def test_embedding_index_csv_schema():
    """embedding_index.csv has all required columns per INTERFACE."""
    idx_path = OUTPUTS / "embedding_index.csv"
    assert idx_path.exists(), f"embedding_index.csv not found at {idx_path}"
    df = pd.read_csv(idx_path)
    missing = REQUIRED_INDEX_COLS - set(df.columns)
    assert not missing, f"Missing columns in embedding_index.csv: {missing}"
    assert len(df) > 0, "embedding_index.csv is empty"


# ── MANIFEST.csv schema ──────────────────────────────────────────────────────

REQUIRED_MANIFEST_COLS = {
    "filename", "sha256", "bytes", "n_records", "created_utc",
    "source_acc", "source_module", "notes"
}


def test_manifest_csv_schema():
    """MANIFEST.csv has all required columns per INTERFACE universal conventions."""
    manifest_path = OUTPUTS / "MANIFEST.csv"
    assert manifest_path.exists(), f"MANIFEST.csv not found at {manifest_path}"
    df = pd.read_csv(manifest_path)
    missing = REQUIRED_MANIFEST_COLS - set(df.columns)
    assert not missing, f"Missing columns in MANIFEST.csv: {missing}"
    assert len(df) > 0, "MANIFEST.csv is empty"
