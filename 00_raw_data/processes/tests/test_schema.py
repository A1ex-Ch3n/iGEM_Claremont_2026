"""
test_schema.py — Validate that MANIFEST.csv exists and has the correct columns
per INTERFACE.md §Universal conventions.

验证 MANIFEST.csv 存在且列名与 INTERFACE.md §Universal conventions 完全一致。
"""
from pathlib import Path
import pandas as pd
import pytest

MODULE_ROOT = Path(__file__).resolve().parents[2]
MANIFEST_PATH = MODULE_ROOT / "MANIFEST.csv"

REQUIRED_COLUMNS = ["filename", "sha256", "bytes", "n_records", "created_utc", "source_acc", "source_module", "notes"]


def test_manifest_exists():
    """MANIFEST.csv must exist at the module root. / MANIFEST.csv 必须在模块根目录。"""
    assert MANIFEST_PATH.exists(), f"MANIFEST.csv not found at {MANIFEST_PATH}"


def test_manifest_columns():
    """MANIFEST.csv must have exactly the required columns. / MANIFEST.csv 必须包含规定的列名。"""
    df = pd.read_csv(MANIFEST_PATH)
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    assert not missing, f"Missing columns: {missing}"


def test_manifest_no_empty_filename():
    """No row in MANIFEST.csv should have an empty filename. / 每行的 filename 不能为空。"""
    df = pd.read_csv(MANIFEST_PATH)
    empty = df[df["filename"].isna() | (df["filename"].str.strip() == "")]
    assert len(empty) == 0, f"Found {len(empty)} rows with empty filename"


def test_manifest_sha256_format():
    """sha256 values must be 64-character hex strings. / sha256 值必须是 64 位十六进制字符串。"""
    df = pd.read_csv(MANIFEST_PATH)
    bad = df[~df["sha256"].str.match(r"^[0-9a-f]{64}$", na=False)]
    assert len(bad) == 0, f"Found {len(bad)} rows with invalid sha256"


def test_manifest_bytes_positive():
    """bytes values must be positive integers. / bytes 值必须为正整数。"""
    df = pd.read_csv(MANIFEST_PATH)
    bad = df[df["bytes"] <= 0]
    assert len(bad) == 0, f"Found {len(bad)} rows with non-positive bytes"


def test_manifest_source_module():
    """source_module must be '00_raw_data' for all rows. / source_module 对所有行必须是 '00_raw_data'。"""
    df = pd.read_csv(MANIFEST_PATH)
    bad = df[df["source_module"] != "00_raw_data"]
    assert len(bad) == 0, f"Found {len(bad)} rows with wrong source_module"


def test_manifest_created_utc_format():
    """created_utc must end with 'Z' (UTC ISO 8601). / created_utc 必须以 'Z' 结尾 (UTC ISO 8601)。"""
    df = pd.read_csv(MANIFEST_PATH)
    bad = df[~df["created_utc"].str.endswith("Z", na=False)]
    assert len(bad) == 0, f"Found {len(bad)} rows with non-UTC timestamp"
