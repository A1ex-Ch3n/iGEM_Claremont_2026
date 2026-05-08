"""
test_schema.py — Module 01 output schema validation.
验证 Module 01 输出文件的列结构是否符合 INTERFACE.md 约定。
"""
import pathlib
import csv
import pytest

# Anchor repo root from tests/ (tests live in 01_data_ground_truth/processes/tests/)
# 从 tests/ 目录反推到 repo 根目录
REPO_ROOT = pathlib.Path(__file__).resolve().parents[3]
MOD_OUT = REPO_ROOT / "01_data_ground_truth" / "outputs"

REQUIRED_INTERACTION_COLS = {
    "phage_acc", "host_acc", "host_organism",
    "label", "source", "confidence", "notes",
}

REQUIRED_INDEX_COLS = {
    "label", "accession", "target_dir", "status", "sha256_genome",
}

REQUIRED_MANIFEST_COLS = {
    "filename", "sha256", "bytes", "n_records",
    "created_utc", "source_acc", "source_module", "notes",
}


def _read_csv_header(path: pathlib.Path) -> set:
    """Read first-row header of a CSV and return column name set.
    读取 CSV 第一行作为列名集合。
    """
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.reader(fh)
        return set(next(reader))


def test_interaction_matrix_exists():
    """interaction_matrix.csv must exist in outputs/.
    outputs/ 下必须有 interaction_matrix.csv。
    """
    p = MOD_OUT / "interaction_matrix.csv"
    assert p.exists(), f"Missing: {p}"


def test_interaction_matrix_columns():
    """interaction_matrix.csv must have all required INTERFACE columns.
    interaction_matrix.csv 必须包含 INTERFACE.md 规定的全部列。
    """
    p = MOD_OUT / "interaction_matrix.csv"
    cols = _read_csv_header(p)
    missing = REQUIRED_INTERACTION_COLS - cols
    assert not missing, f"Missing columns: {missing}"


def test_interaction_matrix_label_values():
    """label column must contain only 1, 0, or -1.
    label 列只能包含 1、0 或 -1。
    """
    import pandas as pd
    df = pd.read_csv(MOD_OUT / "interaction_matrix.csv")
    bad = df["label"][~df["label"].isin([1, 0, -1])]
    assert bad.empty, f"Bad label values: {bad.unique()}"


def test_interaction_matrix_confidence_range():
    """confidence must be in [0.0, 1.0].
    confidence 必须在 [0.0, 1.0] 范围内。
    """
    import pandas as pd
    df = pd.read_csv(MOD_OUT / "interaction_matrix.csv")
    assert df["confidence"].between(0.0, 1.0).all(), "confidence out of [0,1]"


def test_reference_genomes_index_exists():
    """reference_genomes_index.csv must exist.
    reference_genomes_index.csv 必须存在。
    """
    p = MOD_OUT / "reference_genomes_index.csv"
    assert p.exists(), f"Missing: {p}"


def test_reference_genomes_index_columns():
    """reference_genomes_index.csv must have all required columns.
    reference_genomes_index.csv 必须包含所有必要列。
    """
    p = MOD_OUT / "reference_genomes_index.csv"
    cols = _read_csv_header(p)
    missing = REQUIRED_INDEX_COLS - cols
    assert not missing, f"Missing columns: {missing}"


def test_manifest_exists():
    """outputs/MANIFEST.csv must exist.
    outputs/MANIFEST.csv 必须存在。
    """
    p = MOD_OUT / "MANIFEST.csv"
    assert p.exists(), f"Missing: {p}"


def test_manifest_columns():
    """MANIFEST.csv must have INTERFACE-required columns.
    MANIFEST.csv 必须包含 INTERFACE.md 规定的列。
    """
    p = MOD_OUT / "MANIFEST.csv"
    cols = _read_csv_header(p)
    missing = REQUIRED_MANIFEST_COLS - cols
    assert not missing, f"Missing columns: {missing}"


def test_download_log_exists():
    """At least one download_log_*.csv must be present.
    至少有一个 download_log_*.csv 文件存在。
    """
    logs = list(MOD_OUT.glob("download_log_*.csv"))
    assert logs, "No download_log_*.csv found in outputs/"
