"""
test_smoke.py — Smoke tests for Module 01 manifest update logic.
Module 01 MANIFEST 更新逻辑冒烟测试 (不需要访问 NCBI)。
"""
import hashlib
import pathlib
import csv
import datetime
import tempfile
import pytest

# Anchor to repo root / 反推到 repo 根目录
REPO_ROOT = pathlib.Path(__file__).resolve().parents[3]


def _sha256(path: pathlib.Path) -> str:
    """Compute SHA-256 hex digest of a file.
    计算文件的 SHA-256 摘要。
    """
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _append_manifest_row(manifest_path: pathlib.Path, row: dict) -> None:
    """Append a single row to a MANIFEST.csv, creating it if absent.
    向 MANIFEST.csv 追加一行；若文件不存在则新建并写入表头。
    """
    required_cols = [
        "filename", "sha256", "bytes", "n_records",
        "created_utc", "source_acc", "source_module", "notes",
    ]
    write_header = not manifest_path.exists()
    with open(manifest_path, "a", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=required_cols)
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def test_manifest_append_creates_file():
    """_append_manifest_row should create a new MANIFEST.csv when it does not exist.
    当 MANIFEST.csv 不存在时，_append_manifest_row 应新建该文件。
    """
    with tempfile.TemporaryDirectory() as tmp:
        manifest = pathlib.Path(tmp) / "MANIFEST.csv"
        assert not manifest.exists()

        _append_manifest_row(manifest, {
            "filename": "genome.fna",
            "sha256": "abc123",
            "bytes": 100,
            "n_records": 1,
            "created_utc": "2026-05-07T00:00:00Z",
            "source_acc": "TEST001",
            "source_module": "01_data_ground_truth",
            "notes": "",
        })

        assert manifest.exists()
        with open(manifest, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == 1
        assert rows[0]["source_acc"] == "TEST001"


def test_manifest_append_idempotent_header():
    """Appending a second row must NOT duplicate the header.
    追加第二行时不应重复写入表头。
    """
    with tempfile.TemporaryDirectory() as tmp:
        manifest = pathlib.Path(tmp) / "MANIFEST.csv"
        row = {
            "filename": "genome.fna", "sha256": "abc", "bytes": 99,
            "n_records": 1, "created_utc": "2026-05-07T00:00:00Z",
            "source_acc": "A1", "source_module": "01_data_ground_truth", "notes": "",
        }
        _append_manifest_row(manifest, row)
        row2 = dict(row, source_acc="A2")
        _append_manifest_row(manifest, row2)

        with open(manifest, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == 2, f"Expected 2 data rows, got {len(rows)}"


def test_sha256_deterministic():
    """SHA-256 of a known file must be deterministic.
    对已知文件计算 SHA-256 必须保持幂等。
    """
    with tempfile.NamedTemporaryFile(suffix=".fna", delete=False) as tf:
        tf.write(b">test\nACGTACGT\n")
        path = pathlib.Path(tf.name)
    try:
        h1 = _sha256(path)
        h2 = _sha256(path)
        assert h1 == h2
        assert len(h1) == 64  # SHA-256 hex is 64 chars
    finally:
        path.unlink(missing_ok=True)


def test_reference_targets_csv_parseable():
    """inputs/reference_targets.csv must be parseable and have correct columns.
    inputs/reference_targets.csv 必须可解析且含有正确列名。
    """
    p = REPO_ROOT / "01_data_ground_truth" / "inputs" / "reference_targets.csv"
    assert p.exists(), f"Missing reference_targets.csv at {p}"

    with open(p, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)

    required_cols = {"category", "assembly_acc", "nucleotide_acc", "label", "subdir", "priority"}
    assert required_cols <= set(rows[0].keys()), f"Missing columns: {required_cols - set(rows[0].keys())}"
    assert len(rows) >= 3, "reference_targets.csv should have at least 3 entries"

    categories = {r["category"] for r in rows}
    assert "phage" in categories
    assert "bacteria" in categories
