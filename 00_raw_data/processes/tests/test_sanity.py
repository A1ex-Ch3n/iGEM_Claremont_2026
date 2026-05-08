"""
test_sanity.py — Module-specific sanity checks:
  1. phiL7 reference genome (EU717894.1) has expected length 44,080 bp ± 100.
  2. Phage directory count matches phage_list.csv, modulo documented missing accessions.

模块专属合理性检查：
  1. phiL7 参考基因组 (EU717894.1) 长度应为 44,080 bp ± 100。
  2. 噬菌体目录数量与 phage_list.csv 一致（扣除已知缺失条目）。
"""
from pathlib import Path
import pandas as pd
import pytest
from Bio import SeqIO

MODULE_ROOT = Path(__file__).resolve().parents[2]
PHAGE_ROOT = MODULE_ROOT / "phage"

# Known accessions listed in phage_list.csv but missing from disk (documented in AGENT_REPORT.md)
# 已知在 phage_list.csv 中列出但磁盘上缺失的条目（记录于 AGENT_REPORT.md）
KNOWN_MISSING_PHAGE_DIRS = {"NC_013971.1", "NZ_CP007800.1", "NZ_CP008698.1"}

# Known items in phage/ that are NOT directories (e.g. pre-existing MANIFEST.csv)
KNOWN_NON_DIR_ITEMS = {"MANIFEST.csv"}


def test_phil7_genome_length():
    """phiL7 genome (EU717894.1) must be 44,080 bp ± 100. / phiL7 基因组长度应为 44,080 ± 100 bp。"""
    genome_path = PHAGE_ROOT / "EU717894.1" / "genome.fna"
    assert genome_path.exists(), "phiL7 genome.fna not found"
    records = list(SeqIO.parse(str(genome_path), "fasta"))
    assert len(records) == 1, f"Expected 1 record in phiL7 genome, got {len(records)}"
    length = len(records[0].seq)
    assert abs(length - 44080) <= 100, (
        f"phiL7 genome length {length} bp outside expected range 44,080 ± 100"
    )


def test_phage_list_vs_dirs():
    """
    Phage directory count should match phage_list.csv minus known-missing accessions.
    噬菌体目录数量应与 phage_list.csv 匹配（减去已知缺失条目）。
    """
    phage_list = pd.read_csv(MODULE_ROOT / "phage_list.csv")
    csv_accs = set(phage_list["accession"])

    actual_dirs = {
        d.name for d in PHAGE_ROOT.iterdir()
        if d.is_dir() and d.name not in KNOWN_NON_DIR_ITEMS
    }

    # Accessions in CSV but not on disk (should only be known-missing ones)
    missing_from_disk = csv_accs - actual_dirs
    unexpected_missing = missing_from_disk - KNOWN_MISSING_PHAGE_DIRS
    assert unexpected_missing == set(), (
        f"Unexpected missing phage dirs (not in documented known-missing list): {unexpected_missing}"
    )

    # Dirs on disk but not in CSV
    extra_on_disk = actual_dirs - csv_accs
    assert extra_on_disk == set(), (
        f"Phage dirs on disk not in phage_list.csv: {extra_on_disk}"
    )


def test_bacteria_dirs_in_list():
    """All bacteria directories must appear in bacteria_list.csv. / 所有细菌目录必须在 bacteria_list.csv 中。"""
    BACTERIA_ROOT = MODULE_ROOT / "bacteria"
    bacteria_list = pd.read_csv(MODULE_ROOT / "bacteria_list.csv")
    csv_accs = set(bacteria_list["accession"])

    actual_dirs = {d.name for d in BACTERIA_ROOT.iterdir() if d.is_dir()}
    extra_on_disk = actual_dirs - csv_accs
    assert extra_on_disk == set(), (
        f"Bacteria dirs on disk but not in bacteria_list.csv: {extra_on_disk}"
    )


def test_manifest_row_count():
    """MANIFEST.csv row count should equal total files in phage/ + bacteria/ subtrees. / MANIFEST.csv 行数应等于 phage/ + bacteria/ 子树的总文件数。"""
    BACTERIA_ROOT = MODULE_ROOT / "bacteria"
    manifest_path = MODULE_ROOT / "MANIFEST.csv"
    if not manifest_path.exists():
        pytest.skip("MANIFEST.csv not yet generated")

    # Count actual files
    actual_files = []
    for root_dir in [PHAGE_ROOT, BACTERIA_ROOT]:
        for acc_dir in root_dir.iterdir():
            if acc_dir.is_dir():
                for f in acc_dir.iterdir():
                    if f.is_file():
                        actual_files.append(f)

    df = pd.read_csv(manifest_path)
    assert len(df) == len(actual_files), (
        f"MANIFEST.csv has {len(df)} rows but found {len(actual_files)} files in phage/ + bacteria/"
    )
