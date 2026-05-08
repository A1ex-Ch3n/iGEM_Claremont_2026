"""
test_sanity.py — Module 01 post-fetch sanity checks.
Module 01 抓取后的生物学合理性验证。
"""
import pathlib
import pytest
from Bio import SeqIO

# Anchor to repo root / 反推到 repo 根目录
REPO_ROOT = pathlib.Path(__file__).resolve().parents[3]
RAW = REPO_ROOT / "00_raw_data"

# Known reference genome sizes from literature
# 文献中已知的参考基因组大小
XCC_MIN_BP = 5_050_000   # Xcc ATCC 33913: da Silva et al. 2002 reported 5,076,188 bp
XCC_MAX_BP = 5_100_000
T7_MIN_BP  = 39_800      # T7 phage: 39,937 bp (NC_001604.1)
T7_MAX_BP  = 40_100
PHIL7_MIN_BP = 44_000    # phiL7: 44,763 bp (EU717894.1)
PHIL7_MAX_BP = 46_000


def _genome_length(path: pathlib.Path) -> int:
    """Sum all record lengths in a multi-record FASTA.
    累加 FASTA 文件中所有记录的序列长度。
    """
    total = 0
    for rec in SeqIO.parse(str(path), "fasta"):
        total += len(rec.seq)
    return total


def test_phiL7_genome_exists():
    """phiL7 genome must exist (pre-existing in raw data).
    phiL7 基因组文件必须存在（已在原始数据中）。
    """
    p = RAW / "phage" / "EU717894.1" / "genome.fna"
    assert p.exists(), f"Missing phiL7 genome: {p}"


def test_phiL7_genome_length():
    """phiL7 genome must be 44–46 kbp (Lee et al. 2009: 44,763 bp).
    phiL7 基因组长度应为 44–46 kbp（Lee et al. 2009 报告 44,763 bp）。
    """
    p = RAW / "phage" / "EU717894.1" / "genome.fna"
    if not p.exists():
        pytest.skip("phiL7 genome not yet present")
    length = _genome_length(p)
    assert PHIL7_MIN_BP <= length <= PHIL7_MAX_BP, (
        f"phiL7 length {length} bp outside expected {PHIL7_MIN_BP}–{PHIL7_MAX_BP} bp"
    )


def test_xcc_genome_exists():
    """Xcc ATCC 33913 genome must exist after notebook 01 run.
    运行 notebook 01 后 Xcc ATCC 33913 基因组文件必须存在。
    """
    p = RAW / "bacteria" / "GCF_000007145.1" / "genome.fna"
    assert p.exists(), f"Missing Xcc genome: {p}"


def test_xcc_genome_length():
    """Xcc genome must be 5.05–5.10 Mbp (da Silva 2002: 5,076,188 bp).
    Xcc 基因组应为 5.05–5.10 Mbp（da Silva 2002 报告 5,076,188 bp）。
    """
    p = RAW / "bacteria" / "GCF_000007145.1" / "genome.fna"
    if not p.exists():
        pytest.skip("Xcc genome not yet downloaded")
    length = _genome_length(p)
    assert XCC_MIN_BP <= length <= XCC_MAX_BP, (
        f"Xcc length {length} bp outside expected {XCC_MIN_BP}–{XCC_MAX_BP} bp"
    )


def test_t7_genome_exists():
    """T7 phage genome must exist after notebook 01 run.
    运行 notebook 01 后 T7 噬菌体基因组文件必须存在。
    """
    p = RAW / "phage" / "NC_001604.1" / "genome.fna"
    assert p.exists(), f"Missing T7 genome: {p}"


def test_t7_genome_length():
    """T7 genome must be 39.8–40.1 kbp (NC_001604.1: 39,937 bp).
    T7 基因组应为 39.8–40.1 kbp（NC_001604.1 为 39,937 bp）。
    """
    p = RAW / "phage" / "NC_001604.1" / "genome.fna"
    if not p.exists():
        pytest.skip("T7 genome not yet downloaded")
    length = _genome_length(p)
    assert T7_MIN_BP <= length <= T7_MAX_BP, (
        f"T7 length {length} bp outside expected {T7_MIN_BP}–{T7_MAX_BP} bp"
    )


def test_xcc_proteins_faa_exists():
    """Xcc proteins.faa must exist and be non-empty.
    Xcc proteins.faa 必须存在且非空。
    """
    p = RAW / "bacteria" / "GCF_000007145.1" / "proteins.faa"
    assert p.exists(), f"Missing Xcc proteins.faa: {p}"
    assert p.stat().st_size > 0, "proteins.faa is empty"


def test_t7_proteins_faa_exists():
    """T7 proteins.faa must exist and be non-empty.
    T7 proteins.faa 必须存在且非空。
    """
    p = RAW / "phage" / "NC_001604.1" / "proteins.faa"
    assert p.exists(), f"Missing T7 proteins.faa: {p}"
    assert p.stat().st_size > 0, "T7 proteins.faa is empty"


def test_interaction_matrix_has_phiL7():
    """interaction_matrix.csv must contain at least one phiL7 row.
    interaction_matrix.csv 必须包含至少一条 phiL7 记录。
    """
    import pandas as pd
    p = REPO_ROOT / "01_data_ground_truth" / "outputs" / "interaction_matrix.csv"
    if not p.exists():
        pytest.skip("interaction_matrix.csv not yet produced")
    df = pd.read_csv(p)
    has_phil7 = df["phage_acc"].str.contains("EU717894", na=False).any()
    assert has_phil7, "No phiL7 (EU717894.1) row in interaction_matrix.csv"
