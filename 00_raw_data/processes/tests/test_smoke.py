"""
test_smoke.py — Parse the 3 smallest genome.fna files with Bio.SeqIO and verify
they each contain at least 1 sequence record.

用 Bio.SeqIO 解析最小的 3 个 genome.fna 文件，验证每个至少包含 1 条序列。
"""
from pathlib import Path
import pytest
from Bio import SeqIO

MODULE_ROOT = Path(__file__).resolve().parents[2]
PHAGE_ROOT = MODULE_ROOT / "phage"

# Three smallest phage genome files (determined by file size at dataset creation time)
# 三个最小的噬菌体基因组文件（按文件大小排序）
SMALLEST_PHAGE_ACCESSIONS = ["PQ067315.1", "PP895210.1", "PP895209.1"]


@pytest.fixture(params=SMALLEST_PHAGE_ACCESSIONS)
def genome_path(request):
    path = PHAGE_ROOT / request.param / "genome.fna"
    if not path.exists():
        pytest.skip(f"genome.fna missing for {request.param}")
    return path


def test_genome_parseable(genome_path):
    """Each genome.fna must parse without errors. / 每个 genome.fna 必须能无错解析。"""
    records = list(SeqIO.parse(str(genome_path), "fasta"))
    assert len(records) >= 1, f"No records found in {genome_path}"


def test_genome_nonzero_length(genome_path):
    """Each genome sequence must have non-zero length. / 每条基因组序列长度必须大于零。"""
    records = list(SeqIO.parse(str(genome_path), "fasta"))
    for rec in records:
        assert len(rec.seq) > 0, f"Zero-length sequence in {genome_path}: {rec.id}"


def test_phage_dir_structure():
    """Every phage directory must contain at least genome.fna. / 每个噬菌体目录必须包含 genome.fna。"""
    missing_genome = []
    for acc_dir in sorted(PHAGE_ROOT.iterdir()):
        if not acc_dir.is_dir():
            continue
        genome_fna = acc_dir / "genome.fna"
        if not genome_fna.exists():
            missing_genome.append(acc_dir.name)
    assert len(missing_genome) == 0, (
        f"{len(missing_genome)} phage dirs missing genome.fna: {missing_genome[:10]}"
    )
