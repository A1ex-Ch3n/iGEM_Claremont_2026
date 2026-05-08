"""
test_smoke.py — Run annotation helpers on minimal synthetic inputs.
对最小合成输入运行注释辅助函数。
Each test must complete in < 60 seconds. / 每个测试必须在 60 秒内完成。
"""

import sys
import tempfile
from pathlib import Path

import pytest

TESTS_DIR   = Path(__file__).resolve().parent
PROCESS_DIR = TESTS_DIR.parent
sys.path.insert(0, str(PROCESS_DIR))

from annotate_lib import run_phanotate, run_prodigal


# ── Synthetic 5 kb phage genome (single ORF region) / 合成 5 kb 噬菌体基因组（单 ORF 区域）
# Contains an obvious ORF: ATG ... TAA in frame
# 包含一个明显的 ORF: ATG ... TAA 在读码框中
_PHAGE_SEQ = (
    "ATGCGTACGCGTACGCGTACGCGTACGCGTACGCGT"
    "ACGCGTACGCGTACGCGTACGCGTACGCGTACGCGT"
    "ACGCGTACGCGTACGCGTACGCGTACGCGTACGCGT"
    "ACGCGTACGCGTACGCGTACGCGTACGCGTACGCGT"
    "ACGCGTACGCGTACGCGTACGCGTACGCGTACGTAA"   # TAA stop
) * 28   # ~5 kb

# Bacteria-like 22 kb synthetic genome (pyrodigal requires ≥20000 bp for single mode)
# 细菌-like 22 kb 合成基因组（pyrodigal 单基因组模式需要 ≥20000 bp）
_BACT_SEQ = (
    "ATGAAACCCGGGTTTTAA"  # small ORF
    + "GGCTAGCTAGCTAGCTA" * 5
    + "ATGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCTAA"
    + "TTAAGGCCTTAAGGCCT" * 5
) * 220   # ~22 kb


def _write_fasta(seq: str, acc: str, tmp_dir: Path) -> Path:
    """Write a FASTA file for testing. / 为测试写入 FASTA 文件。"""
    fna = tmp_dir / f"{acc}" / "genome.fna"
    fna.parent.mkdir(parents=True, exist_ok=True)
    fna.write_text(f">{acc} synthetic test genome\n" + seq + "\n")
    return fna


class TestPhanotateSmoke:
    """PHANOTATE smoke test on synthetic phage. / PHANOTATE 在合成噬菌体上的冒烟测试。"""

    def test_phanotate_finds_at_least_one_orf(self, tmp_path: Path) -> None:
        """run_phanotate on a 5 kb synthetic genome returns ≥1 ORF. / 对 5 kb 合成基因组返回 ≥1 个 ORF。"""
        fna = _write_fasta(_PHAGE_SEQ, "SMOKE001.1", tmp_path)
        result = run_phanotate(fna, tmp_path)
        assert result["n_orfs"] >= 1, f"Expected ≥1 ORF, got {result['n_orfs']}"

    def test_phanotate_returns_faa_and_gff3(self, tmp_path: Path) -> None:
        """run_phanotate writes both .faa and .gff3 files. / 写出 .faa 和 .gff3 文件。"""
        fna = _write_fasta(_PHAGE_SEQ, "SMOKE002.1", tmp_path)
        result = run_phanotate(fna, tmp_path)
        assert Path(result["faa_path"]).exists(),  "FAA not written"
        assert Path(result["gff3_path"]).exists(), "GFF3 not written"

    def test_phanotate_faa_has_valid_headers(self, tmp_path: Path) -> None:
        """PHANOTATE FAA headers match INTERFACE pattern. / FAA 头符合接口格式。"""
        import re
        fna = _write_fasta(_PHAGE_SEQ, "SMOKE003.1", tmp_path)
        result = run_phanotate(fna, tmp_path)
        faa_text = Path(result["faa_path"]).read_text()
        headers  = [l for l in faa_text.splitlines() if l.startswith(">")]
        RE = re.compile(
            r'^>(\S+) \| source=(\S+) \| length=(\d+)'
            r' \| start=(\d+) \| end=(\d+) \| strand=([+-])'
            r' \| tool=(\S+)$'
        )
        bad = [h for h in headers if not RE.match(h)]
        assert not bad, f"Malformed headers: {bad}"


class TestProdigalSmoke:
    """pyrodigal smoke test on synthetic bacteria. / pyrodigal 在合成细菌上的冒烟测试。"""

    def test_prodigal_finds_at_least_one_orf(self, tmp_path: Path) -> None:
        """run_prodigal on a 10 kb synthetic genome returns ≥1 ORF. / 对 10 kb 合成基因组返回 ≥1 个 ORF。"""
        fna = _write_fasta(_BACT_SEQ, "BSMOKE001", tmp_path)
        result = run_prodigal(fna, tmp_path)
        assert result["n_orfs"] >= 1, f"Expected ≥1 ORF, got {result['n_orfs']}"

    def test_prodigal_returns_faa_and_gff3(self, tmp_path: Path) -> None:
        """run_prodigal writes both .faa and .gff3 files. / 写出 .faa 和 .gff3 文件。"""
        fna = _write_fasta(_BACT_SEQ, "BSMOKE002", tmp_path)
        result = run_prodigal(fna, tmp_path)
        assert Path(result["faa_path"]).exists(),  "FAA not written"
        assert Path(result["gff3_path"]).exists(), "GFF3 not written"

    def test_prodigal_meta_dict_keys(self, tmp_path: Path) -> None:
        """run_prodigal returns all expected metadata keys. / 返回所有期望的元数据键。"""
        fna = _write_fasta(_BACT_SEQ, "BSMOKE003", tmp_path)
        result = run_prodigal(fna, tmp_path)
        for key in ("acc", "n_orfs", "mean_orf_len", "tool", "tool_version",
                    "runtime_s", "faa_path", "gff3_path"):
            assert key in result, f"Missing key: {key}"
