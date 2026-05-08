"""
test_sanity.py — Module-specific biological sanity checks.
模块特定的生物学健全性检查。

Checks against known literature values for reference genomes.
对参考基因组与已知文献数值进行核查。
"""

import csv
import re
import sys
from pathlib import Path

import pytest
from Bio import SeqIO

TESTS_DIR   = Path(__file__).resolve().parent
PROCESS_DIR = TESTS_DIR.parent
MODULE_ROOT = PROCESS_DIR.parent
OUTPUTS_DIR = MODULE_ROOT / "outputs"


class TestPhiL7Sanity:
    """
    Sanity checks for phiL7 (EU717894.1) annotation.
    Lee C.N. et al. (2009) Appl Environ Microbiol 75:7828 reports 59 ORFs.
    PHANOTATE 1.6.7 finds 80 ORFs (includes more small/overlapping ORFs than the 2009 annotation).

    phiL7 (EU717894.1) 注释健全性检查。
    Lee 2009 报告 59 个 ORF；PHANOTATE 1.6.7 找到 80 个。
    """

    PHAGE_FAA  = OUTPUTS_DIR / "phage_proteins" / "EU717894.1.faa"
    PHAGE_GFF3 = OUTPUTS_DIR / "phage_orfs"     / "EU717894.1.gff3"

    def test_phage_faa_exists(self) -> None:
        assert self.PHAGE_FAA.exists(), f"EU717894.1.faa not found: {self.PHAGE_FAA}"

    def test_phage_orf_count_in_range(self) -> None:
        """phiL7 has 50-90 ORFs per PHANOTATE (Lee 2009 reports 59). / phiL7 有 50-90 个 ORF。"""
        records = list(SeqIO.parse(str(self.PHAGE_FAA), "fasta"))
        n = len(records)
        assert 50 <= n <= 90, (
            f"phiL7 ORF count {n} outside expected [50, 90]. "
            "Lee 2009 reports 59; PHANOTATE 1.6.7 finds ~80."
        )

    def test_phage_mean_protein_length(self) -> None:
        """phiL7 mean protein length is in plausible range [50, 800] aa. / 平均蛋白质长度在合理范围内。"""
        records = list(SeqIO.parse(str(self.PHAGE_FAA), "fasta"))
        lengths = [len(r.seq) for r in records]
        mean_len = sum(lengths) / len(lengths)
        assert 50 <= mean_len <= 800, f"Mean protein length {mean_len:.1f} aa outside [50, 800]"

    def test_phage_gff3_orf_count_matches_faa(self) -> None:
        """GFF3 CDS line count matches FAA record count. / GFF3 CDS 行数与 FAA 记录数一致。"""
        cds_count = sum(
            1 for l in self.PHAGE_GFF3.read_text().splitlines()
            if l and not l.startswith("#")
        )
        faa_count = sum(
            1 for l in self.PHAGE_FAA.read_text().splitlines()
            if l.startswith(">")
        )
        assert cds_count == faa_count, (
            f"GFF3 has {cds_count} CDS lines but FAA has {faa_count} records"
        )

    def test_phage_orf_ids_sequential(self) -> None:
        """ORF IDs are consecutive 1..N. / ORF ID 为连续 1..N。"""
        headers = [l for l in self.PHAGE_FAA.read_text().splitlines() if l.startswith(">")]
        ids = [int(h.split("_orf_")[1].split()[0]) for h in headers]
        assert ids == list(range(1, len(ids) + 1)), "ORF IDs are not consecutive from 1"

    def test_phage_no_empty_protein_sequences(self) -> None:
        """All phage protein sequences are non-empty. / 所有噬菌体蛋白质序列非空。"""
        records = list(SeqIO.parse(str(self.PHAGE_FAA), "fasta"))
        empty = [r.id for r in records if len(r.seq) == 0]
        assert not empty, f"Empty protein sequences: {empty}"


class TestXccSanity:
    """
    Sanity checks for Xanthomonas campestris host annotation.
    da Silva et al. (2002) Nature 417:459 reports 4181 ORFs for ATCC 33913.
    NZ_CP155948 (str. 8004) expected ~4000-4500 ORFs.

    Xanthomonas campestris 宿主注释健全性检查。
    da Silva 2002：ATCC 33913 有 4181 个 ORF。NZ_CP155948 预计 ~4000-4500 个。
    """

    HOST_DIR = OUTPUTS_DIR / "host_proteins"

    def _get_xcc_faa(self) -> Path:
        """Find the Xcc FAA file in host_proteins/. / 在 host_proteins/ 中找到 Xcc FAA 文件。"""
        candidates = list(self.HOST_DIR.glob("NZ_CP155948*.faa"))
        if not candidates:
            candidates = list(self.HOST_DIR.glob("*.faa"))
        assert candidates, f"No host FAA files in {self.HOST_DIR}"
        return candidates[0]

    def test_host_faa_exists(self) -> None:
        self._get_xcc_faa()  # asserts inside

    def test_host_orf_count_in_range(self) -> None:
        """Xcc has ≥3000 ORFs (da Silva 2002: ATCC 33913 has 4181). / Xcc 有 ≥3000 个 ORF。"""
        faa = self._get_xcc_faa()
        records = list(SeqIO.parse(str(faa), "fasta"))
        n = len(records)
        assert n >= 3000, f"Host ORF count {n} < 3000 minimum"
        assert n <= 6000, f"Host ORF count {n} implausibly high (>6000)"

    def test_host_mean_protein_length(self) -> None:
        """Xcc mean protein length is in plausible range [100, 600] aa. / 平均蛋白质长度在合理范围内。"""
        faa = self._get_xcc_faa()
        records = list(SeqIO.parse(str(faa), "fasta"))
        lengths = [len(r.seq) for r in records]
        mean_len = sum(lengths) / len(lengths)
        assert 100 <= mean_len <= 600, f"Mean protein length {mean_len:.1f} aa outside [100, 600]"


class TestAnnotationSummary:
    """annotation_summary.csv data integrity. / annotation_summary.csv 数据完整性检查。"""

    SUMMARY = OUTPUTS_DIR / "annotation_summary.csv"

    def test_phage_row_has_correct_tool(self) -> None:
        """Phage row uses PHANOTATE tool. / 噬菌体行使用 PHANOTATE 工具。"""
        with open(self.SUMMARY, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        phage_rows = [r for r in rows if r["type"] == "phage"]
        assert phage_rows, "No phage rows in annotation_summary.csv"
        for r in phage_rows:
            assert "PHANOTATE" in r["tool"], f"Phage row uses wrong tool: {r['tool']}"

    def test_host_row_has_correct_tool(self) -> None:
        """Host row uses pyrodigal tool. / 宿主行使用 pyrodigal 工具。"""
        with open(self.SUMMARY, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        host_rows = [r for r in rows if r["type"] == "host"]
        assert host_rows, "No host rows in annotation_summary.csv"
        for r in host_rows:
            assert "prodigal" in r["tool"].lower() or "pyrodigal" in r["tool"].lower(), (
                f"Host row uses wrong tool: {r['tool']}"
            )
