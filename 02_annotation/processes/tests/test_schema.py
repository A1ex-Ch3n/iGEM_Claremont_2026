"""
test_schema.py — Validate output file schemas match INTERFACE.md spec.
验证输出文件模式符合 INTERFACE.md 规范。
"""

import csv
import re
import sys
from pathlib import Path

import pytest

# Locate module and repo roots / 定位模块和仓库根目录
TESTS_DIR   = Path(__file__).resolve().parent
PROCESS_DIR = TESTS_DIR.parent
MODULE_ROOT = PROCESS_DIR.parent
REPO_ROOT   = MODULE_ROOT.parent
OUTPUTS_DIR = MODULE_ROOT / "outputs"

# INTERFACE-compliant FASTA header regex / INTERFACE 规范的 FASTA 头正则
HEADER_RE = re.compile(
    r'^>(\S+) \| source=(\S+) \| length=(\d+)'
    r' \| start=(\d+) \| end=(\d+) \| strand=([+-])'
    r' \| tool=(\S+)$'
)

SUMMARY_COLS = {"acc", "type", "n_orfs", "mean_orf_len", "tool", "tool_version", "runtime_s"}
MANIFEST_COLS = {"filename", "sha256", "bytes", "n_records", "created_utc",
                 "source_acc", "source_module", "notes"}


class TestFastaHeaders:
    """All FASTA headers must match the INTERFACE spec. / 所有 FASTA 头必须符合接口规范。"""

    def _check_faa(self, faa_path: Path) -> None:
        assert faa_path.exists(), f"FAA not found: {faa_path}"
        headers = [l for l in faa_path.read_text().splitlines() if l.startswith(">")]
        assert headers, f"No FASTA records in {faa_path}"
        bad = [h for h in headers if not HEADER_RE.match(h)]
        assert not bad, (
            f"{faa_path.name}: {len(bad)} malformed headers\n  first: {bad[0][:120]}"
        )

    def test_phage_faa_headers(self) -> None:
        """phiL7 FAA headers match INTERFACE regex. / phiL7 FAA 头符合接口正则。"""
        faa = OUTPUTS_DIR / "phage_proteins" / "EU717894.1.faa"
        self._check_faa(faa)

    def test_host_faa_headers(self) -> None:
        """Xcc FAA headers match INTERFACE regex. / Xcc FAA 头符合接口正则。"""
        # Accept any .faa in host_proteins/ (acc suffix may vary)
        # 接受 host_proteins/ 中任意 .faa 文件（acc 后缀可能不同）
        faa_files = list((OUTPUTS_DIR / "host_proteins").glob("*.faa"))
        assert faa_files, "No host FAA files found in outputs/host_proteins/"
        for faa in faa_files:
            self._check_faa(faa)

    def test_orf_id_format(self) -> None:
        """ORF IDs follow <acc>_orf_<5digit> convention. / ORF id 遵循 <acc>_orf_<5位数字> 格式。"""
        faa = OUTPUTS_DIR / "phage_proteins" / "EU717894.1.faa"
        ORF_ID_RE = re.compile(r'^.+_orf_\d{5}$')
        headers = [l for l in faa.read_text().splitlines() if l.startswith(">")]
        for h in headers:
            orf_id = h.split()[0][1:]  # strip leading >
            assert ORF_ID_RE.match(orf_id), f"Bad ORF ID: {orf_id}"


class TestCSVSchemas:
    """CSV output files have required columns. / CSV 输出文件有必需的列。"""

    def test_annotation_summary_columns(self) -> None:
        """annotation_summary.csv has all required columns. / annotation_summary.csv 有所有必需列。"""
        summary = OUTPUTS_DIR / "annotation_summary.csv"
        assert summary.exists(), f"Missing: {summary}"
        with open(summary, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            cols = set(reader.fieldnames or [])
        missing = SUMMARY_COLS - cols
        assert not missing, f"annotation_summary.csv missing columns: {missing}"

    def test_manifest_columns(self) -> None:
        """MANIFEST.csv has all required columns. / MANIFEST.csv 有所有必需列。"""
        manifest = OUTPUTS_DIR / "MANIFEST.csv"
        assert manifest.exists(), f"Missing: {manifest}"
        with open(manifest, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            cols = set(reader.fieldnames or [])
        missing = MANIFEST_COLS - cols
        assert not missing, f"MANIFEST.csv missing columns: {missing}"

    def test_annotation_summary_has_two_types(self) -> None:
        """annotation_summary.csv has both phage and host rows. / 摘要 CSV 有噬菌体和宿主两种行。"""
        summary = OUTPUTS_DIR / "annotation_summary.csv"
        with open(summary, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        types = {r["type"] for r in rows}
        assert "phage" in types, "No phage rows in annotation_summary.csv"
        assert "host" in types, "No host rows in annotation_summary.csv"

    def test_manifest_source_module(self) -> None:
        """All MANIFEST rows have source_module = '02_annotation'. / 所有 MANIFEST 行来源模块为 02_annotation。"""
        manifest = OUTPUTS_DIR / "MANIFEST.csv"
        with open(manifest, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        bad = [r for r in rows if r.get("source_module") != "02_annotation"]
        assert not bad, f"Wrong source_module in {len(bad)} rows"


class TestGFF3Format:
    """GFF3 files have correct structure. / GFF3 文件有正确的结构。"""

    def _check_gff3(self, gff3_path: Path) -> None:
        assert gff3_path.exists(), f"GFF3 not found: {gff3_path}"
        lines = gff3_path.read_text().splitlines()
        assert lines[0] == "##gff-version 3", f"Missing ##gff-version 3 header in {gff3_path.name}"
        cds_lines = [l for l in lines if l and not l.startswith("#")]
        assert cds_lines, f"No CDS lines in {gff3_path.name}"
        for line in cds_lines[:10]:  # spot-check first 10
            cols = line.split("\t")
            assert len(cols) == 9, f"Wrong column count in GFF3: {line[:80]}"
            assert cols[2] == "CDS", f"Feature type is not CDS: {cols[2]}"
            assert cols[6] in ("+", "-"), f"Invalid strand: {cols[6]}"
            start, end = int(cols[3]), int(cols[4])
            assert start <= end, f"GFF3 start > end: {line[:80]}"

    def test_phage_gff3(self) -> None:
        self._check_gff3(OUTPUTS_DIR / "phage_orfs" / "EU717894.1.gff3")

    def test_host_gff3(self) -> None:
        gff3_files = list((OUTPUTS_DIR / "host_orfs").glob("*.gff3"))
        assert gff3_files, "No host GFF3 files found"
        for g in gff3_files:
            self._check_gff3(g)
