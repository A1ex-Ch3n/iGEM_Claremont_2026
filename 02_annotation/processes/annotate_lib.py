"""
annotate_lib.py — Core annotation helpers for Module 02.
核心注释辅助函数 — 模块 02 使用。

PHANOTATE (phage) and pyrodigal (bacteria) wrappers that produce
INTERFACE-compliant FASTA headers and GFF3 files.
PHANOTATE（噬菌体）和 pyrodigal（细菌）包装器，输出符合 INTERFACE 规范的 FASTA 头和 GFF3 文件。

References / 参考文献:
  McNair et al. (2019) Bioinformatics 35(22):4537 — PHANOTATE
  Hyatt et al. (2010) BMC Bioinformatics 11:119 — Prodigal
  Larralde (2022) JOSS 7(72):4296 — pyrodigal
"""

import hashlib
import logging
import re
import subprocess
import sys
import textwrap
import time
from datetime import datetime, timezone
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

# ── version detection / 版本检测 ──────────────────────────────────────────────
def _phanotate_version() -> str:
    """Return PHANOTATE version string. / 返回 PHANOTATE 版本字符串。"""
    try:
        result = subprocess.run(
            ["phanotate.py", "--version"],
            capture_output=True, text=True
        )
        v = result.stdout.strip() or result.stderr.strip()
        v = re.sub(r"[^\d.]", "", v.split()[-1]) if v else "unknown"
        return v
    except Exception:
        return "unknown"


def _pyrodigal_version() -> str:
    """Return pyrodigal version. / 返回 pyrodigal 版本。"""
    try:
        import pyrodigal
        return pyrodigal.__version__
    except Exception:
        return "unknown"


PHANOTATE_VERSION = _phanotate_version()
PYRODIGAL_VERSION = _pyrodigal_version()

log = logging.getLogger(__name__)


# ── FASTA utilities / FASTA 工具函数 ─────────────────────────────────────────

def _wrap60(seq: str) -> str:
    """Wrap sequence at 60 characters per line. / 每行 60 字符折行。"""
    return "\n".join(textwrap.wrap(seq, 60))


def _orf_id(acc: str, idx: int) -> str:
    """Format ORF id per INTERFACE convention. / 按接口约定格式化 ORF id。"""
    return f"{acc}_orf_{idx:05d}"


def _fasta_header(
    orf_id: str,
    acc: str,
    length_aa: int,
    start: int,
    end: int,
    strand: str,
    tool: str,
) -> str:
    """
    Build an INTERFACE-compliant FASTA header.
    构建符合接口规范的 FASTA 头行。

    Format: ><orf_id> | source=<acc> | length=<aa> | start=<n> | end=<n> | strand=<+/-> | tool=<name>
    """
    return (
        f">{orf_id} | source={acc} | length={length_aa}"
        f" | start={start} | end={end} | strand={strand}"
        f" | tool={tool}"
    )


# ── GFF3 utilities / GFF3 工具函数 ───────────────────────────────────────────

_GFF3_HEADER = "##gff-version 3\n"


def _gff3_line(
    seqid: str,
    source: str,
    start: int,
    end: int,
    strand: str,
    orf_id: str,
    score: str = ".",
) -> str:
    """
    Return one GFF3 CDS line.
    返回一行 GFF3 CDS 记录。
    GFF3 uses 1-based, closed intervals. / GFF3 使用从1开始的闭区间。
    """
    # GFF3: start <= end always; strand specified in column 7
    # GFF3: start 必须 <= end；strand 在第7列指定
    gff_start = min(start, end)
    gff_end   = max(start, end)
    attr = f"ID={orf_id};Name={orf_id}"
    return (
        f"{seqid}\t{source}\tCDS\t{gff_start}\t{gff_end}\t"
        f"{score}\t{strand}\t.\t{attr}"
    )


# ── PHANOTATE / 噬菌体注释 ───────────────────────────────────────────────────

def run_phanotate(genome_fna: Path, output_dir: Path) -> dict:
    """
    Annotate a phage genome with PHANOTATE. Writes:
    用 PHANOTATE 注释噬菌体基因组，输出：
      - <output_dir>/phage_proteins/<acc>.faa  (protein FASTA / 蛋白质 FASTA)
      - <output_dir>/phage_orfs/<acc>.gff3     (ORF coordinates / ORF 坐标)

    Returns metadata dict / 返回元数据字典:
      {acc, n_orfs, mean_orf_len_aa, tool, tool_version, runtime_s,
       faa_path, gff3_path}
    """
    genome_fna = Path(genome_fna)
    output_dir = Path(output_dir)

    # Infer accession from genome file header / 从 FASTA 头推断登录号
    genome_rec = next(SeqIO.parse(genome_fna, "fasta"))
    acc        = genome_rec.id            # e.g. "EU717894.1"
    genome_seq = str(genome_rec.seq).upper()

    faa_path  = output_dir / "phage_proteins" / f"{acc}.faa"
    gff3_path = output_dir / "phage_orfs"     / f"{acc}.gff3"
    faa_path.parent.mkdir(parents=True, exist_ok=True)
    gff3_path.parent.mkdir(parents=True, exist_ok=True)

    # ── Run PHANOTATE CLI in tabular mode / 以表格模式运行 PHANOTATE ──
    tmp_tsv = output_dir / f"_tmp_{acc}_phanotate.tsv"
    t0 = time.perf_counter()
    result = subprocess.run(
        ["phanotate.py", "-f", "tabular", "-o", str(tmp_tsv), str(genome_fna)],
        capture_output=True, text=True
    )
    runtime_s = time.perf_counter() - t0

    # If -o flag doesn't write (old phanotate behaviour), capture stdout
    # 若 -o 标志未写文件（旧版行为），则从 stdout 读取
    if tmp_tsv.stat().st_size == 0:
        tmp_tsv.write_text(result.stdout)

    if result.returncode != 0:
        raise RuntimeError(f"PHANOTATE failed on {acc}:\n{result.stderr}")

    # ── Parse tabular output / 解析表格输出 ──────────────────────────
    orfs = []          # list of (start, stop, strand, score) — 1-based
    with open(tmp_tsv) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            start_col, stop_col, frame_col = parts[0], parts[1], parts[2]
            score_col = parts[4] if len(parts) >= 5 else "."
            try:
                s = int(start_col)
                e = int(stop_col)
                strand = "+" if frame_col.strip() == "+" else "-"
                orfs.append((s, e, strand, score_col))
            except ValueError:
                continue
    tmp_tsv.unlink(missing_ok=True)

    # ── Translate ORFs and write outputs / 翻译 ORF 并写入输出 ───────
    faa_lines  = []
    gff3_lines = [_GFF3_HEADER, f"##sequence-region {acc} 1 {len(genome_seq)}\n"]
    tool_label = f"PHANOTATE_{PHANOTATE_VERSION}"

    lengths_aa = []
    for idx, (start, stop, strand, score) in enumerate(orfs, start=1):
        orf_id = _orf_id(acc, idx)

        # Extract and translate nucleotide sequence
        # 提取并翻译核苷酸序列
        # PHANOTATE coords are 1-based; + strand: start < stop; - strand: start > stop
        # 坐标为 1-based；正链 start < stop；负链 start > stop
        if strand == "+":
            nuc = genome_seq[start - 1 : stop]
        else:
            nuc = genome_seq[stop - 1 : start]
            nuc = str(Seq(nuc).reverse_complement())

        # Pad to multiple of 3 if needed (should not happen with valid ORFs)
        # 若长度不是3的倍数则填充（正常情况不应发生）
        remainder = len(nuc) % 3
        if remainder:
            nuc = nuc + "N" * (3 - remainder)

        protein = str(Seq(nuc).translate())
        # Remove trailing stop codon symbol but keep internal ones for debugging
        # 去除末尾终止密码子符号
        if protein.endswith("*"):
            protein = protein[:-1]

        length_aa = len(protein)
        lengths_aa.append(length_aa)

        # FASTA header / FASTA 头
        header = _fasta_header(orf_id, acc, length_aa, start, stop, strand, tool_label)
        faa_lines.append(header)
        faa_lines.append(_wrap60(protein))

        # GFF3 line / GFF3 行
        gff3_lines.append(
            _gff3_line(acc, "PHANOTATE", start, stop, strand, orf_id, score)
        )

    faa_path.write_text("\n".join(faa_lines) + "\n")
    gff3_path.write_text("\n".join(gff3_lines) + "\n")

    n_orfs       = len(orfs)
    mean_orf_len = round(sum(lengths_aa) / n_orfs, 2) if n_orfs else 0.0

    log.info("PHANOTATE: %s → %d ORFs in %.1fs", acc, n_orfs, runtime_s)

    return {
        "acc":           acc,
        "n_orfs":        n_orfs,
        "mean_orf_len":  mean_orf_len,
        "tool":          "PHANOTATE",
        "tool_version":  PHANOTATE_VERSION,
        "runtime_s":     round(runtime_s, 2),
        "faa_path":      str(faa_path),
        "gff3_path":     str(gff3_path),
    }


# ── Prodigal / pyrodigal (bacteria) / 细菌注释 ──────────────────────────────

def run_prodigal(genome_fna: Path, output_dir: Path) -> dict:
    """
    Annotate a bacterial genome with pyrodigal. Writes:
    用 pyrodigal 注释细菌基因组，输出：
      - <output_dir>/host_proteins/<acc>.faa
      - <output_dir>/host_orfs/<acc>.gff3

    Returns metadata dict / 返回元数据字典.
    """
    import pyrodigal  # local import so callers without pyrodigal still load lib

    genome_fna = Path(genome_fna)
    output_dir = Path(output_dir)

    # Load all records (bacteria may have chromosomes + plasmids)
    # 加载所有记录（细菌可能有染色体和质粒）
    records = list(SeqIO.parse(genome_fna, "fasta"))
    if not records:
        raise ValueError(f"No sequences found in {genome_fna}")

    # Use accession of first record as file-level identifier
    # 用第一个记录的 ID 作为文件级标识符
    # Strip version suffix for cleaner directory names if needed
    acc = records[0].id  # e.g. "NZ_CP155948.1"

    faa_path  = output_dir / "host_proteins" / f"{acc}.faa"
    gff3_path = output_dir / "host_orfs"     / f"{acc}.gff3"
    faa_path.parent.mkdir(parents=True, exist_ok=True)
    gff3_path.parent.mkdir(parents=True, exist_ok=True)

    tool_label = f"pyrodigal_{PYRODIGAL_VERSION}"

    faa_lines  = []
    gff3_lines = [_GFF3_HEADER]
    lengths_aa = []
    global_orf_idx = 0   # global ORF counter across all contigs
    t0 = time.perf_counter()

    for rec in records:
        seq_str = str(rec.seq).upper()
        gff3_lines.append(f"##sequence-region {rec.id} 1 {len(seq_str)}")

        # Train Prodigal on this sequence (single mode)
        # 对当前序列运行 Prodigal 单基因组训练模式
        # single mode: train on this genome; meta mode: use pre-trained models
        # single 模式：对该基因组从头训练；meta 模式：用预训练模型
        orf_finder = pyrodigal.GeneFinder(meta=False)
        orf_finder.train(seq_str)
        genes = orf_finder.find_genes(seq_str)

        for gene in genes:
            global_orf_idx += 1
            orf_id = _orf_id(acc, global_orf_idx)

            protein   = gene.translate()
            if isinstance(protein, bytes):
                protein = protein.decode()
            protein   = protein.rstrip("*")
            length_aa = len(protein)
            lengths_aa.append(length_aa)

            strand = "+" if gene.strand == 1 else "-"
            # pyrodigal uses 1-based coordinates for start/end
            # pyrodigal 使用 1-based 坐标
            start = gene.begin   # 1-based
            end   = gene.end     # 1-based

            header = _fasta_header(
                orf_id, acc, length_aa, start, end, strand, tool_label
            )
            faa_lines.append(header)
            faa_lines.append(_wrap60(protein))

            gff3_lines.append(
                _gff3_line(rec.id, "pyrodigal", start, end, strand, orf_id)
            )

    runtime_s = time.perf_counter() - t0

    faa_path.write_text("\n".join(faa_lines) + "\n")
    gff3_path.write_text("\n".join(gff3_lines) + "\n")

    n_orfs       = len(lengths_aa)
    mean_orf_len = round(sum(lengths_aa) / n_orfs, 2) if n_orfs else 0.0

    log.info("pyrodigal: %s → %d ORFs in %.1fs", acc, n_orfs, runtime_s)

    return {
        "acc":           acc,
        "n_orfs":        n_orfs,
        "mean_orf_len":  mean_orf_len,
        "tool":          "pyrodigal",
        "tool_version":  PYRODIGAL_VERSION,
        "runtime_s":     round(runtime_s, 2),
        "faa_path":      str(faa_path),
        "gff3_path":     str(gff3_path),
    }


# ── MANIFEST / 清单 ──────────────────────────────────────────────────────────

def sha256_file(path: Path) -> str:
    """Compute SHA-256 hex digest. / 计算 SHA-256 十六进制摘要。"""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def append_manifest(manifest_path: Path, rows: list[dict]) -> None:
    """
    Append rows to MANIFEST.csv (create with header if missing).
    向 MANIFEST.csv 追加记录（若文件不存在则创建并写入表头）。
    Columns per INTERFACE.md / 按接口文档规定的列：
      filename, sha256, bytes, n_records, created_utc, source_acc, source_module, notes
    """
    import csv
    COLS = ["filename", "sha256", "bytes", "n_records",
            "created_utc", "source_acc", "source_module", "notes"]
    write_header = not manifest_path.exists()
    with open(manifest_path, "a", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLS)
        if write_header:
            writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in COLS})


def manifest_row_for_file(path: Path, source_acc: str, n_records: int,
                           notes: str = "") -> dict:
    """Build one MANIFEST row for a file. / 为文件构建一行 MANIFEST 记录。"""
    path = Path(path)
    return {
        "filename":      path.name,
        "sha256":        sha256_file(path),
        "bytes":         path.stat().st_size,
        "n_records":     n_records,
        "created_utc":   datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "source_acc":    source_acc,
        "source_module": "02_annotation",
        "notes":         notes,
    }
