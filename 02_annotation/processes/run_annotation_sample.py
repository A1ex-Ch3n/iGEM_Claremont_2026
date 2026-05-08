"""
run_annotation_sample.py — Sample run script for tonight's build.
样本运行脚本，用于今晚的构建。

Runs PHANOTATE on phiL7 (EU717894.1) and pyrodigal on the smallest
Xanthomonas campestris available, then writes all outputs and MANIFEST.
对 phiL7 (EU717894.1) 运行 PHANOTATE，对最小的 Xcc 运行 pyrodigal，
然后写入所有输出文件和 MANIFEST。
"""

import csv
import logging
import sys
from pathlib import Path

# ── Path anchoring per INTERFACE.md convention (script lives in processes/)
# ── 按接口规范锚定路径（脚本在 processes/ 中）
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from annotate_lib import (
    run_phanotate,
    run_prodigal,
    append_manifest,
    manifest_row_for_file,
    PHANOTATE_VERSION,
    PYRODIGAL_VERSION,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s"
)
log = logging.getLogger(__name__)

RAW       = REPO_ROOT / "00_raw_data"
OUTDIR    = REPO_ROOT / "02_annotation" / "outputs"
MANIFEST  = OUTDIR / "MANIFEST.csv"
SUMMARY   = OUTDIR / "annotation_summary.csv"

# ── Target accessions / 目标登录号 ─────────────────────────────────────────
PHAGE_ACC     = "EU717894.1"          # phiL7
BACTERIA_ACC  = "NZ_CP155948"         # Xanthomonas campestris pv. campestris str. 8004
# (closest available relative of Xcc ATCC 33913 that phiL7 infects)
# (与 phiL7 感染的 Xcc ATCC 33913 最近的可用菌株)

PHAGE_FNA    = RAW / "phage"   / PHAGE_ACC   / "genome.fna"
BACTERIA_FNA = RAW / "bacteria"/ BACTERIA_ACC/ "genome.fna"


def write_summary(rows: list[dict]) -> None:
    """Write annotation_summary.csv. / 写入注释汇总 CSV。"""
    cols = ["acc", "type", "n_orfs", "mean_orf_len", "tool", "tool_version", "runtime_s"]
    write_header = not SUMMARY.exists()
    with open(SUMMARY, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        if write_header:
            w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in cols})


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    summary_rows  = []
    manifest_rows = []

    # ── 1. Phage: PHANOTATE on phiL7 / 噬菌体：对 phiL7 运行 PHANOTATE ──
    log.info("=== Annotating phage %s with PHANOTATE %s ===", PHAGE_ACC, PHANOTATE_VERSION)
    meta_phage = run_phanotate(PHAGE_FNA, OUTDIR)

    log.info("phiL7 ORFs: %d  mean_len: %.1f aa  time: %.1fs",
             meta_phage["n_orfs"], meta_phage["mean_orf_len"], meta_phage["runtime_s"])

    # Sanity check — Lee 2009 reports 59 ORFs; phanotate 1.6 finds slightly more
    # 健全性检查 — Lee 2009 报告 59 个 ORF；phanotate 1.6 会多找几个
    if not (50 <= meta_phage["n_orfs"] <= 90):
        log.warning("phiL7 ORF count %d outside expected range [50, 90]", meta_phage["n_orfs"])
    else:
        log.info("phiL7 ORF count %d in expected range ✓", meta_phage["n_orfs"])

    summary_rows.append({
        "acc":          meta_phage["acc"],
        "type":         "phage",
        "n_orfs":       meta_phage["n_orfs"],
        "mean_orf_len": meta_phage["mean_orf_len"],
        "tool":         meta_phage["tool"],
        "tool_version": meta_phage["tool_version"],
        "runtime_s":    meta_phage["runtime_s"],
    })
    for fpath, n_rec, notes in [
        (meta_phage["faa_path"],  meta_phage["n_orfs"], "PHANOTATE protein FASTA"),
        (meta_phage["gff3_path"], meta_phage["n_orfs"], "PHANOTATE GFF3 ORF coords"),
    ]:
        manifest_rows.append(
            manifest_row_for_file(fpath, meta_phage["acc"], n_rec, notes)
        )

    # ── 2. Bacterium: pyrodigal on Xcc / 细菌：对 Xcc 运行 pyrodigal ──
    log.info("=== Annotating bacterium %s with pyrodigal %s ===",
             BACTERIA_ACC, PYRODIGAL_VERSION)
    meta_bact = run_prodigal(BACTERIA_FNA, OUTDIR)

    log.info("%s ORFs: %d  mean_len: %.1f aa  time: %.1fs",
             BACTERIA_ACC, meta_bact["n_orfs"], meta_bact["mean_orf_len"], meta_bact["runtime_s"])

    # Sanity: Xcc str 8004 should have ~4000-4500 ORFs (similar to ATCC 33913's 4181)
    # 健全性检查：Xcc str 8004 应有约 4000-4500 个 ORF（与 ATCC 33913 的 4181 个相近）
    if not (3000 <= meta_bact["n_orfs"] <= 5500):
        log.warning("Xcc ORF count %d outside expected range [3000, 5500]", meta_bact["n_orfs"])
    else:
        log.info("Xcc ORF count %d in expected range ✓", meta_bact["n_orfs"])

    summary_rows.append({
        "acc":          meta_bact["acc"],
        "type":         "host",
        "n_orfs":       meta_bact["n_orfs"],
        "mean_orf_len": meta_bact["mean_orf_len"],
        "tool":         meta_bact["tool"],
        "tool_version": meta_bact["tool_version"],
        "runtime_s":    meta_bact["runtime_s"],
    })
    for fpath, n_rec, notes in [
        (meta_bact["faa_path"],  meta_bact["n_orfs"], "pyrodigal protein FASTA"),
        (meta_bact["gff3_path"], meta_bact["n_orfs"], "pyrodigal GFF3 ORF coords"),
    ]:
        manifest_rows.append(
            manifest_row_for_file(fpath, meta_bact["acc"], n_rec, notes)
        )

    # ── 3. Write summary + manifest / 写入汇总和清单 ───────────────────
    write_summary(summary_rows)
    append_manifest(MANIFEST, manifest_rows)

    log.info("annotation_summary.csv written to %s", SUMMARY)
    log.info("MANIFEST.csv written to %s", MANIFEST)
    log.info("=== Sample annotation run complete ===")


if __name__ == "__main__":
    main()
