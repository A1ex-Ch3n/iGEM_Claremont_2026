#!/usr/bin/env python3
"""
Fetch genome data for all phages in phage_list.csv using NCBI datasets CLI
(virus genome endpoint, since phages are virus records).

Downloads per phage accession:
  genome.fna       — whole-genome nucleotide sequence
  proteins.faa     — NCBI-annotated protein sequences
  cds.fna          — CDS nucleotide sequences
  (no genes.gff — not available via virus endpoint; use pharokka output in 02_annotation/)

Note on sources:
  canonical_xanthomonas (334) — Xanthomonas-infecting phages in the interaction matrix.
  negative_pool (443)         — Non-Xanthomonas phages used as presumed-negative ML training
                                samples. Must have explicit 0-fill rows added to the interaction
                                matrix before ML training (see data_needs.md).

Setup (run once):
  conda install -n igem2026 -c conda-forge ncbi-datasets-cli

Usage:
  python 00_raw_data/processes/fetch_phages.py                     # all 777
  python 00_raw_data/processes/fetch_phages.py --source canonical_xanthomonas
  python 00_raw_data/processes/fetch_phages.py --source negative_pool
  python 00_raw_data/processes/fetch_phages.py --dry-run
  python 00_raw_data/processes/fetch_phages.py --accession AB720063.2
"""

import argparse
import csv
import json
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).parents[1]
PHAGE_DIR = ROOT / "phage"
LIST_CSV = ROOT / "phage_list.csv"
TMP_DIR = ROOT / "processes" / ".tmp"

BATCH_SIZE = 50
EXPECTED_FILES = ["genome.fna", "proteins.faa", "cds.fna"]


def check_datasets_installed() -> bool:
    return subprocess.run(["datasets", "--version"], capture_output=True).returncode == 0


def _split_fasta_by_accession(content: str) -> dict[str, str]:
    """Split a multi-FASTA string into {accession: fasta_block} dict."""
    blocks: dict[str, list[str]] = {}
    current_acc = None
    for line in content.splitlines():
        if line.startswith(">"):
            # Accession is the first token of the header
            current_acc = line[1:].split()[0].split(":")[0]
            blocks.setdefault(current_acc, []).append(line)
        elif current_acc:
            blocks[current_acc].append(line)
    return {acc: "\n".join(lines) + "\n" for acc, lines in blocks.items()}


def _split_proteins_by_accession(
    prot_content: str, gene_counts: dict[str, int]
) -> dict[str, str]:
    """
    Split protein FASTA by genome accession using sequential gene counts.
    gene_counts: {accession: n_proteins} in order of download.
    Proteins for genome i appear immediately after proteins for genome i-1.
    """
    all_records: list[tuple[str, list[str]]] = []
    current_header = None
    current_lines: list[str] = []
    for line in prot_content.splitlines():
        if line.startswith(">"):
            if current_header is not None:
                all_records.append((current_header, current_lines))
            current_header = line
            current_lines = []
        else:
            if current_header is not None:
                current_lines.append(line)
    if current_header is not None:
        all_records.append((current_header, current_lines))

    result: dict[str, str] = {}
    idx = 0
    for acc, count in gene_counts.items():
        chunk = all_records[idx : idx + count]
        result[acc] = "\n".join(
            h + "\n" + "\n".join(lines) for h, lines in chunk
        ) + "\n"
        idx += count
    return result


def fetch_batch(
    accessions: list[str], dry_run: bool = False
) -> dict[str, str]:
    """
    Fetch a batch of phage accessions in one datasets virus call.
    Returns: {accession: "ok" | "skip" | "fail:<reason>"}
    """
    to_fetch = []
    results = {}
    for acc in accessions:
        out_dir = PHAGE_DIR / acc
        missing = [n for n in EXPECTED_FILES if not (out_dir / n).exists()]
        if not missing:
            results[acc] = "skip"
        else:
            to_fetch.append(acc)

    if not to_fetch:
        return results

    if dry_run:
        for acc in to_fetch:
            print(f"  [dry-run] {acc}")
            results[acc] = "dry-run"
        return results

    TMP_DIR.mkdir(parents=True, exist_ok=True)
    batch_id = to_fetch[0].replace(".", "_")
    tmp_zip = TMP_DIR / f"phage_{batch_id}.zip"
    acc_file = TMP_DIR / f"phage_{batch_id}_accs.txt"
    acc_file.write_text("\n".join(to_fetch))

    try:
        cmd = [
            "datasets", "download", "virus", "genome", "accession",
            "--inputfile", str(acc_file),
            "--include", "genome,cds,protein",
            "--filename", str(tmp_zip),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            err = (result.stderr.strip().splitlines() or ["unknown error"])[-1]
            return {acc: f"fail:{err}" for acc in to_fetch} | results

        with zipfile.ZipFile(tmp_zip, "r") as z:
            names = z.namelist()

            # Load gene counts per accession from data_report
            gene_counts: dict[str, int] = {}
            if "ncbi_dataset/data/data_report.jsonl" in names:
                for line in z.read("ncbi_dataset/data/data_report.jsonl").decode("utf-8", "replace").splitlines():
                    if not line.strip():
                        continue
                    try:
                        d = json.loads(line)
                        gene_counts[d["accession"]] = d.get("geneCount", 0)
                    except Exception:
                        pass

            # Read bulk files
            genomic = z.read("ncbi_dataset/data/genomic.fna").decode("utf-8", "replace") if "ncbi_dataset/data/genomic.fna" in names else ""
            cds_raw = z.read("ncbi_dataset/data/cds.fna").decode("utf-8", "replace") if "ncbi_dataset/data/cds.fna" in names else ""
            prot_raw = z.read("ncbi_dataset/data/protein.faa").decode("utf-8", "replace") if "ncbi_dataset/data/protein.faa" in names else ""

        # Split by accession
        genome_map = _split_fasta_by_accession(genomic)
        cds_map = _split_fasta_by_accession(cds_raw)
        prot_map = _split_proteins_by_accession(prot_raw, gene_counts) if gene_counts else {}

        for acc in to_fetch:
            out_dir = PHAGE_DIR / acc
            out_dir.mkdir(parents=True, exist_ok=True)
            written = []

            if acc in genome_map:
                (out_dir / "genome.fna").write_text(genome_map[acc])
                written.append("genome.fna")
            if acc in cds_map:
                (out_dir / "cds.fna").write_text(cds_map[acc])
                written.append("cds.fna")
            if acc in prot_map:
                (out_dir / "proteins.faa").write_text(prot_map[acc])
                written.append("proteins.faa")

            missing_after = [n for n in EXPECTED_FILES if not (out_dir / n).exists()]
            results[acc] = "ok" if not missing_after else f"fail:missing {','.join(missing_after)}"

        return results

    except subprocess.TimeoutExpired:
        return {acc: "fail:timeout" for acc in to_fetch if acc not in results} | results
    except Exception as e:
        return {acc: f"fail:{e}" for acc in to_fetch if acc not in results} | results
    finally:
        if tmp_zip.exists():
            tmp_zip.unlink()
        if acc_file.exists():
            acc_file.unlink()


def load_phages(
    source_filter: str | None = None, only_accession: str | None = None
) -> list[dict]:
    with open(LIST_CSV, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if source_filter:
        rows = [r for r in rows if r["source"] == source_filter]
    if only_accession:
        rows = [r for r in rows if r["accession"] == only_accession]
    return rows


def main():
    parser = argparse.ArgumentParser(description="Fetch phage genome data via NCBI datasets CLI")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--source", choices=["canonical_xanthomonas", "negative_pool"])
    parser.add_argument("--accession", help="Fetch a single accession only")
    args = parser.parse_args()

    if not args.dry_run and not check_datasets_installed():
        print("ERROR: 'datasets' CLI not found.")
        print("Install:  conda install -n igem2026 -c conda-forge ncbi-datasets-cli")
        sys.exit(1)

    rows = load_phages(args.source, args.accession)
    accessions = [r["accession"] for r in rows]
    print(f"Processing {len(accessions)} phages {'[DRY RUN]' if args.dry_run else ''}...")

    stats: dict[str, int] = {}
    failures: list[tuple[str, str]] = []

    for i in range(0, len(accessions), BATCH_SIZE):
        batch = accessions[i : i + BATCH_SIZE]
        batch_num = i // BATCH_SIZE + 1
        total_batches = (len(accessions) + BATCH_SIZE - 1) // BATCH_SIZE
        print(f"\nBatch {batch_num}/{total_batches}  ({batch[0]} … {batch[-1]})")

        batch_results = fetch_batch(batch, dry_run=args.dry_run)

        for acc, status in batch_results.items():
            tag = status.split(":")[0]
            stats[tag] = stats.get(tag, 0) + 1
            if tag == "ok":
                print(f"  [OK]   {acc}")
            elif tag == "skip":
                pass  # silent skip to reduce noise
            elif tag == "fail":
                reason = status.split(":", 1)[1] if ":" in status else status
                print(f"  [FAIL] {acc}: {reason}")
                failures.append((acc, reason))

    skipped = stats.get("skip", 0)
    print(f"\nSummary: {stats.get('ok',0)} downloaded, {skipped} already complete, "
          f"{stats.get('fail',0)} failed")
    if failures:
        print("\nFailed accessions (re-run with --accession to retry):")
        for acc, reason in failures:
            print(f"  {acc}: {reason}")


if __name__ == "__main__":
    main()
