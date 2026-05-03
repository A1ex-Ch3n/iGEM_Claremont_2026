#!/usr/bin/env python3
"""
Fetch genome data for all bacteria in bacteria_list.csv using NCBI datasets CLI.

Downloads per accession:
  genome.fna       — whole-genome nucleotide sequence
  proteins.faa     — NCBI-annotated protein sequences
  genes.gff        — gene annotation (coordinates + features)
  cds.fna          — CDS nucleotide sequences (needed for CAI, Factor 6)

Setup (run once):
  conda install -c conda-forge ncbi-datasets-cli
  OR: conda install -n igem2026 -c conda-forge ncbi-datasets-cli

Usage:
  python 00_raw_data/processes/fetch_bacteria.py
  python 00_raw_data/processes/fetch_bacteria.py --dry-run
  python 00_raw_data/processes/fetch_bacteria.py --accession NZ_CP014028
"""

import argparse
import csv
import json
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).parents[1]          # 00_raw_data/
BACTERIA_DIR = ROOT / "bacteria"
LIST_CSV = ROOT / "bacteria_list.csv"
GCF_MAP_JSON = ROOT / "processes" / "accession_to_gcf.json"
TMP_DIR = ROOT / "processes" / ".tmp"

EXPECTED_FILES = {
    "genome.fna": "genomic sequence (.fna)",
    "proteins.faa": "protein sequences (.faa)",
    "genes.gff": "gene annotation (.gff)",
    "cds.fna": "CDS nucleotides (.fna)",
}


def check_datasets_installed() -> bool:
    result = subprocess.run(["datasets", "--version"], capture_output=True)
    return result.returncode == 0


def fetch_one(accession: str, out_dir: Path, dry_run: bool = False) -> str:
    """
    Download all data types for one accession.
    Returns: "ok", "skip", or "fail:<reason>"
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    missing = [name for name in EXPECTED_FILES if not (out_dir / name).exists()]
    if not missing:
        return "skip"

    if dry_run:
        print(f"  [dry-run] would fetch: {accession} → {out_dir.relative_to(ROOT.parent)}")
        return "dry-run"

    TMP_DIR.mkdir(parents=True, exist_ok=True)
    tmp_zip = TMP_DIR / f"{accession}.zip"
    tmp_extract = TMP_DIR / accession

    try:
        cmd = [
            "datasets", "download", "genome", "accession", accession,
            "--include", "genome,protein,gff3,cds",
            "--no-progressbar",
            "--filename", str(tmp_zip),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            err = result.stderr.strip().splitlines()[-1] if result.stderr.strip() else "unknown error"
            return f"fail:{err}"

        with zipfile.ZipFile(tmp_zip, "r") as z:
            z.extractall(tmp_extract)

        data_root = tmp_extract / "ncbi_dataset" / "data"
        if not data_root.exists():
            return "fail:unexpected zip structure — no ncbi_dataset/data/"

        copied = set()
        for subdir in sorted(data_root.iterdir()):
            if not subdir.is_dir():
                continue
            for src in subdir.iterdir():
                dest = _classify_file(src, out_dir)
                if dest and dest.name not in copied:
                    shutil.copy2(src, dest)
                    copied.add(dest.name)

        still_missing = [n for n in EXPECTED_FILES if not (out_dir / n).exists()]
        if still_missing:
            return f"fail:missing after download: {', '.join(still_missing)}"
        return "ok"

    except subprocess.TimeoutExpired:
        return "fail:timeout"
    except Exception as e:
        return f"fail:{e}"
    finally:
        shutil.rmtree(tmp_extract, ignore_errors=True)
        if tmp_zip.exists():
            tmp_zip.unlink()


def _classify_file(src: Path, out_dir: Path) -> Path | None:
    """Map a downloaded file to its canonical destination name."""
    name_lower = src.name.lower()
    if src.suffix not in (".fna", ".faa", ".gff", ".gff3"):
        return None
    if "cds_from_genomic" in name_lower:
        return out_dir / "cds.fna"
    if src.suffix in (".gff", ".gff3"):
        return out_dir / "genes.gff"
    if src.suffix == ".faa":
        return out_dir / "proteins.faa"
    if src.suffix == ".fna":
        return out_dir / "genome.fna"
    return None


def load_bacteria(only_accession: str | None = None) -> list[dict]:
    gcf_map = json.loads(GCF_MAP_JSON.read_text()) if GCF_MAP_JSON.exists() else {}
    with open(LIST_CSV, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    # Skip unresolved TODOs and invalid accessions
    rows = [r for r in rows if r["source"] not in ("unresolved_TODO", "accession_invalid")]
    # Deduplicate by GCF accession (NZ_CP150073 and NZ_JBWJFR000000000 appear twice)
    seen_gcf = set()
    deduped = []
    for r in rows:
        gcf = r.get("gcf_accession") or gcf_map.get(r["accession"], "")
        r["_gcf"] = gcf
        if gcf and gcf not in seen_gcf:
            seen_gcf.add(gcf)
            deduped.append(r)
        elif not gcf:
            deduped.append(r)
    if only_accession:
        deduped = [r for r in deduped if r["accession"] == only_accession]
    return deduped


def main():
    parser = argparse.ArgumentParser(description="Fetch bacteria genome data via NCBI datasets CLI")
    parser.add_argument("--dry-run", action="store_true", help="Print what would be downloaded without downloading")
    parser.add_argument("--accession", help="Fetch a single accession only")
    args = parser.parse_args()

    if not args.dry_run and not check_datasets_installed():
        print("ERROR: 'datasets' CLI not found.")
        print("Install with:  conda install -n igem2026 -c conda-forge ncbi-datasets-cli")
        sys.exit(1)

    rows = load_bacteria(args.accession)
    if not rows:
        print(f"No matching entries in {LIST_CSV}")
        sys.exit(1)

    print(f"Fetching {len(rows)} bacteria {'[DRY RUN]' if args.dry_run else ''}...")

    stats = {"ok": 0, "skip": 0, "fail": 0}
    failures = []

    for row in rows:
        acc = row["accession"]
        gcf = row.get("_gcf") or row.get("gcf_accession", acc)
        download_acc = gcf if gcf else acc   # use GCF if available
        name = row["organism_name"]
        out_dir = BACTERIA_DIR / acc         # folder named by original seq accession
        print(f"  {acc} ({download_acc})  {name[:45]:<45}", end=" ", flush=True)

        status = fetch_one(download_acc, out_dir, dry_run=args.dry_run)
        tag = status.split(":")[0]
        stats[tag if tag in stats else "fail"] = stats.get(tag if tag in stats else "fail", 0) + 1

        if status == "skip":
            print("[SKIP] already complete")
        elif status == "ok":
            print("[OK]")
        elif status == "dry-run":
            pass
        else:
            reason = status.split(":", 1)[1] if ":" in status else status
            print(f"[FAIL] {reason}")
            failures.append((acc, reason))
            stats["fail"] += 1

    print(f"\nSummary: {stats.get('ok',0)} downloaded, {stats.get('skip',0)} skipped, {stats.get('fail',0)} failed")
    if failures:
        print("\nFailed accessions:")
        for acc, reason in failures:
            print(f"  {acc}: {reason}")


if __name__ == "__main__":
    main()
