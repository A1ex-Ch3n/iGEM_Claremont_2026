"""
Factor #2 (Olivia): pI & Acidity per genome, computed from .faa files.

Outputs are join-ready for the team's integrated feature matrix:
- Primary key: `genome_id` (NCBI accession = folder name under ncbi_dataset/data/).
- Feature columns are namespaced `f02_*` so they will not collide with
  Alex/Weitao/Sarah/Angela/Carol's columns when merged.

Usage:
    python compute_pI_acidity.py \
        --phage-dir "ncbi_dataset/data" \
        --host-faa "Bacterial_Protein_Annotation_Project/results/real_output/proteins_annotated.faa" \
        --out f02_pI_acidity_per_genome.csv

Optional pairwise output (when the Step 1 interaction CSV is available):
    python compute_pI_acidity.py --pairs interactions.csv
"""

import argparse
import csv
import statistics
from pathlib import Path

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def clean(seq: str) -> str:
    return "".join(a for a in seq.upper().rstrip("*") if a in VALID_AA)


def summarize_faa(path: Path):
    pis = []
    total_len = 0
    d_e = 0
    k_r_h = 0
    n = 0
    for rec in SeqIO.parse(str(path), "fasta"):
        s = clean(str(rec.seq))
        if len(s) < 5:
            continue
        try:
            pis.append(ProteinAnalysis(s).isoelectric_point())
        except Exception:
            continue
        total_len += len(s)
        d_e += s.count("D") + s.count("E")
        k_r_h += s.count("K") + s.count("R") + s.count("H")
        n += 1
    if n == 0 or total_len == 0:
        return None
    return {
        "n_proteins": n,
        "f02_pI_median": statistics.median(pis),
        "f02_pI_mean": statistics.fmean(pis),
        "f02_acidity_pct": 100 * d_e / total_len,
        "f02_basicity_pct": 100 * k_r_h / total_len,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phage-dir", default="ncbi_dataset/data",
                    help="Root dir; every *.faa under it is treated as a phage genome. "
                         "genome_id = parent folder name.")
    ap.add_argument("--host-faa", action="append", default=[],
                    help="Path to a host .faa; repeatable. genome_id = file stem.")
    ap.add_argument("--out", default="f02_pI_acidity_per_genome.csv")
    ap.add_argument("--pairs", default=None,
                    help="Optional CSV with columns phage_id,host_id[,y]. "
                         "Emits the pairwise feature CSV.")
    ap.add_argument("--pairs-out", default="f02_pI_acidity_pairs.csv")
    args = ap.parse_args()

    rows = []
    phage_root = Path(args.phage_dir)
    if phage_root.is_dir():
        for faa in sorted(phage_root.rglob("*.faa")):
            s = summarize_faa(faa)
            if not s:
                continue
            rows.append({
                "genome_id": faa.parent.name,
                "kind": "phage",
                "source_path": str(faa),
                **s,
            })

    for hp in args.host_faa:
        p = Path(hp)
        s = summarize_faa(p)
        if not s:
            continue
        rows.append({
            "genome_id": p.stem,
            "kind": "host",
            "source_path": str(p),
            **s,
        })

    fields = [
        "genome_id", "kind", "source_path", "n_proteins",
        "f02_pI_median", "f02_pI_mean", "f02_acidity_pct", "f02_basicity_pct",
    ]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} genomes -> {args.out}")

    if args.pairs:
        idx = {r["genome_id"]: r for r in rows}
        n_pairs = 0
        with open(args.pairs) as f, open(args.pairs_out, "w", newline="") as g:
            reader = csv.DictReader(f)
            writer = csv.DictWriter(
                g,
                fieldnames=["phage_id", "host_id", "f02_x_pI", "f02_x_acidity", "y"],
            )
            writer.writeheader()
            for row in reader:
                p, h = row["phage_id"], row["host_id"]
                if p not in idx or h not in idx:
                    continue
                writer.writerow({
                    "phage_id": p,
                    "host_id": h,
                    "f02_x_pI": abs(idx[p]["f02_pI_median"] - idx[h]["f02_pI_median"]),
                    "f02_x_acidity": abs(idx[p]["f02_acidity_pct"] - idx[h]["f02_acidity_pct"]),
                    "y": row.get("y", ""),
                })
                n_pairs += 1
        print(f"wrote {n_pairs} pairs -> {args.pairs_out}")


if __name__ == "__main__":
    main()
