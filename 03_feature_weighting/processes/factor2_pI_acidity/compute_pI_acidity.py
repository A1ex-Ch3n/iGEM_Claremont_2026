"""
Factor #2 (Olivia): pI & Acidity, computed from .faa files.

Inputs (new repo layout):
  - Phage proteins: 02_annotation/outputs/pharokka_runs/<acc>/phanotate.faa
  - Host  proteins: 02_annotation/outputs/host_proteins/<name>/proteins.faa
  - Pairs (with y): 02_annotation/outputs/host_proteins/pairs.csv

Outputs (under 03_feature_weighting/outputs/per_factor/factor2/):
  - f02_pI_acidity_per_genome.csv   one row per phage / host genome
  - f02_pI_acidity_pairs.csv        one row per phage-host pair, with
                                    x = |Acidity_phage - Acidity_host| and pI diff
"""

import argparse
import csv
import statistics
from pathlib import Path

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
REPO = Path(__file__).resolve().parents[3]


def clean(seq: str) -> str:
    return "".join(a for a in seq.upper().rstrip("*") if a in VALID_AA)


def summarize_faa(path: Path):
    pis, total_len, d_e, k_r_h, n = [], 0, 0, 0, 0
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
        "f02_pI_median": round(statistics.median(pis), 6),
        "f02_pI_mean": round(statistics.fmean(pis), 6),
        "f02_acidity_pct": round(100 * d_e / total_len, 6),
        "f02_basicity_pct": round(100 * k_r_h / total_len, 6),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phage-runs",
                    default=str(REPO / "02_annotation/outputs/pharokka_runs"))
    ap.add_argument("--host-dir",
                    default=str(REPO / "02_annotation/outputs/host_proteins"))
    ap.add_argument("--pairs",
                    default=str(REPO / "02_annotation/outputs/host_proteins/pairs.csv"))
    ap.add_argument("--out-dir",
                    default=str(REPO / "03_feature_weighting/outputs/per_factor/factor2"))
    args = ap.parse_args()

    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    rows = []

    phage_root = Path(args.phage_runs)
    for faa in sorted(phage_root.glob("*/phanotate.faa")):
        s = summarize_faa(faa)
        if not s: continue
        rows.append({"genome_id": faa.parent.name, "kind": "phage",
                     "source_path": str(faa.relative_to(REPO)), **s})

    host_root = Path(args.host_dir)
    for faa in sorted(host_root.glob("*/proteins.faa")):
        s = summarize_faa(faa)
        if not s: continue
        rows.append({"genome_id": faa.parent.name, "kind": "host",
                     "source_path": str(faa.relative_to(REPO)), **s})

    fields = ["genome_id", "kind", "source_path", "n_proteins",
              "f02_pI_median", "f02_pI_mean", "f02_acidity_pct", "f02_basicity_pct"]
    per_genome = out_dir / "f02_pI_acidity_per_genome.csv"
    with open(per_genome, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
    print(f"wrote {len(rows)} genomes -> {per_genome}")

    idx = {r["genome_id"]: r for r in rows}
    pairs_path = Path(args.pairs)
    pairs_out = out_dir / "f02_pI_acidity_pairs.csv"
    n_pairs = 0
    with open(pairs_path) as f, open(pairs_out, "w", newline="") as g:
        reader = csv.DictReader(f)
        writer = csv.DictWriter(g, fieldnames=[
            "phage_id", "host_id",
            "f02_phage_pI_median", "f02_host_pI_median",
            "f02_phage_acidity_pct", "f02_host_acidity_pct",
            "f02_x_pI", "f02_x_acidity", "y",
        ])
        writer.writeheader()
        missing = set()
        for row in reader:
            p, h = row["phage_id"], row["host_id"]
            if p not in idx or h not in idx:
                missing.add(p if p not in idx else h); continue
            rp, rh = idx[p], idx[h]
            writer.writerow({
                "phage_id": p, "host_id": h,
                "f02_phage_pI_median": rp["f02_pI_median"],
                "f02_host_pI_median":  rh["f02_pI_median"],
                "f02_phage_acidity_pct": rp["f02_acidity_pct"],
                "f02_host_acidity_pct":  rh["f02_acidity_pct"],
                "f02_x_pI":      round(abs(rp["f02_pI_median"] - rh["f02_pI_median"]), 6),
                "f02_x_acidity": round(abs(rp["f02_acidity_pct"] - rh["f02_acidity_pct"]), 6),
                "y": row.get("y", ""),
            })
            n_pairs += 1
    print(f"wrote {n_pairs} pairs -> {pairs_out}")
    if missing:
        print(f"  (skipped pairs missing genome rows for: {sorted(missing)})")


if __name__ == "__main__":
    main()
