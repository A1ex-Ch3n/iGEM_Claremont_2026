"""
Per-protein density plot for Factor #2 (Olivia).

Per project spec §4: 'Generate a density plot of phage vs host values to
inspect peak shifts before regression.' Plots the *per-protein* pI and
acidity% distributions across all phage proteins vs all host proteins.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
REPO = Path(__file__).resolve().parents[3]


def clean(seq):
    return "".join(a for a in seq.upper().rstrip("*") if a in VALID_AA)


def per_protein(path):
    pis, acidity = [], []
    for rec in SeqIO.parse(str(path), "fasta"):
        s = clean(str(rec.seq))
        if len(s) < 5:
            continue
        try:
            pis.append(ProteinAnalysis(s).isoelectric_point())
        except Exception:
            continue
        acidity.append(100 * (s.count("D") + s.count("E")) / len(s))
    return pis, acidity


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phage-runs",
                    default=str(REPO / "02_annotation/outputs/pharokka_runs"))
    ap.add_argument("--host-dir",
                    default=str(REPO / "02_annotation/outputs/host_proteins"))
    ap.add_argument("--out",
                    default=str(REPO / "03_feature_weighting/outputs/per_factor/factor2/f02_pI_acidity_density.png"))
    args = ap.parse_args()

    phage_pi, phage_ac = [], []
    for faa in sorted(Path(args.phage_runs).glob("*/phanotate.faa")):
        pis, ac = per_protein(faa)
        phage_pi += pis; phage_ac += ac

    host_pi, host_ac = [], []
    for faa in sorted(Path(args.host_dir).glob("*/proteins.faa")):
        pis, ac = per_protein(faa)
        host_pi += pis; host_ac += ac

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, (pv, hv), title, xl in [
        (axes[0], (phage_pi, host_pi), "Per-protein isoelectric point (pI)", "pI"),
        (axes[1], (phage_ac, host_ac), "Per-protein acidity % (D+E)",       "acidity %"),
    ]:
        pd.Series(pv).plot(kind="kde", ax=ax, color="#d6604d",
                           label=f"phage proteins (n={len(pv)})")
        pd.Series(hv).plot(kind="kde", ax=ax, color="#4393c3",
                           label=f"host proteins (n={len(hv)})")
        ax.set_title(title); ax.set_xlabel(xl); ax.set_ylabel("density")
        ax.legend(); ax.grid(alpha=0.3)

    fig.suptitle("Factor #2: per-protein pI & Acidity (Olivia)", y=1.02)
    fig.tight_layout()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"wrote {args.out}")
    print(f"phage proteins: {len(phage_pi)}  |  host proteins: {len(host_pi)}")


if __name__ == "__main__":
    main()
