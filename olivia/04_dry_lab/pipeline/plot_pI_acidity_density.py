"""
Density plots for Olivia's Factor #2 (pI & Acidity), per project spec §4:
"Generate a density plot of phage vs host values to inspect peak shifts
before regression."

Reads f02_pI_acidity_per_genome.csv. If host rows are absent (currently true —
no bacterial host .faa exists in the workspace yet), only the phage
distribution is shown so the plot can be regenerated once hosts are added.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def kde_or_hist(ax, values, label, color):
    values = values.dropna()
    if len(values) >= 2:
        try:
            values.plot(kind="kde", ax=ax, label=label, color=color)
            return
        except Exception:
            pass
    ax.hist(values, bins=20, alpha=0.5, label=label, color=color, density=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default="f02_pI_acidity_per_genome.csv")
    ap.add_argument("--out", default="f02_pI_acidity_density.png")
    args = ap.parse_args()

    df = pd.read_csv(args.inp)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    for ax, col, title in [
        (axes[0], "f02_pI_median", "Median pI per genome"),
        (axes[1], "f02_acidity_pct", "Acidity % (D+E) per genome"),
    ]:
        for kind, color in [("phage", "#d6604d"), ("host", "#4393c3")]:
            sub = df.loc[df["kind"] == kind, col]
            if len(sub) == 0:
                continue
            kde_or_hist(ax, sub, f"{kind} (n={len(sub)})", color)
        ax.set_title(title)
        ax.set_xlabel(col)
        ax.set_ylabel("density")
        ax.legend()
        ax.grid(alpha=0.3)

    fig.suptitle("Factor #2: pI & Acidity distribution (Olivia)", y=1.02)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"wrote {args.out}")

    # Print quick numeric summary so the team can sanity-check without opening the PNG
    summary = df.groupby("kind")[["f02_pI_median", "f02_acidity_pct",
                                   "f02_basicity_pct"]].agg(["mean", "std", "count"])
    print(summary.round(3))


if __name__ == "__main__":
    main()
