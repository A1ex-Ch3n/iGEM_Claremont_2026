"""
Scatter check for Factor #2 (Olivia), per project spec §4 step 2:
'Plot each factor's x against y to confirm the linear relationship before
running regression.'

Reads f02_pI_acidity_pairs.csv and plots x_acidity vs y (and x_pI vs y).
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp",
                    default=str(REPO / "03_feature_weighting/outputs/per_factor/factor2/f02_pI_acidity_pairs.csv"))
    ap.add_argument("--out",
                    default=str(REPO / "03_feature_weighting/outputs/per_factor/factor2/f02_acidity_vs_y_scatter.png"))
    args = ap.parse_args()

    df = pd.read_csv(args.inp)
    df["y"] = pd.to_numeric(df["y"], errors="coerce")
    df = df.dropna(subset=["y"])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, xcol, title in [
        (axes[0], "f02_x_acidity", "x_acidity = |Acid%_phage - Acid%_host|"),
        (axes[1], "f02_x_pI",      "x_pI = |pI_phage - pI_host|"),
    ]:
        x, y = df[xcol].values, df["y"].values
        colors = ["#d6604d" if v == 1 else "#4393c3" for v in y]
        ax.scatter(x, y, c=colors, s=55, alpha=0.75, edgecolors="white")
        if len(x) >= 2 and np.std(x) > 0:
            a, b = np.polyfit(x, y, 1)
            xx = np.linspace(x.min(), x.max(), 50)
            ax.plot(xx, a * xx + b, "k--", lw=1, label=f"fit: y = {a:.3f}x + {b:.3f}")
            r = np.corrcoef(x, y)[0, 1]
            ax.text(0.04, 0.93, f"Pearson r = {r:.3f}", transform=ax.transAxes,
                    fontsize=9, va="top",
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
            ax.legend(loc="upper right", fontsize=8)
        ax.set_title(title); ax.set_xlabel(xcol); ax.set_ylabel("y (infection)")
        ax.set_ylim(-0.1, 1.1); ax.grid(alpha=0.3)

    fig.suptitle(f"Factor #2 linearity check  (n_pairs = {len(df)})", y=1.02)
    fig.tight_layout()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"wrote {args.out}  ({len(df)} pairs)")


if __name__ == "__main__":
    main()
