import os
import csv
from pathlib import Path
import matplotlib.pyplot as plt
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
HOST_PROTEINS_DIR = REPO_ROOT / "02_annotation/outputs/host_proteins"
PHAROKKA_RUNS_DIR = REPO_ROOT / "02_annotation/outputs/pharokka_runs"
DATASET_CSV = Path(__file__).parent / "phage_host_dataset.csv"


def parse_faa_avg_length(faa_path):
    lengths = []
    current_len = 0
    with open(faa_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_len > 0:
                    lengths.append(current_len)
                current_len = 0
            else:
                current_len += len(line)
        if current_len > 0:
            lengths.append(current_len)
    return sum(lengths) / len(lengths) if lengths else None


def build_bacteria_table():
    table = {}
    for folder in HOST_PROTEINS_DIR.iterdir():
        if not folder.is_dir():
            continue
        faa = folder / "proteins.faa"
        if faa.exists():
            avg = parse_faa_avg_length(faa)
            if avg is not None:
                table[folder.name] = avg
    return table


def build_phage_table():
    table = {}
    for folder in PHAROKKA_RUNS_DIR.iterdir():
        if not folder.is_dir():
            continue
        faa = folder / "phanotate.faa"
        if faa.exists():
            avg = parse_faa_avg_length(faa)
            if avg is not None:
                table[folder.name] = avg
    return table


def main():
    bacteria_avg = build_bacteria_table()
    phage_avg = build_phage_table()

    xs, ys = [], []
    skipped = []

    with open(DATASET_CSV) as f:
        reader = csv.DictReader(f)
        for row in reader:
            phage_id = row["phage_id"].strip()
            host_id = row["host_id"].strip()
            y = float(row["y"].strip())

            p_len = phage_avg.get(phage_id)
            b_len = bacteria_avg.get(host_id)

            if p_len is None or b_len is None:
                skipped.append((phage_id, host_id))
                continue

            x = min(p_len, b_len) / max(p_len, b_len)
            xs.append(x)
            ys.append(y)

    if skipped:
        print(f"Skipped {len(skipped)} rows (missing data): {skipped}")

    slope, intercept, r_value, p_value, std_err = stats.linregress(xs, ys)

    print(f"r = {r_value:.4f}")
    print(f"p = {p_value:.4e}")
    print(f"slope = {slope:.4f}, intercept = {intercept:.4f}")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(xs, ys, alpha=0.7, edgecolors="k", linewidths=0.5, zorder=3)

    x_line = sorted(xs)
    y_line = [slope * x + intercept for x in x_line]
    ax.plot(x_line, y_line, color="red", linewidth=1.5, label=f"r = {r_value:.4f}")

    ax.set_xlabel("min/max average protein length ratio (phage vs host)")
    ax.set_ylabel("y (interaction score)")
    ax.set_title("Phage-Host Average Protein Length Ratio vs Interaction Score")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.4)

    out_path = Path(__file__).parent / "linear_regression_plot.png"
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.show()
    print(f"Plot saved to {out_path}")


if __name__ == "__main__":
    main()
