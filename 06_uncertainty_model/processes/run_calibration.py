"""
run_calibration.py — Calibration check + reliability diagram for Module 06.
校准检查 + reliability diagram 生成脚本（模块 06）。

Generates:
  outputs/cycle_0/calibration.png   — reliability diagram (300 DPI)

Reference / 参考:
  Guo et al. (2017) "On Calibration of Modern Neural Networks." ICML 2017.
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path
from typing import Tuple

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for script use / 脚本模式使用非交互后端

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

SCRIPT_DIR = Path(__file__).resolve().parent
MODULE_DIR = SCRIPT_DIR.parent
OUTPUTS_DIR = MODULE_DIR / "outputs" / "cycle_0"

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from ensemble import DeepEnsemble, EnsembleMember
from synthetic_data import generate_synthetic_dataset

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


def compute_regression_ece(
    y_true: np.ndarray,
    mean_pred: np.ndarray,
    std_pred: np.ndarray,
    n_bins: int = 10,
) -> Tuple[float, np.ndarray, np.ndarray]:
    """
    Regression ECE: weighted MAE between predicted std and actual RMSE of residuals.
    回归 ECE：预测 std 与残差 RMSE 的加权 MAE（按 predicted std 分箱）。

    Returns:
        ece:           scalar ECE value
        bin_pred_std:  mean predicted std per bin (x-axis of reliability diagram)
        bin_actual_std: actual RMSE per bin (y-axis of reliability diagram)
    """
    bins = np.percentile(std_pred, np.linspace(0, 100, n_bins + 1))
    bins = np.unique(bins)
    if len(bins) < 2:
        return 0.0, np.array([]), np.array([])

    bin_pred_std = []
    bin_actual_std = []
    total_weight = 0.0
    ece = 0.0

    for i in range(len(bins) - 1):
        lo, hi = bins[i], bins[i + 1]
        mask = (std_pred >= lo) & (std_pred < hi)
        if i == len(bins) - 2:
            mask = (std_pred >= lo) & (std_pred <= hi)
        if mask.sum() == 0:
            continue
        residuals = y_true[mask] - mean_pred[mask]
        actual_rmse = float(np.sqrt(np.mean(residuals ** 2)))
        pred_std_bin = float(std_pred[mask].mean())
        weight = mask.sum()
        ece += weight * abs(pred_std_bin - actual_rmse)
        total_weight += weight
        bin_pred_std.append(pred_std_bin)
        bin_actual_std.append(actual_rmse)

    ece_val = ece / total_weight if total_weight > 0 else 0.0
    return ece_val, np.array(bin_pred_std), np.array(bin_actual_std)


def compute_coverage_curve(
    y_true: np.ndarray,
    mean_pred: np.ndarray,
    std_pred: np.ndarray,
    levels: np.ndarray | None = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Coverage calibration: for each confidence level, measure actual coverage.
    覆盖率校准：对每个置信水平，测量实际覆盖率。

    Returns:
        expected_coverage: array of nominal CI levels (x-axis)
        actual_coverage:   array of actual fraction of y_true in CI (y-axis)
    """
    if levels is None:
        levels = np.linspace(0.05, 0.95, 19)
    actual = []
    for p in levels:
        z = float(np.abs(np.percentile(np.random.normal(0, 1, 10000), 100 * (1 - (1 - p) / 2))))
        # Compute z-score for two-sided CI at level p
        # 计算 p 水平双侧置信区间的 z 分数
        from scipy import stats as sp_stats
        z = sp_stats.norm.ppf(0.5 + p / 2)
        lo = mean_pred - z * std_pred
        hi = mean_pred + z * std_pred
        coverage = float(np.mean((y_true >= lo) & (y_true <= hi)))
        actual.append(coverage)
    return levels, np.array(actual)


def plot_reliability_diagram(
    y_true: np.ndarray,
    mean_pred: np.ndarray,
    std_pred: np.ndarray,
    ece: float,
    save_path: Path,
    cycle: int = 0,
) -> None:
    """
    Two-panel reliability diagram:
    双面板 reliability diagram：
      Left:  scatter of predicted std vs |residual| per sample
             左图：每样本的预测 std vs 绝对残差散点图
      Right: coverage calibration curve (expected vs actual)
             右图：覆盖率校准曲线（期望 vs 实际）
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(
        f"Cycle {cycle} Uncertainty Calibration\n"
        f"(Synthetic data, ECE = {ece:.4f})",
        fontsize=13,
    )

    # --- Left panel: predicted std vs |residual| ---
    # 左图：预测 std vs 绝对残差
    residuals = np.abs(y_true - mean_pred)
    ax = axes[0]
    ax.scatter(std_pred, residuals, alpha=0.6, s=20, color="#3A86FF", edgecolors="none")

    # Diagonal = perfect calibration / 对角线 = 完美校准
    lim = max(std_pred.max(), residuals.max()) * 1.05
    ax.plot([0, lim], [0, lim], "k--", lw=1.5, label="Perfect calibration")

    # Bin-level overlay / 箱级叠加线
    ece_val, bin_pred, bin_actual = compute_regression_ece(y_true, mean_pred, std_pred, n_bins=8)
    if len(bin_pred):
        ax.plot(bin_pred, bin_actual, "ro-", lw=2, markersize=5, label="Bin avg")

    ax.set_xlabel("Predicted std (total uncertainty)")
    ax.set_ylabel("|Residual|")
    ax.set_title(f"Std vs |Residual|  ECE={ece:.3f}")
    ax.legend(fontsize=9)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, max(residuals.max(), lim) * 1.05)

    # --- Right panel: coverage calibration ---
    # 右图：覆盖率校准
    try:
        levels, actuals = compute_coverage_curve(y_true, mean_pred, std_pred)
        ax2 = axes[1]
        ax2.plot(levels, actuals, "b-o", markersize=4, label="Actual coverage")
        ax2.plot([0, 1], [0, 1], "k--", lw=1.5, label="Perfect calibration")
        ax2.fill_between(levels, actuals, levels,
                         where=(actuals < levels), alpha=0.15, color="red",
                         label="Under-confident")
        ax2.fill_between(levels, actuals, levels,
                         where=(actuals > levels), alpha=0.15, color="green",
                         label="Over-confident")
        ax2.set_xlabel("Expected coverage (confidence level)")
        ax2.set_ylabel("Actual coverage")
        ax2.set_title("Coverage Calibration")
        ax2.legend(fontsize=8)
        ax2.set_xlim(0, 1)
        ax2.set_ylim(0, 1)
    except Exception as e:
        axes[1].text(0.5, 0.5, f"scipy not available:\n{e}",
                     ha="center", va="center", transform=axes[1].transAxes)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    logger.info("Saved calibration plot to %s", save_path)


def main() -> None:
    """
    Load cycle_0 outputs and produce calibration.png + log ECE.
    加载 cycle_0 输出，生成 calibration.png 并记录 ECE。
    """
    logger.info("=== Module 06 Calibration Check (Cycle 0) ===")

    # Load model metadata to reconstruct ensemble
    # 读取模型元数据以重建集成
    meta_path = OUTPUTS_DIR / "model_meta.json"
    if not meta_path.exists():
        raise FileNotFoundError(
            f"model_meta.json not found at {meta_path}. "
            "Run run_cycle0.py first."
        )
    with open(meta_path) as f:
        meta = json.load(f)

    input_dim  = meta["input_dim"]
    hidden_dim = meta["hidden_dim"]
    n_members  = meta["n_members"]
    embed_dim  = meta["embed_dim_per_protein"]
    seed       = 42

    logger.info(
        "Reconstructing ensemble: input_dim=%d, hidden_dim=%d, members=%d",
        input_dim, hidden_dim, n_members,
    )

    # Reconstruct and reload members
    # 重建并重载每个成员
    ens = DeepEnsemble(input_dim=input_dim, hidden_dim=hidden_dim, n_members=n_members)
    members = []
    for i in range(n_members):
        pt_path = OUTPUTS_DIR / f"ensemble_member_{i}.pt"
        if not pt_path.exists():
            raise FileNotFoundError(f"Missing {pt_path}")
        m = EnsembleMember(input_dim=input_dim, hidden_dim=hidden_dim)
        m.load_state_dict(torch.load(pt_path, map_location="cpu", weights_only=True))
        m.eval()
        members.append(m)
    ens.members = members

    # Re-generate the same synthetic dataset (same seed → same pairs, same y)
    # 用相同 seed 重新生成合成数据（相同 pairs 和标签）
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        X, y, data_meta = generate_synthetic_dataset(
            seed=seed,
            embed_dim=embed_dim,
        )

    # Hold-out: same 80/20 split
    # 留出集：相同的 80/20 拆分
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(X))
    split = int(0.8 * len(X))
    val_idx = idx[split:]
    X_val = X[val_idx]
    y_val = y[val_idx]
    val_meta = data_meta.iloc[val_idx].reset_index(drop=True)

    logger.info("Running predictions on %d val samples…", len(X_val))
    mean_pred, std_pred = ens.predict(X_val)

    # ECE
    ece, bin_pred, bin_actual = compute_regression_ece(y_val, mean_pred, std_pred, n_bins=8)
    logger.info("ECE (regression, binned) = %.4f", ece)

    if ece > 0.1:
        logger.warning(
            "ECE = %.4f > 0.1. Consider temperature scaling after real ELISA "
            "data arrives. See AGENT_REPORT.md §calibration.", ece
        )

    # Plot
    plot_reliability_diagram(
        y_true=y_val,
        mean_pred=mean_pred,
        std_pred=std_pred,
        ece=ece,
        save_path=OUTPUTS_DIR / "calibration.png",
        cycle=0,
    )

    # Save calibration metrics into model_meta.json
    # 将校准指标写入 model_meta.json
    meta["calibration_ece_synthetic"] = round(float(ece), 6)
    meta["calibration_n_val"] = len(X_val)
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    logger.info("Updated model_meta.json with calibration metrics")

    logger.info("=== Calibration check complete ===")


if __name__ == "__main__":
    main()
