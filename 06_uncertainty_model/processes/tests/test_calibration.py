"""
test_calibration.py — Calibration quality test.
校准质量测试。

Strategy / 策略:
  Generate well-calibrated synthetic data (label = linear fn + Gaussian noise
  with KNOWN std). Train ensemble. Compute ECE. Assert ECE < 0.2.

  生成校准良好的合成数据（标签 = 线性函数 + 已知标准差的高斯噪声）。
  训练集成。计算 ECE。断言 ECE < 0.2。

ECE definition for regression / 回归 ECE 定义:
  Bin predictions by predicted std; in each bin compare the predicted std to
  the actual RMSE of residuals. ECE = weighted mean |predicted_std - actual_std|.
  按预测 std 分箱；在每箱内比较预测 std 与残差 RMSE。
  ECE = 加权平均 |predicted_std - actual_std|。

  Reference / 参考: Guo et al. (2017) "On Calibration of Modern Neural Networks."
  ICML 2017. (adapted to regression setting)
"""

import sys
from pathlib import Path

import numpy as np

PROCESSES_DIR = Path(__file__).resolve().parents[1]
if str(PROCESSES_DIR) not in sys.path:
    sys.path.insert(0, str(PROCESSES_DIR))

from ensemble import DeepEnsemble


def compute_regression_ece(
    y_true: np.ndarray,
    mean_pred: np.ndarray,
    std_pred: np.ndarray,
    n_bins: int = 5,
) -> float:
    """
    Regression ECE: weighted mean absolute difference between predicted std
    and actual RMSE of residuals, binned by predicted std.
    回归 ECE：按预测 std 分箱后，预测 std 与残差 RMSE 的加权平均绝对差。

    Uses RMSE (not std of absolute errors) to match the prediction std semantics.
    使用 RMSE 而非绝对误差的 std，与预测 std 语义一致。

    Reference: Guo et al. 2017 (adapted for regression).
    """
    bins = np.percentile(std_pred, np.linspace(0, 100, n_bins + 1))
    bins = np.unique(bins)
    if len(bins) < 2:
        return 0.0

    total_weight = 0.0
    ece = 0.0
    for i in range(len(bins) - 1):
        lo, hi = bins[i], bins[i + 1]
        mask = (std_pred >= lo) & (std_pred < hi)
        if i == len(bins) - 2:
            mask = (std_pred >= lo) & (std_pred <= hi)
        if mask.sum() == 0:
            continue
        # Actual RMSE of residuals in this bin (correct denominator)
        # 本箱内残差的 RMSE（正确分母）
        residuals = y_true[mask] - mean_pred[mask]
        actual_rmse = float(np.sqrt(np.mean(residuals ** 2)))
        pred_std = float(std_pred[mask].mean())
        weight = mask.sum()
        ece += weight * abs(pred_std - actual_rmse)
        total_weight += weight

    return ece / total_weight if total_weight > 0 else 0.0


def test_calibration_ece_below_threshold():
    """
    On well-calibrated synthetic data, ensemble ECE must be < 0.2.
    在校准良好的合成数据上，集成 ECE 必须 < 0.2。
    """
    rng = np.random.default_rng(123)
    # Use small dim=8 for fast convergence (model learns the linear signal quickly)
    # 使用 dim=8 以便模型快速收敛（线性信号简单，噪声结构清晰）
    # Small dim=8 for fast convergence; n_train=200 keeps training under 60s
    # 使用 dim=8 快速收敛；n_train=200 保证测试在 60s 内完成
    n_train, n_val, n_test, dim = 200, 40, 80, 8
    TRUE_NOISE_STD = 0.5  # known true noise / 已知真实噪声标准差

    # Well-calibrated data: y = linear projection + known Gaussian noise
    # 校准良好数据：y = 线性投影 + 已知高斯噪声
    W = rng.standard_normal(dim).astype(np.float32) * 0.3
    X_all = rng.standard_normal((n_train + n_val + n_test, dim)).astype(np.float32)
    y_all = (X_all @ W + rng.normal(0, TRUE_NOISE_STD, n_train + n_val + n_test)).astype(np.float32)

    X_tr, y_tr = X_all[:n_train], y_all[:n_train]
    X_val, y_val = X_all[n_train:n_train + n_val], y_all[n_train:n_train + n_val]
    X_test, y_test = X_all[-n_test:], y_all[-n_test:]

    ens = DeepEnsemble(input_dim=dim, hidden_dim=64, n_members=5)
    ens.train(X_tr, y_tr, X_val, y_val, max_epochs=100, patience=15, batch_size=64)

    mean_pred, std_pred = ens.predict(X_test)  # total_std (epistemic + aleatoric)

    # Guard: predictions must be finite
    # 保护：预测值必须有限
    assert np.all(np.isfinite(mean_pred)), "Non-finite mean predictions"
    assert np.all(np.isfinite(std_pred)), "Non-finite std predictions"
    assert np.all(std_pred >= 0), "Negative std"

    ece = compute_regression_ece(y_test, mean_pred, std_pred, n_bins=5)
    assert ece < 0.2, (
        f"ECE = {ece:.4f} ≥ 0.2. Ensemble may be poorly calibrated on "
        f"this synthetic dataset. Consider adjusting noise_std or "
        f"temperature scaling."
    )
