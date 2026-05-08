"""
test_smoke.py — Smoke test: tiny random data, 2-member ensemble, 10 epochs.
冒烟测试：用极小随机数据、2 成员集成、10 个 epoch 验证 pipeline 能跑通。

DOES NOT execute a notebook — imports functions directly from ensemble.py.
不执行 notebook，直接导入 ensemble.py 中的函数，避免超过 60s 时限。
"""

import sys
from pathlib import Path

import numpy as np
import pytest

# Make processes/ importable regardless of how pytest is invoked
# 无论 pytest 从哪里调用都能 import processes/ 下的模块
PROCESSES_DIR = Path(__file__).resolve().parents[1]
if str(PROCESSES_DIR) not in sys.path:
    sys.path.insert(0, str(PROCESSES_DIR))

from ensemble import DeepEnsemble, train_member
from synthetic_data import generate_synthetic_dataset


def make_tiny_dataset(n: int = 50, dim: int = 640, seed: int = 0):
    """50 samples, input_dim=640 (2x320 concat), random labels."""
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, dim)).astype(np.float32)
    y = rng.standard_normal(n).astype(np.float32)
    split = int(0.8 * n)
    return X[:split], y[:split], X[split:], y[split:]


def test_train_member_smoke():
    """
    train_member runs without error and returns a model + log.
    train_member 能无错运行，返回模型和训练日志。
    """
    X_tr, y_tr, X_val, y_val = make_tiny_dataset()
    model, log = train_member(
        X_tr, y_tr, X_val, y_val,
        seed=0, hidden_dim=32, max_epochs=10, patience=5,
    )
    assert model is not None
    assert "train_nll" in log and len(log["train_nll"]) > 0
    assert log["final_val_nll"] is not None
    assert np.isfinite(log["final_val_nll"]), "val NLL is not finite"


def test_deep_ensemble_smoke():
    """
    2-member DeepEnsemble trains and produces finite predictions.
    2 成员集成能训练并输出有限值预测。
    """
    X_tr, y_tr, X_val, y_val = make_tiny_dataset()

    ens = DeepEnsemble(input_dim=640, hidden_dim=32, n_members=2)
    ens.train(X_tr, y_tr, X_val, y_val, max_epochs=10, patience=5)

    assert len(ens.members) == 2, "Expected 2 members"

    X_test = X_val
    mean, std = ens.predict(X_test)  # total_std returned by default

    assert mean.shape == (len(X_test),), "mean shape mismatch"
    assert std.shape == (len(X_test),), "std shape mismatch"
    assert np.all(np.isfinite(mean)), "Non-finite values in mean"
    assert np.all(np.isfinite(std)), "Non-finite values in std"
    assert np.all(std >= 0), "Negative std values (total_std must be ≥ 0)"


def test_synthetic_data_generator():
    """
    generate_synthetic_dataset returns correct shapes and non-NaN values.
    合成数据生成器返回形状正确、无 NaN 的数据。
    """
    X, y, meta = generate_synthetic_dataset(
        n_rbp=5, n_receptor=3, embed_dim=64, seed=0
    )
    n_pairs = 5 * 3
    assert X.shape == (n_pairs, 128), f"Expected X shape ({n_pairs}, 128), got {X.shape}"
    assert y.shape == (n_pairs,), f"Expected y shape ({n_pairs},), got {y.shape}"
    assert len(meta) == n_pairs
    assert np.all(np.isfinite(X)), "Non-finite values in X"
    assert np.all(np.isfinite(y)), "Non-finite values in y"
