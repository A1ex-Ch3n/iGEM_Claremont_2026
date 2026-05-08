"""
test_diversity.py — Hard check: ensemble members MUST disagree on held-out data.
硬性检查：集成成员必须对留出数据产生分歧（std > 0）。

Motivation / 动机:
  If all 5 members predict identically, the uncertainty signal is zero everywhere,
  and the BALD acquisition function (Module 07) becomes useless.
  如果 5 个成员预测完全相同，不确定性信号处处为零，
  Module 07 的 BALD acquisition function 就毫无意义。

Hard requirement (from AGENT_TODO):
  std > 0 for ≥ 50% of held-out points.
  至少 50% 的留出样本 std > 0。
"""

import sys
from pathlib import Path

import numpy as np

PROCESSES_DIR = Path(__file__).resolve().parents[1]
if str(PROCESSES_DIR) not in sys.path:
    sys.path.insert(0, str(PROCESSES_DIR))

from ensemble import DeepEnsemble


def test_ensemble_diversity():
    """
    Train 5 members with seeds 0-4; ≥50% of held-out predictions must have std>0.
    用 seed 0-4 训练 5 个成员；至少 50% 留出样本的 std 必须 > 0。
    """
    rng = np.random.default_rng(99)
    n_train, n_test, dim = 80, 20, 64  # small but non-trivial
    X_tr = rng.standard_normal((n_train, dim)).astype(np.float32)
    y_tr = rng.standard_normal(n_train).astype(np.float32)
    X_test = rng.standard_normal((n_test, dim)).astype(np.float32)
    X_val = rng.standard_normal((10, dim)).astype(np.float32)
    y_val = rng.standard_normal(10).astype(np.float32)

    ens = DeepEnsemble(input_dim=dim, hidden_dim=64, n_members=5)
    ens.train(X_tr, y_tr, X_val, y_val, max_epochs=50, patience=10, batch_size=16)

    assert len(ens.members) == 5, "Expected exactly 5 ensemble members"

    mean, total_std, epistemic_std = ens.predict(X_test, return_epistemic=True)

    assert mean.shape == (n_test,)
    assert total_std.shape == (n_test,)
    assert epistemic_std.shape == (n_test,)
    assert np.all(np.isfinite(mean)), "Non-finite mean predictions"
    assert np.all(np.isfinite(total_std)), "Non-finite total_std predictions"
    assert np.all(total_std >= 0), "Negative total_std — impossible"

    # The diversity check uses epistemic_std (pure disagreement between means)
    # 多样性检查使用 epistemic_std（成员间均值的纯分歧）
    fraction_diverse = (epistemic_std > 1e-8).mean()
    assert fraction_diverse >= 0.50, (
        f"Only {fraction_diverse:.1%} of predictions have epistemic_std > 0 "
        f"(need ≥ 50%). Check that per-member seeds differ and "
        f"weight init is truly random."
    )


def test_members_have_different_weights():
    """
    Directly check that two members trained with different seeds differ in params.
    直接检查：用不同 seed 训练的两个成员参数不同。
    """
    import torch
    from ensemble import train_member

    rng = np.random.default_rng(7)
    dim = 32
    X = rng.standard_normal((40, dim)).astype(np.float32)
    y = rng.standard_normal(40).astype(np.float32)
    X_v = rng.standard_normal((10, dim)).astype(np.float32)
    y_v = rng.standard_normal(10).astype(np.float32)

    m0, _ = train_member(X, y, X_v, y_v, seed=0, hidden_dim=32, max_epochs=5)
    m1, _ = train_member(X, y, X_v, y_v, seed=1, hidden_dim=32, max_epochs=5)

    # At least one layer's weights must differ
    # 至少一层的权重必须不同
    params_0 = list(m0.parameters())
    params_1 = list(m1.parameters())
    any_different = any(
        not torch.equal(p0, p1) for p0, p1 in zip(params_0, params_1)
    )
    assert any_different, (
        "Members with seed=0 and seed=1 have identical weights — "
        "seeding is broken."
    )
