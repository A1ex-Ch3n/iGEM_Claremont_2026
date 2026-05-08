"""
synthetic_data.py — Synthetic dataset generator for Module 06.
合成数据集生成器：在真实 ELISA 数据到来前用于占位和测试。

When real ELISA data arrives (~2026-06-01) replace the label column in the
returned DataFrame with measured log_EC50 values.
真实 ELISA 数据到来（约 2026-06-01）后，用实测 log_EC50 替换合成标签列。
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
import pandas as pd


def generate_synthetic_dataset(
    rbp_embeddings: Optional[np.ndarray] = None,
    receptor_embeddings: Optional[np.ndarray] = None,
    boltz2_priors: Optional[pd.DataFrame] = None,
    noise_std: float = 0.5,
    seed: int = 42,
    n_rbp: int = 20,
    n_receptor: int = 4,
    embed_dim: int = 320,
) -> Tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """
    Build a synthetic (X, y, meta) dataset that mirrors the real pipeline shape.
    生成形状与真实数据一致的合成 (X, y, meta) 数据集。

    Priority / 优先级:
      1. If rbp_embeddings + receptor_embeddings provided: use them.
         提供真实 embedding 时优先使用。
      2. Else: generate random Gaussian embeddings as documented fallback.
         否则用随机高斯 embedding（文档化的 fallback）。

    Label strategy / 标签策略:
      - Pairs present in boltz2_priors → y = predicted_pKd + Gaussian noise.
        在 Boltz-2 先验里的 pair：y = predicted_pKd + 高斯噪声。
      - Pairs NOT in boltz2_priors (or priors absent) → y ~ N(5.0, 1.0).
        不在先验里的 pair（或没有先验时）：y ~ N(5.0, 1.0)，模拟典型 pKd 范围。

    Returns:
        X:    (n_pairs, 2*embed_dim) — concat(rbp_emb, receptor_emb)
        y:    (n_pairs,) — synthetic binding score (proxy for pKd)
        meta: DataFrame with columns rbp_id, receptor_id, true_y, source
    """
    rng = np.random.default_rng(seed)

    # -----------------------------------------------------------------------
    # Step 1: Embeddings / 获取 embedding
    # -----------------------------------------------------------------------
    using_fallback_embeddings = False
    if rbp_embeddings is None or receptor_embeddings is None:
        # Documented fallback: random embeddings (pipeline shape is what matters)
        # 文档化 fallback：随机 embedding（今晚只验证 pipeline 形状，不追求生物精度）
        using_fallback_embeddings = True
        rbp_embeddings = rng.standard_normal((n_rbp, embed_dim)).astype(np.float32)
        receptor_embeddings = rng.standard_normal((n_receptor, embed_dim)).astype(np.float32)

    n_rbp_actual = rbp_embeddings.shape[0]
    n_rec_actual = receptor_embeddings.shape[0]
    actual_embed_dim = rbp_embeddings.shape[1]

    # Build rbp/receptor id lists / 构建 id 列表
    rbp_ids = [f"rbp_{i:02d}" for i in range(n_rbp_actual)]
    rec_ids = [f"rec_{j:02d}" for j in range(n_rec_actual)]

    # -----------------------------------------------------------------------
    # Step 2: Build all (rbp, receptor) pairs / 构建所有 pair
    # -----------------------------------------------------------------------
    rows_rbp, rows_rec = [], []
    for i in range(n_rbp_actual):
        for j in range(n_rec_actual):
            rows_rbp.append(i)
            rows_rec.append(j)

    n_pairs = len(rows_rbp)
    X = np.concatenate(
        [rbp_embeddings[rows_rbp], receptor_embeddings[rows_rec]], axis=1
    ).astype(np.float32)  # (n_pairs, 2*embed_dim)

    # -----------------------------------------------------------------------
    # Step 3: Labels / 标签生成
    # -----------------------------------------------------------------------
    y = np.full(n_pairs, 5.0, dtype=np.float32)
    sources = ["random_prior"] * n_pairs

    if boltz2_priors is not None and not boltz2_priors.empty:
        # Build a lookup: (rbp_id, receptor_id) → predicted_pKd
        # 构建查找表
        prior_lookup = {}
        for _, row in boltz2_priors.iterrows():
            prior_lookup[(str(row["rbp_id"]), str(row["receptor_id"]))] = float(
                row.get("predicted_pKd", 5.0)
            )
        for k in range(n_pairs):
            key = (rbp_ids[rows_rbp[k]], rec_ids[rows_rec[k]])
            if key in prior_lookup:
                y[k] = prior_lookup[key]
                sources[k] = "boltz2_prior"

    # Add Gaussian noise / 加高斯噪声
    noise = rng.normal(0.0, noise_std, size=n_pairs).astype(np.float32)
    y = y + noise

    # -----------------------------------------------------------------------
    # Step 4: Meta DataFrame / 元数据 DataFrame
    # -----------------------------------------------------------------------
    meta = pd.DataFrame(
        {
            "rbp_id": [rbp_ids[rows_rbp[k]] for k in range(n_pairs)],
            "receptor_id": [rec_ids[rows_rec[k]] for k in range(n_pairs)],
            "true_y": y.tolist(),
            "source": sources,
        }
    )
    if using_fallback_embeddings:
        meta["fallback_embeddings"] = True
        import warnings
        warnings.warn(
            "generate_synthetic_dataset: Using RANDOM fallback embeddings "
            "(upstream Module 04/05 outputs not found). "
            "Pipeline shape is correct; biological signal is absent.",
            UserWarning,
            stacklevel=2,
        )

    return X, y, meta
