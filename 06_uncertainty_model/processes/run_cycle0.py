"""
run_cycle0.py — CLI runner for Cycle 0 deep ensemble training on synthetic data.
Cycle 0 深度集成训练的命令行执行脚本（合成数据占位）。

This script is called by the notebook to actually produce the cycle_0 outputs.
Notebook imports from ensemble.py and synthetic_data.py for the same logic.
此脚本由 notebook 调用，用于实际生成 cycle_0 输出文件。
Notebook 引用 ensemble.py 和 synthetic_data.py 中相同的逻辑。

Usage:
    python run_cycle0.py [--hidden_dim 256] [--n_members 5] [--max_epochs 200]

Output:
    06_uncertainty_model/outputs/cycle_0/
        ensemble_member_{0..4}.pt
        predictions.csv
        training_log.json
        model_meta.json
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import logging
import subprocess
import sys
import warnings
from pathlib import Path
from typing import Any, Dict

import numpy as np
import pandas as pd
import torch

# Anchor to repo root / 锚定到仓库根目录
SCRIPT_DIR = Path(__file__).resolve().parent
MODULE_DIR = SCRIPT_DIR.parent          # 06_uncertainty_model/
OUTPUTS_DIR = MODULE_DIR / "outputs" / "cycle_0"
OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)

# Add processes/ to path for imports / 将 processes/ 加入 path
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from ensemble import DeepEnsemble, train_member
from synthetic_data import generate_synthetic_dataset

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


def _get_repo_commit() -> str:
    """Get current git commit sha (short). 获取当前 git commit sha（短）。"""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, cwd=MODULE_DIR
        )
        return result.stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def _sha256_file(path: Path) -> str:
    """Compute SHA-256 of a file. 计算文件 SHA-256。"""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def main(
    hidden_dim: int = 256,
    n_members: int = 5,
    max_epochs: int = 200,
    lr: float = 1e-4,
    patience: int = 10,
    batch_size: int = 32,
    noise_std: float = 0.5,
    embed_dim: int = 320,
    seed: int = 42,
) -> None:
    """
    Train deep ensemble on synthetic data, save all cycle_0 outputs.
    在合成数据上训练深度集成，保存所有 cycle_0 输出文件。
    """
    logger.info("=== Module 06 Cycle 0 Training (Synthetic Data) ===")
    logger.info("Hidden dim=%d, members=%d, max_epochs=%d", hidden_dim, n_members, max_epochs)

    # -----------------------------------------------------------------------
    # Step 1: Generate synthetic data (documented fallback path)
    # 步骤 1: 生成合成数据（文档化的 fallback 路径）
    # -----------------------------------------------------------------------
    logger.info("Generating synthetic dataset (n_rbp=20, n_rec=4, embed_dim=%d)…", embed_dim)

    # Try to load upstream embeddings; fall back to random if unavailable
    # 尝试加载上游 embedding；不可用时退回到随机 embedding
    REPO_ROOT = MODULE_DIR.parent
    emb_path = REPO_ROOT / "04_protein_embedding" / "outputs"

    rbp_emb: np.ndarray | None = None
    rec_emb: np.ndarray | None = None

    if emb_path.exists():
        npz_files = list(emb_path.glob("embeddings_*.npz"))
        if npz_files:
            logger.info("Found %d npz file(s) in Module 04 outputs — using first.", len(npz_files))
            try:
                data = np.load(npz_files[0], allow_pickle=True)
                rbp_emb = data["array"][:20].astype(np.float32)
                embed_dim = rbp_emb.shape[1]
                # Receptor embeddings: use second half or fallback
                if len(data["array"]) > 20:
                    rec_emb = data["array"][20:24].astype(np.float32)
                else:
                    rec_emb = None
            except Exception as e:
                logger.warning("Failed to load Module 04 embeddings: %s. Using fallback.", e)

    with warnings.catch_warnings(record=True) as w_list:
        warnings.simplefilter("always")
        X, y, meta = generate_synthetic_dataset(
            rbp_embeddings=rbp_emb,
            receptor_embeddings=rec_emb,
            noise_std=noise_std,
            seed=seed,
            embed_dim=embed_dim,
        )
    if w_list:
        for w in w_list:
            logger.warning(str(w.message))

    logger.info(
        "Dataset: %d pairs, input_dim=%d, y range=[%.2f, %.2f]",
        len(X), X.shape[1], float(y.min()), float(y.max()),
    )

    # -----------------------------------------------------------------------
    # Step 2: Train / val split (80/20)
    # 步骤 2: 训练/验证集拆分（80/20）
    # -----------------------------------------------------------------------
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(X))
    split = int(0.8 * len(X))
    train_idx, val_idx = idx[:split], idx[split:]

    X_train, y_train = X[train_idx], y[train_idx]
    X_val,   y_val   = X[val_idx],   y[val_idx]
    train_meta = meta.iloc[train_idx].reset_index(drop=True)
    val_meta   = meta.iloc[val_idx].reset_index(drop=True)

    logger.info("Train size=%d, val size=%d", len(X_train), len(X_val))

    # -----------------------------------------------------------------------
    # Step 3: Train DeepEnsemble
    # 步骤 3: 训练深度集成
    # -----------------------------------------------------------------------
    input_dim = X_train.shape[1]
    ens = DeepEnsemble(input_dim=input_dim, hidden_dim=hidden_dim, n_members=n_members)
    ens.train(
        X_train, y_train, X_val, y_val,
        max_epochs=max_epochs,
        lr=lr,
        patience=patience,
        batch_size=batch_size,
    )

    # -----------------------------------------------------------------------
    # Step 4: Save ensemble member state_dicts
    # 步骤 4: 保存每个成员的 state_dict
    # -----------------------------------------------------------------------
    for i, model in enumerate(ens.members):
        path = OUTPUTS_DIR / f"ensemble_member_{i}.pt"
        torch.save(model.state_dict(), path)
        logger.info("Saved %s", path)

    # -----------------------------------------------------------------------
    # Step 5: Sanity assertions (these will raise if something is wrong)
    # 步骤 5: 合理性断言（出错时直接 raise）
    # -----------------------------------------------------------------------
    assert len(ens.members) == n_members, f"Expected {n_members} members, got {len(ens.members)}"

    pred_mean, pred_std = ens.predict(X_val)
    assert np.all(np.isfinite(pred_mean)), "Non-finite mean predictions on val set"
    assert np.all(np.isfinite(pred_std)), "Non-finite std predictions on val set"
    _, _total_std, pred_std_epistemic = ens.predict(X_val, return_epistemic=True)
    frac_diverse = float((pred_std_epistemic > 1e-8).mean())
    logger.info("Val std > 0: %.1f%% of samples (diversity check)", frac_diverse * 100)
    assert frac_diverse > 0.50, (
        f"Only {frac_diverse:.1%} of val predictions have std > 0. "
        "Ensemble members may be converging to same solution."
    )

    # -----------------------------------------------------------------------
    # Step 6: Generate predictions CSV on ALL pairs (train + val)
    # 步骤 6: 在全部 pair 上生成 predictions.csv（训练集 + 验证集）
    # -----------------------------------------------------------------------
    all_mean, all_std, all_epistemic_std = ens.predict(X, return_epistemic=True)
    repo_sha = _get_repo_commit()
    model_version = f"{repo_sha}_cycle_0"

    predictions_df = pd.DataFrame(
        {
            "rbp_id":          meta["rbp_id"].tolist(),
            "receptor_id":     meta["receptor_id"].tolist(),
            "predicted_score": np.round(all_mean, 6).tolist(),
            "std":             np.round(all_std, 6).tolist(),
            # epistemic_std = Std_k[mu_k] = BALD score for Module 07
            # epistemic_std = 各成员均值标准差 = Module 07 的 BALD 打分依据
            "epistemic_std":   np.round(all_epistemic_std, 6).tolist(),
            "lower_95":        np.round(all_mean - 1.96 * all_std, 6).tolist(),
            "upper_95":        np.round(all_mean + 1.96 * all_std, 6).tolist(),
            "model_version":   model_version,
        }
    )
    preds_path = OUTPUTS_DIR / "predictions.csv"
    predictions_df.to_csv(preds_path, index=False)
    logger.info("Saved predictions.csv (%d rows)", len(predictions_df))

    # -----------------------------------------------------------------------
    # Step 7: Save training_log.json
    # 步骤 7: 保存 training_log.json
    # -----------------------------------------------------------------------
    log_path = OUTPUTS_DIR / "training_log.json"
    with open(log_path, "w") as f:
        json.dump(ens.logs, f, indent=2)
    logger.info("Saved training_log.json")

    # -----------------------------------------------------------------------
    # Step 8: Save model_meta.json
    # 步骤 8: 保存 model_meta.json
    # -----------------------------------------------------------------------
    meta_dict: Dict[str, Any] = {
        "arch":                "DeepEnsemble_MLP3",
        "hidden_dim":          hidden_dim,
        "n_members":           n_members,
        "input_dim":           input_dim,
        "train_size":          len(X_train),
        "val_size":            len(X_val),
        "max_epochs":          max_epochs,
        "lr":                  lr,
        "patience":            patience,
        "batch_size":          batch_size,
        "noise_std":           noise_std,
        "embed_dim_per_protein": embed_dim,
        "data_source":         "synthetic_fallback_random",
        "repo_commit_sha":     repo_sha,
        "timestamp_utc":       datetime.datetime.now(datetime.timezone.utc).isoformat().replace("+00:00", "Z"),
        "train_time_s":        ens.meta.get("train_time_s", -1),
        "torch_version":       torch.__version__,
        "numpy_version":       np.__version__,
        "python_version":      sys.version.split()[0],
    }
    meta_path = OUTPUTS_DIR / "model_meta.json"
    with open(meta_path, "w") as f:
        json.dump(meta_dict, f, indent=2)
    logger.info("Saved model_meta.json")

    # -----------------------------------------------------------------------
    # Step 9: Update MANIFEST.csv
    # 步骤 9: 更新 MANIFEST.csv
    # -----------------------------------------------------------------------
    _update_manifest(OUTPUTS_DIR.parent, repo_sha)

    logger.info("=== Cycle 0 training complete. All outputs in %s ===", OUTPUTS_DIR)


def _update_manifest(outputs_dir: Path, repo_sha: str) -> None:
    """
    (Re-)write outputs/MANIFEST.csv with all generated files.
    写出 outputs/MANIFEST.csv，包含所有生成文件。
    """
    rows = []
    for path in sorted(outputs_dir.rglob("*")):
        if path.is_dir() or path.name == "MANIFEST.csv":
            continue
        sha = _sha256_file(path)
        size = path.stat().st_size
        # Estimate n_records based on file type
        # 根据文件类型估算 n_records
        n_records: float | str = ""
        if path.suffix == ".csv":
            try:
                n_records = len(pd.read_csv(path))
            except Exception:
                n_records = ""
        elif path.suffix == ".pt":
            n_records = 1  # one model state_dict
        elif path.suffix in (".json", ".png"):
            n_records = 1

        rows.append(
            {
                "filename":      str(path.relative_to(outputs_dir)),
                "sha256":        sha,
                "bytes":         size,
                "n_records":     n_records,
                "created_utc":   datetime.datetime.now(datetime.timezone.utc).isoformat().replace("+00:00", "Z"),
                "source_acc":    "",
                "source_module": "06_uncertainty_model",
                "notes":         f"cycle_0 synthetic training; repo={repo_sha}",
            }
        )
    manifest_df = pd.DataFrame(rows)
    manifest_path = outputs_dir / "MANIFEST.csv"
    manifest_df.to_csv(manifest_path, index=False)
    logger.info("Wrote MANIFEST.csv (%d files)", len(rows))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train Cycle 0 deep ensemble")
    parser.add_argument("--hidden_dim", type=int, default=256)
    parser.add_argument("--n_members", type=int, default=5)
    parser.add_argument("--max_epochs", type=int, default=200)
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--patience", type=int, default=10)
    parser.add_argument("--batch_size", type=int, default=32)
    parser.add_argument("--noise_std", type=float, default=0.5)
    parser.add_argument("--embed_dim", type=int, default=320)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    main(**vars(args))
