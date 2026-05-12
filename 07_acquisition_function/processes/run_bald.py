"""
run_bald.py — Cycle N+1 acquisition: reads Module 06 predictions, outputs BALD recommendations.
Cycle N+1 采集步骤：读取 Module 06 预测结果，输出 BALD 推荐变体列表。

48-hour SLA pipeline:
  ELISA data (08_cycle_data) → retrain ensemble (06) → this script → next-cycle CSV.

Usage:
    python run_bald.py --cycle 1 [options]

Key options:
    --predictions_csv   Path to Module 06 predictions.csv (auto-detected if omitted).
    --measured_csv      Path to CSV of already-measured (rbp_id, receptor_id) pairs.
                        Omit for Cycle 0 (no measurements yet).
    --cycle             Target cycle number for output directory (e.g. 1 = first BALD run).
    --n_bald            Number of top-BALD picks per cycle (default 4).
    --n_random          Number of random-control picks per cycle (default 1).
    --seed              Random seed (default 42).

Input (from Module 06):
    06_uncertainty_model/outputs/cycle_<N>/predictions.csv
    Required columns: rbp_id, receptor_id, predicted_score, std, epistemic_std

Output:
    07_acquisition_function/outputs/cycle_<N>/
        recommendations.csv  — top-BALD + random-control picks
        random_replay.csv    — what pure-random would have picked
        safe_pick_backup.csv — extended BALD pool for manual PI selection
        run_meta.json        — run parameters and provenance
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import Optional, Set, Tuple

import numpy as np
import pandas as pd

# Anchor repo root / 锚定仓库根目录
SCRIPT_DIR = Path(__file__).resolve().parent
MODULE_DIR = SCRIPT_DIR.parent           # 07_acquisition_function/
REPO_ROOT = MODULE_DIR.parent

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from bald import bald_score, select_batch

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers / 工具函数
# ---------------------------------------------------------------------------

def _get_repo_commit() -> str:
    """Short git SHA for provenance. 获取短 git SHA 用于溯源。"""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, cwd=REPO_ROOT,
        )
        return result.stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def _sha256_file(path: Path) -> str:
    """SHA-256 of a file. 计算文件 SHA-256。"""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_predictions(predictions_csv: Optional[Path], cycle_in: int) -> pd.DataFrame:
    """
    Load predictions CSV from Module 06.
    加载 Module 06 的 predictions.csv。
    Auto-detects path if not provided.
    """
    if predictions_csv is not None:
        path = Path(predictions_csv)
    else:
        # Auto-detect: look for cycle_{cycle_in} first, fall back to cycle_0
        # 自动检测路径：优先使用 cycle_{cycle_in}，没有则退到 cycle_0
        path = REPO_ROOT / "06_uncertainty_model" / "outputs" / f"cycle_{cycle_in}" / "predictions.csv"
        if not path.exists():
            path = REPO_ROOT / "06_uncertainty_model" / "outputs" / "cycle_0" / "predictions.csv"

    if not path.exists():
        raise FileNotFoundError(
            f"predictions.csv not found at {path}. "
            "Run 06_uncertainty_model/processes/run_cycle0.py first."
        )

    df = pd.read_csv(path)

    # Validate required columns / 验证必须的列
    required = {"rbp_id", "receptor_id", "predicted_score", "std"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"predictions.csv is missing required columns: {missing}")

    # Warn if epistemic_std is absent (fallback to total std)
    # 若缺少 epistemic_std 则回退到 total std 并发出警告
    if "epistemic_std" not in df.columns:
        logger.warning(
            "epistemic_std column not found in predictions.csv — "
            "falling back to total std as BALD proxy. "
            "Re-run 06_uncertainty_model/processes/run_cycle0.py to add it."
        )
        df["epistemic_std"] = df["std"]

    logger.info("Loaded predictions.csv: %d pairs from %s", len(df), path)
    return df


def _load_measured_pairs(measured_csv: Optional[Path]) -> Set[Tuple[str, str]]:
    """
    Load already-measured (rbp_id, receptor_id) pairs to exclude.
    加载已测量的 pair 集合，用于排除（Cycle 0 时传 None）。
    """
    if measured_csv is None:
        logger.info("No measured_csv provided — Cycle 0 mode (no exclusions).")
        return set()

    path = Path(measured_csv)
    if not path.exists():
        raise FileNotFoundError(f"measured_csv not found: {path}")

    mdf = pd.read_csv(path)
    if "rbp_id" not in mdf.columns or "receptor_id" not in mdf.columns:
        raise ValueError("measured_csv must have 'rbp_id' and 'receptor_id' columns.")

    pairs = set(zip(mdf["rbp_id"], mdf["receptor_id"]))
    logger.info("Loaded %d measured pairs to exclude.", len(pairs))
    return pairs


# ---------------------------------------------------------------------------
# Main pipeline
# 主流程
# ---------------------------------------------------------------------------

def main(
    cycle: int = 1,
    predictions_csv: Optional[str] = None,
    measured_csv: Optional[str] = None,
    n_bald: int = 4,
    n_random: int = 1,
    seed: int = 42,
) -> None:
    """
    Run BALD acquisition for one cycle. Saves all outputs under
    07_acquisition_function/outputs/cycle_<cycle>/.
    运行一个 Cycle 的 BALD 采集。所有输出保存到
    07_acquisition_function/outputs/cycle_<cycle>/。
    """
    logger.info("=== Module 07 BALD Acquisition — Cycle %d ===", cycle)
    logger.info("n_bald=%d, n_random=%d, seed=%d", n_bald, n_random, seed)

    out_dir = MODULE_DIR / "outputs" / f"cycle_{cycle}"
    out_dir.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------------------------
    # Step 1: Load predictions / 步骤 1: 加载预测结果
    # -----------------------------------------------------------------------
    # Input cycle for predictions = cycle - 1 (we predict for next cycle)
    # 输入的预测对应的是上一个 cycle 训练的模型
    cycle_in = max(cycle - 1, 0)
    df = _load_predictions(
        Path(predictions_csv) if predictions_csv else None,
        cycle_in=cycle_in,
    )

    # -----------------------------------------------------------------------
    # Step 2: Compute BALD scores / 步骤 2: 计算 BALD 分数
    # -----------------------------------------------------------------------
    df["bald_score"] = bald_score(df["epistemic_std"].to_numpy())
    logger.info(
        "BALD score stats — min=%.4f  max=%.4f  mean=%.4f",
        df["bald_score"].min(), df["bald_score"].max(), df["bald_score"].mean(),
    )

    # -----------------------------------------------------------------------
    # Step 3: Load measured pairs (empty for Cycle 0)
    # 步骤 3: 加载已测量 pair（Cycle 0 时为空）
    # -----------------------------------------------------------------------
    measured_pairs = _load_measured_pairs(
        Path(measured_csv) if measured_csv else None
    )
    n_unmeasured = len(df) - sum(
        1 for _, row in df.iterrows()
        if (row["rbp_id"], row["receptor_id"]) in measured_pairs
    )
    logger.info(
        "Candidate pool: %d total pairs, %d unmeasured, %d measured (excluded).",
        len(df), n_unmeasured, len(measured_pairs),
    )

    # -----------------------------------------------------------------------
    # Step 4: Select batch / 步骤 4: 选批次
    # -----------------------------------------------------------------------
    recommendations, random_replay, safe_pick_backup = select_batch(
        df,
        bald_col="bald_score",
        n_bald=n_bald,
        n_random=n_random,
        measured_pairs=measured_pairs,
        seed=seed,
    )

    logger.info(
        "Selected: %d BALD picks + %d random control = %d total.",
        n_bald, n_random, len(recommendations),
    )

    # -----------------------------------------------------------------------
    # Step 5: Save outputs / 步骤 5: 保存输出文件
    # -----------------------------------------------------------------------
    rec_path = out_dir / "recommendations.csv"
    recommendations.to_csv(rec_path, index=False)
    logger.info("Saved recommendations.csv (%d rows)", len(recommendations))

    replay_path = out_dir / "random_replay.csv"
    random_replay.to_csv(replay_path, index=False)
    logger.info("Saved random_replay.csv (%d rows)", len(random_replay))

    backup_path = out_dir / "safe_pick_backup.csv"
    safe_pick_backup.to_csv(backup_path, index=False)
    logger.info("Saved safe_pick_backup.csv (%d rows)", len(safe_pick_backup))

    # -----------------------------------------------------------------------
    # Step 6: Save run_meta.json for provenance / 步骤 6: 保存运行元数据
    # -----------------------------------------------------------------------
    repo_sha = _get_repo_commit()
    meta = {
        "cycle":            cycle,
        "n_bald":           n_bald,
        "n_random":         n_random,
        "seed":             seed,
        "n_candidates":     len(df),
        "n_unmeasured":     n_unmeasured,
        "n_measured_excl":  len(measured_pairs),
        "bald_score_col":   "epistemic_std",
        "bald_score_max":   round(float(df["bald_score"].max()), 6),
        "bald_score_min":   round(float(df["bald_score"].min()), 6),
        "bald_score_mean":  round(float(df["bald_score"].mean()), 6),
        "top_bald_picks":   recommendations[recommendations["selection_reason"].str.startswith("bald_top")][
                                ["rbp_id", "receptor_id", "bald_score", "predicted_score"]
                            ].to_dict(orient="records"),
        "repo_commit_sha":  repo_sha,
        "timestamp_utc":    datetime.datetime.now(datetime.timezone.utc).isoformat().replace("+00:00", "Z"),
        "python_version":   sys.version.split()[0],
        "pandas_version":   pd.__version__,
        "numpy_version":    np.__version__,
    }
    meta_path = out_dir / "run_meta.json"
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    logger.info("Saved run_meta.json")

    # -----------------------------------------------------------------------
    # Step 7: Print summary for wet lab / 步骤 7: 打印湿实验团队摘要
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print(f"CYCLE {cycle} RECOMMENDATIONS (order = priority)")
    print("=" * 60)
    for _, row in recommendations.iterrows():
        reason = row["selection_reason"]
        print(
            f"  [{reason:20s}]  {row['rbp_id']:12s} × {row['receptor_id']:12s}"
            f"  BALD={row['bald_score']:.4f}  pred={row['predicted_score']:.3f}"
        )
    print("=" * 60)
    print(f"Safe-pick backup (top {len(safe_pick_backup)} BALD) → safe_pick_backup.csv")
    print(f"Random-replay baseline → random_replay.csv")
    print(f"All outputs in: {out_dir}")
    print()

    logger.info("=== Cycle %d BALD acquisition complete ===", cycle)


# ---------------------------------------------------------------------------
# CLI entry point / 命令行入口
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Module 07: BALD acquisition function for one cycle."
    )
    parser.add_argument(
        "--cycle", type=int, default=1,
        help="Output cycle number (e.g. 1 for first BALD run after Cycle 0 training).",
    )
    parser.add_argument(
        "--predictions_csv", type=str, default=None,
        help="Path to Module 06 predictions.csv (auto-detected if omitted).",
    )
    parser.add_argument(
        "--measured_csv", type=str, default=None,
        help="Path to CSV of already-measured (rbp_id, receptor_id) pairs. Omit for Cycle 0.",
    )
    parser.add_argument("--n_bald",   type=int, default=4,  help="Number of BALD picks.")
    parser.add_argument("--n_random", type=int, default=1,  help="Number of random-control picks.")
    parser.add_argument("--seed",     type=int, default=42, help="Random seed.")
    args = parser.parse_args()
    main(**vars(args))
