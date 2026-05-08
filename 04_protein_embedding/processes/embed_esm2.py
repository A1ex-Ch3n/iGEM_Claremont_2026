"""
embed_esm2.py — CLI wrapper for ESM-2 protein embedding.
ESM-2 蛋白质嵌入命令行脚本（Laguna GPU 节点使用）。

Usage / 用法:
    python embed_esm2.py \
        --input  <path_to.faa> \
        --output <path_to.npz> \
        --model  esm2_t33_650M_UR50D \
        --batch-size 4 \
        --pooling mean

Reference / 参考文献:
    Lin et al. (2023) Science 379(6637):1123-1130. DOI: 10.1126/science.ade2574
"""
import argparse
import logging
import random
import sys
from pathlib import Path

import numpy as np
import torch

# Ensure esm2_utils is importable / 确保 esm2_utils 可导入
sys.path.insert(0, str(Path(__file__).resolve().parent))
from esm2_utils import embed_sequences, load_fasta, load_model, save_npz

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Embed protein sequences with ESM-2.")
    p.add_argument(
        "--input", required=True, type=Path,
        help="Input FASTA (.faa) file / 输入 FASTA 文件"
    )
    p.add_argument(
        "--output", required=True, type=Path,
        help="Output .npz file path / 输出 .npz 文件路径"
    )
    p.add_argument(
        "--model", default="esm2_t6_8M_UR50D",
        choices=[
            "esm2_t6_8M_UR50D",
            "esm2_t12_35M_UR50D",
            "esm2_t30_150M_UR50D",
            "esm2_t33_650M_UR50D",
            "esm2_t36_3B_UR50D",
        ],
        help="ESM-2 model name / 模型名称 (default: 8M CPU-safe)"
    )
    p.add_argument(
        "--pooling", default="mean", choices=["mean", "per_residue"],
        help="Pooling strategy / 池化策略 (default: mean)"
    )
    p.add_argument(
        "--batch-size", type=int, default=8,
        help="Sequences per forward pass / 每次前向传播的序列数 (default: 8)"
    )
    p.add_argument(
        "--device", default=None,
        help="'cpu' | 'cuda' | None=auto / 运行设备（默认自动选择）"
    )
    p.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducibility / 随机种子（默认 42）"
    )
    return p.parse_args()


def set_seeds(seed: int) -> None:
    """Fix all random seeds for reproducibility / 固定所有随机种子"""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def main() -> None:
    args = parse_args()
    set_seeds(args.seed)

    logger.info("Input:  %s", args.input)
    logger.info("Output: %s", args.output)
    logger.info("Model:  %s | pooling=%s | batch_size=%d", args.model, args.pooling, args.batch_size)

    # Load sequences / 加载序列
    seqs = load_fasta(args.input)
    logger.info("Loaded %d sequences", len(seqs))

    # Load model / 加载模型
    model, alphabet = load_model(args.model, device=args.device)

    # Embed / 生成嵌入向量
    emb_dict = embed_sequences(
        model, alphabet, seqs,
        pooling=args.pooling,
        batch_size=args.batch_size,
        model_name=args.model,
    )
    logger.info("Embeddings shape: %s  dtype=%s", emb_dict["array"].shape, emb_dict["array"].dtype)

    # Save / 保存
    save_npz(args.output, emb_dict)
    logger.info("Done.")


if __name__ == "__main__":
    main()
