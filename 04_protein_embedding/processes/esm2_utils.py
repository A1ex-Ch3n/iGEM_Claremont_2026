"""
esm2_utils.py — ESM-2 embedding utilities for Module 04.
ESM-2 蛋白质嵌入工具函数（Module 04 专用）。

Reference / 参考文献:
    Lin et al. (2023) Science 379(6637):1123-1130. DOI: 10.1126/science.ade2574
"""
import json
import subprocess
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import ssl
import esm
import numpy as np
import torch
from Bio import SeqIO

# Some Python/macOS setups fail SSL cert verification when downloading weights.
# Patch globally so fair-esm's torch.hub download can proceed.
# 部分 Python/macOS 环境下载模型权重时 SSL 证书验证失败，此处全局修复。
ssl._create_default_https_context = ssl._create_unverified_context  # noqa: S323

# REPO_ROOT: two levels up from processes/ / 向上两级到 repo 根目录
REPO_ROOT = Path(__file__).resolve().parents[2]

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# ── git utilities ────────────────────────────────────────────────────────────

def _get_commit_sha() -> str:
    """Return short git HEAD SHA for provenance metadata / 返回短 SHA 用于出处记录"""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, cwd=REPO_ROOT
        )
        return result.stdout.strip() if result.returncode == 0 else "unknown"
    except Exception:
        return "unknown"


# ── model loading ────────────────────────────────────────────────────────────

def load_model(
    model_name: str,
    device: Optional[str] = None,
) -> Tuple[Any, Any]:
    """
    Load an ESM-2 pretrained model and its alphabet.
    加载 ESM-2 预训练模型及其字母表。

    Args:
        model_name: e.g. 'esm2_t6_8M_UR50D', 'esm2_t33_650M_UR50D'
        device:     'cpu' | 'cuda' | None (auto-detect / 自动选择)

    Returns:
        (model, alphabet) — model is in eval() mode with grad disabled
    """
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"

    logger.info("Loading %s on %s …", model_name, device)
    loader = getattr(esm.pretrained, model_name)
    model, alphabet = loader()
    model = model.to(device)
    model.eval()  # Disable dropout → deterministic output / 关闭 dropout 保证确定性
    logger.info("Loaded %s  embed_dim=%d", model_name, model.embed_dim)
    return model, alphabet


# ── sequence embedding ───────────────────────────────────────────────────────

def embed_sequences(
    model: Any,
    alphabet: Any,
    seqs: List[Tuple[str, str]],
    pooling: str = "mean",
    batch_size: int = 8,
    model_name: str = "",
) -> Dict[str, Any]:
    """
    Generate embeddings for a list of protein sequences using ESM-2.
    使用 ESM-2 为蛋白质序列列表生成嵌入向量。

    Token layout for ESM-2:
      position 0   → BOS (<cls>) token  — excluded from pooling
      position 1…L → amino-acid tokens  — pooled over (L = sequence length)
      position L+1 → EOS (<eos>) token  — excluded from pooling
      position L+2… → PAD tokens       — excluded by length mask

    Args:
        model:      ESM-2 model (must be in eval mode / 必须在 eval 模式下)
        alphabet:   ESM-2 alphabet
        seqs:       list of (seq_id, sequence) tuples
        pooling:    'mean' (default) or 'per_residue'
        batch_size: sequences per forward pass
        model_name: model identifier for metadata (optional)

    Returns:
        dict with keys:
            seq_ids  — np.ndarray of str, shape (N,)
            array    — np.float32, shape (N, D) for mean-pooled;
                       (N, max_len, D) for per_residue
            lengths  — np.int32, shape (N,) — residue counts (no BOS/EOS)
            meta     — JSON string per INTERFACE §Module 04
    """
    device = next(model.parameters()).device
    batch_converter = alphabet.get_batch_converter()

    all_ids: List[str] = []
    all_embeddings: List[np.ndarray] = []
    all_lengths: List[int] = []

    for batch_start in range(0, len(seqs), batch_size):
        batch = seqs[batch_start : batch_start + batch_size]
        _, batch_strs, batch_tokens = batch_converter(batch)
        batch_tokens = batch_tokens.to(device)

        with torch.no_grad():
            results = model(
                batch_tokens,
                repr_layers=[model.num_layers],
                return_contacts=False,
            )

        # token_reps: (batch, max_token_len, embed_dim)
        # token_len = max_seq_len_in_batch + 2  (BOS + residues + EOS + padding)
        token_reps = results["representations"][model.num_layers]

        for j, (seq_id, seq_str) in enumerate(batch):
            L = len(seq_str)  # residue count (no BOS/EOS) / 氨基酸长度（不含 BOS/EOS）
            all_ids.append(seq_id)
            all_lengths.append(L)

            if pooling == "mean":
                # Slice positions 1..L (residues only, drop BOS at 0 and EOS at L+1)
                # 只取氨基酸位置（1 到 L），排除 BOS（位置 0）和 EOS（位置 L+1）
                aa_repr = token_reps[j, 1 : L + 1, :]  # (L, D)
                pooled = aa_repr.mean(dim=0)            # (D,)
                all_embeddings.append(pooled.cpu().float().numpy())
            else:  # per_residue
                aa_repr = token_reps[j, 1 : L + 1, :]  # (L, D)
                all_embeddings.append(aa_repr.cpu().float().numpy())

        logger.debug(
            "Embedded batch %d-%d / %d",
            batch_start, batch_start + len(batch) - 1, len(seqs)
        )

    # Stack embeddings / 拼接嵌入向量
    if pooling == "mean":
        array = np.stack(all_embeddings, axis=0).astype(np.float32)  # (N, D)
    else:
        # Zero-pad to max_len for per_residue mode / 零填充至最大长度
        max_len = max(e.shape[0] for e in all_embeddings)
        embed_dim = all_embeddings[0].shape[1]
        array = np.zeros(
            (len(all_embeddings), max_len, embed_dim), dtype=np.float32
        )
        for k, emb in enumerate(all_embeddings):
            array[k, : emb.shape[0], :] = emb

    # Build metadata JSON / 构建元数据 JSON
    _name = model_name or getattr(model, "_model_name", type(model).__name__)
    meta = json.dumps(
        {
            "model": _name,
            "model_revision": "fair-esm-2.0.0",
            "pooling": pooling,
            "date_utc": datetime.now(timezone.utc).isoformat(),
            "repo_commit_sha": _get_commit_sha(),
        }
    )

    return {
        "seq_ids": np.array(all_ids, dtype=object),
        "array": array,
        "lengths": np.array(all_lengths, dtype=np.int32),
        "meta": meta,
    }


# ── I/O helpers ──────────────────────────────────────────────────────────────

def save_npz(out_path: Path, embedding_dict: Dict[str, Any]) -> None:
    """
    Save embedding dict to .npz per INTERFACE §Module 04 §NumPy conventions.
    按 INTERFACE §Module 04 §NumPy 规范将嵌入字典保存为 .npz 文件。

    Keys: seq_ids (object), array (float32), lengths (int32), meta (0-d str)
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        out_path,
        seq_ids=embedding_dict["seq_ids"],
        array=embedding_dict["array"].astype(np.float32),
        lengths=embedding_dict["lengths"].astype(np.int32),
        # Store JSON string as 0-d object array so np.load gives it back as str
        # 将 JSON 字符串存为 0 维对象数组，np.load 时可直接还原为字符串
        meta=np.array(embedding_dict["meta"]),
    )
    logger.info("Saved embeddings to %s", out_path)


def load_fasta(fasta_path: Path) -> List[Tuple[str, str]]:
    """
    Parse a FASTA file into (seq_id, sequence) pairs.
    解析 FASTA 文件，返回 (seq_id, sequence) 列表。

    Strips stop codons (*) and gaps (-); uppercases sequence.
    去除终止密码子（*）和比对间隔（-），序列转大写。
    """
    seqs = []
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        seq = str(record.seq).upper().replace("*", "").replace("-", "")
        seqs.append((record.id, seq))
    return seqs
