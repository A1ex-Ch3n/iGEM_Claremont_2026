"""
test_smoke.py — end-to-end smoke: embed 2 random sequences with ESM-2 8M.
端到端冒烟测试：用 ESM-2 8M 对 2 条随机序列生成嵌入向量。
"""
import sys
import random
import numpy as np
import pytest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from esm2_utils import load_model, embed_sequences

MODEL_NAME = "esm2_t6_8M_UR50D"
EMBED_DIM = 320  # ESM-2 8M embedding dimension


@pytest.fixture(scope="module")
def model_and_alphabet():
    """Load once per test module / 每个测试模块仅加载一次模型"""
    return load_model(MODEL_NAME)


def _random_seqs(n: int, length: int, seed: int = 42):
    rng = random.Random(seed)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    return [(f"smoke_seq_{i:02d}", "".join(rng.choices(aa, k=length))) for i in range(n)]


def test_embed_two_sequences_shape(model_and_alphabet):
    """Output array shape is (2, 320) for two 50-aa sequences."""
    model, alphabet = model_and_alphabet
    seqs = _random_seqs(n=2, length=50)
    result = embed_sequences(model, alphabet, seqs, pooling="mean", batch_size=2)
    arr = result["array"]
    assert arr.shape == (2, EMBED_DIM), f"Expected (2, {EMBED_DIM}), got {arr.shape}"


def test_embed_two_sequences_no_nan(model_and_alphabet):
    """Embeddings contain no NaN / 嵌入向量不含 NaN"""
    model, alphabet = model_and_alphabet
    seqs = _random_seqs(n=2, length=50)
    result = embed_sequences(model, alphabet, seqs, pooling="mean", batch_size=2)
    assert not np.any(np.isnan(result["array"])), "NaN values in embedding"


def test_embed_two_sequences_dtype(model_and_alphabet):
    """Embedding array dtype is float32 / 嵌入数组数据类型为 float32"""
    model, alphabet = model_and_alphabet
    seqs = _random_seqs(n=2, length=50)
    result = embed_sequences(model, alphabet, seqs, pooling="mean", batch_size=2)
    assert result["array"].dtype == np.float32


def test_lengths_match_sequences(model_and_alphabet):
    """lengths[i] == len(seq_i), not len(seq_i) + 2 (no BOS/EOS counted)."""
    model, alphabet = model_and_alphabet
    seqs = _random_seqs(n=3, length=50)
    result = embed_sequences(model, alphabet, seqs, pooling="mean", batch_size=8)
    for i, (_, seq) in enumerate(seqs):
        assert result["lengths"][i] == len(seq), \
            f"lengths[{i}] = {result['lengths'][i]}, expected {len(seq)}"


def test_result_has_all_keys(model_and_alphabet):
    """embed_sequences returns dict with seq_ids, array, lengths, meta."""
    model, alphabet = model_and_alphabet
    seqs = _random_seqs(n=1, length=30)
    result = embed_sequences(model, alphabet, seqs)
    for key in ("seq_ids", "array", "lengths", "meta"):
        assert key in result, f"Missing key: {key}"


def test_meta_is_valid_json(model_and_alphabet):
    """meta field parses as JSON with required keys / meta 字段可解析为 JSON"""
    import json
    model, alphabet = model_and_alphabet
    seqs = _random_seqs(n=1, length=30)
    result = embed_sequences(model, alphabet, seqs)
    meta = json.loads(result["meta"])
    for key in ("model", "pooling", "date_utc", "repo_commit_sha"):
        assert key in meta, f"meta missing key: {key}"
