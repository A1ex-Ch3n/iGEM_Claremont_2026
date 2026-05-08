"""
test_determinism.py — ESM-2 in eval mode must be deterministic.
ESM-2 在 eval 模式下必须保证确定性（相同输入→相同输出）。
"""
import sys
import numpy as np
import pytest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from esm2_utils import load_model, embed_sequences

MODEL_NAME = "esm2_t6_8M_UR50D"


@pytest.fixture(scope="module")
def model_and_alphabet():
    return load_model(MODEL_NAME)


def test_deterministic_mean_pooled(model_and_alphabet):
    """
    Same sequences embedded twice → identical float32 arrays.
    相同序列嵌入两次 → 数组完全相同（eval 模式禁用 dropout）。
    """
    model, alphabet = model_and_alphabet
    seqs = [
        ("det_seq_A", "MKTAYIAKQRQISFVKSHFSRQLEE"),
        ("det_seq_B", "ACDEFGHIKLMNPQRSTVWYACDEF"),
        ("det_seq_C", "WVTSRQPNMLKIHGFEDCAWVTSRQ"),
    ]
    run1 = embed_sequences(model, alphabet, seqs, pooling="mean")
    run2 = embed_sequences(model, alphabet, seqs, pooling="mean")

    assert np.allclose(run1["array"], run2["array"], atol=1e-6), \
        "Embeddings differ between two identical runs — model not in eval() mode?"


def test_different_sequences_different_embeddings(model_and_alphabet):
    """
    Different sequences should produce different embeddings.
    不同序列应生成不同的嵌入向量（基础可区分性检查）。
    """
    model, alphabet = model_and_alphabet
    seq_a = [("seqA", "ACDEFGHIKLMNPQRSTVWY" * 2)]
    seq_b = [("seqB", "WVTSRQPNMLKIHGFEDCAW" * 2)]
    emb_a = embed_sequences(model, alphabet, seq_a)["array"]
    emb_b = embed_sequences(model, alphabet, seq_b)["array"]
    assert not np.allclose(emb_a, emb_b), "Different sequences produced identical embeddings"
