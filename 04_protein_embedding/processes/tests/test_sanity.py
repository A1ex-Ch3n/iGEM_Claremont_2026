"""
test_sanity.py — module-specific sanity checks for embeddings.
模块级完整性检查：验证 phiL7 RBP 嵌入向量的合理性。
"""
import sys
import json
import numpy as np
import pytest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from esm2_utils import load_model, embed_sequences, load_fasta

MODULE_DIR = Path(__file__).resolve().parents[2]
INPUTS = MODULE_DIR / "inputs"
OUTPUTS = MODULE_DIR / "outputs"

MODEL_NAME = "esm2_t6_8M_UR50D"
EMBED_DIM = 320


@pytest.fixture(scope="module")
def model_and_alphabet():
    return load_model(MODEL_NAME)


def test_phiL7_rbp01_embedding_finite_and_normed(model_and_alphabet):
    """
    phiL7 RBP candidate 1 (~412 aa): embedding must be finite with L2 norm ∈ [0.1, 100].
    phiL7 RBP 候选 1（约 412 aa）的嵌入必须有限，且 L2 范数在 [0.1, 100] 之间。
    """
    model, alphabet = model_and_alphabet
    seqs = load_fasta(INPUTS / "mock_phiL7_rbp_sequences.faa")
    # rbp_01 is the first record (412 aa) / rbp_01 是第一条记录（412 aa）
    rbp_01 = [seqs[0]]
    result = embed_sequences(model, alphabet, rbp_01, pooling="mean")
    arr = result["array"]

    assert np.all(np.isfinite(arr)), "Embedding contains non-finite values (NaN/Inf)"
    l2_norm = float(np.linalg.norm(arr[0]))
    assert 0.1 <= l2_norm <= 100.0, \
        f"L2 norm {l2_norm:.4f} outside expected range [0.1, 100] — embedding may be degenerate"


def test_short_sequence_not_zero_padded(model_and_alphabet):
    """
    Short sequences (95 aa) must not produce near-zero mean from pad tokens bleeding in.
    短序列（95 aa）的均值不应因 padding token 拉低而趋近于零（正确 mask 验证）。
    """
    model, alphabet = model_and_alphabet
    seqs = load_fasta(INPUTS / "mock_phiL7_rbp_sequences.faa")
    # rbp_04 is 95 aa — shortest candidate / rbp_04 是最短候选（95 aa）
    short_seq = [s for s in seqs if "rbp_04" in s[0]]
    assert short_seq, "Could not find rbp_04 in mock FASTA"
    result = embed_sequences(model, alphabet, short_seq, pooling="mean")
    arr = result["array"]
    l2_norm = float(np.linalg.norm(arr[0]))
    # If padding leaks in, norm would collapse toward zero
    # 如果 padding 漏入均值，范数会向零坍缩
    assert l2_norm > 0.5, \
        f"Short-sequence L2 norm {l2_norm:.4f} suspiciously small — padding mask may be wrong"


def test_all_five_rbps_embedded_with_correct_dim(model_and_alphabet):
    """
    All 5 mock phiL7 RBPs embed to (5, 320) with distinct embeddings.
    5 条 mock RBP 序列应嵌入为 (5, 320)，且彼此不同。
    """
    model, alphabet = model_and_alphabet
    seqs = load_fasta(INPUTS / "mock_phiL7_rbp_sequences.faa")
    assert len(seqs) == 5, f"Expected 5 mock RBPs, got {len(seqs)}"
    result = embed_sequences(model, alphabet, seqs, pooling="mean", batch_size=4)
    arr = result["array"]
    assert arr.shape == (5, EMBED_DIM), f"Expected (5, {EMBED_DIM}), got {arr.shape}"
    # All rows should be distinct / 所有行应彼此不同
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            assert not np.allclose(arr[i], arr[j]), \
                f"Embeddings for seq {i} and {j} are identical (suspicious)"


def test_xcc_receptors_embedded(model_and_alphabet):
    """
    4 Xcc receptor sequences embed to (4, 320) with correct receptor IDs.
    4 条 Xcc 受体序列嵌入为 (4, 320)，且 seq_ids 正确。
    """
    model, alphabet = model_and_alphabet
    seqs = load_fasta(INPUTS / "xcc_receptors.faa")
    assert len(seqs) == 4, f"Expected 4 Xcc receptors, got {len(seqs)}"
    result = embed_sequences(model, alphabet, seqs, pooling="mean")
    arr = result["array"]
    assert arr.shape == (4, EMBED_DIM), f"Expected (4, {EMBED_DIM}), got {arr.shape}"
    # Check receptor IDs / 检查受体 ID
    ids = list(result["seq_ids"])
    for expected_id in ("GCF_000007145.1_tonB", "GCF_000007145.1_exbB"):
        assert expected_id in ids, f"Receptor ID '{expected_id}' missing from seq_ids"


def test_saved_npz_mean_matches_recomputed(model_and_alphabet):
    """
    Saved .npz embeddings match a freshly recomputed embedding.
    保存的 .npz 嵌入向量与重新计算的结果应匹配。
    """
    npz_path = OUTPUTS / "embeddings_esm2_t6_8M_phiL7_rbps.npz"
    if not npz_path.exists():
        pytest.skip("Output .npz not yet generated — run embedding first")
    data = np.load(npz_path, allow_pickle=True)
    saved_arr = data["array"]

    model, alphabet = model_and_alphabet
    seqs = load_fasta(INPUTS / "mock_phiL7_rbp_sequences.faa")
    fresh = embed_sequences(model, alphabet, seqs, pooling="mean")
    assert np.allclose(saved_arr, fresh["array"], atol=1e-5), \
        "Saved embeddings differ from freshly computed — model or input may have changed"
