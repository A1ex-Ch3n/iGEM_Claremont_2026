"""
Smoke tests: run the pipeline on a synthetic 3-protein FASTA.
冒烟测试：在合成的三蛋白质 FASTA 上运行流程。

Tests that the pipeline completes without crashing.
Do NOT assert on scores — random sequences produce arbitrary scores.
测试流程能无错完成。不对评分做断言——随机序列会产生任意评分。
"""

import random
import tempfile
from pathlib import Path

import pandas as pd
import pytest

PROCESSES_DIR = Path(__file__).resolve().parents[1]
MODULE_DIR = PROCESSES_DIR.parent

import sys
sys.path.insert(0, str(PROCESSES_DIR))

from rbp_pipeline import run_hmmscan, classify_hmm_hits, build_candidates_df, run_phage

HMM_DB = MODULE_DIR / 'inputs' / 'phagerbpdetect_data' / 'RBPdetect_phageRBPs.hmm'
AA_LETTERS = 'ACDEFGHIKLMNPQRSTVWY'


def _synthetic_faa(n_proteins: int = 3, seed: int = 42) -> Path:
    """
    Create a temporary FASTA with n random protein sequences (~100 aa each).
    创建包含 n 条随机蛋白质序列（每条约 100 aa）的临时 FASTA 文件。
    """
    rng = random.Random(seed)
    tmp = tempfile.NamedTemporaryFile(
        suffix='.faa', mode='w', delete=False, prefix='smoke_test_'
    )
    for i in range(1, n_proteins + 1):
        length = rng.randint(80, 120)
        seq = ''.join(rng.choice(AA_LETTERS) for _ in range(length))
        tmp.write(f'>synthetic_phage_orf_{i:05d} | synthetic protein {i}\n')
        for j in range(0, len(seq), 60):
            tmp.write(seq[j:j+60] + '\n')
    tmp.close()
    return Path(tmp.name)


class TestSmoke:
    """End-to-end smoke tests on synthetic data. 合成数据端到端冒烟测试。"""

    def test_hmmscan_runs(self, tmp_path):
        """hmmscan completes on synthetic FASTA without error.
        hmmscan 在合成 FASTA 上无误完成。"""
        faa = _synthetic_faa()
        try:
            scan_df = run_hmmscan(faa, HMM_DB)
            # Must be a DataFrame (possibly empty for random sequences)
            # 必须为 DataFrame（随机序列可能为空）
            assert isinstance(scan_df, pd.DataFrame)
        finally:
            faa.unlink(missing_ok=True)

    def test_classify_hmm_hits_empty(self):
        """classify_hmm_hits handles empty input gracefully.
        classify_hmm_hits 对空输入能优雅处理。"""
        empty = pd.DataFrame(columns=['seq_id', 'domain_id', 'domain_pfam_acc',
                                      'bitscore', 'bias', 'evalue'])
        result = classify_hmm_hits(empty)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0

    def test_build_candidates_df_smoke(self, tmp_path):
        """build_candidates_df produces correct columns for synthetic input.
        build_candidates_df 对合成输入生成正确的列。"""
        faa = _synthetic_faa(n_proteins=3)
        try:
            no_hits = pd.DataFrame(columns=['seq_id', 'n_blocks', 'c_blocks',
                                            'best_domain', 'best_bitscore', 'best_pfam_acc'])
            df = build_candidates_df(faa, 'TEST_ACC', no_hits)
            assert len(df) == 3, f'Expected 3 rows, got {len(df)}'
            required = [
                'orf_id', 'phage_acc', 'length_aa',
                'hmm_score', 'hmm_match_pfam', 'ml_score',
                'combined_score', 'evidence_track', 'rank', 'passes_threshold',
            ]
            missing = set(required) - set(df.columns)
            assert not missing, f'Missing columns: {missing}'
        finally:
            faa.unlink(missing_ok=True)

    def test_full_pipeline_smoke(self, tmp_path):
        """Full run_phage call completes without exception.
        完整 run_phage 调用无异常完成。"""
        faa = _synthetic_faa(n_proteins=3)
        try:
            df = run_phage(
                faa_path=faa,
                phage_acc='SMOKE_TEST',
                hmm_db=HMM_DB,
                output_dir=tmp_path,
                top_k=5,
            )
            assert isinstance(df, pd.DataFrame)
            assert len(df) == 3
        finally:
            faa.unlink(missing_ok=True)
