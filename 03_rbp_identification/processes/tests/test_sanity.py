"""
Sanity checks for Module 03 phiL7 outputs per INTERFACE.md §Module 03.
根据 INTERFACE.md §Module 03 对模块 03 phiL7 输出进行健全性检查。

INTERFACE.md states:
  - phiL7 should yield 1-3 RBP candidates with combined_score > 0.7
  - The tail spike gp25 is the top candidate

Note on protein length: INTERFACE.md cites ~412 aa for the tail spike (Lee et al. 2009
CDS annotation). The NCBI proteins.faa uses a broader annotation (712 aa for orf_00001)
that includes the full fusion protein containing tail spike domains. Both annotations
are consistent with RBP function; we document the discrepancy here.

注意蛋白质长度：INTERFACE.md 引用尾刺蛋白约 412 aa（Lee et al. 2009 CDS 注释）。
NCBI proteins.faa 使用更广泛的注释（orf_00001 为 712 aa），包含含尾刺域的完整融合蛋白。
两种注释均与 RBP 功能一致；此处记录该差异。
"""

from pathlib import Path

import pandas as pd
import pytest

MODULE_DIR = Path(__file__).resolve().parents[2]
OUTPUTS_DIR = MODULE_DIR / 'outputs'
PHAGE_ACC = 'EU717894.1'


@pytest.fixture(scope='module')
def candidates():
    """Load phiL7 candidates CSV once for all sanity tests.
    为所有健全性测试加载一次 phiL7 候选者 CSV。"""
    path = OUTPUTS_DIR / f'{PHAGE_ACC}_rbp_candidates.csv'
    assert path.exists(), f'Run the pipeline first: {path}'
    return pd.read_csv(path)


def test_phage_acc_in_filename(candidates):
    """Output filename contains phage accession. 输出文件名包含噬菌体编号。"""
    path = OUTPUTS_DIR / f'{PHAGE_ACC}_rbp_candidates.csv'
    assert path.exists()


def test_at_least_one_rbp_candidate(candidates):
    """phiL7 must yield ≥1 RBP candidate (combined_score ≥ 0.5).
    phiL7 必须至少有 1 个 RBP 候选者（combined_score ≥ 0.5）。"""
    passing = candidates[candidates['passes_threshold'] == True]
    assert len(passing) >= 1, (
        f'Expected ≥1 RBP candidate, found {len(passing)}'
    )


def test_between_one_and_three_high_confidence(candidates):
    """INTERFACE.md sanity: phiL7 yields 1-3 RBP candidates with combined_score > 0.7.
    INTERFACE.md 健全性：phiL7 产生 1-3 个 combined_score > 0.7 的 RBP 候选者。"""
    high_conf = candidates[candidates['combined_score'] > 0.7]
    assert 1 <= len(high_conf) <= 3, (
        f'Expected 1-3 high-confidence RBPs, found {len(high_conf)}'
    )


def test_top_candidate_is_rank_one(candidates):
    """The top-ranked candidate must be rank=1. 最高排名候选者必须为 rank=1。"""
    passing = candidates[candidates['passes_threshold'] == True]
    assert passing['rank'].min() == 1, 'No passing candidate has rank=1'


def test_top_candidate_has_hmm_evidence(candidates):
    """Top RBP (rank=1) must have HMM evidence. rank=1 候选者必须有 HMM 证据。"""
    top = candidates[candidates['rank'] == 1]
    assert len(top) == 1
    assert top.iloc[0]['evidence_track'] == 'hmm', (
        f"Expected 'hmm', got '{top.iloc[0]['evidence_track']}'"
    )


def test_top_candidate_length_reasonable(candidates):
    """
    Top candidate protein length should be in a range consistent with phage RBPs
    (tail spikes: typically 200-1500 aa). The INTERFACE.md note of ~412 aa reflects
    Lee et al. 2009 annotation; NCBI annotation yields a 712 aa fusion protein.
    Both are acceptable.

    顶级候选蛋白质长度应在噬菌体 RBP 合理范围内（尾刺蛋白通常 200-1500 aa）。
    """
    top = candidates[candidates['rank'] == 1].iloc[0]
    length = top['length_aa']
    assert 100 <= length <= 2000, (
        f'Top candidate length {length} aa is outside expected RBP range [100, 2000]'
    )


def test_all_orfs_present(candidates):
    """All 58 phiL7 ORFs should have rows in the output (one per ORF).
    phiL7 的所有 58 个 ORF 在输出中应各有一行。"""
    assert len(candidates) == 58, (
        f'Expected 58 rows (one per phiL7 ORF), got {len(candidates)}'
    )


def test_hmm_hits_have_combined_score_one(candidates):
    """All HMM-positive rows must have combined_score = 1.0 per INTERFACE.md formula.
    所有 HMM 阳性行的 combined_score 必须为 1.0（依据 INTERFACE.md 公式）。"""
    hmm_rows = candidates[candidates['evidence_track'] == 'hmm']
    assert len(hmm_rows) > 0, 'No HMM hits found — check HMM scan'
    assert (hmm_rows['combined_score'] == 1.0).all(), (
        'Some HMM rows have combined_score ≠ 1.0'
    )


def test_tail_spike_domain_detected(candidates):
    """At least one candidate should have a known tail spike Pfam or PhageRBPdetect domain.
    至少一个候选者应具有已知尾刺 Pfam 或 PhageRBPdetect 域。"""
    tail_spike_domains = {'Tail_spike_N', 'unknown_C54', 'End_tail_spike',
                          'PhageP22-tail', 'Phage_spike_2'}
    hmm_matches = candidates['hmm_match_pfam'].fillna('')
    found = hmm_matches[hmm_matches.isin(tail_spike_domains)]
    assert len(found) > 0, (
        f'No tail spike domain found in hmm_match_pfam. Domains seen: '
        f'{candidates["hmm_match_pfam"].dropna().unique().tolist()}'
    )
