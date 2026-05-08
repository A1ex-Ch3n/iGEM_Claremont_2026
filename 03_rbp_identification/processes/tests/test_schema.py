"""
Schema validation tests for Module 03 outputs.
模块 03 输出的模式验证测试。

Checks that all output files conform to INTERFACE.md §Module 03.
验证所有输出文件是否符合 INTERFACE.md §Module 03。
"""

import re
from pathlib import Path

import pandas as pd
import pytest

# Output directory relative to test location
# 相对于测试位置的输出目录
MODULE_DIR = Path(__file__).resolve().parents[2]
OUTPUTS_DIR = MODULE_DIR / 'outputs'

PHAGE_ACC = 'EU717894.1'


# ── helpers ────────────────────────────────────────────────────────────────

def _candidates_path() -> Path:
    return OUTPUTS_DIR / f'{PHAGE_ACC}_rbp_candidates.csv'


def _sequences_path() -> Path:
    return OUTPUTS_DIR / f'{PHAGE_ACC}_rbp_sequences.faa'


def _manifest_path() -> Path:
    return OUTPUTS_DIR / 'MANIFEST.csv'


# ── candidates CSV ─────────────────────────────────────────────────────────

class TestCandidatesCSV:
    """Validate EU717894.1_rbp_candidates.csv against INTERFACE.md schema.
    验证 EU717894.1_rbp_candidates.csv 是否符合 INTERFACE.md 模式。"""

    REQUIRED_COLUMNS = [
        'orf_id', 'phage_acc', 'length_aa',
        'hmm_score', 'hmm_match_pfam', 'ml_score',
        'combined_score', 'evidence_track', 'rank', 'passes_threshold',
    ]

    def test_file_exists(self):
        """Output CSV must exist. 输出 CSV 必须存在。"""
        assert _candidates_path().exists(), (
            f'Missing output: {_candidates_path()}'
        )

    def test_required_columns(self):
        """All INTERFACE.md columns must be present. 所有 INTERFACE.md 列必须存在。"""
        df = pd.read_csv(_candidates_path())
        missing = set(self.REQUIRED_COLUMNS) - set(df.columns)
        assert not missing, f'Missing columns: {missing}'

    def test_orf_id_format(self):
        """orf_id must match <acc>_orf_<5-digit>. orf_id 必须符合格式 <acc>_orf_<5位数字>。"""
        df = pd.read_csv(_candidates_path())
        pattern = re.compile(r'^[A-Z0-9_.]+_orf_\d{5}$')
        bad = df['orf_id'][~df['orf_id'].str.match(pattern)]
        assert bad.empty, f'Malformed orf_ids: {bad.tolist()[:5]}'

    def test_phage_acc_consistent(self):
        """phage_acc column must be uniform. phage_acc 列必须一致。"""
        df = pd.read_csv(_candidates_path())
        unique = df['phage_acc'].unique()
        assert len(unique) == 1, f'Multiple phage_acc values: {unique}'
        assert unique[0] == PHAGE_ACC

    def test_rank_is_contiguous(self):
        """rank column must be 1..N with no gaps. rank 列必须连续（1..N，无间隙）。"""
        df = pd.read_csv(_candidates_path())
        expected = list(range(1, len(df) + 1))
        actual = sorted(df['rank'].tolist())
        assert actual == expected, f'rank column is not contiguous 1..N'

    def test_combined_score_range(self):
        """combined_score must be in [0,1] or NaN. combined_score 必须在 [0,1] 或 NaN 范围内。"""
        df = pd.read_csv(_candidates_path())
        scores = df['combined_score'].dropna()
        assert ((scores >= 0.0) & (scores <= 1.0)).all(), (
            f'combined_score out of [0,1]: {scores[~((scores>=0)&(scores<=1))].tolist()}'
        )

    def test_hmm_hit_combined_score_is_one(self):
        """Per INTERFACE.md: if HMM hit then combined_score = 1.0.
        根据 INTERFACE.md：HMM 命中时 combined_score = 1.0。"""
        df = pd.read_csv(_candidates_path())
        hmm_rows = df[df['evidence_track'] == 'hmm']
        if not hmm_rows.empty:
            assert (hmm_rows['combined_score'] == 1.0).all(), (
                'HMM-positive rows must have combined_score = 1.0'
            )

    def test_passes_threshold_consistent(self):
        """passes_threshold = (combined_score >= 0.5). 通过阈值判断须与 combined_score 一致。"""
        df = pd.read_csv(_candidates_path())
        expected = (df['combined_score'] >= 0.5).fillna(False)
        actual = df['passes_threshold'].astype(bool)
        mismatches = (expected != actual).sum()
        assert mismatches == 0, f'{mismatches} rows have inconsistent passes_threshold'


# ── RBP sequences FASTA ────────────────────────────────────────────────────

class TestRBPSequencesFAA:
    """Validate EU717894.1_rbp_sequences.faa header format.
    验证 EU717894.1_rbp_sequences.faa 头部格式。"""

    # INTERFACE.md header: >EU717894.1_rbp_01 | source_orf=... | length=... | combined_score=... | evidence=...
    HEADER_PATTERN = re.compile(
        r'^>[A-Z0-9_.]+_rbp_\d{2} \| source_orf=[A-Z0-9_.]+_orf_\d{5} '
        r'\| length=\d+ \| combined_score=[\d.nan]+ \| evidence=\w+$'
    )

    def test_file_exists(self):
        """RBP FASTA must exist. RBP FASTA 文件必须存在。"""
        assert _sequences_path().exists(), f'Missing: {_sequences_path()}'

    def test_header_format(self):
        """Every FASTA header must match INTERFACE.md format.
        每条 FASTA 头部必须符合 INTERFACE.md 格式。"""
        bad_headers = []
        with open(_sequences_path()) as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith('>') and not self.HEADER_PATTERN.match(line):
                    bad_headers.append(line)
        assert not bad_headers, (
            f'Malformed FASTA headers:\n' + '\n'.join(bad_headers)
        )

    def test_rbp_ids_sequential(self):
        """RBP ids must be sequential: rbp_01, rbp_02, …
        RBP 编号必须连续：rbp_01, rbp_02, …"""
        ids = []
        with open(_sequences_path()) as fh:
            for line in fh:
                if line.startswith('>'):
                    ids.append(line.split()[0][1:])  # strip '>'
        for i, rbp_id in enumerate(ids, start=1):
            expected_suffix = f'rbp_{i:02d}'
            assert rbp_id.endswith(expected_suffix), (
                f'Expected {expected_suffix}, got {rbp_id}'
            )


# ── MANIFEST ───────────────────────────────────────────────────────────────

class TestManifest:
    """Validate MANIFEST.csv per INTERFACE.md universal conventions.
    验证 MANIFEST.csv 是否符合 INTERFACE.md 通用规范。"""

    REQUIRED_COLUMNS = [
        'filename', 'sha256', 'bytes', 'n_records',
        'created_utc', 'source_acc', 'source_module', 'notes',
    ]

    def test_manifest_exists(self):
        assert _manifest_path().exists(), f'Missing: {_manifest_path()}'

    def test_required_columns(self):
        df = pd.read_csv(_manifest_path())
        missing = set(self.REQUIRED_COLUMNS) - set(df.columns)
        assert not missing, f'Missing MANIFEST columns: {missing}'

    def test_source_module_correct(self):
        df = pd.read_csv(_manifest_path())
        assert (df['source_module'] == '03_rbp_identification').all()
