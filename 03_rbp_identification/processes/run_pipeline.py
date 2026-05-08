"""
Runner: generate RBP identification outputs for all available phages.
运行器：为所有可用噬菌体生成 RBP 鉴定输出。

Usage: python run_pipeline.py
"""

import logging
import sys
from pathlib import Path

import pandas as pd

# Path anchoring / 路径锚定
# notebooks live in processes/, so parents[1] = module root, parents[2] = repo root
PROCESSES_DIR = Path(__file__).resolve().parent
MODULE_DIR = PROCESSES_DIR.parent
REPO_ROOT = MODULE_DIR.parents[1]

sys.path.insert(0, str(PROCESSES_DIR))
from rbp_pipeline import run_phage, write_manifest, sha256_file

logging.basicConfig(level=logging.INFO, format='%(levelname)s  %(message)s')

HMM_DB = MODULE_DIR / 'inputs' / 'phagerbpdetect_data' / 'RBPdetect_phageRBPs.hmm'
OUTPUT_DIR = MODULE_DIR / 'outputs'

# Input: prefer 02_annotation output; fallback to 00_raw_data
# 输入：优先使用 02_annotation 输出；回退到 00_raw_data
ANNOTATION_DIR = REPO_ROOT / '02_annotation' / 'outputs' / 'phage_proteins'
FALLBACK_DIR = REPO_ROOT / '00_raw_data' / 'phage'

PHAGE_ACC = 'EU717894.1'


def resolve_faa(acc: str) -> Path:
    """
    Resolve phage protein FASTA path. Priority order:
    1. 02_annotation/outputs/phage_proteins/<acc>.faa (canonical)
    2. 00_raw_data/phage/<acc>/proteins.faa (raw data fallback)
    3. inputs/<acc>_fallback_proteins.faa (pre-copied worktree fallback)

    解析噬菌体蛋白质 FASTA 路径。优先级：
    1. 02_annotation 正式输出（最优）
    2. 00_raw_data 原始数据回退
    3. inputs/ 预复制的工作树回退
    """
    annotation_faa = ANNOTATION_DIR / f'{acc}.faa'
    if annotation_faa.exists():
        logging.info('Using 02_annotation output: %s', annotation_faa)
        return annotation_faa

    fallback_faa = FALLBACK_DIR / acc / 'proteins.faa'
    if fallback_faa.exists():
        logging.warning(
            'Module 02 output not found for %s; using 00_raw_data fallback %s. '
            'ORF IDs assigned in FASTA order (EU717894.1_orf_00001…). '
            'Headers are NCBI-format and differ from PHANOTATE output.',
            acc, fallback_faa,
        )
        return fallback_faa

    # Worktree-local fallback: pre-copied into inputs/ for isolated worktree builds
    # 工作树本地回退：已预复制到 inputs/ 以支持隔离工作树构建
    worktree_fallback = MODULE_DIR / 'inputs' / f'{acc}_fallback_proteins.faa'
    if worktree_fallback.exists():
        logging.warning(
            'Using worktree-local fallback %s (copied from 00_raw_data). '
            'ORF IDs assigned in FASTA order (EU717894.1_orf_00001…).',
            worktree_fallback,
        )
        return worktree_fallback

    raise FileNotFoundError(f'No protein FASTA found for {acc}')


def main():
    all_dfs = []

    faa = resolve_faa(PHAGE_ACC)
    candidates = run_phage(
        faa_path=faa,
        phage_acc=PHAGE_ACC,
        hmm_db=HMM_DB,
        output_dir=OUTPUT_DIR,
        top_k=5,
    )
    all_dfs.append(candidates)

    # Unified all_rbp_candidates.csv (only when >1 phage processed)
    # 多个噬菌体时才写 all_rbp_candidates.csv
    if len(all_dfs) > 1:
        all_df = pd.concat(all_dfs, ignore_index=True)
        all_path = OUTPUT_DIR / 'all_rbp_candidates.csv'
        all_df.to_csv(all_path, index=False, float_format='%.6f')
        logging.info('Wrote all_rbp_candidates.csv with %d rows', len(all_df))
    else:
        # Single phage: write per-phage CSV only; document in notes
        # 单个噬菌体：只写每噬菌体 CSV；文档中注明
        logging.info('Only one phage processed; skipping all_rbp_candidates.csv per AGENT_TODO anti-goals')

    # Print sanity summary / 打印健全性摘要
    passing = candidates[candidates['passes_threshold']]
    print(f'\n=== phiL7 ({PHAGE_ACC}) RBP Summary ===')
    print(f'Total ORFs: {len(candidates)}')
    print(f'HMM-positive candidates: {len(passing)}')
    for _, row in passing.head(5).iterrows():
        print(
            f'  rank={row["rank"]:2d}  {row["orf_id"]}  '
            f'{row["length_aa"]:4d} aa  '
            f'hmm_score={row["hmm_score"]:.1f}  '
            f'domain={row["hmm_match_pfam"]}  '
            f'evidence={row["evidence_track"]}'
        )


if __name__ == '__main__':
    main()
