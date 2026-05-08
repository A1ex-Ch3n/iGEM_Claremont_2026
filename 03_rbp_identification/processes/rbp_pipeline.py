"""
RBP identification pipeline for iGEM Claremont 2026.
RBP（受体结合蛋白）鉴定流程 — 2026年 iGEM Claremont 项目

Implements HMM-based detection using PhageRBPdetect (Boeckaerts 2022) curated HMMs.
ML track (XGBoost-on-ESM) is documented as blocked: v2 model requires 1024-dim ProtTrans
embeddings (bio_embeddings install fails on Python 3.13+); v4 requires Zenodo model download.
ml_score is NaN throughout; combined_score = 1.0 for HMM hits per INTERFACE.md.

Reference: Boeckaerts et al. (2022) Viruses 14(6):1329. DOI: 10.3390/v14061329
"""

import os
import math
import subprocess
import hashlib
import logging
import tempfile
from pathlib import Path
from datetime import datetime, timezone
from typing import Optional

import pandas as pd
import numpy as np
from Bio import SeqIO

logger = logging.getLogger(__name__)

# ────────────────────────────────────────────────────────────────────────────
# HMM-track constants (from PhageRBPdetect v2 standalone.py)
# HMM 轨道常量（来自 PhageRBPdetect v2 standalone.py）
# ────────────────────────────────────────────────────────────────────────────

# N-terminal structural domains — presence marks a candidate
# N 端结构域 — 存在即标记为候选者
N_BLOCKS = {
    'Phage_T7_tail', 'Tail_spike_N', 'Prophage_tail', 'BppU_N', 'Mtd_N',
    'Head_binding', 'DUF3751', 'End_N_terminal', 'phage_tail_N', 'Prophage_tailD1',
    'DUF2163', 'Phage_fiber_2',
    'unknown_N0', 'unknown_N1', 'unknown_N2', 'unknown_N3', 'unknown_N4',
    'unknown_N6', 'unknown_N10', 'unknown_N11', 'unknown_N12', 'unknown_N13',
    'unknown_N17', 'unknown_N19', 'unknown_N23', 'unknown_N24', 'unknown_N26',
    'unknown_N29', 'unknown_N36', 'unknown_N45', 'unknown_N48', 'unknown_N49',
    'unknown_N53', 'unknown_N57', 'unknown_N60', 'unknown_N61', 'unknown_N65',
    'unknown_N73', 'unknown_N82', 'unknown_N83', 'unknown_N101', 'unknown_N114',
    'unknown_N119', 'unknown_N122', 'unknown_N163', 'unknown_N174', 'unknown_N192',
    'unknown_N200', 'unknown_N206', 'unknown_N208',
}

# C-terminal enzymatic/binding domains — stricter bitscore threshold (25 bits)
# C 端酶/结合域 — 更严格的 bitscore 阈值（25 bits）
C_BLOCKS = {
    'Lipase_GDSL_2', 'Pectate_lyase_3', 'gp37_C', 'Beta_helix', 'Gp58',
    'End_beta_propel', 'End_tail_spike', 'End_beta_barrel', 'PhageP22-tail',
    'Phage_spike_2', 'gp12-short_mid', 'Collar',
    'unknown_C2', 'unknown_C3', 'unknown_C8', 'unknown_C15', 'unknown_C35',
    'unknown_C54', 'unknown_C76', 'unknown_C100', 'unknown_C105', 'unknown_C112',
    'unknown_C123', 'unknown_C179', 'unknown_C201', 'unknown_C203', 'unknown_C228',
    'unknown_C234', 'unknown_C242', 'unknown_C258', 'unknown_C262', 'unknown_C267',
    'unknown_C268', 'unknown_C274', 'unknown_C286', 'unknown_C292', 'unknown_C294',
    'Peptidase_S74', 'Phage_fiber_C', 'S_tail_recep_bd', 'CBM_4_9',
    'DUF1983', 'DUF3672',
}

# Minimum bitscore thresholds (from v2 source; checked with OM_score > OM_bias filter)
# 最低 bitscore 阈值（来自 v2 源码；通过 OM_score > OM_bias 过滤器验证）
BITSCORE_MIN_ANY = 18       # global minimum for any domain to be considered
BITSCORE_MIN_C = 25         # C-block domains must also clear this bar
BITSCORE_MIN_N = 18         # N-block domains need only global minimum


def _hmmer_bin() -> Path:
    """Find hmmscan binary. 查找 hmmscan 二进制文件。"""
    candidates = [
        Path('/opt/homebrew/opt/hmmer/bin/hmmscan'),
        Path('/usr/local/bin/hmmscan'),
        Path('/usr/bin/hmmscan'),
    ]
    for c in candidates:
        if c.exists():
            return c
    # fallback: hope it's on PATH
    return Path('hmmscan')


def press_hmm_db(hmm_path: Path) -> None:
    """
    Press HMM database if not already pressed. Idempotent.
    如果尚未 press，则 press HMM 数据库（幂等操作）。
    """
    index_file = hmm_path.with_suffix(hmm_path.suffix + '.h3i')
    if index_file.exists():
        logger.info('HMM database already pressed: %s', hmm_path)
        return
    hmmscan = _hmmer_bin().parent / 'hmmpress'
    result = subprocess.run(
        [str(hmmscan), str(hmm_path)],
        capture_output=True, text=True, check=True,
    )
    logger.info('hmmpress output: %s', result.stdout.strip())


def run_hmmscan(faa_path: Path, hmm_db: Path) -> pd.DataFrame:
    """
    Run hmmscan and return a tidy DataFrame of per-sequence best hits.
    运行 hmmscan，返回每条序列最佳命中的整洁 DataFrame。

    Returns columns: seq_id, domain_id, domain_pfam_acc, bitscore, bias, evalue
    """
    hmmscan_bin = _hmmer_bin()
    with tempfile.NamedTemporaryFile(suffix='.tbl', delete=False) as tbl_f:
        tbl_path = tbl_f.name

    try:
        subprocess.run(
            [str(hmmscan_bin), '--tblout', tbl_path, str(hmm_db), str(faa_path)],
            capture_output=True, text=True, check=True,
        )
        rows = []
        with open(tbl_path) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) < 13:
                    continue
                # tblout format: target_name acc query_name acc E-value score bias ...
                domain_id = parts[0]
                domain_pfam_acc = parts[1] if parts[1] != '-' else ''
                seq_id = parts[2]
                evalue = float(parts[4])
                bitscore = float(parts[5])
                bias = float(parts[6])
                rows.append({
                    'seq_id': seq_id,
                    'domain_id': domain_id,
                    'domain_pfam_acc': domain_pfam_acc,
                    'bitscore': bitscore,
                    'bias': bias,
                    'evalue': evalue,
                })
        return pd.DataFrame(rows)
    finally:
        os.unlink(tbl_path)


def _passes_om_filter(bitscore: float, bias: float) -> bool:
    """
    Order-of-magnitude filter from PhageRBPdetect v2: log(score) > log(bias).
    使用自然对数（ln），与 PhageRBPdetect v2 源码一致。
    """
    om_score = math.floor(math.log(max(bitscore, 1e-9)))
    om_bias = math.floor(math.log(max(bias, 1e-9)))
    return om_score > om_bias


def classify_hmm_hits(scan_df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply PhageRBPdetect v2 domain classification rules.
    应用 PhageRBPdetect v2 域分类规则。

    Returns one row per sequence that passes HMM filter:
    columns: seq_id, n_blocks, c_blocks, best_domain, best_bitscore, best_pfam_acc
    """
    if scan_df.empty:
        return pd.DataFrame(columns=['seq_id', 'n_blocks', 'c_blocks',
                                     'best_domain', 'best_bitscore', 'best_pfam_acc'])

    result_rows = []
    for seq_id, group in scan_df.groupby('seq_id'):
        n_hits, c_hits = [], []
        best_score, best_domain, best_pfam = 0.0, '', ''
        for _, row in group.iterrows():
            dom = row['domain_id']
            score = row['bitscore']
            bias = row['bias']
            pfam = row['domain_pfam_acc']

            if score < BITSCORE_MIN_ANY:
                continue
            if not _passes_om_filter(score, bias):
                continue

            if dom in N_BLOCKS and score >= BITSCORE_MIN_N:
                n_hits.append(dom)
            elif dom in C_BLOCKS and score >= BITSCORE_MIN_C:
                c_hits.append(dom)

            if score > best_score:
                best_score, best_domain, best_pfam = score, dom, pfam

        if n_hits or c_hits:
            result_rows.append({
                'seq_id': seq_id,
                'n_blocks': n_hits,
                'c_blocks': c_hits,
                'best_domain': best_domain,
                'best_bitscore': best_score,
                'best_pfam_acc': best_pfam,
            })

    return pd.DataFrame(result_rows)


def build_candidates_df(
    faa_path: Path,
    phage_acc: str,
    hmm_hits: pd.DataFrame,
) -> pd.DataFrame:
    """
    Build the output DataFrame matching INTERFACE.md §Module 03 schema.
    构建符合 INTERFACE.md §Module 03 模式的输出 DataFrame。

    NOTE: ml_score is NaN throughout because bio_embeddings (required by v2 XGBoost)
    fails to install on Python 3.13+ (gensim build error) and the v4 Zenodo model
    was not downloaded in this sprint. combined_score = 1.0 for HMM hits.
    注意：整个 ml_score 为 NaN，原因是 bio_embeddings（v2 XGBoost 所需）在 Python 3.13+
    上安装失败，且 v4 Zenodo 模型本次未下载。HMM 命中时 combined_score = 1.0。
    """
    records = list(SeqIO.parse(faa_path, 'fasta'))
    hmm_hit_ids = set(hmm_hits['seq_id']) if not hmm_hits.empty else set()

    rows = []
    for idx, rec in enumerate(records, start=1):
        orf_id = f'{phage_acc}_orf_{idx:05d}'
        seq_id = rec.id
        length_aa = len(rec.seq)

        if seq_id in hmm_hit_ids:
            hit_row = hmm_hits[hmm_hits['seq_id'] == seq_id].iloc[0]
            hmm_score = hit_row['best_bitscore']
            hmm_match = hit_row['best_domain']
            ml_score = float('nan')
            combined_score = 1.0
            evidence_track = 'hmm'
        else:
            hmm_score = float('nan')
            hmm_match = ''
            ml_score = float('nan')
            combined_score = float('nan')
            evidence_track = ''

        passes = bool(not math.isnan(combined_score) and combined_score >= 0.5)

        rows.append({
            'orf_id': orf_id,
            'phage_acc': phage_acc,
            'length_aa': length_aa,
            'hmm_score': hmm_score,
            'hmm_match_pfam': hmm_match,
            'ml_score': ml_score,
            'combined_score': combined_score,
            'evidence_track': evidence_track,
            'rank': 0,          # filled below
            'passes_threshold': passes,
            '_seq_id': seq_id,  # internal; dropped before save
        })

    df = pd.DataFrame(rows)

    # Rank: HMM-hits first (by hmm_score desc), then non-hits
    # 排名：HMM 命中者优先（按 hmm_score 降序），其余排后
    hmm_mask = df['passes_threshold']
    hmm_part = df[hmm_mask].sort_values('hmm_score', ascending=False)
    non_hmm_part = df[~hmm_mask]
    ranked = pd.concat([hmm_part, non_hmm_part], ignore_index=True)
    ranked['rank'] = range(1, len(ranked) + 1)

    # Drop internal helper column
    ranked = ranked.drop(columns=['_seq_id'])
    return ranked


def write_candidates_csv(df: pd.DataFrame, out_path: Path) -> None:
    """
    Write candidates CSV per INTERFACE.md conventions.
    按 INTERFACE.md 规范写入候选者 CSV。
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False, float_format='%.6f')
    logger.info('Wrote %d candidate rows to %s', len(df), out_path)


def write_rbp_sequences_faa(
    faa_path: Path,
    candidates_df: pd.DataFrame,
    phage_acc: str,
    out_path: Path,
    top_k: int = 5,
) -> None:
    """
    Write top-K RBP candidates as FASTA per INTERFACE.md header format.
    按 INTERFACE.md 头部格式写入 top-K RBP 候选序列的 FASTA 文件。

    Header: >EU717894.1_rbp_01 | source_orf=... | length=... | combined_score=... | evidence=...
    """
    records = {rec.id: rec for rec in SeqIO.parse(faa_path, 'fasta')}

    # Use internal _seq_id mapping: reconstruct from faa order
    # 通过 faa 顺序重建内部 _seq_id 映射
    faa_records = list(SeqIO.parse(faa_path, 'fasta'))
    orf_to_seqid = {}
    for idx, rec in enumerate(faa_records, start=1):
        orf_id = f'{phage_acc}_orf_{idx:05d}'
        orf_to_seqid[orf_id] = rec.id

    top_candidates = candidates_df[candidates_df['passes_threshold']].head(top_k)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as fh:
        for rbp_rank, (_, row) in enumerate(top_candidates.iterrows(), start=1):
            orf_id = row['orf_id']
            seq_id = orf_to_seqid.get(orf_id, '')
            if not seq_id or seq_id not in records:
                logger.warning('Sequence not found for %s', orf_id)
                continue

            rbp_id = f'{phage_acc}_rbp_{rbp_rank:02d}'
            combined = row['combined_score']
            combined_str = f'{combined:.6f}' if not math.isnan(combined) else 'nan'
            evidence = row['evidence_track']
            length = row['length_aa']
            header = (
                f'>{rbp_id} | source_orf={orf_id} | length={length} '
                f'| combined_score={combined_str} | evidence={evidence}'
            )
            seq = str(records[seq_id].seq)
            # Wrap at 60 chars / 每行 60 个字符
            wrapped = '\n'.join(seq[i:i+60] for i in range(0, len(seq), 60))
            fh.write(f'{header}\n{wrapped}\n')

    logger.info('Wrote %d RBP sequences to %s', len(top_candidates), out_path)


def sha256_file(path: Path) -> str:
    """Compute SHA-256 hex digest of a file. 计算文件的 SHA-256 十六进制摘要。"""
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(output_dir: Path, entries: list[dict], phage_acc: str) -> None:
    """
    Write MANIFEST.csv per INTERFACE.md universal conventions.
    按 INTERFACE.md 通用规范写入 MANIFEST.csv。

    Required columns: filename, sha256, bytes, n_records, created_utc, source_acc,
                      source_module, notes
    """
    now_utc = datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    rows = []
    for e in entries:
        path = Path(e['path'])
        if not path.exists():
            continue
        rows.append({
            'filename': path.name,
            'sha256': sha256_file(path),
            'bytes': path.stat().st_size,
            'n_records': e.get('n_records', ''),
            'created_utc': now_utc,
            'source_acc': phage_acc,
            'source_module': '03_rbp_identification',
            'notes': e.get('notes', ''),
        })
    manifest_path = output_dir / 'MANIFEST.csv'
    pd.DataFrame(rows).to_csv(manifest_path, index=False)
    logger.info('Wrote MANIFEST.csv with %d entries', len(rows))


def run_phage(
    faa_path: Path,
    phage_acc: str,
    hmm_db: Path,
    output_dir: Path,
    top_k: int = 5,
) -> pd.DataFrame:
    """
    Full pipeline for one phage: HMM scan → score → write outputs.
    单个噬菌体的完整流程：HMM 扫描 → 评分 → 写入输出。

    Returns the candidates DataFrame.
    返回候选者 DataFrame。
    """
    logger.info('Processing %s from %s', phage_acc, faa_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    # HMM scan / HMM 扫描
    press_hmm_db(hmm_db)
    raw_hits = run_hmmscan(faa_path, hmm_db)
    classified = classify_hmm_hits(raw_hits)
    logger.info('HMM classified %d protein(s) as candidates for %s', len(classified), phage_acc)

    # Build candidates table / 构建候选者表格
    candidates = build_candidates_df(faa_path, phage_acc, classified)

    # Write outputs / 写入输出
    csv_path = output_dir / f'{phage_acc}_rbp_candidates.csv'
    faa_out = output_dir / f'{phage_acc}_rbp_sequences.faa'
    write_candidates_csv(candidates, csv_path)
    write_rbp_sequences_faa(faa_path, candidates, phage_acc, faa_out, top_k=top_k)

    # MANIFEST / 清单
    n_rbp = int(candidates['passes_threshold'].sum())
    write_manifest(output_dir, [
        {'path': csv_path, 'n_records': len(candidates), 'notes': f'{n_rbp} HMM-positive RBPs'},
        {'path': faa_out, 'n_records': n_rbp, 'notes': f'top-{top_k} RBP sequences'},
    ], phage_acc)

    return candidates
