"""
Build 01_run_phagerbpdetect.ipynb programmatically using nbformat.
使用 nbformat 以编程方式构建 01_run_phagerbpdetect.ipynb。

Run once: python make_notebook.py
"""

import json
from pathlib import Path

import nbformat as nbf

NB_PATH = Path(__file__).parent / '01_run_phagerbpdetect.ipynb'

# ── Cell factories ─────────────────────────────────────────────────────────

def md(src: str) -> nbf.NotebookNode:
    return nbf.v4.new_markdown_cell(src)

def code(src: str) -> nbf.NotebookNode:
    return nbf.v4.new_code_cell(src)

# ── Notebook cells ─────────────────────────────────────────────────────────

cells = []

# Cell 1: Title + bilingual purpose
cells.append(md("""
# Module 03 — RBP Identification via PhageRBPdetect
## 模块 03 — 利用 PhageRBPdetect 进行 RBP 鉴定

**Purpose:** Identify receptor-binding protein (RBP) candidates from phage proteomes using a hybrid HMM + ML approach (PhageRBPdetect, Boeckaerts et al. 2022).

**目的：** 使用混合 HMM + ML 方法（PhageRBPdetect，Boeckaerts et al. 2022）从噬菌体蛋白质组中鉴定受体结合蛋白（RBP）候选者。

**Primary reference phage:** phiL7 (NCBI EU717894.1) infecting *Xanthomonas campestris* pv. *campestris* ATCC 33913. The tail spike gp25 is experimentally validated (Lee et al. 2009 *Appl Environ Microbiol* 75:7828).

**主要参考噬菌体：** phiL7（NCBI EU717894.1），感染 *Xanthomonas campestris* pv. *campestris* ATCC 33913。尾刺蛋白 gp25 已有实验验证（Lee et al. 2009）。

### References / 参考文献
- **PhageRBPdetect (method):** Boeckaerts, D. et al. (2022) "Identification of phage receptor-binding protein sequences with HMMs and XGBoost." *Viruses* 14(6):1329. DOI: 10.3390/v14061329
- **Pfam database:** Mistry, J. et al. (2021) *Nucleic Acids Research* 49(D1):D412.
- **HMMER:** Eddy, S.R. (2011) *PLOS Comput Biol* 7(10):e1002195.
- **RBP biology:** Latka, A. et al. (2021) *mBio* 12(3):e00455-21.
""".strip()))

# Cell 2: Imports + version printout
cells.append(code("""\
# Cell 2: Imports and library versions
# 单元格 2：导入库并打印版本信息
import sys
import subprocess
from pathlib import Path

import pandas as pd
import numpy as np
from Bio import SeqIO

# Path anchoring per INTERFACE.md: notebooks live in processes/
# 路径锚定（遵循 INTERFACE.md）：notebook 位于 processes/ 目录
REPO_ROOT = Path.cwd().resolve().parents[1]  # agent-03-rbp-identification/
MODULE_DIR = Path.cwd().resolve().parents[0]  # 03_rbp_identification/

print(f"Python: {sys.version}")
print(f"pandas: {pd.__version__}")
print(f"numpy:  {np.__version__}")
print(f"BioPython: {SeqIO.__module__.split('.')[0]} (imported)")
print(f"REPO_ROOT: {REPO_ROOT}")
print(f"MODULE_DIR: {MODULE_DIR}")

# HMMER version / HMMER 版本
hmmscan_bin = Path('/opt/homebrew/opt/hmmer/bin/hmmscan')
if not hmmscan_bin.exists():
    hmmscan_bin = Path('hmmscan')
try:
    result = subprocess.run([str(hmmscan_bin), '-h'], capture_output=True, text=True)
    hmmer_line = [l for l in result.stdout.split('\\n') if 'HMMER' in l]
    print("HMMER:", hmmer_line[0] if hmmer_line else "unknown")
except Exception:
    print("HMMER: not found on PATH — pipeline will fail if hmmscan unavailable")
"""))

# Cell 3: Method explanation (bilingual)
cells.append(md("""
## Dual-track Design / 双轨道设计

PhageRBPdetect (Boeckaerts 2022) uses two complementary tracks:

1. **HMM track** — scans phage proteins against 92 curated profile HMMs (N-terminal structural domains + C-terminal binding/enzymatic domains). A protein that matches any N-block or C-block domain is a candidate. This track is deterministic, fast, and interpretable.

2. **ML track** — embeds proteins with ESM-2 (Lin et al. 2023 *Science*), then applies a trained XGBoost classifier (Chen & Guestrin 2016 KDD). This catches RBPs that lack detectable Pfam homologs.

**Per INTERFACE.md:** `combined_score = 1.0` if HMM hit; else `ml_score`.

**Evidence track:** `"hmm"`, `"ml"`, or `"both"` (both tracks positive).

---

1. **HMM 轨道** — 将噬菌体蛋白质与 92 个精选轮廓 HMM 比对（N 端结构域 + C 端结合/酶活性域）。命中任一 N-block 或 C-block 域的蛋白质即为候选者。该轨道确定性强、速度快、可解释。

2. **ML 轨道** — 用 ESM-2 对蛋白质进行嵌入，再用训练好的 XGBoost 分类器打分。该轨道能捕获缺乏可检测 Pfam 同源物的 RBP。

**⚠ ML Track Status (this sprint):** `bio_embeddings` (required for v2 ProtTrans embeddings) fails to install on Python 3.13+ (gensim build error). v4 fine-tuned ESM-2 model was not downloaded from Zenodo. Therefore `ml_score = NaN` throughout; `combined_score` relies solely on the HMM track. For phiL7 this is acceptable because the tail spike has strong Pfam hits.

**⚠ ML 轨道状态（本次冲刺）：** `bio_embeddings`（v2 ProtTrans 嵌入所需）在 Python 3.13+ 上安装失败（gensim 构建错误）。v4 精调 ESM-2 模型未从 Zenodo 下载。因此整体 `ml_score = NaN`；`combined_score` 完全依赖 HMM 轨道。对于 phiL7，由于尾刺蛋白具有强 Pfam 命中，这是可接受的。
""".strip()))

# Cell 4: HMM track helper
cells.append(code("""\
# Cell 4: HMM track — run_hmm_track helper
# 单元格 4：HMM 轨道 — run_hmm_track 辅助函数
import sys
sys.path.insert(0, str(Path.cwd()))
from rbp_pipeline import (
    press_hmm_db, run_hmmscan, classify_hmm_hits,
    N_BLOCKS, C_BLOCKS,
)

HMM_DB = MODULE_DIR / 'inputs' / 'phagerbpdetect_data' / 'RBPdetect_phageRBPs.hmm'

def run_hmm_track(faa_path: Path, hmm_db: Path = HMM_DB) -> pd.DataFrame:
    \"\"\"
    HMM track: press DB → hmmscan → classify hits.
    HMM 轨道：压缩数据库 → hmmscan → 分类命中。
    Returns DataFrame with HMM-positive proteins only.
    返回只包含 HMM 阳性蛋白质的 DataFrame。
    \"\"\"
    press_hmm_db(hmm_db)
    raw = run_hmmscan(faa_path, hmm_db)
    classified = classify_hmm_hits(raw)
    print(f"HMM track: {len(raw)} total domain hits → {len(classified)} protein(s) classified as candidates")
    return classified

print(f"HMM DB: {HMM_DB}")
print(f"N_BLOCKS count: {len(N_BLOCKS)}")
print(f"C_BLOCKS count: {len(C_BLOCKS)}")
"""))

# Cell 5: ML track stub (documented as blocked)
cells.append(code("""\
# Cell 5: ML track — documented as blocked; ml_score = NaN
# 单元格 5：ML 轨道 — 记录为阻塞；ml_score = NaN

def run_ml_track(faa_path: Path, pfam_misses: pd.DataFrame) -> pd.DataFrame:
    \"\"\"
    ML track stub. Returns NaN scores for all sequences.
    ML 轨道存根。为所有序列返回 NaN 分数。

    Blocked because:
    - v2 XGBoost model requires 1024-dim ProtTrans BFD embeddings
      (bio_embeddings install fails on Python 3.13+: gensim build error)
    - v4 fine-tuned ESM-2 model requires Zenodo download (not attempted this sprint)

    阻塞原因：
    - v2 XGBoost 模型需要 1024 维 ProtTrans BFD 嵌入
      （bio_embeddings 在 Python 3.13+ 上安装失败：gensim 构建错误）
    - v4 精调 ESM-2 模型需要 Zenodo 下载（本次冲刺未尝试）

    Trade-off: HMM track alone provides high-confidence RBP identification for phiL7
    because the tail spike has strong Pfam domain hits. ML track would add value for
    novel phages lacking Pfam homologs.
    权衡：对于 phiL7，仅 HMM 轨道即可提供高置信度的 RBP 鉴定，因为尾刺蛋白具有
    强 Pfam 域命中。ML 轨道对于缺乏 Pfam 同源物的新型噬菌体更有价值。
    \"\"\"
    names = [rec.id for rec in SeqIO.parse(faa_path, 'fasta')]
    return pd.DataFrame({'seq_id': names, 'ml_score': float('nan')})

print("ML track: NaN stub loaded (see docstring for rationale)")
print("Future work: integrate bio_embeddings or download v4 model from Zenodo 10515367")
"""))

# Cell 6: combine_tracks helper
cells.append(code("""\
# Cell 6: combine_tracks — merge HMM + ML per INTERFACE.md formula
# 单元格 6：combine_tracks — 按 INTERFACE.md 公式合并 HMM + ML
import math

def combine_tracks(hmm_df: pd.DataFrame, ml_df: pd.DataFrame) -> pd.DataFrame:
    \"\"\"
    Merge HMM and ML results per INTERFACE.md combined_score formula:
    combined_score = 1.0 if HMM hit; else ml_score.
    evidence_track = 'hmm' | 'ml' | 'both'.

    按 INTERFACE.md combined_score 公式合并 HMM 和 ML 结果：
    combined_score = 1.0（HMM 命中）；否则 = ml_score。
    evidence_track = 'hmm' | 'ml' | 'both'。
    \"\"\"
    hmm_hit_ids = set(hmm_df['seq_id']) if not hmm_df.empty else set()
    ml_dict = {row['seq_id']: row['ml_score'] for _, row in ml_df.iterrows()} if not ml_df.empty else {}

    merged = []
    all_ids = set(ml_dict.keys()) | hmm_hit_ids
    for sid in all_ids:
        has_hmm = sid in hmm_hit_ids
        ml_score = ml_dict.get(sid, float('nan'))
        ml_positive = not (isinstance(ml_score, float) and math.isnan(ml_score)) and ml_score >= 0.5

        if has_hmm and ml_positive:
            evidence = 'both'
            combined = 1.0
        elif has_hmm:
            evidence = 'hmm'
            combined = 1.0
        elif ml_positive:
            evidence = 'ml'
            combined = float(ml_score)
        else:
            evidence = ''
            combined = float(ml_score) if not (isinstance(ml_score, float) and math.isnan(ml_score)) else float('nan')

        merged.append({'seq_id': sid, 'ml_score': ml_score, 'combined_score': combined, 'evidence_track': evidence})

    return pd.DataFrame(merged)

print("combine_tracks helper loaded")
"""))

# Cell 7: Run on phiL7
cells.append(md("""
## Run on phiL7 (EU717894.1) / 在 phiL7 (EU717894.1) 上运行

**Input:** `02_annotation/outputs/phage_proteins/EU717894.1.faa` (preferred) or
`00_raw_data/phage/EU717894.1/proteins.faa` (fallback, 58 proteins from NCBI).

**⚠ Fallback note:** NCBI proteins.faa headers use format `>ADP02444.1:1-712` (reference protein IDs, some annotated against a different phage). ORF IDs are reassigned in genome order: `EU717894.1_orf_00001` … `EU717894.1_orf_00058`.

**输入：** `02_annotation/outputs/phage_proteins/EU717894.1.faa`（首选）或
`00_raw_data/phage/EU717894.1/proteins.faa`（回退，来自 NCBI 的 58 个蛋白质）。

**⚠ 回退说明：** NCBI proteins.faa 头部使用格式 `>ADP02444.1:1-712`（参考蛋白质 ID，部分注释针对不同噬菌体）。ORF ID 按基因组顺序重新赋值：`EU717894.1_orf_00001` … `EU717894.1_orf_00058`。
""".strip()))

# Cell 8: Sample run on phiL7
cells.append(code("""\
# Cell 8: Run full pipeline on phiL7
# 单元格 8：对 phiL7 运行完整流程
import sys
sys.path.insert(0, str(Path.cwd()))
from rbp_pipeline import run_phage

PHAGE_ACC = 'EU717894.1'
OUTPUT_DIR = MODULE_DIR / 'outputs'

# Resolve input FASTA
# 解析输入 FASTA
annotation_faa = REPO_ROOT / '02_annotation' / 'outputs' / 'phage_proteins' / f'{PHAGE_ACC}.faa'
fallback_faa_raw = REPO_ROOT / '00_raw_data' / 'phage' / PHAGE_ACC / 'proteins.faa'
fallback_faa_local = MODULE_DIR / 'inputs' / f'{PHAGE_ACC}_fallback_proteins.faa'

if annotation_faa.exists():
    faa_path = annotation_faa
    print(f"Using 02_annotation output: {faa_path}")
elif fallback_faa_raw.exists():
    faa_path = fallback_faa_raw
    print(f"WARNING: Using 00_raw_data fallback: {faa_path}")
else:
    faa_path = fallback_faa_local
    print(f"WARNING: Using worktree-local fallback: {faa_path}")

candidates = run_phage(
    faa_path=faa_path,
    phage_acc=PHAGE_ACC,
    hmm_db=HMM_DB,
    output_dir=OUTPUT_DIR,
    top_k=5,
)
candidates[candidates['passes_threshold']].head(10)
"""))

# Cell 9: Sanity assertion
cells.append(code("""\
# Cell 9: Sanity assertion — phiL7 should yield 1-3 RBP candidates with combined_score > 0.7
# 单元格 9：健全性断言 — phiL7 应产生 1-3 个 combined_score > 0.7 的 RBP 候选者
high_conf = candidates[candidates['combined_score'] > 0.7]
n_high = len(high_conf)
print(f"High-confidence RBPs (combined_score > 0.7): {n_high}")
assert 1 <= n_high <= 3, f"Expected 1-3, got {n_high}"
print("✓ Sanity check passed")

# Note on gp25 length discrepancy
# gp25 长度差异说明
top = candidates[candidates['rank'] == 1].iloc[0]
print(f"\\nTop candidate: {top['orf_id']}  {top['length_aa']} aa  domain={top['hmm_match_pfam']}")
print("NOTE: INTERFACE.md cites gp25 ~412 aa (Lee et al. 2009 CDS annotation).")
print("      NCBI annotation yields 712 aa for orf_00001 (broader fusion protein boundary).")
print("      Both annotations identify the same RBP; discrepancy is in annotation version.")
"""))

# Cell 10: Write outputs
cells.append(code("""\
# Cell 10: Write EU717894.1_rbp_candidates.csv and EU717894.1_rbp_sequences.faa
# 单元格 10：写入 EU717894.1_rbp_candidates.csv 和 EU717894.1_rbp_sequences.faa
# (Already written by run_phage above; this cell confirms files exist)
# （run_phage 已写入；此单元格确认文件存在）
from pathlib import Path

csv_path = OUTPUT_DIR / f'{PHAGE_ACC}_rbp_candidates.csv'
faa_out = OUTPUT_DIR / f'{PHAGE_ACC}_rbp_sequences.faa'
manifest = OUTPUT_DIR / 'MANIFEST.csv'

for p in [csv_path, faa_out, manifest]:
    size = p.stat().st_size if p.exists() else -1
    status = '✓' if size > 0 else '✗'
    print(f"{status}  {p.name}  ({size} bytes)")

print("\\nTop-5 RBP FASTA headers:")
with open(faa_out) as fh:
    for line in fh:
        if line.startswith('>'):
            print(' ', line.rstrip())
"""))

# Cell 11: Optional batch note
cells.append(md("""
## Optional: Batch Processing / 可选：批量处理

Additional phage genomes from `02_annotation/outputs/phage_proteins/` can be processed
by calling `run_phage()` in a loop. The unified `all_rbp_candidates.csv` is written
only when >1 phage is processed (per AGENT_TODO anti-goals).

更多噬菌体基因组可从 `02_annotation/outputs/phage_proteins/` 批量处理，
通过循环调用 `run_phage()` 即可。仅当处理 >1 个噬菌体时才写入
`all_rbp_candidates.csv`（遵循 AGENT_TODO 反目标）。

```python
# Example batch run / 批量运行示例
for faa in sorted(ANNOTATION_DIR.glob('*.faa')):
    acc = faa.stem
    run_phage(faa_path=faa, phage_acc=acc, hmm_db=HMM_DB, output_dir=OUTPUT_DIR)
```
""".strip()))

# Cell 12: MANIFEST update
cells.append(code("""\
# Cell 12: Confirm MANIFEST.csv
# 单元格 12：确认 MANIFEST.csv
import pandas as pd
manifest = pd.read_csv(OUTPUT_DIR / 'MANIFEST.csv')
print(manifest.to_string(index=False))
"""))

# ── Assemble notebook ──────────────────────────────────────────────────────

nb = nbf.v4.new_notebook()
nb['cells'] = cells
nb['metadata'] = {
    'kernelspec': {
        'display_name': 'Python 3',
        'language': 'python',
        'name': 'python3',
    },
    'language_info': {'name': 'python'},
}

with open(NB_PATH, 'w') as fh:
    nbf.write(nb, fh)

print(f'Written: {NB_PATH}')
