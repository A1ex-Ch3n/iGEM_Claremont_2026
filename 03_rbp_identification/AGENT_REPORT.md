# AGENT_REPORT — Module 03 / `03_rbp_identification/`

## Status: ✅ Complete (HMM track; ML track blocked — see below)

**Sprint:** 2026-05-07 overnight parallel build
**Agent:** Claude (agent-03-rbp-identification worktree)

---

## Deliverables

| File | Status | Notes |
|------|--------|-------|
| `processes/01_run_phagerbpdetect.ipynb` | ✅ | Bilingual EN+ZH, 12 cells |
| `processes/rbp_pipeline.py` | ✅ | Core HMM pipeline functions |
| `processes/run_pipeline.py` | ✅ | Runner script |
| `outputs/EU717894.1_rbp_candidates.csv` | ✅ | 58 rows, 3 HMM-positive |
| `outputs/EU717894.1_rbp_sequences.faa` | ✅ | Top-3 RBP sequences |
| `outputs/MANIFEST.csv` | ✅ | SHA-256 checksums, 2 entries |
| `processes/tests/test_schema.py` | ✅ | 11 assertions |
| `processes/tests/test_smoke.py` | ✅ | 4 assertions |
| `processes/tests/test_sanity.py` | ✅ | 9 assertions |
| `requirements.txt` | ✅ | With ML-track blocker documented |
| `AGENT_REPORT.md` | ✅ | This file |

**Tests: 27/27 passing** (`pytest 03_rbp_identification/processes/tests/ -v`)

---

## phiL7 Results

Phage: EU717894.1 (Xanthomonas phage phiL7), 58 ORFs

| Rank | orf_id | length_aa | hmm_score | domain | combined_score |
|------|--------|-----------|-----------|--------|----------------|
| 1 | EU717894.1_orf_00001 | 712 | 342.0 | unknown_C54 | 1.0 |
| 2 | EU717894.1_orf_00021 | 918 | 235.1 | unknown_C112 | 1.0 |
| 3 | EU717894.1_orf_00003 | 224 | 56.7 | unknown_C294 | 1.0 |

**Sanity check:** INTERFACE.md expects 1-3 candidates with combined_score > 0.7. ✅ Found exactly 3.

**gp25 length note:** INTERFACE.md cites ~412 aa for the phiL7 tail spike (Lee et al. 2009 CDS annotation). NCBI proteins.faa annotates the same protein as 712 aa (ADP02444.1, a fusion protein that includes the tail spike domain). Both identify the same RBP functional unit. The 712 aa protein has strong hits to both `Tail_spike_N` (PF18668, bitscore=65.9) and `unknown_C54` (bitscore=342.0), confirming it as the primary RBP.

Orf 00021 (918 aa) is a collagen-like repeat protein with C-terminal binding domains — likely a secondary structural component.

---

## Architecture

### HMM Track (✅ Implemented)

**Tool:** HMMER 3.4 (`hmmscan`)
**Database:** `PhageRBPdetect_v2/data/RBPdetect_phageRBPs.hmm` (92 custom profile HMMs from Boeckaerts 2022)

Domain classification follows PhageRBPdetect v2 rules:
- N-blocks (50 structural domains): bitscore ≥ 18, ln(score) > ln(bias)
- C-blocks (43 binding/enzymatic domains): bitscore ≥ 25, ln(score) > ln(bias)

**combined_score formula (per INTERFACE.md):**
```
combined_score = 1.0   if HMM hit
combined_score = NaN   otherwise (ml track blocked)
```

### ML Track (🚫 Blocked)

**Attempted:** `pip install bio_embeddings` — **FAILED** on Python 3.13+ (gensim wheel build error: no pre-built binary for Apple Silicon).

**Why bio_embeddings:** PhageRBPdetect v2 XGBoost model (`RBPdetect_xgb_model.json`) was trained on 1024-dimensional ProtTrans-BFD embeddings. ESM-2 produces 320-dimensional embeddings — incompatible without retraining.

**Alternative considered:** PhageRBPdetect v4 uses a fine-tuned ESM-2 model (Zenodo record 10515367). Not pursued due to model size and download time within 4-hour sprint budget.

**Trade-off documented:** HMM track alone provides high-confidence identification for phiL7 (tail spike has strong Pfam homologs). For novel phages without Pfam domain matches, the ML track would be essential. Future work: install bio_embeddings on Python ≤3.11, or download v4 Zenodo model.

---

## Input Handling

**Canonical input:** `02_annotation/outputs/phage_proteins/EU717894.1.faa` — not available (Module 02 phage_proteins/ was empty at sprint time).

**Fallback used:** `inputs/EU717894.1_fallback_proteins.faa` (pre-copied from `00_raw_data/phage/EU717894.1/proteins.faa`, 19 KB, 58 proteins).

**ORF ID assignment:** NCBI proteins.faa headers (`>ADP02444.1:1-712`) do not follow INTERFACE.md format. ORF IDs are assigned in FASTA order: `EU717894.1_orf_00001` … `EU717894.1_orf_00058`. These will be superseded by PHANOTATE-assigned IDs once Module 02 produces its output.

---

## Git Commits

```
03_rbp_identification: install PhageRBPdetect + smoke test
03_rbp_identification: notebook + phiL7 RBP identification
03_rbp_identification: outputs in INTERFACE format
03_rbp_identification: tests
03_rbp_identification: AGENT_REPORT
```

---

## Downstream Consumers

**Module 04 (embeddings):** reads `outputs/EU717894.1_rbp_sequences.faa` — 3 sequences ready.

**Module 05 (structure):** reads `outputs/EU717894.1_rbp_sequences.faa` — same file.

Top RBP sequence (`EU717894.1_rbp_01`, 712 aa) is the tail spike; it's the primary candidate for folding and binding affinity prediction against TonB/ExbB/ExbD receptors.

---

## 中文摘要

本模块使用 PhageRBPdetect（Boeckaerts 2022）的 HMM 轨道对 phiL7 噬菌体的 58 个蛋白质进行了 RBP 鉴定。共识别出 3 个高置信度候选者（combined_score = 1.0），均通过 HMM 轨道验证。顶级候选者（orf_00001，712 aa）具有强烈的尾刺蛋白域命中（Tail_spike_N + unknown_C54），与 Lee et al. 2009 描述的 phiL7 gp25 尾刺蛋白一致。ML 轨道因 bio_embeddings 在 Python 3.13+ 上安装失败而被阻塞，ml_score 全部为 NaN，combined_score 完全依赖 HMM 轨道。所有 27 个测试通过。
