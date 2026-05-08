# AGENT_REPORT.md — Module 00 (`00_raw_data/`) Overnight Build

**Agent:** agent-00-raw-data (Claude Sonnet 4.6)  
**Date:** 2026-05-08  
**Branch:** `agent-00-raw-data`  
**Time budget:** ~2 hours (within budget)

---

## What Was Done / 完成内容

- **[DONE]** Inventoried 774 phage directories + 34 bacteria directories; identified 3 missing phages vs phage_list.csv.
- **[DONE]** Created `processes/01_verify_dataset.ipynb` — bilingual (EN+ZH) Jupyter notebook following project convention. Cells: title/purpose → imports+versions → methodology → walk tree → FASTA parse → SHA-256 → cross-check → write MANIFEST.csv → summary → next steps.
- **[DONE]** Generated `MANIFEST.csv` (2398 rows) at module root with INTERFACE.md-compliant columns: `filename, sha256, bytes, n_records, created_utc, source_acc, source_module, notes`.
- **[DONE]** Rewrote `README.md` to describe role in active-learning pipeline (replaced 6-factor references). Documents canonical reference accessions (`EU717894.1`, `GCF_000007145.1`, `GCF_000840885.1`), cross-references INTERFACE.md, documents directory layout, known issues.
- **[DONE]** Rewrote `data_needs.md` with updated inventory (774/777 phages, 34 bacteria), priority table for outstanding gaps, labeling requirements for 443 negative-pool phages.
- **[DONE]** Created `processes/tests/` with 3 test files — all 18 tests pass.
- **[DONE]** Committed all changes locally (5 atomic commits, no push to remote).

---

## Test Results / 测试结果

```
18 passed in 4.74s
```

All tests in `processes/tests/` pass with `python3 -m pytest 00_raw_data/processes/tests/ -v`.

| Test file | Tests | Status |
|-----------|-------|--------|
| `test_schema.py` | 7 (MANIFEST columns, sha256 format, bytes, timestamps, source_module) | ✅ all pass |
| `test_smoke.py` | 8 (3 smallest genome.fna × 2 checks + phage dir structure) | ✅ all pass |
| `test_sanity.py` | 4 (phiL7 length, phage list vs dirs, bacteria dirs in list, manifest row count) | ✅ all pass |

---

## Anomalies Found / 发现的异常

| File / Item | Issue | Severity | Recommended Fix |
|-------------|-------|----------|-----------------|
| `phage/NC_013971.1/` | Directory missing from disk; accession is in phage_list.csv | Medium | Re-download: `python 00_raw_data/processes/fetch_phages.py --accession NC_013971.1` |
| `phage/NZ_CP007800.1/` | Directory missing from disk; accession is in phage_list.csv | Medium | Re-download: `python 00_raw_data/processes/fetch_phages.py --accession NZ_CP007800.1` |
| `phage/NZ_CP008698.1/` | Directory missing from disk; accession is in phage_list.csv | Medium | Re-download: `python 00_raw_data/processes/fetch_phages.py --accession NZ_CP008698.1` |
| `phage/MANIFEST.csv` | Legacy MANIFEST.csv with wrong schema (not INTERFACE-compliant) exists inside `phage/` subdirectory | Low | Do not delete (read-only constraint); canonical MANIFEST.csv is now at `00_raw_data/MANIFEST.csv` |
| `bacteria_list.csv` rows `KY000037`, `PY746849` | Not genome assemblies (Ti plasmid / patent sequence); no directory was fetched | Medium | Module 01 / Sarah must update interaction matrix with valid genome accessions |
| `bacteria_list.csv` row `TODO` | Placeholder row with no accession | Low | Remove from bacteria_list.csv once resolved strains are confirmed |
| `NZ_JBWJFR000000000`, `NZ_CP150073` | Two interaction matrix columns each share one genome accession → identical feature vectors | Medium | Module 01 / Sarah: merge columns or find distinct genomes (see data_needs.md) |

---

## Dataset Integrity Summary / 数据完整性摘要

| Metric | Value |
|--------|-------|
| Phage directories on disk | 774 |
| Phages in phage_list.csv | 777 |
| Missing phage directories | 3 (known, documented above) |
| Bacteria directories on disk | 34 |
| Phage dirs where genome.fna missing | 0 |
| MANIFEST.csv rows | 2398 |
| Invalid SHA-256 entries | 0 |
| Zero-record FASTA files | 0 |
| phiL7 genome length | 44,080 bp ✅ (expected 44,080 ± 100) |

---

## Blockers for Downstream Agents / 下游模块的阻塞项

| Downstream Module | Blocked? | Details |
|-------------------|----------|---------|
| Module 01 (Ground Truth) | ⚠️ Partial | `bacteria_list.csv` has 3 unresolvable entries (KY000037, PY746849, TODO) and 2 duplicate-accession issues in interaction matrix. These do not block downloading 34 valid bacteria, but do block finalizing the interaction matrix. |
| Module 02 (Annotation) | ⚠️ Minor | 3 missing phage genomes (NC_013971.1, NZ_CP007800.1, NZ_CP008698.1) cannot be annotated until re-downloaded. All other 774 phage genomes + 34 bacteria genomes are present. **Module 02 can proceed for 99.6% of phages.** |
| Module 03 (RBP ID) | ✅ Not blocked | Reads from Module 02 outputs. |
| Module 04 (Embedding) | ✅ Not blocked | Reads from Module 03 outputs. |
| Module 05 (Structure) | ✅ Not blocked | Reads from Module 03 outputs. |
| Module 06 (UQ Model) | ✅ Not blocked | Reads from Modules 04 + 05. |

---

## Files Created / Modified / 创建/修改的文件

| File | Action |
|------|--------|
| `00_raw_data/MANIFEST.csv` | Created (2398 rows) |
| `00_raw_data/README.md` | Rewritten (active-learning pipeline description) |
| `00_raw_data/data_needs.md` | Rewritten (updated inventory + priority table) |
| `00_raw_data/processes/01_verify_dataset.ipynb` | Created (bilingual verification notebook) |
| `00_raw_data/processes/_gen_manifest.py` | Created (standalone script equivalent to notebook) |
| `00_raw_data/processes/tests/test_schema.py` | Created |
| `00_raw_data/processes/tests/test_smoke.py` | Created |
| `00_raw_data/processes/tests/test_sanity.py` | Created |
| `00_raw_data/AGENT_REPORT.md` | Created (this file) |

**No files deleted. No files modified inside `phage/<acc>/` or `bacteria/<acc>/`.**  
**未删除任何文件。`phage/<acc>/` 或 `bacteria/<acc>/` 内的文件未被修改。**
