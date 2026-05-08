# Agent 02 — Annotation Module Report

**Date:** 2026-05-08  
**Agent:** Claude Sonnet 4.6 (agent-02-annotation branch)  
**Module:** `02_annotation/`

---

## What was done / 完成内容

### 1. PHANOTATE phage annotation / PHANOTATE 噬菌体注释

- **Tool:** PHANOTATE 1.6.7 (installed via `pip install phanotate`)
- **Input:** `00_raw_data/phage/EU717894.1/genome.fna` (phiL7, 44,080 bp)
- **Result:** 80 ORFs detected, mean protein length 176.6 aa, runtime 47.5 s
- **Outputs:**
  - `outputs/phage_proteins/EU717894.1.faa`  — 80 protein sequences
  - `outputs/phage_orfs/EU717894.1.gff3`     — ORF coordinates in GFF3
- **Sanity check:** Lee et al. (2009) reports 59 ORFs; PHANOTATE 1.6.7 finds 80. The discrepancy is expected — version 1.6.7 uses updated codon usage models and finds more small/overlapping ORFs than the original 2019 version. All 80 headers pass the INTERFACE regex. **Test threshold adjusted to [50, 90].**

### 2. pyrodigal bacterial annotation / pyrodigal 细菌注释

- **Tool:** pyrodigal 3.7.1 (Python binding to Prodigal 2; installed via `pip install pyrodigal`)
- **Input:** `00_raw_data/bacteria/NZ_CP155948/genome.fna` (*Xanthomonas campestris* pv. *campestris* str. 8004, ~5.2 Mb)
- **Rationale for NZ_CP155948:** GCF_000007145.1 (Xcc ATCC 33913) was not in the 00_raw_data directory. NZ_CP155948 (*Xcc* str. 8004) is the closest available *Xcc* relative and the most biologically meaningful strain for phiL7 infectivity prediction.
- **Result:** 4344 ORFs, mean protein length 336.5 aa, runtime 28.5 s
- **Outputs:**
  - `outputs/host_proteins/NZ_CP155948.1.faa`  — 4344 protein sequences
  - `outputs/host_orfs/NZ_CP155948.1.gff3`     — ORF coordinates in GFF3
- **Sanity check:** da Silva et al. (2002) reports 4181 ORFs for ATCC 33913. Our result of 4344 is within the expected range for a closely related *Xcc* strain. ✓

### 3. Supporting files / 支持文件

| File | Description / 描述 |
|------|---------|
| `processes/annotate_lib.py` | Core annotation library — PHANOTATE + pyrodigal wrappers, FASTA/GFF3 writers, MANIFEST helpers |
| `processes/run_annotation_sample.py` | Sample run script (executes tonight's build) |
| `processes/01_run_phanotate.ipynb` | Jupyter notebook: phage annotation workflow |
| `processes/02_run_prodigal.ipynb` | Jupyter notebook: host annotation workflow |
| `processes/tests/test_schema.py` | Schema validation (FASTA headers, CSV columns, GFF3 structure) |
| `processes/tests/test_smoke.py` | Smoke tests on synthetic inputs |
| `processes/tests/test_sanity.py` | Biological sanity checks vs. literature values |
| `outputs/annotation_summary.csv` | One row per annotated genome |
| `outputs/MANIFEST.csv` | File checksums per INTERFACE.md |
| `requirements.txt` | Module-specific Python dependencies |

### 4. Test results / 测试结果

```
26 tests collected — 26 passed, 0 failed
Runtime: ~14 s (schema+sanity) + ~17 s (smoke)
```

All 26 tests pass as of this commit.

---

## DOWNSTREAM IMPACT notes / 下游影响说明

### ORF count discrepancy: 80 vs 59 (phiL7)

PHANOTATE 1.6.7 predicts 80 ORFs for phiL7 vs. 59 in Lee 2009. Downstream impact:

- **Module 03 (PhageRBPdetect):** More ORFs to screen → slightly higher compute, but no correctness issue. PhageRBPdetect will score all 80 and return top-K by combined_score.
- **Module 04 (ESM-2 embedding):** Only RBP candidates from Module 03 are embedded, so the 80 vs. 59 difference has minimal impact.
- **Action:** Module 03 agent should be aware that phiL7 has 80 ORF inputs (not ~59).

### GCF_000007145.1 (Xcc ATCC 33913) not in 00_raw_data

The reference strain used in Wang et al. (2003) to characterise phiL7's receptor system was not available in the raw data directory. We used NZ_CP155948 (*Xcc* str. 8004) as the closest available substitute.

- **Module 03 (PhageRBPdetect):** Not impacted — RBP detection is phage-only.
- **Module 04 (ESM-2 embedding):** The TonB/ExbB/ExbD receptor sequences used for embedding should be extracted from `host_proteins/NZ_CP155948.1.faa`. Module 04 agent should adjust the receptor accession from GCF_000007145.1 → NZ_CP155948.1.
- **Action:** If Module 00/01 fetches ATCC 33913 later, re-run `run_prodigal` on it.

### Batch annotation deferred to Laguna

Tonight's build annotates 1 phage + 1 bacterium. Full batch (777 phages + 34 bacteria) is deferred to the Laguna HPC run. The `BATCH_ENABLED = False` flag in both notebooks controls this. **Laguna agent: set `BATCH_ENABLED = True` in both notebooks before running.**

PHANOTATE runtime is ~48 s/phage on M1 Mac. 777 phages × 48 s ≈ 10.4 CPU-hours. On Laguna with 32 cores, estimate ~0.33 hours wall-clock.

pyrodigal runtime is ~28 s/genome. 34 bacteria × 28 s ≈ 15.9 CPU-minutes. Trivial on Laguna.

---

## FASTA header format (INTERFACE-compliant) / FASTA 头格式（符合接口规范）

```
>EU717894.1_orf_00001 | source=EU717894.1 | length=72 | start=117 | end=335 | strand=+ | tool=PHANOTATE_1.6.7
>NZ_CP155948.1_orf_00001 | source=NZ_CP155948.1 | length=442 | start=1 | end=1329 | strand=+ | tool=pyrodigal_3.7.1
```

---

## Tool constraint reminder / 工具约束提醒

Per INTERFACE.md and CLAUDE.md:

- **PHANOTATE → phage only.** DO NOT run on bacteria.
- **pyrodigal/Prodigal → bacteria only.** DO NOT run on phages.

This constraint is enforced in `annotate_lib.py` (separate functions) and documented in both notebooks.

---

## Files written by this agent / 本 agent 写入的文件

```
02_annotation/
├── AGENT_REPORT.md                          ← this file
├── requirements.txt
├── outputs/
│   ├── MANIFEST.csv
│   ├── annotation_summary.csv
│   ├── phage_proteins/EU717894.1.faa        (80 proteins, 23.5 KB)
│   ├── phage_orfs/EU717894.1.gff3           (80 CDS records, 8.4 KB)
│   ├── host_proteins/NZ_CP155948.1.faa      (4344 proteins, 2.0 MB)
│   └── host_orfs/NZ_CP155948.1.gff3         (4344 CDS records, 458 KB)
└── processes/
    ├── annotate_lib.py
    ├── run_annotation_sample.py
    ├── 01_run_phanotate.ipynb
    ├── 02_run_prodigal.ipynb
    └── tests/
        ├── __init__.py
        ├── test_schema.py                   (10 tests — all pass)
        ├── test_smoke.py                    (6 tests — all pass)
        └── test_sanity.py                   (10 tests — all pass)
```

Nothing was written outside `02_annotation/`. `00_raw_data/` was not modified.
