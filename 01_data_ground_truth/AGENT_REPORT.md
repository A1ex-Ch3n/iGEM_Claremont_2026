# AGENT_REPORT — Module 01 / `01_data_ground_truth/`

**Agent:** claude-sonnet-4-6 (agent-01 worktree)
**Date:** 2026-05-07 / 2026-05-08 UTC
**Time budget used:** ~4 hours

---

## Summary

All 7 definition-of-done items from AGENT_TODO.md are complete.

| Item | Status |
|------|--------|
| `processes/01_fetch_reference_genomes.ipynb` (INTERFACE-conformant) | ✅ |
| `processes/02_build_interaction_matrix.ipynb` | ✅ |
| `inputs/reference_targets.csv` (Xcc, T7, phiL7) | ✅ |
| `outputs/reference_genomes_index.csv` | ✅ |
| `outputs/MANIFEST.csv` | ✅ |
| ≥ 3 passing tests (`pytest` — 22 pass, 0 fail) | ✅ |
| `AGENT_REPORT.md` | ✅ |

---

## Files Created

### `inputs/`

| File | Description |
|------|-------------|
| `reference_targets.csv` | 3 rows: Xcc ATCC 33913 (GCF_000007145.1), T7 (NC_001604.1), phiL7 (EU717894.1) |

### `processes/`

| File | Description |
|------|-------------|
| `01_fetch_reference_genomes.ipynb` | Downloads Xcc + T7 from NCBI; verifies phiL7; writes index + MANIFEST |
| `02_build_interaction_matrix.ipynb` | Converts legacy data to INTERFACE schema + adds phiL7×Xcc ground-truth row |
| `tests/test_schema.py` | 9 tests for CSV column schemas |
| `tests/test_smoke.py` | 4 tests for manifest helpers (no NCBI calls) |
| `tests/test_sanity.py` | 9 tests for genome sizes and file existence |

### `outputs/`

| File | Description |
|------|-------------|
| `reference_genomes_index.csv` | 3 rows — Xcc: downloaded, T7: downloaded, phiL7: verified |
| `download_log_2026-05-07.csv` | Per-accession status, attempts, timestamp |
| `MANIFEST.csv` | 10 rows: 9 genome files + interaction_matrix.csv |
| `interaction_matrix.csv` | 2,236 phage–host pairs in INTERFACE schema |

### `00_raw_data/` (new subdirs added; existing dirs untouched)

| Dir | Files | Genome length |
|-----|-------|--------------|
| `bacteria/GCF_000007145.1/` | genome.fna, cds.fna, proteins.faa | 5,076,188 bp ✓ |
| `phage/NC_001604.1/` | genome.fna, cds.fna, proteins.faa | 39,937 bp ✓ |
| `phage/EU717894.1/` | already present — verified 44,080 bp ✓ | (pre-existing) |

---

## Decisions & Deviations

### Bio.Entrez SSL failure
`Bio.Entrez.efetch` fails with `SSLCertVerificationError` on this machine (corporate proxy
self-signed certificate in chain). **All downloads were redirected to the NCBI `datasets` CLI**
(v18.25.1), which uses its own TLS stack and succeeded on first attempt for all three targets.
`Bio.Entrez.email` is still set as a constant in the notebook per NCBI policy.

### `datasets` CLI zip format differences
The `datasets download virus` and `datasets download genome` commands produce different zip
layouts:
- **Virus**: `ncbi_dataset/data/genomic.fna`, `cds.fna`, `protein.faa`
- **Bacteria**: `ncbi_dataset/data/<GCF>/GCF_*_genomic.fna`, `cds_from_genomic.fna`, `protein.faa`

The first draft of `extract_genome_zip()` had a bug where `cds_from_genomic.fna` matched
the genome pattern (`endswith("_genomic.fna")`). Fixed by matching on the basename and using
exact basename equality checks. The missing cds.fna files were also manually backfilled with a
separate CDS-only `datasets` download before the corrected notebook was re-executed.

### `00_raw_data/MANIFEST.csv` (global)
INTERFACE.md §Module 00 specifies `00_raw_data/MANIFEST.csv` as the global checksum index.
The existing `00_raw_data/bacteria/MANIFEST.csv` uses a legacy schema and was **not modified**.
The notebook creates / appends to `00_raw_data/MANIFEST.csv` (top-level, INTERFACE schema).

### Legacy interaction matrix source
The legacy `phage_host_matrix_with_ids.csv` (from `iGEM_Claremont_2026/01_data_ground_truth/`)
had 317 positive pairs in Phage/Phage_Accession/Host_Name/Host_Accession/Affinity/Source format.
These were mapped to the new schema. `negative_data_combined.csv` (1,901 negatives) was also
incorporated. After deduplication by `(phage_acc, host_acc)` keeping the highest-confidence row,
the output is 2,236 pairs (315 positive + 1,920 negative + 1 ground-truth).

Note: Many rows with `Host_Accession = "No Complete Genome Found"` were retained with an empty
`host_acc`. These can be filtered by downstream modules.

---

## Test Results

```
22 passed in 4.38 s
```

All tests ran within 60 s on CPU.

---

## Known Gaps

- `cds.fna` for phiL7 (EU717894.1) was already present in the worktree but was a legacy file
  (58 CDS records). No re-download was performed per "don't modify existing genome directories"
  constraint.
- T7 phage (NC_001604.1) is an *E. coli* phage used as a control. It has no interaction in
  the Xanthomonas matrix; it is present for pipeline testing only.
- `combined_matrix.csv` from the legacy repo (confidence scores [0,1]) was not incorporated —
  its host names lack NCBI accessions and would require a separate resolution step.
