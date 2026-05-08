# AGENT TODO — Module 01 / `01_data_ground_truth/`

## Read first (mandatory)

1. `/INTERFACE.md` — your outputs must conform exactly to §Module 01.
2. `/CLAUDE.md` — notebook-first, bilingual.
3. The existing `processes/fetch_reference_genomes.ipynb` — refactor, don't rewrite from scratch.
4. The existing `processes/{download_genomes,fetch_positive_pairs,fetch_negative_pairs}.py` scripts — these built the legacy interaction matrix. Read them to understand the data, then replace with a clean notebook.
5. Invoke `superpowers:test-driven-development` and `superpowers:verification-before-completion`.

## Your scope

You own NCBI fetching + the canonical interaction matrix. Reference genomes you download go into `00_raw_data/{phage,bacteria}/<acc>/` (not into your own `outputs/`); only the *index* of what you fetched lives in your `outputs/`.

You can write to `00_raw_data/phage/<new_acc>/` and `00_raw_data/bacteria/<new_acc>/` for newly-fetched genomes ONLY. Do not modify existing genome directories.

## Goal (definition of done)

By morning, `01_data_ground_truth/` contains:

1. ✅ Refactored `processes/01_fetch_reference_genomes.ipynb` matching INTERFACE §Module 01 exactly.
2. ✅ A new `processes/02_build_interaction_matrix.ipynb` that reads the legacy CSV in `outputs/interaction_matrix/` (if any) and the existing legacy scripts to produce `outputs/interaction_matrix.csv` per INTERFACE schema.
3. ✅ `inputs/reference_targets.csv` listing what we plan to fetch (Xcc, T7, phiL7).
4. ✅ `outputs/reference_genomes_index.csv` after a successful run.
5. ✅ `outputs/MANIFEST.csv`.
6. ✅ ≥3 passing tests.
7. ✅ `AGENT_REPORT.md`.

The two reference notebooks (01 + 02) must run end-to-end without manual intervention.

## Setup

```bash
cd /path/to/agent-01-data-ground-truth
conda activate igem2026

# Verify NCBI Datasets CLI is on PATH (already in env)
datasets --version
```

If `datasets` isn't found: `conda install -c conda-forge ncbi-datasets-cli` and document in SETUP.md.

## Step-by-step plan

### Step 1 — Audit existing legacy scripts (20 min)

Read `processes/download_genomes.py`, `fetch_positive_pairs.py`, `fetch_negative_pairs.py`. Take notes on:
- What CSVs they produce
- What columns those CSVs have
- Where the interaction labels come from (positive: Millard lab? PhageHostLearn dataset? negative: random sampling?)

This audit informs Step 4 (interaction matrix construction).

### Step 2 — Create `inputs/reference_targets.csv` (5 min)

Single source of truth for what Module 01 fetches. Schema:

| Column | Example |
|---|---|
| category | `bacteria` or `phage` |
| assembly_acc | `GCF_000007145.1` (or empty if Entrez-only) |
| nucleotide_acc | `AE008922` |
| label | `Xcc ATCC 33913 (host reference)` |
| subdir | `GCF_000007145.1` (target dir name under `00_raw_data/{cat}/`) |
| priority | `1` (1=must, 2=nice-to-have) |

Three rows for tonight: Xcc, T7, phiL7 (phiL7 already present — verify-only).

### Step 3 — Refactor `01_fetch_reference_genomes.ipynb` (60 min)

The current notebook is a starting point but written before INTERFACE.md was locked. Refactor to:
- Read `inputs/reference_targets.csv` (no hard-coded list of references).
- Write fetched files into `00_raw_data/{phage,bacteria}/<subdir>/{genome.fna, cds.fna, proteins.faa}`.
- Write `outputs/reference_genomes_index.csv` with status per accession.
- Write `outputs/download_log_<YYYY-MM-DD>.csv`.
- Update both `00_raw_data/MANIFEST.csv` AND `01_data_ground_truth/outputs/MANIFEST.csv` (your additions to the global manifest must merge cleanly — don't overwrite).
- Email `[email protected]` if a fetch fails 3 times in a row (use `Bio.Entrez.email = "[email protected]"` per Entrez requirement; the email is for NCBI rate-limit reporting, not actual notifications — the field is required by NCBI).

Bilingual cells throughout. Pattern: existing notebook is the canonical example.

### Step 4 — Build `02_build_interaction_matrix.ipynb` (75 min)

Convert the legacy interaction data into the standardized format from INTERFACE §Module 01.

If the legacy `outputs/interaction_matrix/` exists:
- Read the legacy CSVs.
- Map columns: legacy `phage_id, host_id, lyses?` → new `phage_acc, host_acc, label, source, confidence, notes`.
- Set `source = "literature_curated"` for legacy positive pairs.
- Set `source = "inferred_taxonomy"` for legacy negative pairs (since they were inferred, not directly tested).
- Set `confidence = 0.8` for literature-curated, `0.5` for inferred.

If legacy data is missing or unusable: write the matrix from scratch using just the 5-10 well-known phiL7 / Xcc / Xanthomonas records you can confirm from PubMed (Lee 2009, Wang 2003, etc.). Document choices in the notebook.

Output: `outputs/interaction_matrix.csv` per INTERFACE schema.

### Step 5 — Sample-then-batch test of fetcher (30 min)

Per project convention (sample first, then batch):
1. Run the fetch notebook with only the T7 phage entry first (smallest of the 3, ~40 kb).
2. Verify the output structure matches INTERFACE.
3. Then run all 3 references.

If T7 fetch fails, fix the notebook before attempting Xcc.

### Step 6 — Tests (45 min)

`processes/tests/`:
- `test_schema.py` — Load `outputs/interaction_matrix.csv`, assert all required columns. Load `outputs/reference_genomes_index.csv`, assert columns.
- `test_smoke.py` — In a temp dir, simulate fetching a tiny genome (use a dummy 100 bp FASTA file via mock; skip actual NCBI call to keep test < 30s). Assert the manifest update logic works.
- `test_sanity.py` — After running notebook 01, `00_raw_data/bacteria/GCF_000007145.1/genome.fna` exists and has total length ~5,076,187 bp (Xcc reference). For T7 phage, total length ~39,937 bp.

Run: `pytest 01_data_ground_truth/processes/tests/ -v`

### Step 7 — Commit + report

Commits:
- `01_data_ground_truth: refactor fetch notebook to match INTERFACE.md`
- `01_data_ground_truth: build interaction_matrix.csv per new schema`
- `01_data_ground_truth: add tests`
- `01_data_ground_truth: AGENT_REPORT`

Don't push.

## References (cite in notebook markdown cells)

- **NCBI Datasets CLI**: O'Leary, N.A. et al. (2024) "Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets." *Scientific Data* 11:732.
- **Bio.Entrez**: Cock, P.J.A. et al. (2009) "Biopython..." *Bioinformatics* 25:1422.
- **Reference Xcc genome**: da Silva, A.C.R. et al. (2002) "Comparison of the genomes of two *Xanthomonas* pathogens with differing host specificities." *Nature* 417:459. NCBI: AE008922 / GCF_000007145.1.
- **Reference phiL7 genome**: Lee, C.N. et al. (2009) "Genomic Characterization of the Intron-Containing T7-Like Phage phiL7 of *Xanthomonas campestris*." *Appl Environ Microbiol* 75:7828. NCBI: EU717894.
- **Phage host range datasets** (for interaction matrix curation): Boeckaerts, D. et al. (2024) "Predicting phage-host interactions in *Klebsiella* with PhageHostLearn." *Nat Commun* 15:4768. (Methodology reference even though species differs.)

## Anti-goals

- ❌ Don't re-download phiL7 — verify it.
- ❌ Don't modify the existing 777 phage / 34 bacteria directories.
- ❌ Don't write into `00_raw_data/` outside of new fetched subdirectories.
- ❌ Don't push to remote.
- ❌ Don't break the 6 PHANOTATE / Prodigal / RBP downstream agents by changing schemas — they're reading from your INTERFACE-locked outputs.

## Time budget

~4 hours. If exceeding, prioritize: notebook 01 working > notebook 02 working > tests > docs.

## If stuck

Most likely failure: NCBI rate-limit or the Datasets CLI bundle path differing between versions. If so:
1. Add explicit retries with exponential backoff (60s → 120s → 240s).
2. Fall back to `Bio.Entrez.efetch` for nucleotide-only download (no CDS/protein bundle), and document the limitation.
3. Log in `AGENT_REPORT.md`.
