# AGENT TODO — Module 00 / `00_raw_data/`

## Read first (mandatory)

1. `/INTERFACE.md` — data contract for **all** outputs you write.
2. `/CLAUDE.md` — general project conventions (notebook-first, bilingual comments).
3. The existing files in `00_raw_data/` — do **not** delete; verify and document.
4. Invoke `superpowers:test-driven-development` and `superpowers:verification-before-completion`.

## Your scope

You are responsible for `00_raw_data/` only. Do NOT touch any other module folder, do NOT touch `docs/`, `shared/`, `archive/`. The 777 phage + 34 bacteria genomes already in this folder are **read-only** from your perspective; you verify them and document them, you don't re-download.

If you find a corrupted file, log it in `AGENT_REPORT.md` and create a stub fetch notebook to re-download just that one. Don't mass re-download.

## Goal (definition of done)

By morning, `00_raw_data/` contains:

1. ✅ A runnable notebook `processes/01_verify_dataset.ipynb` that audits every file, computes SHA-256, validates FASTA parseability, and writes `MANIFEST.csv`.
2. ✅ A `MANIFEST.csv` at module root conforming to INTERFACE §Universal conventions.
3. ✅ Updated `README.md` describing the dataset's role in the new active-learning pipeline (replace any references to the old 6-factor approach).
4. ✅ Updated `data_needs.md` listing what still needs to be fetched (likely just T7 phage if Module 01 hasn't grabbed it).
5. ✅ At least 3 passing tests in `processes/tests/`.
6. ✅ `AGENT_REPORT.md` summarizing what was done, any anomalies found, and any blockers.

## Setup

```bash
cd /path/to/agent-00-raw-data       # your worktree
conda activate igem2026              # env from shared/env/environment.yml

# No extra packages needed for this module — biopython + pandas + tqdm
# are already in the conda env. If you need anything else, document why
# in 00_raw_data/requirements.txt.
```

## Step-by-step plan

### Step 1 — Inventory the existing dataset (15 min)

Walk `00_raw_data/phage/` and `00_raw_data/bacteria/`. For each subdirectory:
- Confirm presence of `genome.fna`, `cds.fna`, `proteins.faa` (note any missing — some phages may only have genome.fna).
- Use `Bio.SeqIO` to parse each FASTA, count records, capture sequence lengths.
- Capture file sizes and SHA-256.

Output an in-memory DataFrame; will be saved to MANIFEST.csv in step 4.

### Step 2 — Cross-check against `phage_list.csv` and `bacteria_list.csv` (10 min)

For each accession listed in the CSV:
- Confirm a corresponding directory exists.
- Flag mismatches in an "audit_issues" column (missing dir, wrong record count, sequence length differs from CSV's `length_bp` by > 1%).

Don't fix mismatches yourself — log them.

### Step 3 — Build the verification notebook (30 min)

Create `processes/01_verify_dataset.ipynb` following the bilingual convention from `01_data_ground_truth/processes/fetch_reference_genomes.ipynb`. Cells:

1. **Markdown** — Title + bilingual purpose statement.
2. **Code** — Imports, paths, version printouts.
3. **Markdown** — Methodology (one English paragraph + one Chinese paragraph).
4. **Code** — Walk file tree, build inventory DataFrame.
5. **Code** — FASTA parse + record count.
6. **Code** — SHA-256 computation (use `hashlib`, chunked read for large files).
7. **Code** — Cross-check against `phage_list.csv` / `bacteria_list.csv`.
8. **Code** — Write `MANIFEST.csv` per INTERFACE schema.
9. **Markdown** — Summary table + any anomalies.
10. **Markdown** — Next steps / hand-off notes.

Notebook must produce `MANIFEST.csv` with columns exactly as in INTERFACE §Universal conventions:
`filename, sha256, bytes, n_records, created_utc, source_acc, source_module, notes`

### Step 4 — Update `README.md` (15 min)

The current `README.md` describes the dataset for the old 6-factor pipeline. Rewrite to:
- Describe role in active-learning pipeline (reference dataset + Layer 3 contrastive prior corpus).
- Document the directory layout.
- Cross-reference `INTERFACE.md` for the schema.
- List the canonical reference accessions used by Module 01 (`EU717894.1`, `GCF_000007145.1`, `GCF_000840885.1`).

### Step 5 — Update `data_needs.md` (10 min)

Should be a short doc listing:
- Currently present (777 phages, 34 bacteria, 1 phage reference at full quality)
- Currently missing (fill in based on your audit — likely Xcc full assembly + T7 reference)
- Pointer to `01_data_ground_truth/` for the actual fetching

### Step 6 — Write tests (45 min)

Create `processes/tests/`:
- `test_schema.py` — Open `MANIFEST.csv`, assert columns match INTERFACE.
- `test_smoke.py` — Pick the smallest 3 `genome.fna` files, parse with `SeqIO`, assert ≥1 record each.
- `test_sanity.py` — Assert `phage/EU717894.1/genome.fna` has expected length 44,080 bp ± 100. Assert `phage_list.csv` row count matches number of subdirectories in `phage/`.

Run: `pytest 00_raw_data/processes/tests/ -v`

All 3 tests must pass.

### Step 7 — Commit + report (15 min)

Commit messages format (one per logical chunk):
```
00_raw_data: <verb> <what>

<details>
```

Final commit creates `00_raw_data/AGENT_REPORT.md` with:
- What was done (bullet list)
- Anomalies found (table: file, issue, recommended fix)
- Total runtime
- Any blockers for downstream agents

Push policy: do NOT push to remote. Alex reviews + pushes in the morning.

## References (cite in notebook markdown cells)

- **NCBI Datasets CLI documentation**: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/
- **FASTA format**: https://en.wikipedia.org/wiki/FASTA_format (de facto standard)
- **BioPython `SeqIO`**: Cock, P.J.A. et al. (2009) "Biopython: freely available Python tools for computational molecular biology and bioinformatics." *Bioinformatics* 25:1422.

## Anti-goals — DO NOT DO

- ❌ Re-download all 777 phages — they're already there.
- ❌ Modify any file inside `phage/<acc>/` or `bacteria/<acc>/`.
- ❌ Touch other modules' folders.
- ❌ Push to git remote.
- ❌ Use absolute paths or `~`.
- ❌ Skip the bilingual comment requirement.

## Time budget

~2 hours total. If you exceed 3 hours, log a blocker and stop.

## If you get stuck

Document the blocker in `AGENT_REPORT.md` under a "BLOCKERS" section with:
- What you tried
- The exact error message
- Your guess at the cause
- Whether downstream agents can proceed regardless

Then move on to the next non-blocked task.
