# AGENT TODO — Module 03 / `03_rbp_identification/`

## Read first (mandatory)

1. `/INTERFACE.md` §Module 03 — your output spec.
2. `/CLAUDE.md` — notebook-first, bilingual.
3. `processes/README.md` — already drafted in last commit, review and extend.
4. Invoke `superpowers:test-driven-development` and `superpowers:verification-before-completion`.

## Your scope

Identify which proteins in a phage proteome are receptor-binding proteins (RBPs / tail spikes / tail fibers). You read protein FASTAs from `02_annotation/outputs/phage_proteins/`, you produce ranked candidate CSVs and FASTA subsets in `03_rbp_identification/outputs/`.

The downstream consumers are Module 04 (embeddings) and Module 05 (structure prediction), both of which need RBP sequences only — not the whole phage proteome.

## Goal (definition of done)

By morning, `03_rbp_identification/` contains:

1. ✅ `processes/01_run_phagerbpdetect.ipynb` — runs PhageRBPdetect on phage proteomes.
2. ✅ `outputs/EU717894.1_rbp_candidates.csv` (per INTERFACE schema).
3. ✅ `outputs/EU717894.1_rbp_sequences.faa` (top-K candidates as FASTA).
4. ✅ `outputs/all_rbp_candidates.csv` (unified across whatever phages you process).
5. ✅ `outputs/MANIFEST.csv`.
6. ✅ ≥3 passing tests.
7. ✅ `AGENT_REPORT.md`.

## Setup

```bash
cd /path/to/agent-03-rbp-identification
conda activate igem2026

# PhageRBPdetect — author's GitHub repo
pip install git+https://github.com/dimiboeckaerts/PhageRBPdetection.git
# OR clone + dev install if pip fails:
git clone https://github.com/dimiboeckaerts/PhageRBPdetection.git /tmp/phagerbpdetect
pip install -e /tmp/phagerbpdetect

# Dependencies pulled in:
# - hmmer (for Pfam HMM scan) — verify with: which hmmscan
# - xgboost
# - fair-esm (for ESM embedding the ML track)

# If hmmscan missing:
conda install -c bioconda hmmer

# Pfam-A HMM database (~ 1.5 GB)
mkdir -p inputs/pfam
cd inputs/pfam
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
cd -
```

If Pfam download fails or is too slow tonight, document and fall back to ML-only track.

## Step-by-step plan

### Step 1 — Verify PhageRBPdetect installation (30 min)

Quick smoke test before doing real work:
- Import the package; print version.
- Run on the package's bundled example FASTA (if any) to confirm pipeline works.
- Confirm Pfam HMM DB is loaded.
- Confirm ESM-2 weights download (PhageRBPdetect uses ESM-2 8M internally — small enough for CPU).

If anything fails, log to SETUP.md and SETUP_BLOCKERS.md, then proceed with what works.

### Step 2 — Build the notebook (`01_run_phagerbpdetect.ipynb`) (90 min)

Cells:

1. **Markdown** — Title + bilingual purpose. Cite Boeckaerts 2022.
2. **Code** — Imports, paths, version printouts.
3. **Markdown** — Method explanation (bilingual): the dual-track design (HMM + ML) and why both are needed. Cite Latka 2021 mBio for RBP modularity context.
4. **Code** — Helper: `run_hmm_track(faa_path, pfam_db_path) → DataFrame`. Returns ORFs that hit RBP-related Pfams. Use the curated Pfam list from PhageRBPdetect's source code.
5. **Code** — Helper: `run_ml_track(faa_path, pfam_misses) → DataFrame`. Embeds Pfam-missed ORFs with ESM-2, runs XGBoost classifier, returns probabilities.
6. **Code** — Helper: `combine_tracks(hmm_df, ml_df) → DataFrame`. Merges per INTERFACE schema (combined_score = 1.0 if HMM hit; else ml_score).
7. **Markdown** — Sample-then-batch.
8. **Code** — Sample run on phiL7 (`02_annotation/outputs/phage_proteins/EU717894.1.faa`).
9. **Code** — Sanity assertion: phiL7 should yield 1-3 RBP candidates with combined_score > 0.7. The known tail spike gp25 (~412 aa) should be the top candidate.
10. **Code** — Write `EU717894.1_rbp_candidates.csv` and `EU717894.1_rbp_sequences.faa` (top-5 by rank).
11. **Code** — Optional batch: process additional phages from `02_annotation/outputs/phage_proteins/`. Skip if Module 02 only annotated 1.
12. **Code** — Write `all_rbp_candidates.csv` (unified) and update MANIFEST.

### Step 3 — Top-K RBP FASTA subset (30 min)

After ranking, write `<phage_acc>_rbp_sequences.faa` containing top-K sequences (default K=5). Header format per INTERFACE:
```
>EU717894.1_rbp_01 | source_orf=EU717894.1_orf_00031 | length=412 | combined_score=0.94 | evidence=both
```

These are the proteins Module 04 will embed and Module 05 will fold.

### Step 4 — Tests (45 min)

`processes/tests/`:
- `test_schema.py` — Open `EU717894.1_rbp_candidates.csv`, assert columns match INTERFACE. Open `EU717894.1_rbp_sequences.faa`, assert headers parse correctly.
- `test_smoke.py` — Construct a tiny synthetic FASTA with 3 fake "proteins" (random ~100 aa each), run the combined pipeline, assert it returns a DataFrame without crashing. Don't assert on scores (they'll be junk for random sequences).
- `test_sanity.py` — `EU717894.1_rbp_candidates.csv` has ≥1 row with `combined_score ≥ 0.7`. The top candidate's `length_aa` should be ∈ [200, 600] (tail spike size range).

Run: `pytest 03_rbp_identification/processes/tests/ -v`

### Step 5 — Commit + report

Commits:
- `03_rbp_identification: install PhageRBPdetect + smoke test`
- `03_rbp_identification: notebook + phiL7 RBP identification`
- `03_rbp_identification: outputs in INTERFACE format`
- `03_rbp_identification: tests`
- `03_rbp_identification: AGENT_REPORT`

## References (cite in notebook markdown cells)

- **PhageRBPdetect (the tool)**: Boeckaerts, D. et al. (2022) "Identification of phage receptor-binding protein sequences with HMMs and XGBoost." *Viruses* 14(6):1329. DOI: 10.3390/v14061329. **THIS IS THE PRIMARY METHODOLOGY REFERENCE.**
- **Pfam database**: Mistry, J. et al. (2021) "Pfam: The protein families database in 2021." *Nucleic Acids Research* 49(D1):D412-D419.
- **HMMER (used internally)**: Eddy, S.R. (2011) "Accelerated profile HMM searches." *PLOS Computational Biology* 7(10):e1002195.
- **ESM-2 (used in ML track)**: Lin, Z. et al. (2023) "Evolutionary-scale prediction of atomic-level protein structure." *Science* 379(6637):1123-1130.
- **XGBoost (used in ML track)**: Chen, T. & Guestrin, C. (2016) "XGBoost: A scalable tree boosting system." *KDD '16*.
- **RBP biology context**:
  - Latka, A. et al. (2021) "Engineering the modular receptor-binding proteins of *Klebsiella* phages switches their capsule serotype specificity." *mBio* 12(3):e00455-21. DOI: 10.1128/mBio.00455-21
  - Yehl, K. et al. (2019) "Engineering Phage Host-Range and Suppressing Bacterial Resistance through Phage Tail Fiber Mutagenesis." *Cell* 179(2):459-469.
  - Yang, J. et al. (2024) "SpikeHunter: identifying tail spike proteins from phage genomes." *GigaScience*.
- **phiL7 RBP biology**: Lee, C.N. et al. (2009) *Appl Environ Microbiol* 75:7828. The tail spike gene is in their Table 1 (gp25 region).

## Anti-goals

- ❌ Don't try to identify RBPs from bacterial proteomes (RBPs are phage-only).
- ❌ Don't trust ML-track scores < 0.5 — these are noise.
- ❌ Don't write the unified `all_rbp_candidates.csv` if you only processed 1 phage; just have the per-phage file plus the FASTA subset.
- ❌ Don't push to remote.

## Time budget

~4 hours. PhageRBPdetect install + Pfam DB download is the main risk. If install drags on past hour 2, escalate.

## If stuck

- Pfam DB download too slow → use the smaller `Pfam-A.full` selection PhageRBPdetect bundles, or skip HMM track entirely (ml_score-only); document trade-off.
- ESM weights download blocked → set `TORCH_HOME` to a writable path; pre-download with `transformers` library if `fair-esm` blocked.
- Combined pipeline gives 0 candidates on phiL7 → likely a header parsing bug in your `02_annotation` input. Inspect the FASTA headers, fix parsing in your helper, retest.
