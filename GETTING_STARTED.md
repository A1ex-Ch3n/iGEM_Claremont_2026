# Getting Started — iGEM Claremont 2026
# 快速入門指南

> One-stop guide: what each module does, which files to open, and how to know when you're done.
> 一站式指南：每個模組做什麼、從哪個檔案開始、怎麼知道跑完了。

---

## First-time setup / 第一次使用

```bash
# 1. Create the conda environment (one-time, ~5 min)
conda env create -f shared/env/environment.yml
conda activate igem2026

# 2. Download the minimum genomes needed for local development
python 00_raw_data/processes/fetch_phages.py --accession EU717894.1      # phiL7
python 00_raw_data/processes/fetch_phages.py --accession NC_001604.1     # T7 (tests)
python 00_raw_data/processes/fetch_bacteria.py --accession GCF_000007145.1  # Xcc

# 3. Launch JupyterLab
jupyter lab
```

> **Why so few genomes?** The full 777 phage + 34 bacteria dataset (~630 MB) is gitignored.
> You only need phiL7 + Xcc for local development. Everything else is downloaded on Laguna
> before batch jobs. See Module 00 Cell 7 for a full stage-by-stage breakdown.

---

## Module map / 模組地圖

```
00_raw_data/          Verify data integrity
01_data_ground_truth/ Known phage-host interaction pairs
02_annotation/        DNA → protein sequences (gene calling)
03_rbp_identification/ Find receptor-binding proteins in phage
04_protein_embedding/ Protein sequences → numbers (ESM-2)
05_structure_prediction/ Predict 3D binding complex + affinity
06_uncertainty_model/ Train model that knows what it doesn't know
07_acquisition_function/ Pick the next experiment (BALD) ✅ Done
08_cycle_data/        Store wet-lab results each cycle ← starts ~June 1
```

---

## Module 00 — Raw Data Verification / 原始數據驗證

**What it does:** Audits all genome files on disk, computes SHA-256 checksums, and writes a
reproducibility manifest. Tells you what's missing and what each pipeline stage needs.

**Open:** `00_raw_data/processes/01_verify_dataset.ipynb`
**Produces:** `00_raw_data/MANIFEST.csv`

**Run tests:**
```bash
pytest 00_raw_data/processes/tests/ -v
```

**Checklist:**
- [ ] Cell 7 shows `[OK] All priority genomes present`
- [ ] `MANIFEST.csv` exists at `00_raw_data/MANIFEST.csv`
- [ ] `EU717894.1` and `GCF_000007145.1` are on disk
- [ ] 15+ tests pass (3 expected failures — GCF/T7 not in original phage_list)

---

## Module 01 — Ground Truth Data / 基準真實數據

**What it does:** Downloads reference genomes (phiL7, Xcc, T7) from NCBI and builds the
phage-host interaction matrix (2,236 pairs: 315 positive + 1,920 negative + 1 ground-truth).

**Before you run:** Download the 3 prerequisite genomes (see top of notebook).

**Open in order:**
1. `01_data_ground_truth/processes/01_fetch_reference_genomes.ipynb`
2. `01_data_ground_truth/processes/02_build_interaction_matrix.ipynb`

**Produces:**
- `01_data_ground_truth/outputs/reference_genomes_index.csv`
- `01_data_ground_truth/outputs/interaction_matrix.csv`
- `01_data_ground_truth/outputs/MANIFEST.csv`

**Run tests:**
```bash
pytest 01_data_ground_truth/processes/tests/ -v
```

**Checklist:**
- [ ] `interaction_matrix.csv` has 2,236 rows
- [ ] phiL7 × Xcc row exists with `label=1`
- [ ] `reference_genomes_index.csv` shows all 3 genomes as `downloaded`
- [ ] 22/22 tests pass

---

## Module 02 — Genome Annotation / 基因組注釋

**What it does:** Calls ORFs (protein-coding genes) from raw DNA.
Uses PHANOTATE for phage (handles overlapping genes) and pyrodigal for bacteria.
**Never swap the tools — wrong results if you do.**

**Open in order:**
1. `02_annotation/processes/01_run_phanotate.ipynb` — phiL7
2. `02_annotation/processes/02_run_prodigal.ipynb` — Xcc

**Produces:**
- `02_annotation/outputs/phage_proteins/EU717894.1.faa` (80 proteins)
- `02_annotation/outputs/phage_orfs/EU717894.1.gff3`
- `02_annotation/outputs/host_proteins/NZ_CP155948.1.faa` (4,344 proteins)
- `02_annotation/outputs/host_orfs/NZ_CP155948.1.gff3`
- `02_annotation/outputs/annotation_summary.csv`

**Run tests:**
```bash
pytest 02_annotation/processes/tests/ -v
```

**Checklist:**
- [ ] `EU717894.1.faa` has 80 records
- [ ] FASTA headers match format: `>EU717894.1_orf_00001 | source=... | length=...`
- [ ] `annotation_summary.csv` has 2 rows (phiL7 + Xcc)
- [ ] 26/26 tests pass

---

## Module 03 — RBP Identification / 受體結合蛋白鑑定

**What it does:** Scans all phage proteins to find receptor-binding proteins (RBPs) —
the "keys" phage uses to attach to host cells. Uses HMM domain search (PhageRBPdetect).

**Prerequisite:** Module 02 output (`EU717894.1.faa`) must exist.

**Open:** `03_rbp_identification/processes/01_run_phagerbpdetect.ipynb`

**Produces:**
- `03_rbp_identification/outputs/EU717894.1_rbp_candidates.csv` (all 58 proteins scored)
- `03_rbp_identification/outputs/EU717894.1_rbp_sequences.faa` (top-3 RBP candidates)
- `03_rbp_identification/outputs/MANIFEST.csv`

**Run tests:**
```bash
pytest 03_rbp_identification/processes/tests/ -v
```

**Checklist:**
- [ ] At least 1 candidate with `combined_score ≥ 0.7`
- [ ] Top candidate (`rbp_01`) is ~712 aa — this is the phiL7 tail spike
- [ ] `EU717894.1_rbp_sequences.faa` has 3 sequences
- [ ] 25+ tests pass (2 expected failures — HMM binary file needs local setup)

> **HMM binary setup** (if tests fail with hmmscan error):
> ```bash
> cd 03_rbp_identification/inputs/phagerbpdetect_data
> hmmpress RBPdetect_phageRBPs.hmm
> ```

---

## Module 04 — Protein Embedding / 蛋白質嵌入

**What it does:** Converts protein sequences into 320-dimensional number vectors using
ESM-2 (a protein language model). These vectors are the input to the ML model.

**Prerequisite (important):** Run with Module 03's **real** RBP sequences, not the mock file.
Check that `03_rbp_identification/outputs/EU717894.1_rbp_sequences.faa` exists first.

**Open in order:**
1. `04_protein_embedding/processes/01_embed_esm2.ipynb` — embed phiL7 RBPs
2. `04_protein_embedding/processes/02_extract_receptors.ipynb` — embed Xcc receptors

**Produces:**
- `04_protein_embedding/outputs/embeddings_esm2_t6_8M_phiL7_rbps.npz` — shape (3, 320)
- `04_protein_embedding/outputs/embeddings_esm2_t6_8M_xcc_receptors.npz` — shape (4, 320)
- `04_protein_embedding/outputs/embedding_index.csv`

**Run tests:**
```bash
pytest 04_protein_embedding/processes/tests/ -v
```

**Checklist:**
- [ ] RBP npz: `seq_ids` matches `EU717894.1_rbp_01/02/03`, shape `(3, 320)`
- [ ] Receptor npz: 4 entries for TonB, ExbB, ExbD1, ExbD2
- [ ] No NaN values in either array
- [ ] 17/17 tests pass

> **Laguna upgrade:** Swap `esm2_t6_8M_UR50D` → `esm2_t33_650M_UR50D` for production
> (320-dim → 1280-dim). See `04_protein_embedding/AGENT_REPORT.md` for sbatch command.

---

## Module 05 — Structure Prediction / 結構預測

**What it does:** Predicts the 3D structure of phage RBP bound to Xcc receptor using
Boltz-2, and produces a binding confidence score (ipTM). This becomes the synthetic
prior for Module 06 before real ELISA data arrives.

**Important:** Full structure prediction requires Laguna GPU (~15 min/pair on A100).
The CPU smoke test times out after 30 min — use Laguna for real results.

**Open in order:**
1. `05_structure_prediction/processes/01_run_boltz2.ipynb` — Boltz-2 predict
2. `05_structure_prediction/processes/02_run_af3.ipynb` — AF3 docs (needs Google approval)

**Produces:**
- `05_structure_prediction/outputs/affinity_priors.csv`
- `05_structure_prediction/outputs/boltz2/<rbp>__<receptor>/` (PDB + confidence JSON)

**Run tests:**
```bash
pytest 05_structure_prediction/processes/tests/ -v
```

**Checklist:**
- [ ] Input FASTAs exist: `inputs/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB.fasta`
- [ ] `affinity_priors.csv` has correct columns (even if `run_success=False` for now)
- [ ] 28/28 tests pass (1 expected skip — PDB sanity test needs completed GPU run)
- [ ] Laguna sbatch command ready (see `05_structure_prediction/AGENT_REPORT.md`)

> **Boltz-2 note:** `predicted_dG_kcal_mol` is always NaN for protein-protein pairs.
> Use `confidence` (ipTM) as the binding signal instead.

---

## Module 06 — Uncertainty Model / 不確定性模型

**What it does:** Trains a 5-member deep ensemble that predicts binding scores AND
quantifies how confident it is. Low confidence = experiment is worth running.
Currently uses synthetic data; swaps to real ELISA data ~June 1.

**Prerequisite:** Module 04 embeddings must be re-generated with real RBP sequences
(not the mock sequences from the overnight build).

**Open in order:**
1. `06_uncertainty_model/processes/01_train_deep_ensemble_synthetic.ipynb`
2. `06_uncertainty_model/processes/02_calibration_check.ipynb`

**Produces:**
- `06_uncertainty_model/outputs/cycle_0/ensemble_member_{0..4}.pt`
- `06_uncertainty_model/outputs/cycle_0/predictions.csv`
- `06_uncertainty_model/outputs/cycle_0/calibration.png`
- `06_uncertainty_model/outputs/cycle_0/model_meta.json`

**Run tests:**
```bash
pytest 06_uncertainty_model/processes/tests/ -v
```

**Checklist:**
- [ ] 5 `.pt` model files saved in `outputs/cycle_0/`
- [ ] `predictions.csv` has `std > 0` and `epistemic_std > 0` for all rows
- [ ] `calibration.png` generated (ECE reported in `model_meta.json`)
- [ ] 9/9 tests pass
- [ ] `model_meta.json` has `data_source: synthetic_fallback_random` ← will change to `elisa_cycle_0` when real data arrives

> **When ELISA data arrives (~June 1):** See `06_uncertainty_model/AGENT_REPORT.md`
> "When ELISA Arrives" section — only 2 lines of notebook change.

---

## Module 07 — BALD Acquisition Function / BALD 採集函數

**What it does:** Reads Module 06 predictions, ranks all unmeasured (variant, receptor) pairs
by epistemic uncertainty (BALD score = Std of ensemble member means), and outputs the
next batch for wet lab. Includes a mandatory random-control pick for retrospective comparison.

**Prerequisite:** Module 06 `predictions.csv` must exist (with `epistemic_std` column).

**Run (Cycle 0 — no ELISA data yet):**
```bash
python 07_acquisition_function/processes/run_bald.py --cycle 1 --n_bald 4 --n_random 1
```

**Run (Cycle N — with ELISA data):**
```bash
python 07_acquisition_function/processes/run_bald.py \
  --cycle N \
  --measured_csv 08_cycle_data/outputs/cycle_<N-1>/elisa_processed.csv
```

**Produces:**
- `07_acquisition_function/outputs/cycle_<N>/recommendations.csv` — variants to test next
- `07_acquisition_function/outputs/cycle_<N>/safe_pick_backup.csv` — manual fallback
- `07_acquisition_function/outputs/cycle_<N>/random_replay.csv` — retrospective baseline
- `07_acquisition_function/outputs/cycle_<N>/run_meta.json`

**Run tests:**
```bash
pytest 07_acquisition_function/processes/tests/ -v
```

**Checklist:**
- [ ] 18/18 tests pass
- [ ] `recommendations.csv` has `n_bald + n_random` rows with `selection_reason` column
- [ ] Top pick has highest `bald_score` (= `epistemic_std`) in the pool
- [ ] `run_meta.json` contains `repo_commit_sha` and `timestamp_utc`

---

## Laguna HPC / 高性能計算

For batch jobs (all 777 phages, ESM-2 650M, Boltz-2 GPU):

```bash
ssh <username>@laguna.hpc.institution.edu
cd $SCRATCH && git clone https://github.com/A1ex-Ch3n/iGEM_Claremont_2026.git
git checkout active-learning-pipeline   # ← all pipeline code is on this branch
conda env create -f shared/env/environment.yml
conda activate igem2026

# Download all genomes
python 00_raw_data/processes/fetch_phages.py
python 00_raw_data/processes/fetch_bacteria.py

# Submit batch jobs (see LAGUNA.md for full templates)
sbatch scripts/slurm/embed_esm2_650m.slurm
sbatch scripts/slurm/boltz2_screen.slurm
```

Full SLURM templates with correct partition/account flags: see **`LAGUNA.md`**.

---

## Key reference files / 重要參考文件

| File | Purpose |
|------|---------|
| `CLAUDE.md` | Project conventions, team roles, data flow rules |
| `INTERFACE.md` | Data contracts between modules (column schemas, file formats) |
| `LAGUNA.md` | HPC setup, SLURM job templates |
| `docs/planning/PI_briefing_2026-05-11.md` | Current status briefing — all outputs, work log, attachments (bilingual) |
| `docs/planning/iGEM_2026_Project_Plan.md` | Full project plan (English, PI-facing) |
| `docs/planning/iGEM_2026_项目大纲_中文版.md` | Full project plan (Chinese, team-facing) |
| `docs/reference/papers.md` | Per-module reading guide (19 papers, annotation protocol) |

---

## When something breaks / 出問題時

1. Check the module's `AGENT_REPORT.md` — known issues and workarounds are documented there.
2. Re-run the module's tests: `pytest <module>/processes/tests/ -v`
3. Check `INTERFACE.md` to confirm input/output file formats match expectations.
4. Ask Claude Code — context for this entire build is preserved in the conversation history.
