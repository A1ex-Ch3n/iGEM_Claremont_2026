# iGEM Claremont 2026 — Active-Learning Phage Engineering for *Xanthomonas* Biocontrol

## Project overview

We are building a **closed-loop active-learning pipeline** that integrates a machine-learning model of phage receptor-binding-protein (RBP) — host receptor interactions with iterative wet-lab validation. The system uses Bayesian Optimal Experimental Design (BOED) to let the model nominate the next-most-informative experiment, addressing the central pain point of phage-host prediction: data scarcity.

**Reference scaffold (dry lab):** *Xanthomonas campestris* pv. *campestris* (Xcc) ATCC 33913 + phage phiL7 (NCBI EU717894) — phiL7's receptor system (TonB-ExbB-ExbD1 is experimentally identified as essential (ExbD2 NOT required per Hung et al. 2003, BBRC 302:878–884).

**Wet lab strategy:** Self-isolate Xanthomonas + lytic phage from California symptomatic crops (per PI consultation, 2026-05-07; bypasses USDA APHIS PPQ-526 permit).

**Targeting:** iGEM Best Agriculture Project · Best Model · Best Composite Part.

**Core engineer:** Alex Chen. **Wet lab:** Sarah, Olivia, Weitao. **Other contributors:** Ryan, Leah. **PI:** Prof. J. Cesar Ignacio-Espinoza. **Faculty advisor:** Prof. Ran Libeskind-Hadas.

Project planning documents (see `docs/README.md` for the full index):
- `docs/planning/iGEM_2026_Project_Plan.md` — full PI-facing English plan
- `docs/planning/iGEM_2026_项目大纲_中文版.md` — Chinese version with deep dry-lab module mechanics + 6-layer data-scarcity strategy
- `docs/planning/iGEM_2026_Project_Flow.png` — full pipeline flow chart
- `docs/reference/glossary.md` — technical vocabulary reference
- `docs/protocols/` — wet lab Benchling SOPs (cultivation, transformation, plaque, infection curves, lysate amp)

---

## Pipeline at a glance

| Step | Folder | Description | Owner(s) | Status |
|------|--------|-------------|----------|--------|
| 0 | `00_raw_data/` | Raw genome dump (phages + bacteria, gitignored binaries) | — | ✅ Done |
| 1 | `01_data_ground_truth/` | NCBI download scripts + reference interaction matrix | Sarah, Weitao | ✅ Partial |
| 2 | `02_annotation/` | PHANOTATE (phages) + Prodigal (hosts) + pharokka | Weitao, Olivia | 🟡 Partial |
| 3 | `03_rbp_identification/` | PhageRBPdetect (HMM + XGBoost) on phage proteomes | Alex | ✅ 3 RBP candidates (rbp_01 = 712aa tail spike) |
| 4 | `04_protein_embedding/` | ESM-2 650M / 3B embeddings (+ optional PLM-interact transfer) | ESM-experienced teammate | 🔄 In progress |
| 5 | `05_structure_prediction/` | AlphaFold 3 (structures) + Boltz-2 (ipTM structural confidence) | Alex | ✅ rbp_01 712aa × TonB done; ipTM=0.365 |
| 6 | `06_uncertainty_model/` | Deep ensemble (5 MLPs) over ESM-2 → calibrated UQ | Alex | ✅ 5-member MLP, calibrated; predictions.csv includes epistemic_std |
| 7 | `07_acquisition_function/` | BALD acquisition + 48-hour cycle infrastructure | Alex | ✅ bald.py + run_bald.py; 18 tests pass; first cycle run done |
| 8 | `08_cycle_data/` | Per-cycle wet-lab data + model checkpoints + retrospective replay | All | ⏳ Cycle 0 starts ~6/1 |

---

## Repository layout

```
iGEM_Claremont_2026/
├── CLAUDE.md                         ← this file
├── docs/                             ← all planning documents (PI plan, brief, glossary, flow chart, Benchling protocols)
├── Proposal_General.pdf              ← original 2026 iGEM proposal
├── phage_host_prediction_project_brief.md  ← pivot starting point (May 2026)
│
├── 00_raw_data/                      ← raw NCBI dump: 777 phage + 34 bacteria genomes
├── 01_data_ground_truth/             ← NCBI download scripts + canonical interaction matrix
│   ├── inputs/                       ← manual seed lists
│   ├── processes/                    ← fetch_positive_pairs.py, fetch_negative_pairs.py, download_genomes.py
│   └── outputs/                      ← interaction matrix CSVs + FASTA pools
├── 02_annotation/                    ← PHANOTATE + Prodigal + pharokka
│   ├── processes/phage_phanotate/    ← PHANOTATE binary + batch_phanotate.py
│   ├── processes/host_prodigal/      ← batch_prodigal.py
│   └── outputs/                      ← per-accession .faa files
├── 03_rbp_identification/            ← PhageRBPdetect outputs
├── 04_protein_embedding/             ← ESM-2 embeddings (active development)
├── 05_structure_prediction/          ← AlphaFold 3 + Boltz-2 outputs
├── 06_uncertainty_model/             ← deep ensemble checkpoints + calibration
├── 07_acquisition_function/          ← BALD recommendations + cycle infra
├── 08_cycle_data/                    ← per-cycle ELISA data + model snapshots
│
├── shared/                           ← cross-step utilities, conda env
│   └── env/environment.yml
└── archive/                          ← superseded artifacts
    ├── 2026-05-pivot/                ← pre-pivot folders (6-factor weighting, KNN baseline, fastISM, digital phagogram, weitao/, TASK/)
    ├── legacy_master_pipeline.py
    └── tasks/                        ← old planning docs (e.g. 0428 6-factor regression spec)
```

`docs/` is organized into 4 sub-folders (see `docs/README.md` for full index):

```
docs/
├── README.md                ← navigation index
├── planning/                ← project plans, pivot summary, flow chart
├── reference/               ← glossary, tool manuals
├── protocols/               ← Benchling wet lab SOPs (5 PDFs)
└── archive/                 ← pre-pivot artifacts (workflow chart, integrated guide)
```

Each step folder has its own `README.md` describing inputs, processes, outputs, and current status.

---

## Data flow contract

1. **`inputs/`** — read-only references to upstream step `outputs/` (or external seeds). Never write generated data here.
2. **`processes/`** — the only place with code. All scripts read from `inputs/` pointers or upstream `outputs/`; write to their step's `outputs/`.
3. **`outputs/`** — canonical products consumed by the next step. Large artifact trees are gitignored; commit `MANIFEST.csv` files instead.

**Canonical files:**
- Reference phage genome: `00_raw_data/phage_genomes/<accession>.fna` (phiL7 = EU717894)
- Reference host genome: `00_raw_data/bacteria/<accession>.fna` (Xcc ATCC 33913 = GCF_000007145.1)
- Per-cycle ELISA data: `08_cycle_data/outputs/cycle_<N>/elisa_processed.csv`
- Per-cycle BALD recommendations: `07_acquisition_function/outputs/cycle_<N>/recommendations.csv`

---

## Team ownership map

| Step | Module | Owner | Key deliverable |
|------|--------|-------|----------------|
| 1 — Data | NCBI download + interaction matrix | Sarah, Weitao | `01_data_ground_truth/outputs/` |
| 2 — Annotation (phage) | PHANOTATE | Weitao | `02_annotation/outputs/phage_proteins/<acc>.faa` |
| 2 — Annotation (host) | Prodigal | Olivia | `02_annotation/outputs/host_proteins/<acc>.faa` |
| 3 — RBP ID | PhageRBPdetect | Alex | `03_rbp_identification/outputs/<phage>_rbp_candidates.csv` |
| 4 — Embedding | ESM-2 (650M / 3B on Laguna) | ESM-experienced teammate | `04_protein_embedding/outputs/*.npy` |
| 5 — Structure | AlphaFold 3 + Boltz-2 | Laguna-trained teammate | `05_structure_prediction/outputs/` |
| 6 — UQ Model | Deep ensemble | Alex | `06_uncertainty_model/outputs/cycle_<N>/predictions.csv` |
| 7 — Acquisition | BALD + cycle infra | Alex | `07_acquisition_function/outputs/cycle_<N+1>/recommendations.csv` |
| 8 — Wet/dry sync | Cycle data ingestion | All | `08_cycle_data/outputs/cycle_<N>/` |
| Wet — Strain isolation | Self-isolate from CA brassicas | Sarah, Olivia, Weitao | In-house Xanthomonas isolates + WGS |
| Wet — Phage isolation | Co-isolation + plaque purification | Sarah, Olivia, Weitao | In-house lytic phage + WGS |
| Wet — Cloning + ELISA | Gibson + BL21 + Ni-NTA + ELISA | Sarah, Olivia, Weitao | Per-cycle binding curves |
| Wet — Receptor KO | pK18mobsacB markerless deletion | Sarah, Olivia, Weitao | ΔtonB / ΔexbB / ΔexbD strains |

---

## Conventions

### Code authoring (notebook-first workflow / 笔记本优先开发)

All new code is written as **Jupyter notebooks first** for exploratory development, then frozen as `.py` modules once the logic is stable. Comments are **bilingual** (English + 中文) so the entire team can read them.

- **Notebook location:** `<step_folder>/processes/<NN>_<short_name>.ipynb` (numbered sequentially within each module's `processes/`).
- **Bilingual comment style:**
  - Markdown cells: one English paragraph followed by a Chinese paragraph (or vice versa) — see `01_data_ground_truth/processes/fetch_reference_genomes.ipynb` for the canonical example.
  - Inline code comments: `# English / 中文` on the same line for short comments; separate lines for longer ones.
- **Freeze trigger:** Once a notebook (a) runs end-to-end without manual intervention and (b) its outputs match a verification cell, freeze it as `<NN>_<short_name>.py` in the same folder. Mark the original notebook in its filename: `<NN>_<short_name>__frozen.ipynb` and add a top-cell pointer to the `.py`.
- **Module 07 exception:** The cycle-orchestration code in `07_acquisition_function/` is a production pipeline, not exploratory analysis — write directly as `.py` from day 1.
- **Notebook git hygiene:** Install `nbstripout` (`pip install nbstripout && nbstripout --install`) so cell outputs don't pollute git diffs. Already added to `shared/env/environment.yml`.

### General

- **No hard-coded absolute paths.** Use `pathlib.Path(__file__).parents[N]` (in `.py`) or `Path.cwd().resolve().parents[N]` (in notebooks) to anchor at the repo root or step folder.
- **Output filenames:** use descriptive prefixes (`cycle_<N>_*`, `<phage_acc>_*`, etc.) so files self-document their lineage.
- **Docs:** bilingual EN + ZH (Simplified or Traditional Chinese acceptable depending on intended audience).
- **Large artifacts:** gitignored. Add accession/filename to a `MANIFEST.csv` in the same `outputs/` folder so the set is reproducible.
- **Tool split:** PHANOTATE for phage ORF calling; Prodigal / pyrodigal for bacterial hosts. Never swap.
- **Cycle versioning:** every model checkpoint, prediction CSV, and ELISA dataset is tagged with `cycle_<N>` and tracked via MLflow.
- **Genome binaries gitignored:** `00_raw_data/phage/*/` and `00_raw_data/bacteria/*/` are not in git (630MB). Re-download: `python 00_raw_data/processes/fetch_phages.py && python 00_raw_data/processes/fetch_bacteria.py`
- **Boltz-2 reference pair:** Use `EU717894.1_rbp_01` (712aa tail spike, not P25 85aa) × `GCF_000007145.1_tonB`. Script: `scripts/boltz2_phiL7_tonB.slurm`. HPC config: see `LAGUNA.md`.
- **Boltz-2 affinity head = small molecule only:** For protein–protein pairs (e.g., RBP × TonB receptor), the affinity head outputs NaN. Use **ipTM** (0–1 interface confidence score) as a structural confidence proxy — it is NOT a quantitative binding affinity. Do not describe Boltz-2 as providing a "zero-shot affinity prior" for protein–protein interactions.
- **HPC jobs (Boltz-2, ESM-2 3B):** See `LAGUNA.md` for full cluster setup, CUDA compatibility, and SLURM templates.
- **Merging agent branches:** Use `git checkout <branch> -- <dir>` (additive) not `git merge` — merge creates large deletion diffs when branch histories diverge.
- **setuptools platform split:** macOS: `pip install "setuptools<70"` for phanotate (pkg_resources). Linux/Laguna: keep default (≥75) for trifast. Do NOT put the macOS pin in `environment.yml`.

---

## Quick start

```bash
# 1. Create environment + install notebook hygiene tool
# 创建环境 + 安装 notebook diff 净化工具
conda env create -f shared/env/environment.yml
conda activate igem2026
pip install nbstripout && nbstripout --install   # one-time per clone

# 2. Launch JupyterLab and open the reference fetch notebook
# 启动 JupyterLab 并打开参考基因组抓取 notebook
jupyter lab
# Then open: 01_data_ground_truth/processes/fetch_reference_genomes.ipynb
# 打开: 01_data_ground_truth/processes/fetch_reference_genomes.ipynb

# Pipeline progression (each step has a notebook in its processes/ folder):
# Pipeline 流程 (每个 step 在自己的 processes/ 里有对应 notebook):
#   01 → fetch_reference_genomes.ipynb              (DONE — example notebook)
#   02 → phage_phanotate/  +  host_prodigal/        (existing scripts; will rewrite as notebooks)
#   03 → run_phagerbpdetect.ipynb                   (DONE — 3 RBP candidates)
#   04 → embed_esm2.ipynb                           (8M done locally; 650M pending Laguna)
#   05 → run_af3.ipynb  +  run_boltz2.ipynb         (DONE — rbp_01 712aa × TonB, ipTM=0.365)
#   06 → train_ensemble.ipynb                       (DONE — 5-member MLP, epistemic_std exported)
#   07 → run_bald.py                                (DONE — bald.py + run_bald.py + 18 tests)
```

---

## Key files to read first (for new contributors)

| File | Why |
|------|-----|
| `docs/iGEM_2026_Project_Plan.md` | Full project plan in English (PI version) |
| `docs/iGEM_2026_项目大纲_中文版.md` | Chinese version with extra detail on dry-lab modules + data-scarcity strategy |
| `docs/iGEM_2026_Project_Flow.png` | Visual flow chart of the entire pipeline |
| `docs/glossary.md` | Technical vocabulary (ELISA, BALD, ESM-2, etc.) |
| `docs/project_pivot_summary_for_team.md` | Chinese team-sync explaining the May 2026 pivot |
| `Proposal_General.pdf` | Original iGEM 2026 proposal |

---

## Changelog

| Date | Who | What |
|------|-----|------|
| 2026-05-10 | Alex (Claude session) | 7 overnight agent branches cherry-picked into `active-learning-pipeline`; cleaned up worktrees |
| 2026-05-10 | Alex (Claude session) | Fixed 6 `00_raw_data` issues: contaminated phiL7 proteins.faa, invalid bacteria accessions (KY000037→GCF_000092025.1, PY746849→GCF_000006765.1), gitignored 630MB genome binaries |
| 2026-05-10 | Alex (Claude session) | Set up CARC Laguna: OnDemand, Code Server, `boltz2` conda env (torch 2.5.1+cu121 — 8 attempts to resolve CUDA/trifast/optree version conflicts) |
| 2026-05-10 | Alex (Claude session) | First Boltz-2 3D structure (phiL7 P25 85aa × Xcc TonB, ipTM=0.345) — wrong protein; re-run with rbp_01 712aa pending |
| 2026-05-10 | Alex (Claude session) | Generated `GETTING_STARTED.md`, `docs/alex_project_guide.md`, `docs/pipeline_build_report.md` |
| 2026-05-11 | Alex (Claude session) | CLAUDE.md: fixed Xcc accession (AE008922→GCF_000007145.1), updated pipeline table, added HPC/git/setuptools gotchas |
| 2026-05-11 | Alex (Claude session) | Paper reading pass: verified 19 papers, fixed ExbD2/Greenman/Boltz-2/Hie2024 errors; see docs/reference/paper_reading_notes.md |
| 2026-05-11 | Alex (Claude session) | Correct Boltz-2 run (job 59986): rbp_01 712aa × TonB, ipTM=0.365, chain_A_ptm=0.808 |
| 2026-05-12 | Alex (Claude session) | Module 07 BALD complete: bald.py + run_bald.py + 18 tests; Module 06 patched to export epistemic_std; full pipeline (00–07) done |
| 2026-05-12 | Alex (Claude session) | PI_briefing updated to reflect full completion; progress_report merged into PI_briefing |
