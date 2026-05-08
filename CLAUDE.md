# iGEM Claremont 2026 — Active-Learning Phage Engineering for *Xanthomonas* Biocontrol

## Project overview

We are building a **closed-loop active-learning pipeline** that integrates a machine-learning model of phage receptor-binding-protein (RBP) — host receptor interactions with iterative wet-lab validation. The system uses Bayesian Optimal Experimental Design (BOED) to let the model nominate the next-most-informative experiment, addressing the central pain point of phage-host prediction: data scarcity.

**Reference scaffold (dry lab):** *Xanthomonas campestris* pv. *campestris* (Xcc) ATCC 33913 + phage phiL7 (NCBI EU717894) — phiL7's receptor system (TonB-ExbB-ExbD1D2) is experimentally identified (Wang et al. 2003).

**Wet lab strategy:** Self-isolate Xanthomonas + lytic phage from California symptomatic crops (per PI consultation, 2026-05-07; bypasses USDA APHIS PPQ-526 permit).

**Targeting:** iGEM Best Agriculture Project · Best Model · Best Composite Part.

**Core engineer:** Alex Chen. **Wet lab:** Sarah, Olivia, Weitao. **Other contributors:** Ryan, Leah. **PI:** Prof. J. Cesar Ignacio-Espinoza. **Faculty advisor:** Prof. Ran Libeskind-Hadas.

Project planning documents:
- `docs/iGEM_2026_Project_Plan.md` — full PI-facing English plan
- `docs/iGEM_2026_项目大纲_中文版.md` — Chinese version with deep dry-lab module mechanics + 6-layer data-scarcity strategy
- `docs/iGEM_2026_Project_Flow.png` — full pipeline flow chart
- `docs/glossary.md` — technical vocabulary reference

---

## Pipeline at a glance

| Step | Folder | Description | Owner(s) | Status |
|------|--------|-------------|----------|--------|
| 0 | `00_raw_data/` | Raw genome dump (777 phages + 34 bacteria) | — | ✅ Done |
| 1 | `01_data_ground_truth/` | NCBI download scripts + reference interaction matrix | Sarah, Weitao | ✅ Partial |
| 2 | `02_annotation/` | PHANOTATE (phages) + Prodigal (hosts) + pharokka | Weitao, Olivia | 🟡 Partial |
| 3 | `03_rbp_identification/` | PhageRBPdetect (HMM + XGBoost) on phage proteomes | TBD | ⏳ Sprint deliverable |
| 4 | `04_protein_embedding/` | ESM-2 650M / 3B embeddings (+ optional PLM-interact transfer) | ESM-experienced teammate | 🔄 In progress |
| 5 | `05_structure_prediction/` | AlphaFold 3 (structures) + Boltz-2 (zero-shot affinity prior) | Laguna-trained teammate | ⏳ Sprint deliverable |
| 6 | `06_uncertainty_model/` | Deep ensemble (5 MLPs) over ESM-2 → calibrated UQ | Alex | ⏳ Sprint deliverable |
| 7 | `07_acquisition_function/` | BALD acquisition + 48-hour cycle infrastructure | Alex | ⏳ Sprint deliverable |
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

Each step folder has its own `README.md` describing inputs, processes, outputs, and current status.

---

## Data flow contract

1. **`inputs/`** — read-only references to upstream step `outputs/` (or external seeds). Never write generated data here.
2. **`processes/`** — the only place with code. All scripts read from `inputs/` pointers or upstream `outputs/`; write to their step's `outputs/`.
3. **`outputs/`** — canonical products consumed by the next step. Large artifact trees are gitignored; commit `MANIFEST.csv` files instead.

**Canonical files:**
- Reference phage genome: `00_raw_data/phage_genomes/<accession>.fna` (phiL7 = EU717894)
- Reference host genome: `00_raw_data/bacteria/<accession>.fna` (Xcc = AE008922)
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

- **No hard-coded absolute paths.** Use `pathlib.Path(__file__).parents[N]` to anchor at the repo root or step folder.
- **Output filenames:** use descriptive prefixes (`cycle_<N>_*`, `<phage_acc>_*`, etc.) so files self-document their lineage.
- **Code comments:** English only.
- **Docs:** bilingual EN + ZH (Simplified or Traditional Chinese acceptable depending on intended audience).
- **Large artifacts:** gitignored. Add accession/filename to a `MANIFEST.csv` in the same `outputs/` folder so the set is reproducible.
- **Tool split:** PHANOTATE for phage ORF calling; Prodigal / pyrodigal for bacterial hosts. Never swap.
- **Cycle versioning:** every model checkpoint, prediction CSV, and ELISA dataset is tagged with `cycle_<N>` and tracked via MLflow.

---

## Quick start

```bash
# 1. Create environment
conda env create -f shared/env/environment.yml
conda activate igem2026

# 2. Step 1 — build interaction matrix (positive pairs)  [legacy reference data]
python 01_data_ground_truth/processes/fetch_positive_pairs.py

# 3. Step 2 — annotate reference genomes (phiL7 + Xcc)
python 02_annotation/processes/phage_phanotate/batch_phanotate.py
python 02_annotation/processes/host_prodigal/batch_prodigal.py

# 4. Step 3 — identify RBP candidates with PhageRBPdetect  [TBD: Alex sprint]
python 03_rbp_identification/processes/run_phagerbpdetect.py \
    --input 02_annotation/outputs/phage_proteins/EU717894.faa \
    --out 03_rbp_identification/outputs/

# 5. Step 4 — ESM-2 embeddings  [in development]
python 04_protein_embedding/processes/embed_esm2.py \
    --input 03_rbp_identification/outputs/EU717894_rbp_sequences.faa \
    --model esm2_t33_650M_UR50D \
    --out 04_protein_embedding/outputs/

# 6. Step 5 — Structure + affinity prior
python 05_structure_prediction/processes/run_af3.py    # Laguna HPC
python 05_structure_prediction/processes/run_boltz2.py # Laguna HPC

# 7. Step 6 — Train deep ensemble  [TBD: Alex sprint]
python 06_uncertainty_model/processes/train_ensemble.py --cycle 0

# 8. Step 7 — BALD recommendations for next wet-lab batch  [TBD: Alex sprint]
python 07_acquisition_function/processes/run_bald.py --cycle 0
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
