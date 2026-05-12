# iGEM Claremont 2026 — Active-Learning Phage Engineering for *Xanthomonas* Biocontrol

A closed-loop machine learning pipeline that uses Bayesian Active Learning (BALD) to engineer phage receptor-binding proteins (RBPs) with improved host-range specificity against *Xanthomonas campestris* pv. *campestris*.

**Targeting:** iGEM Best Agriculture Project · Best Model · Best Composite Part

---

## Getting started

All pipeline code and outputs live on the `active-learning-pipeline` branch:

```bash
git clone https://github.com/A1ex-Ch3n/iGEM_Claremont_2026.git
cd iGEM_Claremont_2026
git checkout active-learning-pipeline
```

Then open **[`GETTING_STARTED.md`](GETTING_STARTED.md)** for a step-by-step module guide.

---

## Pipeline at a glance

| Module | What it does | Status |
|--------|-------------|--------|
| 00 Raw Data | 777 phage + 34 bacteria genome library | ✅ |
| 01 Ground Truth | 2,236 phage-host interaction pairs | ✅ |
| 02 Annotation | DNA → protein sequences (PHANOTATE + Prodigal) | ✅ |
| 03 RBP ID | Find receptor-binding proteins (PhageRBPdetect HMM) | ✅ |
| 04 Embedding | Protein sequences → ESM-2 vectors | ✅ |
| 05 Structure | 3D binding complex (Boltz-2 on Laguna GPU) | ✅ |
| 06 Ensemble | 5-member deep ensemble with epistemic uncertainty | ✅ |
| 07 BALD | Acquisition function — picks the next experiment | ✅ |
| 08 Cycle Data | Wet-lab ELISA ingestion + model retraining | ⏳ Cycle 0 ~June 1 |

---

## Key documents

| File | Purpose |
|------|---------|
| [`GETTING_STARTED.md`](GETTING_STARTED.md) | Step-by-step module guide with checklists |
| [`CLAUDE.md`](CLAUDE.md) | Project conventions, team roles, data contracts |
| [`LAGUNA.md`](LAGUNA.md) | CARC Laguna HPC setup + SLURM templates |
| [`docs/planning/PI_briefing_2026-05-11.md`](docs/planning/PI_briefing_2026-05-11.md) | Current status briefing (bilingual EN/ZH) — all outputs listed |
| [`docs/planning/iGEM_2026_Project_Plan.md`](docs/planning/iGEM_2026_Project_Plan.md) | Full project plan (English) |
| [`docs/planning/iGEM_2026_项目大纲_中文版.md`](docs/planning/iGEM_2026_项目大纲_中文版.md) | Full project plan (Chinese) |
| [`docs/reference/papers.md`](docs/reference/papers.md) | Per-module reading guide (19 papers) |

---

## Team

**Core engineer (dry lab):** Alex Chen  
**Wet lab:** Sarah, Olivia, Weitao, Carol  
**Contributors:** Ryan, Leah  
**PI:** Prof. J. Cesar Ignacio-Espinoza  
**Faculty advisor:** Prof. Ran Libeskind-Hadas
