# Live Demo Runbook — iGEM Claremont 2026 Onboarding

**~15 minutes**, runs Modules 03 → 04 → 06 → 07 on screen-share. Pairs with `slides_en.pdf` (the deck pauses around slide 48 for the demo, or run at the end).

> All commands are **copy-paste** from the repo root. Every step has runtime, expected output, and a "if it fails" line. Do a dry-run within 24 h of the meeting.

---

## Before the meeting (5-minute checklist)

```bash
# Anchor the working directory
cd "/Users/alexy/Desktop/Claude Workspace/iGEM_Claremont_2026"

# 1) Right branch — main is empty
git checkout active-learning-pipeline
git pull --rebase

# 2) Right environment
conda activate igem2026
python -c "import torch, numpy, pandas, sklearn; print('ok')"

# 3) Genomes present (~5 MB total)
#    Genomes live in 00_raw_data/phage/<acc>/ and 00_raw_data/bacteria/<acc>/
ls 00_raw_data/phage/EU717894.1/ 2>/dev/null \
   00_raw_data/bacteria/GCF_000007145.1/ 2>/dev/null
# If either is missing:
#   python 00_raw_data/processes/fetch_phages.py   --accession EU717894.1
#   python 00_raw_data/processes/fetch_bacteria.py --accession GCF_000007145.1

# 4) HMM pressed (Module 03) — the .hmm file is fetched, not committed
ls 03_rbp_identification/inputs/phagerbpdetect_data/RBPdetect_phageRBPs.hmm 2>/dev/null
# If missing, run the one-time setup script (clones PhageRBPdetect repo + hmmpress):
#   bash 03_rbp_identification/inputs/setup_inputs.sh

# 5) Cached outputs present (we'll mostly look at these; live runs regenerate them)
ls 06_uncertainty_model/outputs/cycle_0/predictions.csv \
   07_acquisition_function/outputs/cycle_1/recommendations.csv \
   03_rbp_identification/outputs/EU717894.1_rbp_candidates.csv

# 6) Browser tabs open:
#    - docs/planning/PI_briefing_2026-05-11.md
#    - this file (DEMO.md)
#    - pae_heatmap.png   (open via 'open pae_heatmap.png')
```

If any of (1)–(5) fail, **skip live execution and walk through committed outputs only**. The cached `outputs/` files are the source of truth on the slides.

---

## Demo 1 — Module 03 RBP identification (≈2 min)

**What this shows.** PhageRBPdetect finds candidate receptor-binding proteins in phage proteomes via HMM domain search. For phiL7 it returns 3 candidates; the top is `rbp_01` at 712 aa.

```bash
pytest 03_rbp_identification/processes/tests/ -v --tb=line | tail -15
```

| Field | Expected |
|-------|----------|
| Runtime | ~15 s |
| Result | `25+ passed, 2 skipped (HMM binary)` |
| If it fails | Press the HMM (one-liner above) and retry; otherwise skip and open the committed CSV. |

Then show the output:

```bash
head 03_rbp_identification/outputs/EU717894.1_rbp_candidates.csv | column -t -s,
```

**What to say.** "Top row, rank 1: `orf_00001`, 712 aa, HMM score 342 — that's `rbp_01`. This is the tail spike Lee 2009 *searched for and couldn't find by sequence similarity*. HMMs detect proteins too diverged for BLAST."

---

## Demo 2 — Module 04 protein embedding (≈1 min, cached)

**What this shows.** ESM-2 embeddings turn each RBP sequence into a 320-dim vector (8M local proof-of-concept; production is 1280-dim on Laguna). We don't run ESM-2 live — the model download is too slow. Instead we inspect the cached `.npz`, which is what every downstream module actually consumes.

```bash
python -c "
import numpy as np
d = np.load('04_protein_embedding/outputs/embeddings_esm2_t6_8M_phiL7_rbps.npz', allow_pickle=True)
print('files :', d.files)
print('seq_ids:', d['seq_ids'])
print('shape  :', d['array'].shape, d['array'].dtype)
print('lengths:', d['lengths'])
"
```

| Field | Expected |
|-------|----------|
| Runtime | <1 s |
| `files` | `['seq_ids', 'array', 'lengths', 'meta']` |
| `seq_ids` | `['EU717894.1_rbp_01' 'EU717894.1_rbp_02' 'EU717894.1_rbp_03']` |
| `array.shape` | `(3, 320) float32` |
| `lengths` | `[712 918 224]` (matches Module 03's rank-1, 2, 3 candidates) |
| If it fails | The `.npz` is committed — re-clone or `git checkout 04_protein_embedding/outputs/` to restore. |

**What to say.** "Three sequences from Module 03 — `rbp_01` at 712 aa, `rbp_02` at 918 aa, `rbp_03` at 224 aa — go in, three 320-dim vectors come out. ESM-2 was trained on ~65 M protein sequences (Lin 2023). On Laguna we'd upgrade to the 650M parameter model — 1280-dim vectors. The contract is the same."

---

## Demo 3 — Module 06 deep ensemble (≈3 min)

**What this shows.** 5 MLPs train independently in seconds on synthetic embeddings, then emit `predicted_score`, `std`, and `epistemic_std` per (RBP, receptor) pair. This is the input BALD consumes.

```bash
python 06_uncertainty_model/processes/run_cycle0.py 2>&1 | tail -10
```

| Field | Expected |
|-------|----------|
| Runtime | ~3 seconds on CPU (verified 2026-05-17) |
| Result | `=== Cycle 0 training complete. All outputs in .../cycle_0 ===` |
| If it fails | `git checkout -- 06_uncertainty_model/outputs/cycle_0` and demo the committed files instead. |

```bash
column -t -s, < 06_uncertainty_model/outputs/cycle_0/predictions.csv | head -6
```

Then open the calibration plot:

```bash
open 06_uncertainty_model/outputs/cycle_0/calibration.png
```

**What to say.** "Look at `epistemic_std` — that's the *only* column BALD reads. Notice the variance: 0.04 for `rbp_01 × rec_03` (model is confident), 0.22 for `rbp_07 × rec_02` (members disagree). The calibration plot is near-diagonal on synthetic data; the post-Cycle-0 version will be the one we actually trust."

---

## Demo 4 — Module 07 BALD acquisition (≈2 min)

**What this shows.** BALD ranks every unmeasured pair by epistemic uncertainty, takes the top 4, adds 1 random control, writes `recommendations.csv` + provenance.

```bash
python 07_acquisition_function/processes/run_bald.py --cycle 1 --n_bald 4 --n_random 1
```

| Field | Expected |
|-------|----------|
| Runtime | <1 second |
| Result | Writes 4 files under `07_acquisition_function/outputs/cycle_1/` |
| If it fails | Run again with `--seed 42`; if persistent, fall back to the committed CSV. |

```bash
column -t -s, < 07_acquisition_function/outputs/cycle_1/recommendations.csv
```

**What to say.** "rbp_07 × rec_02 is BALD top-1 — that's the variant the model is least confident about, so measuring it shrinks epistemic uncertainty the most. The last row is `random_control` with `bald_rank = 19` — far from the top — that's our retrospective control arm for the AL-vs-random comparison at project end."

Show the provenance:

```bash
cat 07_acquisition_function/outputs/cycle_1/run_meta.json | head -25
```

**What to say.** "Every recommendation is tagged with `repo_commit_sha` and `timestamp_utc` — fully reproducible. The same input → same picks."

---

## Demo 5 — Module 05 Boltz-2 result (read-only, ≈3 min)

**What this shows.** First protein-protein complex prediction for the project. We don't re-run Boltz-2 live (it's a ~15-min Laguna GPU job).

Open the PAE heatmap:

```bash
open pae_heatmap.png
```

Then `affinity.json`:

```bash
cat 05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/affinity.json
```

| Field | Expected |
|-------|----------|
| `interface_ipTM` | 0.365 |
| `chain_A_ptm` | 0.808 |
| `confidence_score` | 0.683 |
| `predicted_dG` | `null` (small-molecule head — protein-protein returns NaN) |

**What to say.**
- *"Low ipTM (0.365) is not failure — it defines the experiment. ELISA + active learning shrinks that uncertainty."*
- *"High chain-A pTM (0.808) means rbp_01 monomer is well-predicted — reliable basis for variant design."*
- *"`predicted_dG` is null because Boltz-2's affinity head was trained on small-molecule × protein only (PubChem / ChEMBL / BindingDB). For RBP × receptor we use **ipTM**, not affinity. This was one of the May 11 audit corrections."*

If you have PyMOL or ChimeraX on the laptop:

```bash
open 05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/boltz_results_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/predictions/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB_model_0.pdb
```

Chain A (rbp_01) and Chain B (TonB) render as ribbons; B-factor column encodes pLDDT.

---

## Demo 6 — Round-trip (optional, ≈3 min)

**What this shows.** Two commands end-to-end: retrain the ensemble, then call BALD. This is the entire dry-lab loop, minus ELISA ingestion.

```bash
python 06_uncertainty_model/processes/run_cycle0.py
python 07_acquisition_function/processes/run_bald.py --cycle 1 --n_bald 4 --n_random 1
head 07_acquisition_function/outputs/cycle_1/recommendations.csv
```

| Field | Expected |
|-------|----------|
| Runtime | ~4 seconds total |
| Result | Fresh `recommendations.csv` identical to the committed one (deterministic seed). |

**What to say.** *"That's the 48-hour SLA in 4 seconds. Module 08 is the part with the 10–14-day SDM cycle behind it."*

---

## After the demo — pointers for the team

| You want to … | Open this |
|---------------|-----------|
| Run your first module locally | `GETTING_STARTED.md` |
| Understand the project end to end | `docs/planning/PI_briefing_2026-05-11.md` |
| Find a particular output file | "Where outputs live" in `guide_en.md` (and the PI briefing) |
| Author new code | `CLAUDE.md` — notebook-first workflow, bilingual comments |
| Run on Laguna | `LAGUNA.md` |
| Read a paper | `docs/reference/papers.md` (annotated, with priority tags) |
| Check what we fixed in the literature audit | `docs/reference/paper_reading_notes.md` |

Questions back to: Alex (CChen29@cmc.edu).
