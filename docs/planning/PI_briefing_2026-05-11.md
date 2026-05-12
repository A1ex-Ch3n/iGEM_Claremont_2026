# iGEM Claremont 2026 — PI Briefing & Project Status
**Date:** 2026-05-12 (created 2026-05-11, updated 2026-05-12) | **By:** Alex Chen | **For:** Prof. J. Cesar Ignacio-Espinoza

---

> **To view all outputs referenced in this document, you must be on the correct branch.**
>
> ```bash
> git clone https://github.com/A1ex-Ch3n/iGEM_Claremont_2026.git
> cd iGEM_Claremont_2026
> git checkout active-learning-pipeline
> ```
>
> All pipeline outputs (PDB structure, predictions, BALD recommendations, calibration plots) are committed on `active-learning-pipeline`. The `main` branch does not contain the latest dry-lab work.

---

# ENGLISH VERSION

## TL;DR

The **entire computational pipeline (Modules 00–07) is fully built and tested.** We have the first Boltz-2 3D structure of phiL7 rbp_01 × Xcc TonB (ipTM = 0.365), a calibrated 5-member deep ensemble, and a working BALD acquisition function that produces ranked variant recommendations in under 1 second. The system will begin active learning as soon as wet lab delivers the first ELISA measurements (~June 1). **No critical dry-lab work remains before May 17 wet lab launch.**

Five factual corrections were made after a full literature audit (19 papers). All documented in `docs/reference/paper_reading_notes.md`.

---

## Pipeline Status

| Module | Status | Key fact |
|--------|--------|---------|
| 00 Raw Data | ✅ | 777 phage + 34 bacteria genomes |
| 01 Interaction Matrix | ✅ | 2,236 pairs; 1 confirmed (phiL7 × Xcc) |
| 02 Annotation | ✅ | phiL7: 80 ORFs; Xcc: 4,344 ORFs |
| 03 RBP Identification | ✅ | 3 candidates; rbp_01 (712 aa) primary |
| 04 Embeddings | ✅ | ESM-2 vectors ready (650M version pending Laguna) |
| 05 Structure | ✅ | rbp_01 × TonB PDB + ipTM = 0.365 |
| 06 Deep Ensemble | ✅ | 5-member MLP, calibrated, outputs epistemic_std |
| 07 BALD | ✅ | Scorer + CLI orchestrator; 18 tests pass; first cycle run |
| 08 Cycle Data | ⏳ Cycle 0 starts ~6/1 | Waiting for ELISA measurements |

---

## Computational Results

### 1. phiL7 rbp_01 × Xcc TonB — First 3D Complex Prediction (Boltz-2, job 59986)

| Metric | Value | Interpretation |
|--------|-------|----------------|
| `interface_ipTM` | **0.365** | Low — model uncertain about HOW they dock. Expected for a novel system with no PDB template. ELISA will resolve this. |
| `chain A ptm` | **0.808** | High — rbp_01 monomer structure is well-predicted. Reliable basis for variant design. |
| `confidence_score` | **0.683** | Overall complex quality — moderate. |

The low ipTM is not a failure — it defines the experiment. That uncertainty IS the question the ELISA + active learning loop is designed to answer. The high chain A ptm (0.808) means rbp_01 is structurally well-constrained — good news for expression and stability.

**File:** `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/` (view in PyMOL or ChimeraX)

### 2. RBP Candidates from phiL7

| Candidate | Length | Domain hit | Priority |
|-----------|--------|-----------|---------|
| rbp_01 | **712 aa** | Tail_spike_N + C-terminal binding | **Primary — Cycle 0 target** |
| rbp_02 | 918 aa | Collagen-like repeat | Backup / chimera source |
| rbp_03 | 224 aa | Short C-terminal domain | Low priority |

Note: rbp_01 is computationally identified. Lee et al. 2009 (phiL7 genome) suggests p20 (1105 aa) for host range but does not name a tail spike. Our identification is based on the Tail_spike_N HMM hit from PhageRBPdetect.

### 3. BALD Acquisition — First Cycle Recommendations (synthetic prior)

Before any ELISA data arrives, BALD ranks by epistemic uncertainty (Var_k[μ_k], variance of ensemble member means). On the synthetic Cycle 0 predictions:

| Priority | Variant | Receptor | BALD score | Rationale |
|----------|---------|----------|------------|-----------|
| 1 (BALD top) | rbp_07 | rec_02 | 0.218 | Highest ensemble disagreement |
| 2 | rbp_03 | rec_01 | 0.197 | — |
| 3 | rbp_05 | rec_02 | 0.197 | — |
| 4 | rbp_01 | rec_02 | 0.190 | — |
| 5 (random) | rbp_03 | rec_03 | 0.127 | Control arm for retrospective comparison |

These are synthetic placeholders — actual Cycle 1 recommendations will be based on real rbp_01 variants × TonB after Cycle 0 ELISA data is in.

---

## What Each Module Does

| Module | Role | Key output |
|--------|------|-----------|
| 00 Raw Data | Genome library for training data | 777 phage + 34 bacteria genomes |
| 01 Ground Truth | Labeled phage-host pairs (the only real labels in the system) | interaction_matrix.csv: 2,236 pairs |
| 02 Annotation | Translates raw DNA → protein sequences | phiL7: 80 ORFs; Xcc: 4,344 ORFs |
| 03 RBP ID | Finds the "key" protein determining host range | rbp_01 (712 aa): primary tail spike candidate |
| 04 Embedding | Converts protein sequences → neural network input vectors | 1280-dim ESM-2 embeddings |
| 05 Structure | 3D structure of RBP × receptor complex; structural confidence prior | PDB + ipTM score |
| 06 Ensemble | 5 independent neural networks; their disagreement = what the model doesn't know | (predicted score, epistemic_std) per variant |
| 07 BALD | Uses uncertainty to pick the next experiment: select the variant the model is most uncertain about | Ranked variant CSV for wet lab |
| 08 Cycle Data | Ingests ELISA results, triggers retraining, closes the loop | Per-cycle data + model checkpoints |

---

## Validation Strategy

Three practical tiers, all using the same workflow (ELISA-first, then add layers):

| Tier | What you do | Scientific story | Feasibility |
|------|------------|-----------------|-------------|
| 1 — ELISA only | Binding curves on WT Xcc, 4–6 variants/cycle | "We found variants that bind better." | ✅ Guaranteed — this is Cycle 0–2 |
| 2 — + Plaque assay | +2 days: spot assay on WT confirms lytic infection | "Binding → infection confirmed." | ✅ Near-zero incremental cost |
| 3 — + ΔtonB/ΔexbB/ΔexbD1 | Markerless deletions (pK18mobsacB, 4–6 weeks) | "Binding is causally receptor-specific." This is a **mechanistic paper-quality claim**. | 🟡 If knockouts start May 17 |

**Built-in negative control:** Hung 2003 confirms ΔexbD2 does NOT affect phiL7 infection — generating this strain validates the entire knockout system for free. If ΔexbD2 still allows infection and ΔtonB does not, the receptor specificity claim is experimentally supported.

**Recommendation:** Commit to Tier 3 now. Start pK18mobsacB knockouts May 17. If knockouts succeed by early July, the result is a five-star iGEM story and paper-grade mechanistic claim.

---

## Corrected Biological Facts (Literature Audit)

After reading Hung et al. 2003 (BBRC 302:878–884, PMID 12646254) directly:

**ExbD2 is NOT required for phiL7 infection.** Only TonB, ExbB, ExbD1 are essential. The exbD2 knockout (CH620) retained full sensitivity. All planning documents have been corrected.

Four additional corrections from the 19-paper audit:
1. Boltz-2 affinity head = small molecule–protein only. Protein-protein pairs output NaN. We use ipTM as structural confidence proxy — NOT a binding affinity.
2. Greenman 2025 journal = *PLoS Comput Biol* 21(1):e1012639 (not NAR Genomics); conclusion: "no single best UQ method."
3. Hie 2024 uses ESM-1b/1v (not ESM-2); ~20 variants per antibody (not ~50).
4. Lee 2009 explicitly searched for and was **unable to find** a homolog of OP1's tail fiber (ORF25) in phiL7. This is stronger than "didn't name a tail spike" — they actively looked and found nothing by sequence homology. Our rbp_01 identification via Tail_spike_N HMM is consistent: HMMs detect structurally similar proteins too diverged for BLAST to find.

---

## Key Decisions Needed from PI

### 1. Receptor knockout system — pK18mobsacB or CRISPRi? 🔴 Before May 17
- **pK18mobsacB (our plan):** Markerless deletion, proven in *Xanthomonas* (Hung 2003). 4–6 weeks, permanent clean knockouts. Plasmid on Addgene (#87097, ~$95).
- **CRISPRi:** Knockdown, reversible, faster (~2 weeks). Less established in Xcc; guides need design.
- **Recommendation:** pK18mobsacB for tonB/exbB/exbD1 in parallel. Start May 17.

### 2. AlphaFold 3 model weights application 🔴 This week
- AF3 gives higher-quality structures than Boltz-2, including trimer predictions.
- Requires Google form (academic use). 1–7 day approval.
- Should Alex or PI submit? (Institutional email preferred.)

### 3. Validation tier commitment 🟡 Before Cycle 0 (June 1)
See Validation Strategy section above. Recommend Tier 3 if knockouts start May 17.

### 4. Phage enrichment source 🟡 Before June 1
Co-isolation requires enrichment substrate broader than crop tissue alone. What agricultural runoff or sewage sources does the lab have access to in LA County?

### 5. Manuscript ambition 🟡 Before Cycle 0
Aim for concurrent submission (*Bioinformatics* / *NAR Genomics & Bioinformatics*) alongside iGEM wiki? This sets data documentation granularity from the first ELISA measurement.

---

## Discussion-Ready Outputs (for PI meeting)

| # | Item | File | Key talking point |
|---|------|------|-------------------|
| 1 | rbp_01 × TonB 3D structure | `05_structure_prediction/outputs/boltz2/.../predictions/*.pdb` | View in PyMOL/ChimeraX. Does predicted interface match known TonB-ExbB face from *E. coli* literature? |
| 2 | Structural confidence numbers | `...affinity.json`: ipTM=0.365, chain_A_ptm=0.808 | Low ipTM = experiment needed, not model failure. High chain A = reliable for variant design. |
| 3 | RBP candidate list | `03_rbp_identification/outputs/EU717894.1_rbp_candidates.csv` | rbp_01 primary. Should we express rbp_02 in Cycle 0 as chimera source? |
| 4 | Interaction matrix | `01_data_ground_truth/outputs/interaction_matrix.csv` | 2,236 pairs; 1 confirmed. Thin baseline compensated by Boltz-2 structural prior. |
| 5 | ESM-2 embeddings | `04_protein_embedding/outputs/*.npz` | 8M proof-of-concept; 650M needed on Laguna before Cycle 0. |
| 6 | BALD recommendations | `07_acquisition_function/outputs/cycle_1/recommendations.csv` | 4 BALD picks + 1 random control. Synthetic now; real picks after Cycle 0 ELISA. |
| 7 | Validation strategy comparison | This document, "Validation Strategy" section | Recommend Tier 3. Hung 2003 used the same approach in Xcc — proven feasible. |
| 8 | Literature audit corrections | `docs/reference/paper_reading_notes.md` | 5 errors corrected. Most important: ExbD2 is NOT required → free negative control. |

---

## Dry-Lab → Wet-Lab Handoff Protocol

Each cycle has two handoffs:

**Dry lab → Wet lab (48-hour SLA after ELISA data arrives):**
| File | Content |
|------|---------|
| `recommendations.csv` | Ranked variants: ID, mutation, BALD score, predicted EC50 ± uncertainty |
| `primer_sequences.txt` | NEB Q5-compatible primers for SDM variants |
| `uncertainty_bands.png` | Calibration plot: predicted vs measured from previous cycle |
| `safe_pick_backup.csv` | Pre-selected fallback list for PI/team if 48-h SLA is missed |

**Wet lab → Dry lab (end of each cycle):**
| File | Required columns |
|------|-----------------|
| `elisa_processed.csv` | variant_id, receptor_id, ec50_nM, hill_slope, r2, plate_id, date |
| `plaque_results.csv` | variant_id, strain_id, pfu_per_ml, date |
| `qc_report.md` | SDS-PAGE image path, concentration, expression notes |

Minimum for retraining: ≥3 valid EC50 measurements (R² > 0.9) per new variant. Failed variants marked `ec50_nM = NaN` with `failed_reason` — the model handles missing data.

**Cycle 0 exception:** No ELISA yet — variant design is structure-based (Boltz-2 interface + expert picks) + gene synthesis order (not SDM).

---

## How We Got Here — Work Log (May 7–12)

A detailed record of every action taken and why, for full reproducibility and PI context.

### May 7 — Project Pivot

Decided to abandon the original 6-factor biophysical scoring pipeline (GC content, pI, tail fiber phylogeny, etc.) in favor of a **closed-loop active learning system**. The core argument: the 6-factor approach can only predict whether a phage *might* infect a host based on sequence similarity; it cannot learn from wet-lab feedback or improve with each cycle. A BALD-driven deep ensemble, by contrast, directly models uncertainty and shrinks it with each ELISA measurement.

### May 7–8 — Overnight Parallel Build

Seven AI agents worked in parallel overnight (May 7→8), each building one pipeline module from scratch (Modules 00–06). By morning: working code, test suites, and READMEs for all 7 modules.

### May 10 — Post-Build Fixes + Laguna Setup + First Boltz-2 Run (wrong protein)

Post-build audit found 6 data issues requiring manual correction:
1. phiL7 `proteins.faa` in `00_raw_data/` had contaminated headers — re-downloaded raw genome from NCBI and regenerated protein sequences via PHANOTATE.
2. Two invalid bacteria accessions in the genome manifest (`KY000037` → `GCF_000092025.1`; `PY746849` → `GCF_000006765.1`) — corrected via NCBI Assembly search.
3. 630 MB of genome binaries committed to git — moved to `.gitignore`; re-downloadable via `fetch_phages.py` + `fetch_bacteria.py`.
4. Set up CARC Laguna HPC: OnDemand portal, Code Server, conda `boltz2` environment. Required 8 debugging iterations to resolve CUDA/trifast/optree version conflicts (final solution: torch 2.5.1+cu121, NVIDIA L40S driver 550.90.07 = CUDA 12.4; setuptools<70 on macOS only).
5. First Boltz-2 GPU run on Laguna (job 59949, 47 seconds) — **wrong protein**: FASTA named `rbp_01` but contained P25 (85 aa), not rbp_01 (712 aa). ipTM = 0.345 result discarded.
6. Diagnosed cache bug: Boltz-2 `--override` flag only re-runs predictions, NOT the `processed/` cache. Fix: `rm -rf` the entire output directory before each run.

### May 11 — Literature Audit + Correct Boltz-2 Run

Read all 19 core papers and cross-checked every quantitative claim in the project documents. Found and corrected 5 factual errors (see "Corrected Biological Facts" section above). The most consequential: Wang 2003 (the citation supporting ExbD2 as essential) does not exist — it is a hallucinated citation. The real paper is Hung 2003, which shows ExbD2 is NOT essential.

Fixed the FASTA file by overwriting Chain A with the real rbp_01 712aa sequence from `03_rbp_identification/outputs/EU717894.1_rbp_sequences.faa`. Cleared the full Boltz-2 output directory and re-ran on Laguna (job 59986).

Result: **ipTM = 0.365, chain A ptm = 0.808.** The PAE matrix (1316×1316) and per-residue pLDDT are saved as `.npz` files.

Module 04 re-run locally with the real 3 RBP sequences (rbp_01, rbp_02, rbp_03) — removed the 2 mock placeholder sequences (rbp_04, rbp_05).

### May 12 — Module 07 BALD + Full Pipeline Completion

Identified that Module 06's `predictions.csv` only exported `std` (total uncertainty = epistemic + aleatoric), but BALD requires the epistemic component alone (Var_k[μ_k]). Patched `run_cycle0.py` to also export `epistemic_std` via `ens.predict(return_epistemic=True)`. Re-ran Module 06 (~25 seconds).

Built Module 07:
- `bald.py`: BALD score = epistemic_std; `select_batch()` with `selection_reason` audit trail, random control arm, random replay, safe-pick backup
- `run_bald.py`: CLI 48-h SLA orchestrator; auto-detects predictions path; handles Cycle 0 (no measured pairs) and Cycle 1+ (exclusion filter); saves all outputs + `run_meta.json`
- 18 tests covering BALD math (red-green: highest epistemic ranked first AND lowest not in top-K), measured-pair exclusion, Cycle 0 empty path, schema, fallback warning

First live BALD run on synthetic data: top pick is `rbp_07 × rec_02` (BALD=0.218). All 5 outputs committed.

All planning documents updated to reflect full pipeline completion. Two commits pushed to `active-learning-pipeline` (`6422a6b`, `c541726`).

---

## Key Outputs & Attachments

All files below exist on branch `active-learning-pipeline`. All paths relative to repo root.

> **To view the 3D structure:** open the `.pdb` file in PyMOL (free download: pymol.org) or ChimeraX (free: rbvi.ucsf.edu). The two chains will appear as ribbons; Chain A = rbp_01 (purple/blue), Chain B = TonB (orange).

---

### Module 01 — Interaction Matrix (Training Labels)

| File | Description |
|------|-------------|
| `01_data_ground_truth/outputs/interaction_matrix.csv` | **Main output.** 2,236 phage-host pairs. 315 positive (literature-curated), 1,921 negative. Columns: phage_acc, host_acc, label, confidence, source, notes. |
| `01_data_ground_truth/outputs/negative_samples/negative_data_combined.csv` | Full negative set with method labels (cross-genus + pathovar inference). |

---

### Module 03 — RBP Candidates

| File | Description |
|------|-------------|
| `03_rbp_identification/outputs/EU717894.1_rbp_candidates.csv` | **3 RBP candidates** from phiL7's 80 ORFs. Columns: orf_id, length_aa, hmm_score, hmm_match_pfam, combined_score, rank. rbp_01 (orf_00001, 712 aa, HMM score 342) is the primary target. |

---

### Module 04 — Protein Embeddings

| File | Description |
|------|-------------|
| `04_protein_embedding/outputs/embeddings_esm2_t6_8M_phiL7_rbps.npz` | 3 × 320-dim ESM-2 (8M) vectors for rbp_01, rbp_02, rbp_03. Load: `np.load(path)['array']`. Production version (1280-dim, 650M) pending Laguna. |
| `04_protein_embedding/outputs/embeddings_esm2_t6_8M_xcc_receptors.npz` | ESM-2 embeddings for Xcc receptor proteins (TonB, ExbB, ExbD1). |

---

### Module 05 — 3D Structure (Boltz-2, job 59986) ⭐ Primary visual output

> **Note on path depth:** Boltz-2 creates a deeply nested output directory. All key files are inside:
> `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/`
> The sub-path `boltz_results_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/predictions/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/` is Boltz-2's internal naming convention. Full paths are given below.

| File (full path from repo root) | Description |
|------|-------------|
| `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/boltz_results_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/predictions/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB_model_0.pdb` | **3D atomic structure.** Chain A = rbp_01 (712 aa), Chain B = TonB (604 aa). Open in PyMOL or ChimeraX. |
| `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/boltz_results_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/predictions/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/pae_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB_model_0.npz` | **PAE matrix (heatmap data).** 1316×1316 float32. Entry [i,j] = predicted aligned error (Å) between residue i and j. Low = high confidence. Interface block (rows 712–1315 × cols 0–711) = rbp_01 vs TonB. Load: `np.load(path)['pae']`. |
| `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/boltz_results_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/predictions/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/plddt_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB_model_0.npz` | Per-residue pLDDT. 1316 values, range 0.30–0.98, mean 0.76. Load: `np.load(path)['plddt']`. |
| `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/boltz_results_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/predictions/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/confidence_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB_model_0.json` | Per-chain confidence. Key: `iptm=0.365`, `chains_ptm={'0': 0.808, '1': 0.494}`. |
| `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/affinity.json` | Summary scores. `interface_ipTM=0.365`, `chain_A_ptm=0.808`, `confidence_score=0.683`. `predicted_dG=null` — affinity head is small-molecule only. |

**How to generate the PAE heatmap image (run from repo root):**
```python
import numpy as np, matplotlib.pyplot as plt
BASE = "05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB"
PRED = f"{BASE}/boltz_results_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/predictions/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB"
pae = np.load(f"{PRED}/pae_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB_model_0.npz")['pae']
plt.figure(figsize=(8, 8))
plt.imshow(pae, cmap='bwr_r', vmin=0, vmax=30)
plt.axhline(712, color='k', lw=1); plt.axvline(712, color='k', lw=1)
plt.colorbar(label='PAE (Å)'); plt.title('rbp_01 × TonB — Predicted Aligned Error')
plt.xlabel('Residue (Chain A: 0–711 = rbp_01 / Chain B: 712–1315 = TonB)')
plt.savefig('pae_heatmap.png', dpi=150); print("Saved pae_heatmap.png")
```

---

### Module 06 — Deep Ensemble (Cycle 0, Synthetic)

| File | Description |
|------|-------------|
| `06_uncertainty_model/outputs/cycle_0/calibration.png` | **Calibration plot.** Predicted uncertainty interval vs fraction of true labels captured. Well-calibrated = points on the diagonal. Generated on synthetic data; meaningful version will come after Cycle 0 ELISA. |
| `06_uncertainty_model/outputs/cycle_0/predictions.csv` | 80 RBP×receptor pair predictions. Columns: rbp_id, receptor_id, predicted_score, std, **epistemic_std**, lower_95, upper_95, model_version. |
| `06_uncertainty_model/outputs/cycle_0/training_log.json` | Per-member training curves (train NLL and val NLL per epoch for all 5 members). |
| `06_uncertainty_model/outputs/cycle_0/model_meta.json` | Architecture (hidden_dim=256, n_members=5, input_dim=640), training parameters, calibration ECE. |

---

### Module 07 — BALD Acquisition (Cycle 1 Recommendations)

| File | Description |
|------|-------------|
| `07_acquisition_function/outputs/cycle_1/recommendations.csv` | **5 variants recommended for next wet lab cycle.** 4 BALD-optimal picks + 1 random control. Columns: rbp_id, receptor_id, predicted_score, epistemic_std, bald_score, bald_rank, selection_reason. |
| `07_acquisition_function/outputs/cycle_1/safe_pick_backup.csv` | Top-10 BALD picks for manual PI/team selection if 48-h SLA is missed. |
| `07_acquisition_function/outputs/cycle_1/random_replay.csv` | What pure-random selection would have picked. Saved for retrospective learning-curve comparison at project end. |
| `07_acquisition_function/outputs/cycle_1/run_meta.json` | Provenance: cycle, seed, candidate pool size, BALD score stats, top picks, git SHA, timestamp. |

---

### Archive — Pre-Pivot Baseline (Historical Context Only)

These files document the pre-May-2026 approach. Kept for reference; superseded by the active-learning pipeline.

| File | Description |
|------|-------------|
| `archive/2026-05-pivot/05_predictive_modeling/outputs/baseline_taxonomy_knn/heatmap_combined.png` | Taxonomy-KNN baseline prediction heatmap (pre-pivot). Shows the limitation: predictions based on phylogenetic distance only, no per-protein specificity. |
| `archive/2026-05-pivot/05_predictive_modeling/outputs/baseline_taxonomy_knn/heatmap_original.png` | Ground-truth interaction matrix (same scale as predicted). |
| `archive/2026-05-pivot/03_feature_weighting/outputs/per_factor/factor2/f02_acidity_vs_y_scatter.png` | pI vs host-range scatter plot (pre-pivot factor analysis). |

---

## Timeline

```
2026-05-07   Project pivot decided (6-factor → active learning)
2026-05-07 → 2026-05-08   Overnight parallel build (7 agents, Modules 00–06)
2026-05-10   Post-build fixes; Laguna HPC setup; first Boltz-2 run (wrong protein)
2026-05-11   Literature audit (19 papers, 5 corrections); correct Boltz-2 run ✅
2026-05-12   Module 07 BALD complete; full pipeline (00–07) committed and pushed ✅

2026-05-17   Wet lab launches:
             • Brassica sampling (LA County)
             • Receptor knockout plasmid construction (pK18mobsacB)
             • Cycle 0 variants ordered (gene synthesis, IDT/Twist)

2026-05-17 → 2026-06-01
             Strain + phage isolation; Cycle 0 variants expressed in BL21

2026-06-01 → 2026-06-14
             Cycle 0: ELISA optimization + first binding measurements

2026-06-14 → 2026-06-28
             Cycle 1: Ensemble retrained; BALD recommends next variants (SDM)

2026-06-28 → 2026-07-12
             Cycle 2: Round 2 + receptor knockouts complete + causal validation
```

---

*Full pipeline: `docs/planning/iGEM_2026_Project_Plan.md` (EN) | `docs/planning/iGEM_2026_项目大纲_中文版.md` (ZH)*
*Boltz-2 structure: `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/`*
*Literature audit: `docs/reference/paper_reading_notes.md`*

---
---

# 中文版（PI 简报）

**日期：** 2026-05-12（2026-05-11 撰写，2026-05-12 更新）| **撰写：** Alex Chen | **呈送：** Prof. J. Cesar Ignacio-Espinoza

---

## 三句话总结

**整个计算 pipeline（Module 00–07）已全部构建完成并通过测试。** 我们已得到 phiL7 rbp_01 × Xcc TonB 复合体的首个 Boltz-2 3D 结构（ipTM = 0.365），一个经过校准的 5 成员深度集成模型，以及能在 1 秒内输出 variant 排名的 BALD 采集函数。系统已准备好在 wet lab 交出第一批 ELISA 数据（~6/1）后即刻启动主动学习循环。**5/17 wet lab 启动前，干实验室已无关键待办事项。**

完整阅读 19 篇核心论文后，对规划文档进行了 5 处事实性修正，全部记录于 `docs/reference/paper_reading_notes.md`。

---

## Pipeline 状态

| 模块 | 状态 | 关键事实 |
|------|------|---------|
| 00 原始数据 | ✅ | 777 phage + 34 bacteria 基因组 |
| 01 相互作用矩阵 | ✅ | 2,236 对；1 个已确认（phiL7 × Xcc）|
| 02 注释 | ✅ | phiL7: 80 ORF；Xcc: 4,344 ORF |
| 03 RBP 识别 | ✅ | 3 个候选；rbp_01（712 aa）为主 |
| 04 Embedding | ✅ | ESM-2 向量就绪（650M 版本待 Laguna）|
| 05 结构预测 | ✅ | rbp_01 × TonB PDB + ipTM = 0.365 |
| 06 深度集成 | ✅ | 5 成员 MLP，已校准，输出 epistemic_std |
| 07 BALD | ✅ | 采集函数 + CLI 流程脚本；18 个测试通过；首轮跑通 |
| 08 循环数据 | ⏳ Cycle 0 ~6/1 | 等待 ELISA 测量结果 |

---

## 计算结果

### 1. phiL7 rbp_01 × Xcc TonB — 首个 3D 复合体预测（Boltz-2，job 59986）

| 指标 | 数值 | 解读 |
|------|------|------|
| `interface_ipTM` | **0.365** | 低——模型对两个蛋白如何对接不确定。对于无 PDB 模板的全新系统，这是预期结果。ELISA 数据将填补这个空白。 |
| `chain A ptm` | **0.808** | 高——rbp_01 单体结构预测可信。是 variant design 的可靠基础。 |
| `confidence_score` | **0.683** | 整体复合体质量——中等。 |

低 ipTM 不是失败——它定义了实验目标。这种不确定性正是 ELISA + 主动学习循环要解答的问题。高 chain A ptm（0.808）说明 rbp_01 是结构约束良好的蛋白——对表达和稳定性是好消息。

**文件：** `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/`（用 PyMOL 或 ChimeraX 打开）

### 2. phiL7 的 RBP 候选

| 候选 | 长度 | 结构域 | 优先级 |
|------|------|--------|--------|
| rbp_01 | **712 aa** | Tail_spike_N + C 端结合域 | **主要候选——Cycle 0 目标** |
| rbp_02 | 918 aa | Collagen-like repeat | 备选 / Chimera 来源 |
| rbp_03 | 224 aa | 短 C 端结构域 | 低优先级 |

注意：rbp_01 是计算识别的。Lee et al. 2009 建议 p20（1105 aa）可能参与宿主范围决定，但未明确命名任何蛋白为 tail spike。我们的识别依据是 Tail_spike_N HMM。

### 3. BALD 采集函数——第一轮推荐（合成先验）

当前基于合成 Cycle 0 预测（真实 ELISA 数据到位后将替换）：

| 优先级 | Variant | Receptor | BALD 分数 | 说明 |
|-------|---------|----------|----------|-----|
| 1（BALD 最优） | rbp_07 | rec_02 | 0.218 | 集成成员分歧最大 |
| 2 | rbp_03 | rec_01 | 0.197 | — |
| 3 | rbp_05 | rec_02 | 0.197 | — |
| 4 | rbp_01 | rec_02 | 0.190 | — |
| 5（随机对照） | rbp_03 | rec_03 | 0.127 | 用于事后学习曲线对比 |

---

## 每个 Module 的作用

| 模块 | 作用 | 关键输出 |
|------|------|---------|
| 00 原始数据 | 基因组库——训练数据来源 | 777 phage + 34 bacteria 基因组 |
| 01 Ground Truth | 已知 phage-host 相互作用标签（系统中唯一的真实标签） | interaction_matrix.csv：2,236 对 |
| 02 注释 | 原始 DNA → 蛋白质序列 | phiL7：80 ORF；Xcc：4,344 ORF |
| 03 RBP 识别 | 找到决定宿主范围的「钥匙蛋白」 | rbp_01（712 aa）：主要尾刺蛋白候选 |
| 04 Embedding | 蛋白质序列 → 神经网络输入向量 | 1280 维 ESM-2 embedding |
| 05 结构预测 | RBP × 受体复合体 3D 结构；提供结构先验 | PDB + ipTM 分数 |
| 06 深度集成 | 5 个独立神经网络；分歧 = 不确定性 = 模型不知道的地方 | 每个 variant 的（预测分数，epistemic_std）|
| 07 BALD | 用不确定性选下一个实验：选模型最不确定的 variant | 给 wet lab 的 ranked variant CSV |
| 08 循环数据 | 摄入 ELISA 数据，触发重训练，闭合循环 | 每轮 ELISA 数据 + 模型 checkpoint |

---

## 验证策略

三个实际可行的层次，用同一个工作流（ELISA 为基础，逐层叠加）：

| 层次 | 具体操作 | 科学故事 | 可行性 |
|------|---------|---------|--------|
| 第一层：纯 ELISA | WT Xcc 上测结合曲线，每轮 4–6 个 variant | 「我们找到了结合更好的 variant。」 | ✅ 有保障——这是 Cycle 0–2 实际内容 |
| 第二层：+ 噬斑测定 | +2 天：WT 上做 spot assay 确认裂解性感染 | 「结合 → 感染已确认。」 | ✅ 几乎零边际成本 |
| 第三层：+ ΔtonB/ΔexbB/ΔexbD1 | Markerless deletion（pK18mobsacB，4–6 周） | 「结合信号是受体特异性的。」这是**论文级别的机制性 claim**。 | 🟡 若 5/17 启动敲除 |

**内置阴性对照：** Hung 2003 已确认 ΔexbD2 不影响 phiL7 感染——生成此菌株即可免费验证整个敲除实验体系。若 ΔexbD2 仍允许感染而 ΔtonB 不允许，则受体特异性 claim 得到实验支持。

**建议：** 现在就承诺第三层。5/17 启动 pK18mobsacB 敲除。若 7 月初前成功，你们有五星级 iGEM 故事和论文级机制性 claim。

---

## 文献核查修正的事实

直接阅读 Hung et al. 2003（BBRC 302:878–884，PMID 12646254）后：

**ExbD2 不是 phiL7 感染的必需基因。** 只有 TonB、ExbB、ExbD1 是 essential。exbD2 敲除株（CH620）保留对 phiL7 的完全敏感性。所有规划文档已修正。

另外 4 处修正：
1. Boltz-2 affinity head 只支持小分子-蛋白。蛋白-蛋白对输出 NaN。用 ipTM 作结构信心 proxy——**不是**结合亲和力。
2. Greenman 2025 期刊是 *PLoS Comput Biol* 21(1):e1012639（不是 NAR Genomics）；结论：「没有单一最佳 UQ 方法」。
3. Hie 2024 用的是 ESM-1b/1v（不是 ESM-2）；每个抗体 ~20 个 variant（不是 ~50）。
4. Lee 2009 主动搜索并**明确找不到** OP1 tail fiber（ORF25）在 phiL7 中的同源物。这比「没提到 tail spike」更强——他们找了，但序列相似性方法找不到。我们的 rbp_01 用 Tail_spike_N HMM 识别，与此一致：HMM 能检测序列已分歧到 BLAST 看不见的蛋白。

---

## 需要 PI 做决定的事项

### 1. 受体敲除系统——pK18mobsacB 还是 CRISPRi？🔴 5/17 前
- **pK18mobsacB（当前计划）：** Markerless deletion，已在 *Xanthomonas* 上验证（Hung 2003）。4–6 周，永久干净的敲除株。质粒在 Addgene（#87097，~$95）。
- **CRISPRi：** 基因沉默，可逆，更快（~2 周）。在 Xcc 上不如 pK18mobsacB 成熟。
- **建议：** 对 tonB/exbB/exbD1 并行使用 pK18mobsacB。5/17 启动。

### 2. AlphaFold 3 模型权重申请 🔴 本周
- AF3 提供比 Boltz-2 更高质量的结构，含 trimer 预测。
- 需通过 Google 表单申请（学术用途）。1–7 天审批。
- 由 Alex 还是 PI 提交？（学术机构 email 优先。）

### 3. 验证层次承诺 🟡 6/1 前
见上方「验证策略」。若 5/17 启动敲除，建议承诺第三层。

### 4. 噬菌体富集来源 🟡 6/1 前
共分离需要多样性更丰富的富集底物。实验室在 LA 县有哪些农业污水或灌溉水来源？

### 5. 论文投稿意向 🟡 Cycle 0 前
是否计划与 iGEM wiki 并行投稿（*Bioinformatics* / *NAR Genomics & Bioinformatics*）？这会影响从第一个 ELISA 测量开始的数据文档化粒度。

---

## 过去这段时间做了什么——详细工作记录（5月7–12日）

### 5月7日 — 项目方向转型

决定放弃原有的 6 因子生物物理打分 pipeline（GC 含量、pI、尾纤维系统发育等），转向**闭环主动学习系统**。核心论点：6 因子方法只能根据序列相似性预测噬菌体是否可能感染宿主，无法从 wet lab 数据中学习，也无法随每个实验循环改进。BALD 驱动的深度集成可以直接建模不确定性，并通过每次 ELISA 测量缩小它。

### 5月7–8日 — 隔夜并行构建

7 个 AI agent 在隔夜（5/7→5/8）并行构建，每个负责一个 pipeline 模块（Module 00–06）。早晨起来：所有 7 个模块都有可运行的代码、测试套件和 README。

### 5月10日 — 修复数据问题 + Laguna 设置 + 第一次 Boltz-2 跑（错误蛋白）

构建后审查发现 6 个数据问题：
1. phiL7 `proteins.faa` headers 污染——从 NCBI 重新下载原始基因组，用 PHANOTATE 重新生成蛋白质序列。
2. 基因组清单中两个无效细菌 accession（`KY000037` → `GCF_000092025.1`；`PY746849` → `GCF_000006765.1`）——通过 NCBI Assembly 搜索修正。
3. 630MB 基因组 binary 被提交到 git——移至 `.gitignore`；可通过 `fetch_phages.py` + `fetch_bacteria.py` 重新下载。
4. 设置 CARC Laguna HPC：OnDemand 门户、Code Server、conda `boltz2` 环境。需要 8 次调试才能解决 CUDA/trifast/optree 版本冲突（最终方案：torch 2.5.1+cu121，NVIDIA L40S 驱动 550.90.07 = CUDA 12.4；仅在 macOS 上需要 setuptools<70）。
5. 在 Laguna 上首次 Boltz-2 GPU 运行（job 59949，47 秒）——**错误蛋白**：FASTA 文件命名为 `rbp_01` 但包含 P25（85 aa）而非真正的 rbp_01（712 aa）。ipTM = 0.345 结果作废。
6. 诊断缓存 bug：Boltz-2 的 `--override` 标志只重跑预测，不清除 `processed/` 缓存。修复方案：每次跑之前 `rm -rf` 整个输出目录。

### 5月11日 — 文献审查 + 正确的 Boltz-2 结果

阅读全部 19 篇核心论文，逐条核查 project documents 中的每个定量声明。发现并修正 5 个事实性错误（见上方「文献核查修正的事实」部分）。最严重的：Wang 2003（支持 ExbD2 为必需基因的引用）根本不存在——是一个幻觉引用。真实论文是 Hung 2003，其结论恰好相反：ExbD2 不是必需的。

用 `03_rbp_identification/outputs/EU717894.1_rbp_sequences.faa` 中真正的 rbp_01 712aa 序列覆写了 FASTA 文件的 Chain A。清除全部 Boltz-2 输出目录后在 Laguna 重跑（job 59986）。

结果：**ipTM = 0.365，chain A ptm = 0.808**。PAE 矩阵（1316×1316）和逐残基 pLDDT 以 `.npz` 格式保存。

Module 04 在本地用真实的 3 个 RBP 序列（rbp_01、rbp_02、rbp_03）重跑——移除了 2 个 mock 占位序列（rbp_04、rbp_05）。

### 5月12日 — Module 07 BALD + 全 pipeline 完成

发现 Module 06 的 `predictions.csv` 只导出了 `std`（总不确定性 = epistemic + aleatoric），但 BALD 只需要 epistemic 分量（Var_k[μ_k]）。修改 `run_cycle0.py`，通过 `ens.predict(return_epistemic=True)` 额外导出 `epistemic_std`。重跑 Module 06（约 25 秒）。

构建 Module 07：
- `bald.py`：BALD 分数 = epistemic_std；`select_batch()` 含 `selection_reason` 审计列、随机对照臂、随机重放、安全备选
- `run_bald.py`：CLI 48 小时 SLA 流程脚本；自动检测预测路径；处理 Cycle 0（无测量对）和 Cycle 1+（排除过滤）；保存全部输出 + `run_meta.json`
- 18 个测试，覆盖 BALD 数学正确性（red-green：最高 epistemic 排第一，最低的不在 top-K）、已测量 pair 排除、Cycle 0 空路径、schema 验证、fallback 警告

在合成数据上首次 BALD 运行：top pick 是 `rbp_07 × rec_02`（BALD=0.218）。全部 5 个输出文件已提交。

所有规划文档更新以反映完整 pipeline。两个 commits 推送到 `active-learning-pipeline`（`6422a6b`、`c541726`）。

---

## 关键输出与附件目录

所有文件均位于 `active-learning-pipeline` 分支上。路径相对于仓库根目录。

> **查看 3D 结构：** 用 PyMOL（免费下载：pymol.org）或 ChimeraX（免费：rbvi.ucsf.edu）打开 `.pdb` 文件。两条链以 ribbon 形式显示：Chain A = rbp_01（712 aa），Chain B = TonB（604 aa）。B-factor 列编码逐残基 pLDDT。

---

### Module 01 — 相互作用矩阵（训练标签）

| 文件 | 说明 |
|------|------|
| `01_data_ground_truth/outputs/interaction_matrix.csv` | **主要输出。** 2,236 个 phage-host 对。315 个阳性（文献核查），1,921 个阴性。列：phage_acc、host_acc、label、confidence、source、notes。 |
| `01_data_ground_truth/outputs/negative_samples/negative_data_combined.csv` | 完整阴性集，含方法标签（跨属 + pathovar 推断）。 |

---

### Module 03 — RBP 候选

| 文件 | 说明 |
|------|------|
| `03_rbp_identification/outputs/EU717894.1_rbp_candidates.csv` | phiL7 80 个 ORF 中**3 个 RBP 候选**。列：orf_id、length_aa、hmm_score、hmm_match_pfam、combined_score、rank。rbp_01（orf_00001，712 aa，HMM score 342）是主要靶点。 |

---

### Module 04 — 蛋白质 Embedding

| 文件 | 说明 |
|------|------|
| `04_protein_embedding/outputs/embeddings_esm2_t6_8M_phiL7_rbps.npz` | rbp_01、rbp_02、rbp_03 的 3 × 320 维 ESM-2（8M）向量。加载：`np.load(path)['array']`。生产版（1280 维，650M）待 Laguna 上跑。 |
| `04_protein_embedding/outputs/embeddings_esm2_t6_8M_xcc_receptors.npz` | Xcc 受体蛋白（TonB、ExbB、ExbD1）的 ESM-2 embedding。 |

---

### Module 05 — 3D 结构（Boltz-2，job 59986）⭐ 主要可视化输出

> **路径说明：** Boltz-2 会生成深层嵌套目录。所有关键文件都在：
> `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/`
> 其中子路径 `boltz_results_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/predictions/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/` 是 Boltz-2 内部命名规则自动生成的。完整路径如下。

| 文件（从仓库根目录的完整路径） | 说明 |
|------|------|
| `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/boltz_results_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/predictions/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB_model_0.pdb` | **3D 原子结构。** Chain A = rbp_01（712 aa），Chain B = TonB（604 aa）。用 PyMOL 或 ChimeraX 打开。 |
| `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/boltz_results_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/predictions/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/pae_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB_model_0.npz` | **PAE 矩阵（热图数据）。** 1316×1316 float32。元素 [i,j] = 残基 i 与 j 相对位置的预测对齐误差（Å）。低值 = 高置信度。界面块（行 712–1315 × 列 0–711）= rbp_01 vs TonB。加载：`np.load(path)['pae']`。 |
| `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/boltz_results_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/predictions/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/plddt_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB_model_0.npz` | 逐残基 pLDDT。1316 个值，范围 0.30–0.98，均值 0.76。加载：`np.load(path)['plddt']`。 |
| `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/boltz_results_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/predictions/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/confidence_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB_model_0.json` | 逐链置信度。关键值：`iptm=0.365`，`chains_ptm={'0': 0.808, '1': 0.494}`。 |
| `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/affinity.json` | 汇总分数。`interface_ipTM=0.365`，`chain_A_ptm=0.808`，`confidence_score=0.683`。`predicted_dG=null`——亲和力头仅支持小分子-蛋白。 |

**生成 PAE 热图的代码（从仓库根目录运行）：**
```python
import numpy as np, matplotlib.pyplot as plt
BASE = "05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB"
PRED = f"{BASE}/boltz_results_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/predictions/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB"
pae = np.load(f"{PRED}/pae_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB_model_0.npz")['pae']
plt.figure(figsize=(8, 8))
plt.imshow(pae, cmap='bwr_r', vmin=0, vmax=30)
plt.axhline(712, color='k', lw=1); plt.axvline(712, color='k', lw=1)
plt.colorbar(label='PAE (Å)'); plt.title('rbp_01 × TonB — Predicted Aligned Error')
plt.xlabel('残基（Chain A: 0–711 = rbp_01 / Chain B: 712–1315 = TonB）')
plt.savefig('pae_heatmap.png', dpi=150); print("Saved pae_heatmap.png")
```

---

### Module 06 — 深度集成（Cycle 0，合成数据）

| 文件 | 说明 |
|------|------|
| `06_uncertainty_model/outputs/cycle_0/calibration.png` | **校准图。** 预测不确定性区间 vs 真实标签覆盖率。良好校准 = 点落在对角线上。当前基于合成数据；真实版本将在 Cycle 0 ELISA 后产生。 |
| `06_uncertainty_model/outputs/cycle_0/predictions.csv` | 80 个 RBP×受体对的预测。列：rbp_id、receptor_id、predicted_score、std、**epistemic_std**、lower_95、upper_95、model_version。 |
| `06_uncertainty_model/outputs/cycle_0/training_log.json` | 各成员训练曲线（5 个成员，每个成员的 train NLL 和 val NLL 逐 epoch 记录）。 |
| `06_uncertainty_model/outputs/cycle_0/model_meta.json` | 模型架构（hidden_dim=256、n_members=5、input_dim=640）、训练参数、校准 ECE。 |

---

### Module 07 — BALD 采集（Cycle 1 推荐）

| 文件 | 说明 |
|------|------|
| `07_acquisition_function/outputs/cycle_1/recommendations.csv` | **下一个 wet lab 循环的 5 个推荐 variant。** 4 个 BALD 最优 + 1 个随机对照。列：rbp_id、receptor_id、predicted_score、epistemic_std、bald_score、bald_rank、selection_reason。 |
| `07_acquisition_function/outputs/cycle_1/safe_pick_backup.csv` | Top-10 BALD picks，用于 48 小时 SLA 未达时 PI/团队手动选择。 |
| `07_acquisition_function/outputs/cycle_1/random_replay.csv` | 纯随机策略会选什么——保存用于项目结束时的事后学习曲线对比。 |
| `07_acquisition_function/outputs/cycle_1/run_meta.json` | 溯源记录：cycle、seed、候选池大小、BALD 分数统计、top picks、git SHA、时间戳。 |

---

### 归档 — 转型前基线（仅供历史参考）

| 文件 | 说明 |
|------|------|
| `archive/2026-05-pivot/05_predictive_modeling/outputs/baseline_taxonomy_knn/heatmap_combined.png` | 分类学 KNN 基线预测热图（转型前方法）。展示其局限性：预测仅基于系统发育距离，无蛋白特异性。 |
| `archive/2026-05-pivot/05_predictive_modeling/outputs/baseline_taxonomy_knn/heatmap_original.png` | 原始 ground-truth 相互作用矩阵（与预测热图同尺度）。 |
| `archive/2026-05-pivot/03_feature_weighting/outputs/per_factor/factor2/f02_acidity_vs_y_scatter.png` | pI 与宿主范围散点图（转型前因子分析）。 |

---

## 时间线

```
2026-05-07   项目转型决定（6 因子 → 主动学习）
2026-05-07 → 2026-05-08   隔夜并行构建（7 个 agent，Module 00–06）
2026-05-10   修复数据问题；Laguna HPC 设置；首次 Boltz-2 跑（错误蛋白）
2026-05-11   文献审查（19 篇论文，5 处修正）；正确的 Boltz-2 结果 ✅
2026-05-12   Module 07 BALD 完成；全 pipeline（00–07）提交并推送 ✅

2026-05-17   Wet lab 启动：
             • 十字花科蔬菜采样（LA 县）
             • 受体敲除质粒构建（pK18mobsacB）
             • Cycle 0 variants 下单（基因合成，IDT/Twist）

2026-05-17 → 2026-06-01
             菌株 + 噬菌体分离；Cycle 0 variants 在 BL21 中表达

2026-06-01 → 2026-06-14
             Cycle 0：ELISA 优化 + 第一批结合测量

2026-06-14 → 2026-06-28
             Cycle 1：集成重训，BALD 推荐下一批 variants（SDM）

2026-06-28 → 2026-07-12
             Cycle 2：第二轮 + 受体敲除完成 + 因果验证
```

---

*完整 pipeline：`docs/planning/iGEM_2026_Project_Plan.md`（EN）| `docs/planning/iGEM_2026_项目大纲_中文版.md`（ZH）*
*Boltz-2 结构：`05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/`*
*文献审查笔记：`docs/reference/paper_reading_notes.md`*
