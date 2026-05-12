# iGEM Claremont 2026 — Technical Brief for PI Review

**Companion to the full plan** (`iGEM_2026_Project_Plan.md`). This brief covers Sections 3, 4, 5, 9, and 11 only — the technical core, where we are, and the questions we need PI input on.

**Date:** May 7, 2026 

---

## Status Legend

| Symbol | Meaning |
|---|---|
| ✅ | Complete or available |
| 🔄 | In progress (named team member actively working) |
| ⏳ | Queued — start within 2 weeks |
| 🔧 | Needs revision / rework |
| ⬜ | Not started |

---

## 3. Dry-Lab (Computational) Modules

### 3.1 Genome processing & annotation
**Methods:** Prodigal for bacterial ORFs (Hyatt et al., 2010); PHANOTATE for phage ORFs (McNair et al., 2019); pharokka for functional annotation (Bouras et al., 2023).

**Where we are:** ✅ Genome database built (`00_raw_data/`: 34 bacteria + 777 phages downloaded). ✅ PHANOTATE + Prodigal scripts implemented (`02_annotation/processes/`). 🔄 Annotation runs partial across the dataset.

### 3.2 RBP identification
**Methods:** PhageRBPdetect — HMM + XGBoost ensemble; 93.8% precision-recall AUC for phage RBP identification (Boeckaerts et al., 2022, *Viruses*).

**Where we are:** ⏳ Tool identified, integration pending. Will run on phiL7 (NCBI EU717894) reference genome first, then on in-house phage isolates as they come online.

### 3.3 Structural prediction
**Methods:** AlphaFold 3 for static RBP and receptor structures (Abramson et al., 2024, *Nature*); Boltz-2 for zero-shot binding affinity priors (Passaro et al., 2025, bioRxiv).

**Where we are:** ⏳ Will execute on Laguna HPC. AF3 needed for both phiL7 RBP and Xcc TonB-dependent receptor (no PDB structure exists for either).

### 3.4 Protein language model embeddings
**Methods:** ESM-2 (Lin et al., 2023, *Science*) — 650M variant for routine use, 3B for final benchmarking. Optional: PLM-interact (Liu et al., 2025, *Nat Commun*) as transfer prior from human PPI to phage-bacteria binding.

**Where we are:** 🔄 **ESM embedding pipeline actively in development by team member with prior ESM experience.** Module 3.4 is our most-advanced dry-lab component.

### 3.5 Uncertainty-aware regression model
**Methods:** Deep ensemble of 5 MLPs over ESM-2 embeddings (Lakshminarayanan et al., 2017, *NeurIPS*). Benchmarked across multiple protein fitness landscapes (Greenman et al., 2025, *PLoS Comput Biol* 21(1):e1012639); no single UQ method dominates all scenarios.

**Where we are:** ⬜ Architecture designed; implementation begins after Module 3.4 outputs are stable.

### 3.6 Acquisition function (active learning core)
**Methods:** BALD — Bayesian Active Learning by Disagreement (Houlsby et al., 2011); recommends variants where ensemble members maximally disagree. Random selection is the control baseline (Hie et al., 2022, *Cell*).

**Where we are:** ⬜ To be implemented. Reference: Yang et al., 2025, *Nat Commun* — Active Learning-assisted Directed Evolution (ALDE) — provides the most directly applicable benchmark.

### 3.7 Cycle infrastructure
**Methods:** Versioned data pipeline that ingests new ELISA results, retrains the ensemble, and outputs next-cycle recommendations within 48 hours. PyTorch + MLflow + Docker for reproducibility on Laguna.

**Where we are:** ⬜ To be built during dry-lab sprint (May 7–17).

> **What's deprecated from the original plan:** The 6-feature biophysical scoring pipeline (`03_feature_weighting/`, factors 1–6) is being retired. ESM-2 + PLM-based scoring (Module 3.4–3.5) replaces it with substantially stronger literature support.

---

## 4. Wet-Lab (Experimental) Modules

### 4.1 Strain isolation from California sources ⭐ **Per PI consultation, this replaces commercial culture-collection ordering.**
**Methods:** Surface-sterilize symptomatic brassica tissue (cabbage, broccoli, kale showing V-shaped chlorotic lesions); enrich on YPGA + cycloheximide (250 mg/L) at 28°C; pick yellow mucoid colonies; identify by 16S rRNA + MALDI-TOF + *gyrB*/*rpoD* MLST. Protocol from EPPO PM 7/110 (2013) and Mwangi et al. (2007, *Plant Pathol*); identification per Parkinson et al. (2009, *J Appl Microbiol*).

**Where we are:** ⬜ Not started. Source plan: local farmers' markets + (with grower permission) commercial fields in LA County. Target: 3–5 verified Xanthomonas isolates within 2 weeks of starting.

**Why self-isolate:** USDA APHIS PPQ-526 permits average 127 days to process; in-state environmental isolation does not require this permit and yields registry-quality novel isolates.

### 4.2 Phage isolation
**Methods:** Co-isolation from same source material + sewage / agricultural-runoff supplementation. Standard enrichment per Adams (1959) and Bonilla et al. (2016, *PeerJ*); 3-round single-plaque purification.

**Where we are:** ⬜ Not started. Will run in parallel with Module 4.1.

### 4.3 RBP cloning & expression
**Methods:** N-terminal His6 in pET-28a; Gibson assembly (Gibson et al., 2009, *Nat Methods*) for variant construction; expression in *E. coli* BL21(DE3) at 18°C overnight with 0.5 mM IPTG (Studier & Moffatt, 1986).

**Where we are:** ⬜ Not started. Cycle 0 variant design (4–6 variants: truncations + chimeras + point mutations) is the immediate dry-lab deliverable.

### 4.4 Protein purification (Ni-NTA)
**Methods:** IMAC with imidazole gradient elution (Hochuli et al., 1988); SDS-PAGE + anti-His6 Western for QC.

**Where we are:** ⬜ Not started. Standard column protocol — needs PI confirmation that Ni-NTA resin and FPLC/gravity column setup are available.

### 4.5 ELISA-based binding assay ⭐ **Primary readout for active-learning cycles.**
**Methods:** Plate-based binding of His6-RBP serial dilutions (1 nM – 1 µM) against immobilized Xanthomonas cells (10⁸ CFU/well); HRP-anti-His6 detection; OD450; 4-parameter logistic fit for EC50 / apparent Kd. Adapted from Boeckaerts et al. (2024, *Nat Commun*) and Latka et al. (2021, *mBio*).

**Where we are:** ⬜ Not started. **Highest wet-lab risk:** ELISA optimization (positive control = T7 gp17 with published binding data) must be locked before Cycle 0 launches.

### 4.6 Receptor knockout for causal validation
**Methods:** Markerless deletion of *tonB* / *exbB* / *exbD1* / *exbD2* using pK18mobsacB suicide vector (Schäfer et al., 1994, *Gene*; Addgene #87097). KanR + sucrose counter-selection. Targets: *tonB*, *exbB*, *exbD1* (all essential per Hung et al. 2003) + *exbD2* as negative control (not required for infection per Hung et al., 2003, *BBRC* 302:878–884, PMID 12646254).

**Where we are:** ⬜ Not started. Existing Benchling Xcc transformation protocol (electroporation, 14 kV/cm, 10 µF) provides the required tool. Plasmid order needed.

### 4.7 Plaque assay validation
**Methods:** NYG soft agar (0.6–0.7%); 100 µL log-phase cells + serial-diluted phage; 30°C overnight.

**Where we are:** ✅ **Protocol exists on Benchling, no changes needed.** This and Module 4.1 cultivation/transformation are the two wet-lab capabilities we already have running.

> **Existing Benchling protocols inventory:**
> ✅ Xanthomonas Cultivation · ✅ Xanthomonas Transformation · ✅ Plaque Assay · ✅ Liquid Amplification · 🔧 Phage Infection Curves (current version is bacterial density study — needs rewrite to add phage challenge)

---

## 5. The Active-Learning Cycle

### Cycle structure (~2 weeks per cycle)

```
Week 1                Week 2-3                Week 4
─────────────────────────────────────────────────────────
Model recommends     Wet lab clones,         Results → dry lab
4-6 variants    ──▶  expresses, purifies, ──▶ model retrains,
                     runs ELISA               next cycle begins
                     + 1 random control
```

### Cycle 0 (Weeks 0–2): Cold start
- Variant selection: Boltz-2 zero-shot affinity prior (Module 3.3) flags maximally divergent candidates; team adds 1–2 expert "intuition" picks for sanity-checking.
- 4–6 variants spanning truncation / mutation / chimera classes.
- ELISA optimization runs in parallel (T7 gp17 positive control).

### Cycles 1–2 (Weeks 2–6): Closed loop
1. Dry lab retrains ensemble on cumulative dataset.
2. BALD ranks unmeasured (variant, receptor) combinations.
3. Recommend top 4–5 + 1 random control.
4. Wet lab measures.
5. Predicted vs measured: residuals + calibration check.
6. Document and version.

### Evaluation metrics
- Per-cycle test R² on held-out variants.
- Calibration: predicted uncertainty vs absolute residual.
- **AL vs random comparison** (the key claim): retrospective replay using all final data, plus the 1-per-cycle random control measurements.
- Information gain per experiment (KL divergence of posterior from prior).

### Honest reporting commitment
If by end of Cycle 2 the AL trajectory's test R² is statistically equivalent to random, we will not claim AL advantage in the writeup. Negative results are still scientifically informative for the iGEM community.

---

## 9. Resource Requirements

### Computational
- **Laguna HPC access** for ESM-2 3B inference, AlphaFold 3 batch, Boltz-2 affinity. Estimated ~200 GPU-hours over project duration.
- **In-team Laguna operator** identified ✅
- **In-team ESM expertise** identified ✅ (currently driving Module 3.4)
- Personal workstation for model development and acquisition function.

### Wet-lab consumables (Cycle 0 priority order)
| Item | Vendor / source | Approx. cost |
|---|---|---|
| pET-28a backbone | Addgene #69864-3 | $95 |
| pK18mobsacB | Addgene #87097 | $95 |
| BL21(DE3) competent cells | NEB | $200 |
| Gibson Assembly Master Mix (50 rxn) | NEB | $330 |
| Ni-NTA agarose | Qiagen | $280 |
| 96-well ELISA plates + HRP-anti-His6 | Thermo / Sigma | $400 |
| YPGA + cycloheximide media components | Sigma | $200 |
| Gene synthesis (4–6 RBP variants for Cycle 0) | Twist or IDT | $400–600 |
| **Total Cycle 0 estimated spend** | | **~$2,000–2,500** |

### Personnel allocation (current state)
- Alex — dry-lab core, full-time during sprints.
- Sarah, Olivia, Maggie, Camille — three parallel wet-lab pipelines (cloning / expression+purification / ELISA).
- Ryan (ESM-experienced team member) — currently driving Module 3.4 🔄
- Carol (Laguna-trained team member) — on call for HPC submissions.
- All wet-lab leads cross-train on culture handling and plaque assays so any can substitute.

---

## 11. Open Questions for PI Consultation

We're flagging five concrete decisions where PI input would unblock or accelerate the next two weeks:

### 11.1 Laguna HPC allocation
Confirm the team has access to Laguna with sufficient quota for ~200 GPU-hours and ~500 GB scratch storage. We have an in-team trained operator ready to submit jobs.

### 11.2 Receptor-knockout system selection
We propose pK18mobsacB markerless deletion (Schäfer 1994; Addgene #87097) as primary, with CRISPRi as fallback if double-crossover efficiency is poor in our isolated strain. Does the PI prefer a different system based on lab experience?

### 11.3 USDA / CDFA permit posture
*X. citri* (citrus canker) is a CDFA quarantine concern. If our isolation effort happens to recover a citrus-pathogenic isolate (possible if we sample citrus tissue), what is the lab's documentation and reporting obligation? Does the lab have any standing PPQ-526 permit we could use?

### 11.4 Sewage / agricultural runoff sampling
For phage co-isolation, we plan to enrich from agricultural-runoff sources beyond direct plant tissue. Is there a preferred local source the PI recommends (campus agricultural plots, partner growers)?

### 11.5 Manuscript ambition
Is the team scoping toward a parallel manuscript submission (target: *Bioinformatics* or *NAR Genom Bioinform*), or focusing on the iGEM wiki / promotion video as primary deliverables? This affects how we structure the data documentation from Cycle 0 onward.

---

## Snapshot — where we are right now

| Component                                                       | Status                                                        |
| --------------------------------------------------------------- | ------------------------------------------------------------- |
| Dry-lab genome dataset (777 phages + 34 bacteria)               | ✅ Built                                                       |
| Phage / host annotation pipeline                                | 🔄 Partial runs done                                          |
| ESM-2 embedding pipeline                                        | 🔄 In development                                             |
| RBP identification, structure prediction, AL model, cycle infra | ⏳ Next 10 days                                                |
| Old 6-feature scoring pipeline                                  | 🔧 Being retired                                              |
| Wet-lab cultivation, transformation, plaque assay, lysate prep  | ✅ Benchling protocols ready                                   |
| Phage Infection Curves protocol                                 | 🔧 Needs rewrite (current version is bacterial density study) |
| RBP cloning + expression + Ni-NTA + ELISA + receptor KO         | ⬜ To stand up before Cycle 0                                  |
| Strain + phage isolation from California sources                | ⬜ Starts immediately after PI sign-off                        |

**The two-week critical path from today (May 7):** strain isolation pipeline online + dry-lab core (Modules 3.4–3.7) functional, in time for Cycle 0 to launch after both converge (~late May / early June).
