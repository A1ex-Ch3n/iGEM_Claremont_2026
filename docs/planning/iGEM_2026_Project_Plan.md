# iGEM Claremont 2026 — Project Plan
## Active-Learning–Guided Phage Engineering for *Xanthomonas* Biocontrol

**Date:** May 7, 2026
**Faculty advisors:** Prof. J. Cesar Ignacio-Espinoza (PI), Prof. Ran Libeskind-Hadas
**Awards targeted:** Best Agriculture Project · Best Model · Best Composite Part

---

## Executive Summary

We are building a **closed-loop active-learning pipeline** that integrates a machine-learning model of phage receptor-binding-protein (RBP) – host receptor interactions with iterative wet-lab validation. The central innovation is the use of **Bayesian Optimal Experimental Design** (BOED; Lindley, 1956) to let the model nominate the most informative next experiment, rather than passively training on whatever data is available. This addresses the most-cited bottleneck in computational phage biology — labeled phage-host interaction data is scarce, expensive, and species-bound (Boeckaerts et al., 2024, *Nat Commun*; Mutalik group benchmark, 2025, bioRxiv).

To our knowledge, no published system applies closed-loop active learning to phage RBP – receptor binding prediction, although the underlying methodology is well-established in adjacent domains: directed antibody evolution (Hie et al., 2022, *Cell*), enzyme engineering (Yang et al., 2025, *Nat Commun* — "Active Learning-assisted Directed Evolution", ALDE), and Bayesian optimization for protein design (Romero & Arnold, 2009, *Nat Chem Biol*).

The agricultural target is *Xanthomonas*, a genus causing major losses in brassicas, citrus, and solanaceous crops worldwide. **Per PI consultation, we will self-isolate the host strain and a paired lytic phage from California symptomatic plant tissue**, bypassing dependency on commercial culture collections (which require multi-month USDA APHIS PPQ-526 permits). Dry-lab development proceeds in parallel on public reference genomes — *X. campestris* pv. *campestris* (Xcc) ATCC 33913 (NCBI: AE008922) and phage phiL7 (NCBI: EU717894), whose receptor was experimentally identified as a TonB-dependent system (Hung et al., 2003, *BBRC* 302:878–884, PMID 12646254).

The deliverable is a system, not a single model: a quantitative motif-level RBP–receptor binding atlas, a benchmarked active-learning framework, and a small library of in-house *Xanthomonas* isolates and phages contributed back to the research community via the iGEM Registry.

---

## 0. What We Know, What We Don't Know, and What We Learn

### What is already known

- phiL7 infects *Xcc* through the TonB-ExbB-ExbD1 receptor system (TonB, ExbB, ExbD1 are essential; ExbD2 is co-expressed in the same operon but not required for infection) (Hung et al., 2003, *BBRC* 302:878–884).
- The phiL7 tail spike protein (rbp_01, 712 aa) is the receptor-binding protein responsible for host recognition (Lee et al., 2009, *AEM*; Module 03 HMM result).
- Mutating RBP C-terminal domains can redirect phage host range (Latka et al., 2021, *mBio*; Yehl et al., 2019, *Cell*).
- The best existing phage-host prediction model reaches AUC 0.82 within one genus but degrades to 0.67 across genera — because labeled binding data is extremely scarce (Boeckaerts et al., 2024, *Nat Commun*; PAML benchmark, 2025).

### What is not known — the core research gap

**For rbp_01 (712 amino acids), which positions determine binding strength to TonB? How does changing those positions affect affinity?**

This is a *fitness landscape* problem: the space of possible rbp_01 variants has 20^712 points, and we have essentially zero quantitative measurements. No one has systematically mapped RBP–receptor binding affinity for any *Xanthomonas* phage system.

The deeper open question: **can computational priors (Boltz-2 ipTM structural confidence proxy, ESM-2 protein language model embeddings, PLM-interact PPI transfer learning) meaningfully accelerate learning this landscape from a small number of wet-lab measurements?** This has never been tested in phage biology.

### What we learn through this pipeline

| Question | How we answer it |
|----------|-----------------|
| Which rbp_01 variants bind TonB better or worse? | ELISA binding curves, ~30–60 variants across 2–3 cycles |
| Does active learning (BALD) outperform random selection? | Retrospective replay: same data, two selection policies |
| Do computational priors (Boltz-2, ESM-2) transfer to phage-bacteria? | Compare model performance with vs without each prior layer |
| How much of phage infection is receptor-specific vs defense-dependent? | Layer 2: ΔtonB / ΔexbB knockout causal decomposition |

The core scientific contribution is **methodological**: a validated framework for efficiently exploring phage RBP fitness landscapes with limited wet-lab budget. If active learning provides 2–5× efficiency gain (as shown in antibody and enzyme engineering), Cycle 0–2 with ~30 ELISA measurements is equivalent in information to 60–150 random measurements — at a fraction of the cost.

### Why wet-lab cost stays manageable

Cycle 0 (4–6 variants) requires gene synthesis (~$500). Subsequent cycles primarily use site-directed mutagenesis (~$50/variant, 4 days), not new gene synthesis, because BALD selects point mutations on existing constructs. Total wet-lab spend across all cycles: ~$2,000–2,500.

---

## 1. Background and Motivation

### 1.1 *Xanthomonas* impact and current control limitations

The genus *Xanthomonas* contains pathogens of more than 400 plant species (Ryan et al., 2011, *Nat Rev Microbiol*). California-relevant pathovars include Xcc (black rot of brassicas; widely present in commercial vegetable production), *X. citri* subsp. *citri* (citrus canker), and *X. perforans* / *X. euvesicatoria* (bacterial spot of tomato and pepper). Current management depends heavily on copper-based bactericides; resistance is now widespread (Aiello et al., 2019, *Plant Dis*) and ecological externalities — particularly impact on beneficial rhizosphere microbiota — have prompted regulatory attention (Lamichhane et al., 2018, *Crop Prot*).

Bacteriophage biocontrol offers an attractive alternative: phages are host-specific, self-amplifying, and environmentally degradable. Field trials report plaque-based formulations matching copper hydroxide efficacy for bacterial spot of tomato (Iriarte et al., 2018, *Front Microbiol*) and successful suppression of Xcc in cabbage (Holtappels et al., 2022, *Microb Biotechnol*).

### 1.2 The prediction problem

Despite these advances, *predicting* which phage will infect which bacterial strain remains unreliable. State-of-the-art models (PhageHostLearn, AUC 0.82; Boeckaerts et al., 2024, *Nat Commun*) train on relatively few labeled (phage, strain) pairs (10²–10³) and degrade sharply across host genera (PAML benchmark; Mutalik group, 2025, bioRxiv). The fundamental bottleneck is **data scarcity** — generating new phage-host interaction labels is slow and labor-intensive.

### 1.3 Why active learning is the right framing

Active learning (AL) is a mathematical response to data scarcity (MacKay, 1992, *Neural Comput*; Settles, 2009 *Active Learning Literature Survey*). Rather than training on whatever data happens to exist, an AL system uses its current uncertainty to identify which next experiment, if performed, would most improve the model. Acquisition functions such as Bayesian Active Learning by Disagreement (BALD; Houlsby et al., 2011) and Expected Information Gain (EIG; Lindley, 1956) formalize this.

Recent demonstrations show AL-guided directed evolution outperforms random or expert-guided variant selection by 2–5× in protein engineering (Hie et al., 2022, *Cell*: human antibody evolution from PLMs; Yang et al., 2025, *Nat Commun*: ALDE benchmark; Wittmann et al., 2021, *Cell Syst*: ML for directed evolution). Critically, **none of these have been applied to phage RBPs or phage host-range prediction.** This is the gap our project occupies.

### 1.4 Why this is iGEM-aligned

The project naturally targets three iGEM evaluation dimensions:
- **Best Agriculture Project**: a deployable framework for targeted phage discovery against an economically important plant pathogen.
- **Best Model**: a closed-loop active-learning system with rigorous uncertainty quantification, control comparisons, and causal validation.
- **Best Composite Part**: a registered library of RBP-His6 expression constructs for community use.

---

## 2. System Architecture

The pipeline integrates three layers running in parallel:

```
   ┌─────────────────────────────────────────────────────────────┐
   │  LAYER 0 — Reference Initialization (Pre-cycle, public data)│
   │  • Public phage RBP corpora (PhageRBPdetect; ~5K sequences)│
   │  • Boltz-2 ipTM structural confidence (binding proxy)      │
   │  • PLM-interact transfer-learning prior (human PPI)        │
   └─────────────────────────────────────────────────────────────┘
                              │ informative prior
                              ▼
   ┌─────────────────────────────────────────────────────────────┐
   │  LAYER 1 — Closed-Loop Active Learning (~2-3 cycles)       │
   │                                                             │
   │   ┌──────────────┐    ┌────────────────┐    ┌────────────┐│
   │   │ ESM-2 +      │    │ Acquisition    │    │ Wet lab    ││
   │   │ deep         │───▶│ function       │───▶│ ELISA      ││
   │   │ ensemble     │    │ (BALD / EIG)   │    │ binding    ││
   │   │ regressor    │◀───│ recommends 4-6 │◀───│ on RBP     ││
   │   │              │    │ next variants  │    │ variants   ││
   │   └──────────────┘    └────────────────┘    └────────────┘│
   └─────────────────────────────────────────────────────────────┘
                              │ binding scores
                              ▼
   ┌─────────────────────────────────────────────────────────────┐
   │  LAYER 2 — Causal Validation                                │
   │  • Receptor knockout (pK18mobsacB markerless deletion)      │
   │  • Plaque assay on WT vs ΔReceptor strains                 │
   │  • Decouples binding signal from defense systems           │
   └─────────────────────────────────────────────────────────────┘
```

The two-layer design (binding prediction → infection calibration) explicitly addresses a fundamental confound: phage RBP binding is necessary but not sufficient for productive infection (Farquharson et al., 2021, *J Virol*; Doud et al., 2025, *Front Cell Infect Microbiol* review). Layer 2's receptor knockouts let us quantify how much of the model's binding signal translates to infection across genetic backgrounds.

---

## 3. Dry-Lab (Computational) Modules

### Module 3.1 — Genome processing and annotation

**Inputs:** Public reference (Xcc ATCC 33913 genome, NCBI AE008922; phiL7 genome, NCBI EU717894); subsequently, in-house assembled genomes from isolated strains/phages.

**Methods:**
- **Bacterial gene prediction:** Prodigal (Hyatt et al., 2010, *BMC Bioinformatics*) — the standard for Gram-negative bacterial ORF calling, with documented sensitivity > 99% on closed genomes.
- **Phage gene prediction:** PHANOTATE (McNair et al., 2019, *Bioinformatics*) — purpose-built for phages, handles overlapping ORFs that defeat general-purpose callers.
- **Functional annotation:** pharokka (Bouras et al., 2023, *Bioinformatics*) — integrated phage annotation pipeline with PHROG/VFDB/CARD databases.

**Output:** Curated `.faa` protein files per genome.

### Module 3.2 — RBP identification

**Methods:** PhageRBPdetect (Boeckaerts et al., 2022, *Viruses* 14:1329) — a dual-track HMM + XGBoost pipeline that identifies receptor-binding proteins from phage genomes with 93.8% precision-recall AUC. Combines Pfam-derived HMMs with custom HMMs for RBP-relevant Pfam-absent domains.

**Why this method:** Existing RBP-identification approaches (PhANNs, HMMER alone) miss RBPs that have diverged from canonical Pfam motifs, which is common in *Xanthomonas* phages (Lee et al., 2009, *AEM* — phiL7 has unique features absent in Pfam).

**Output:** Ranked RBP candidates per phage with confidence scores.

### Module 3.3 — Structural prediction

**Methods:**
- **AlphaFold 3** (Abramson et al., 2024, *Nature*) — for predicting both phage RBP trimer structures and host receptor (TonB / ExbB / ExbD1 / ExbD2) structures.
- **Boltz-2** (Passaro et al., 2025, bioRxiv) — for batched zero-shot structural confidence estimation of RBP variants against receptors. Includes an explicit affinity head trained on small molecule–protein binding data (ChEMBL, BindingDB); outputs NaN for protein-protein pairs — we use ipTM as structural confidence proxy, which AlphaFold 3 does not expose per-complex.

**Why both:** AF3 provides higher-quality static structures; Boltz-2 provides ipTM structural confidence scores as a binding quality proxy during the cold-start phase. This dual-model strategy follows the Toronto iGEM 2025 PHORAGER team's design philosophy.

**Output:** PDB-format predictions; per-variant predicted ΔΔG estimates (used as soft labels for Layer 0 prior).

### Module 3.4 — Protein language model embeddings

**Methods:** ESM-2 (Lin et al., 2023, *Science*), specifically the 650M-parameter variant for routine use and the 3B variant for final benchmarking. Embeddings are mean-pooled over residues for sequence-level prediction, or extracted per-residue for motif-level analysis.

**Why ESM-2:** Trained on 65M unique protein sequences with masked-language-modeling objective; demonstrated state-of-the-art on protein function prediction tasks. Critically, protein language model embeddings transfer to phage proteins despite phage proteins being underrepresented in training data (phage protein transfer demonstrated in Hie et al., 2024, *Nat Biotechnol* using ESM-1b).

**Optional layer:** PLM-interact (Liu et al., 2025, *Nat Commun*) — fine-tuned ESM-2 on human PPI data, demonstrated AUPR improvement of 16–28% when transferred to virus-host PPI in mouse, fly, worm, yeast, and *E. coli*. We will benchmark whether this transfer extends to phage-bacteria binding.

**Compute requirements:** ESM-2 650M runs on single 12 GB GPU; ESM-2 3B requires 24 GB+ (will use Laguna HPC cluster).

**Output:** Numerical embeddings (1280-D for 650M; 2560-D for 3B).

### Module 3.5 — Uncertainty-aware regression model

**Methods:** Deep ensemble of 5 multilayer perceptrons (Lakshminarayanan et al., 2017, *NeurIPS*) trained on ESM-2 embeddings to predict ELISA binding scores. Each ensemble member uses a different random seed; predictive uncertainty is the variance across members.

**Why deep ensembles over Gaussian Processes:**
- Better calibration than MC Dropout (Beluch et al., 2018, *CVPR*; Ovadia et al., 2019, *NeurIPS*).
- Handles ESM-2's 1280-D inputs without GP scaling issues.
- Competitive UQ method for protein engineering across multiple fitness landscapes (Greenman et al., 2025, *PLoS Comput Biol*: "Benchmarking uncertainty quantification for protein engineering").
- ALDE (Yang et al., 2025) used deep ensembles + Gaussian Processes; we will start with ensembles for simplicity and add GPs as a comparative.

**Output:** Predicted binding score + uncertainty per (RBP variant, receptor) pair.

### Module 3.6 — Acquisition function

**Methods:**
- **Primary:** BALD (Bayesian Active Learning by Disagreement; Houlsby et al., 2011) — selects variants where ensemble members maximally disagree.
- **Secondary baseline:** Random selection — required to demonstrate AL is outperforming random, the first question any reviewer will ask.

**Why BALD over greedy / UCB:** BALD directly optimizes information gain, not predicted performance. This is critical when training data is small (early cycles), where greedy strategies tend to over-exploit local optima (Yang et al., 2025).

**Control arm design:** Following Hie et al. (2022), each cycle's recommended batch will include 1 randomly-selected variant; we will additionally retrospectively replay the experiment as if random selection had been used throughout, to quantify the AL advantage.

**Output:** Ranked variant list per cycle (top 6 selected; 4-5 enter wet lab, 1-2 held in reserve).

### Module 3.7 — Cycle infrastructure

A reproducible data pipeline that:
1. Ingests new ELISA results from wet lab.
2. Re-trains the ensemble on the cumulative dataset.
3. Runs the acquisition function to recommend next variants.
4. Outputs a wet-lab task list within 48 hours of receiving new data.
5. Versions all model checkpoints, data snapshots, and predictions for retrospective analysis.

**Tools:** Python, PyTorch, scikit-learn, ONNX for model versioning, MLflow for experiment tracking, all containerized with Docker for reproducibility on Laguna HPC.

---

## 4. Wet-Lab (Experimental) Modules

### Module 4.1 — Strain isolation from California sources

**Approach:** *Per PI consultation, we will isolate Xanthomonas strains directly from symptomatic plant tissue collected in California, bypassing commercial culture collections.*

**Source material:** Symptomatic brassica leaves (cabbage, broccoli, kale) showing characteristic black-rot V-shaped chlorotic lesions, sourced from local farmers' markets and (with permission) commercial growers within Los Angeles County.

**Protocol:** Standard semi-selective enrichment following EPPO PM 7/110 guidance (EPPO Bulletin, 2013) and Mwangi et al. (2007, *Plant Pathol*):
1. Surface-sterilize symptomatic tissue (1% NaOCl, 1 min; 70% EtOH rinse).
2. Macerate in sterile saline.
3. Plate on YPGA + cycloheximide (Yeast Peptone Glycerol Agar with 250 mg/L cycloheximide to suppress fungal contaminants).
4. Pick yellow mucoid colonies after 48 hr at 28°C (characteristic of *Xanthomonas*).
5. **Identification:** 16S rRNA Sanger sequencing + MALDI-TOF (HMC core facility); pathovar confirmation by *gyrB*/*rpoD* multilocus sequence typing (Parkinson et al., 2009, *J Appl Microbiol*).

**Why self-isolation is preferred:** USDA APHIS PPQ-526 permits for plant-pathogen import average 127 days processing time; in-state environmental isolation does not require this permit. The approach also yields novel isolates that are scientifically interesting and registrable as iGEM Parts.

**Risk mitigation:** If primary isolation yields non-Xanthomonas yellow contaminants (*Pantoea*, *Erwinia*), 16S/MALDI-TOF screening filters them at the colony stage. If no Xanthomonas is recovered after two rounds, secondary sources (citrus canker fruit, tomato bacterial spot) will be sampled.

**Output:** ~3-5 verified Xanthomonas isolates with whole-genome shotgun sequencing (Illumina MiSeq, HMC core).

### Module 4.2 — Phage isolation

**Approach:** Co-isolation from the same source material, plus a sewage/agricultural-runoff enrichment for broader phage diversity.

**Protocol:** Adams (1959) classical enrichment, modernized per Bonilla et al. (2016, *PeerJ*):
1. Filter (0.22 µm) aqueous extracts of source material.
2. Mix with mid-log Xanthomonas isolate from Module 4.1.
3. Incubate 30°C, 200 rpm, 12-16 hr.
4. Re-filter; spot on Xanthomonas lawn.
5. Single-plaque purification × 3 rounds to ensure clonality.

**Characterization:**
- Genome sequencing (Illumina MiSeq).
- TEM morphology (HMC core, optional).
- Lytic/temperate determination from genome (presence of integrase, repressor genes).
- Host range against the panel of in-house isolates.

**Output:** ~2-3 lytic phage isolates with sequenced genomes, archived as glycerol stocks (≥10⁹ PFU/mL).

### Module 4.3 — RBP cloning and expression

**Strategy:** Express RBPs as N-terminal His6-tagged proteins in *E. coli* BL21(DE3) using the pET system (Studier & Moffatt, 1986, *J Mol Biol*).

**Cloning method:** Gibson assembly (Gibson et al., 2009, *Nat Methods*) — 2 fragments (linearized vector + RBP insert), thermocycler-free, supports library construction for variant series.

**Variant design:** For Cycle 0 seed batch, 4–6 variants spanning three classes:
- **Truncation series:** N-terminal anchor vs C-terminal specificity domain (per Latka et al., 2021, *mBio* — Klebsiella RBP modularity; Yang et al., 2024, *GigaScience* — SpikeHunter analysis showing C-terminal swap drives specificity).
- **Targeted point mutations:** Predicted binding-interface residues from AlphaFold structure (per Yehl et al., 2019, *Cell* — T7 RBP DMS).
- **Chimeras:** N-terminal anchor of one RBP + C-terminal head of another (where >1 RBP is identified across our isolates).

**Expression conditions:** BL21(DE3) + pET-28a backbone; IPTG induction (0.5 mM, 18°C overnight) to favor soluble expression of trimeric tail proteins.

**Output:** ~4-6 soluble His6-RBP proteins per cycle.

### Module 4.4 — Protein purification

**Method:** Ni-NTA immobilized metal affinity chromatography (Hochuli et al., 1988, *Bio/Technology*); standard column-based protocol with imidazole gradient elution.

**Quality control:**
- SDS-PAGE for size and purity (target ≥ 90% by densitometry).
- Western blot with anti-His6 antibody to confirm identity.
- Bradford assay for concentration.
- (Optional) Size-exclusion chromatography to confirm trimerization.

**Output:** Purified RBP variants at ≥ 0.1 mg/mL, stored at 4°C in PBS + 10% glycerol.

### Module 4.5 — ELISA-based binding assay

**Approach:** Plate-based binding assay quantifying RBP attachment to immobilized Xanthomonas cells. This is the **primary readout** for the active-learning cycle.

**Protocol** (adapted from Boeckaerts et al., 2024 *Nat Commun* and Latka et al., 2021 *mBio*):
1. Coat 96-well ELISA plate with heat-inactivated Xanthomonas cells (10⁸ CFU/well).
2. Block (3% BSA, 1 hr).
3. Apply serial dilutions of purified His6-RBP (range 1 nM – 1 µM).
4. Wash; detect with HRP-conjugated anti-His6 antibody.
5. Develop with TMB; read OD450.
6. Fit to 4-parameter logistic to extract EC50 / apparent Kd.

**Controls per plate:**
- BSA-only well (no Xanthomonas) — non-specific binding.
- Wild-type RBP at fixed concentration — internal standard for inter-plate normalization.
- Heat-denatured RBP — confirms binding specificity to folded protein.
- 3 technical replicates per concentration.

**Why ELISA over SPR/BLI:** SPR/BLI are technically preferable but capital-intensive; ELISA is the standard for higher-throughput phage RBP studies (Boeckaerts et al., 2024). Will validate against a single SPR measurement at HMC core if available, as a proof of equivalence.

**Output per cycle:** ~12-30 quantitative binding curves (4–6 RBP variants × 2-5 receptor sources × 3 replicates).

### Module 4.6 — Receptor knockout for causal validation

**Approach:** Markerless gene deletion of candidate receptor genes in our Xanthomonas isolates using the pK18mobsacB suicide-vector system (Schäfer et al., 1994, *Gene*; Addgene #87097).

**Receptor targets:** Based on phiL7 reference biology (Hung et al., 2003, *BBRC* 302:878–884) — *tonB*, *exbB*, *exbD1* (all three essential); *exbD2* is also tested as a negative control (Hung 2003 shows it is NOT required for phiL7 penetration). If our isolated phage targets a different receptor, we will identify it via comparative genomics of phage-resistant escape mutants.

**Protocol:**
1. Construct deletion plasmid: ~500 bp upstream + ~500 bp downstream homology arms flanking the target gene, cloned into pK18mobsacB (kanamycin resistance + sacB sucrose counter-selection).
2. Electroporate into Xanthomonas isolate (per existing Benchling protocol; 14 kV/cm, 10 µF).
3. Single-crossover selection on kanamycin.
4. Double-crossover selection on sucrose-containing medium (sacB-mediated counter-selection; sucrose-tolerant kanamycin-sensitive colonies have completed deletion).
5. PCR-confirm deletion; sequence verify.

**Validation experiment (the core causal test):**
- Plaque assay of in-house phage on (WT, ΔtonB, ΔexbB, ΔexbD) panel.
- ELISA binding of WT-RBP and select variants on (WT, ΔReceptor) cell coats.
- Compare:
  - WT-RBP × WT-strain → baseline binding + infection.
  - WT-RBP × ΔReceptor → binding loss × infection loss = receptor-specific contribution.
  - Variant-RBP × WT-strain vs WT-RBP × WT-strain → variant effect.
  - Variant-RBP × ΔReceptor → epistasis test.

**Output:** Receptor-knockout strains; quantitative decomposition of binding vs receptor-specific infection vs background.

### Module 4.7 — Plaque assay validation (Layer 2 readout)

**Use of existing protocols:** Sarah's Benchling small-scale plaque assay protocol (NYG soft agar, 0.6–0.7%; 100 µL log-phase cells + serial-diluted phage; 30°C overnight) is suitable as-is.

**Output:** Plaque-forming units per mL (PFU/mL) and plaque morphology per (phage, strain) pair.

---

## 5. The Active-Learning Cycle

### Cycle structure

Each cycle is approximately 2 weeks:

```
   Week 1                Week 2-3                Week 4
   ─────────────────────────────────────────────────────────
   Model recommends     Wet lab clones,         Results returned
   4-6 variants    ──▶  expresses, purifies,   ──▶  to dry lab,
                         runs ELISA              model retrains
                         + 1 random control
```

### Cycle 0 (Weeks 0-2): Cold-start

**Goal:** Generate the first batch of training data with maximally informative seed variants.

**Variant selection:** 4-6 variants chosen by:
- **Prior-driven:** Boltz-2 ipTM structural confidence scores (Module 3.3) flag variants with widely-divergent predicted binding geometry (encourages information spread).
- **Expert-driven:** team selects 1-2 "biological intuition" picks for sanity-checking the model.
- **Coverage-driven:** distribute across truncation / mutation / chimera variant classes.

### Cycles 1-2 (Weeks 2-6): Closed loop

**Goal:** Demonstrate improving model performance under AL guidance.

**Per-cycle workflow:**
1. Dry lab retrains ensemble on all data from prior cycles + Cycle 0.
2. BALD acquisition function ranks all *unmeasured* variant-receptor combinations.
3. Recommend top 4-5 + 1 random control.
4. Wet lab measures.
5. Predicted vs measured comparison: per-variant residual + ensemble calibration.
6. Document and version.

### Cycle 2-3 evaluation

**Performance metrics:**
- Per-cycle test-set R² on held-out variants (cross-validated).
- Calibration: predicted uncertainty vs absolute residual.
- AL vs random comparison: retrospective replay of "what if we had used random selection from cycle 0?"
- Information gain per experiment (KL divergence of posterior from prior).

**Decision criterion:** If by end of Cycle 2 the AL trajectory's test R² is statistically equivalent to random, we do not claim AL advantage in the writeup. We commit to honest reporting per iGEM transparency norms.

---

## 6. Critical Risks and Mitigations

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| **Strain isolation fails** (no recoverable Xanthomonas from local sources) | Low | High | Dual-source strategy: brassica primary, citrus/tomato backup; budget for 2 rounds of sampling within first 3 weeks |
| **Phage isolation yields no lytic phage on our isolate** | Moderate | High | Sewage/agricultural runoff as supplementary source; if 4 weeks pass without success, request phage stock from collaborators (Phage Directory community channel) |
| **RBP expression failure** (insolubility, mis-folding in trimer assembly) | Moderate | Medium | Validate on T7 gp17 (positive control with published binding data) before committing variant library; budget GCN4 leucine zipper trimerization tag as backup (Frey et al., 2008, *J Mol Biol*) |
| **ELISA dynamic range insufficient** | Moderate | High | Run optimization phase weeks 2-3 of Cycle 0 (positive control + dilution series) before locking format; backup readout: bio-layer interferometry at HMC core |
| **Active learning underperforms random** | Low | High (but reportable) | Standard control arm; honest reporting; useful negative result documented |
| **Wet/dry sync breaks** (dry lab can't deliver next-cycle recommendation in 48 hr) | Moderate | Medium | Pre-defined "safe-pick" backup variant list — wet lab proceeds without dry lab if SLA missed |
| **Receptor knockout difficulty** (poor electroporation, sacB selection escapes) | Moderate | Medium | Multiple receptor targets in parallel (tonB, exbB, exbD); CRISPRi as alternative knockdown if knockout fails |

---

## 7. Models and Methods Summary

| Module | Method | Reference |
|---|---|---|
| Phage gene calling | PHANOTATE | McNair et al., 2019, *Bioinformatics* |
| Bacterial gene calling | Prodigal | Hyatt et al., 2010, *BMC Bioinf* |
| Phage annotation | pharokka | Bouras et al., 2023, *Bioinformatics* |
| RBP identification | PhageRBPdetect | Boeckaerts et al., 2022, *Viruses* |
| Structure prediction | AlphaFold 3 | Abramson et al., 2024, *Nature* |
| Affinity prior | Boltz-2 | Passaro et al., 2025, bioRxiv |
| Sequence embedding | ESM-2 (650M / 3B) | Lin et al., 2023, *Science* |
| PPI transfer prior | PLM-interact | Liu et al., 2025, *Nat Commun* |
| Regression + UQ | Deep ensemble | Lakshminarayanan et al., 2017, *NeurIPS* |
| Acquisition function | BALD | Houlsby et al., 2011 |
| AL framework reference | ALDE | Yang et al., 2025, *Nat Commun* |
| Closest published analog | PLM-guided antibody evolution | Hie et al., 2022, *Cell* |
| Strain isolation | Semi-selective YPGA + cycloheximide | Mwangi et al., 2007; EPPO PM 7/110, 2013 |
| Strain ID | 16S + MALDI-TOF + *gyrB*/*rpoD* MLST | Parkinson et al., 2009, *J Appl Microbiol* |
| Phage isolation | Enrichment + plaque purification | Adams, 1959; Bonilla et al., 2016, *PeerJ* |
| Cloning | Gibson assembly | Gibson et al., 2009, *Nat Methods* |
| Expression | pET / BL21(DE3) | Studier & Moffatt, 1986, *J Mol Biol* |
| Purification | Ni-NTA IMAC | Hochuli et al., 1988, *Bio/Technology* |
| Binding readout | Cell-based ELISA | Boeckaerts et al., 2024, *Nat Commun* |
| Receptor knockout | pK18mobsacB | Schäfer et al., 1994, *Gene* |

---

## 8. Timeline and Milestones

```
2026-05-07   TODAY
2026-05-07 → 05-17 (10 days)   DRY-LAB SPRINT
              • Baseline model + acquisition function operational
              • Reference genome analysis (Xcc + phiL7)
              • Variant library design
              • Gene synthesis ordered

2026-05-17 → 06-01              STRAIN ISOLATION PHASE
              • Source collection from California crops
              • Enrichment + plaque + identification
              • Whole-genome sequencing of isolates
              • Parallel: Cycle 0 variants cloned + expressed

2026-06-01 → 06-14              CYCLE 0 EXECUTION
              • ELISA optimization (T7 gp17 positive control)
              • First binding measurements on isolated strain
              • Layer 0 model retrained with ground truth

2026-06-14 → 06-28              CYCLE 1 EXECUTION
              • Acquisition function recommends round-1 variants
              • Wet-lab clones, expresses, measures
              • First quantitative AL vs random comparison

2026-06-28 → 07-12              CYCLE 2 EXECUTION
              • Round-2 recommendations
              • Receptor knockouts complete
              • Layer 2 causal validation

2026-07-12 → 08-09              ANALYSIS & EXTENSION
              • Apply trained model to soil-derived phage genomes
              • Encapsulation experiments (Pluronic F-127 hydrogel)
              • Cross-strain validation

2026-08-09 → 09-12              WRITEUP & PACKAGING
              • Wiki, software, parts registration
              • Manuscript draft

2026-09-12 → 10-21              POLISH
              • Wiki freeze 10/21
              • Presentation video

2026-11-13   GRAND JAMBOREE
```

---

## 9. Resource Requirements

### Computational
- Access to Laguna HPC cluster (1 team member trained in operation) — for ESM-2 3B inference, AlphaFold 3 batch jobs, and Boltz-2 affinity prediction.
- Personal workstation for model development.
- Estimated GPU-hours: ~200 over project duration (well within Laguna academic allocation).

### Wet lab consumables (priority for spring order)
- pET-28a plasmid (Addgene #69864-3, ~$95)
- pK18mobsacB plasmid (Addgene #87097 or similar, ~$95)
- BL21(DE3) competent cells (NEB, ~$200)
- Gibson Assembly Master Mix (NEB, ~$330 / 50 rxn)
- Ni-NTA agarose (Qiagen, ~$280)
- 96-well ELISA plates + His6 detection antibody (~$400 total)
- YPGA + cycloheximide media components (~$200)
- Gene synthesis (4-6 variants for Cycle 0; ~$400-600 from IDT/Twist)
- **Total estimated Cycle 0 wet lab spend:** ~$2,000-2,500

### Personnel
- 1 dry-lab core (Alex; full-time during sprints, partial otherwise).
- 3 wet-lab core (Sarah, Olivia, Weitao; parallelized as cloning / expression / ELISA pipelines).
- 1 ESM-experienced advisor (in-team consult).
- 1 Laguna-trained operator (in-team, for HPC jobs).

---

## 10. Expected Deliverables

### For iGEM
- **Wiki sections** documenting Engineering Cycles (DBTL), Modeling, Hardware/Software, and Human Practices.
- **Composite Part submission:** RBP-His6 expression library — at minimum 4-6 characterized RBP-His6 constructs derived from in-house phage isolates.
- **Software repository:** Open-source closed-loop active-learning pipeline.
- **Promotion video** demonstrating live execution of one closed-loop cycle.

### For the broader research community
- **In-house Xanthomonas isolate collection** (deposited with University accession, with contact for sharing).
- **Sequenced phage genomes** (deposited to NCBI).
- **Quantitative motif-level RBP–receptor binding atlas** (the dataset itself, released as supplementary).
- **Benchmark of closed-loop AL on phage RBP engineering** (manuscript, intended submission to *Bioinformatics* or *Nucleic Acids Research*).

---

## 11. Open Questions for PI Consultation

1. **HPC allocation:** Confirm Laguna access and storage quota for the team (estimated 200 GPU-hours, ~500 GB disk).
2. **Receptor knockout system selection:** PI's preference between pK18mobsacB markerless deletion vs Tn-mutant library vs CRISPRi. We have base protocol for the first.
3. **Permit posture:** Does the Ignacio-Espinoza lab hold or is willing to apply for any APHIS permits, in case our isolate identification reveals a quarantine-status pathovar (e.g., *X. citri* citrus canker would trigger CDFA quarantine concerns).
4. **Sewage / agricultural runoff sampling:** Approval and source recommendation for phage enrichment substrate.
5. **Manuscript timeline:** Is the team aiming for a parallel manuscript submission, or focusing on iGEM wiki + promotion video alone?

---

## Appendix A — Annotated Key References

### Core methodological framework
- **Lindley, D.V. (1956).** "On a measure of the information provided by an experiment." *Ann Math Stat* 27:986. — Foundational work for Bayesian experimental design.
- **Settles, B. (2009).** *Active Learning Literature Survey*. CS Tech Report, U. Wisconsin–Madison. — Canonical AL reference.
- **Houlsby, N. et al. (2011).** "Bayesian Active Learning for Classification and Preference Learning." *arXiv:1112.5745*. — BALD acquisition function.

### Active learning in protein engineering (closest analogs)
- **Hie, B.L. et al. (2022).** "Efficient evolution of human antibodies from general protein language models." *Cell* 185:1-15. — PLM-guided antibody evolution; methodologically the closest published work.
- **Yang, J. et al. (2025).** "Active Learning-assisted Directed Evolution." *Nat Commun*. — Recent benchmark of acquisition functions (greedy, UCB, Thompson) and surrogates (ESM-2, GP, DNN ensemble) on standardized landscapes.
- **Wittmann, B.J. et al. (2021).** "Advances in machine learning for directed evolution." *Curr Opin Struct Biol*. — Field overview.

### Phage host-range prediction (the field we extend)
- **Boeckaerts, D. et al. (2024).** "Predicting phage-host interactions in Klebsiella with PhageHostLearn." *Nat Commun* 15:4768. — SOTA strain-level model.
- **Liu, Y. et al. (2025).** "PLM-interact: extending protein language models." *Nat Commun*. — Demonstrates PPI transfer learning to virus-host interactions; we extend to phage-bacteria.
- **Boeckaerts, D. et al. (2022).** "Identification of phage receptor-binding protein sequences with HMMs and XGBoost." *Viruses* 14:1329. — RBP detection tool (PhageRBPdetect).

### *Xanthomonas* and phiL7 (our reference system)
- **da Silva, A.C.R. et al. (2002).** "Comparison of the genomes of two *Xanthomonas* pathogens." *Nature* 417:459. — Xcc ATCC 33913 genome (NCBI AE008922).
- **Lee, C.N. et al. (2009).** "Genomic characterization of the intron-containing T7-like phage phiL7." *Appl Environ Microbiol* 75:7828. — phiL7 genome (NCBI EU717894).
- **Hung, C.-H. et al. (2003).** "Involvement of *tonB-exbBD1D2* operon in infection of *Xanthomonas campestris* phage ϕL7." *Biochem Biophys Res Commun* 302(4):878–884. PMID: 12646254. — Receptor system identification.

### Plant-pathology and isolation methods
- **EPPO Bulletin (2013).** PM 7/110: "*Xanthomonas* spp. causing bacterial spot of tomato and sweet pepper." — Standard isolation protocol.
- **Mwangi, M. et al. (2007).** "Semi-selective medium for *X. campestris* pv. *musacearum*." *Plant Pathol*. — YPGA-cycloheximide variant.
- **Schäfer, A. et al. (1994).** "Small mobilizable multi-purpose cloning vectors derived from the *E. coli* plasmids pK18 and pK19." *Gene* 145:69. — pK18mobsacB system.

### Tools and machine learning
- **Hyatt, D. et al. (2010).** "Prodigal: prokaryotic gene recognition." *BMC Bioinformatics* 11:119.
- **McNair, K. et al. (2019).** "PHANOTATE: a novel approach to gene identification in phage genomes." *Bioinformatics* 35:4537.
- **Lin, Z. et al. (2023).** "Evolutionary-scale prediction of atomic-level protein structure." *Science* 379:1123. — ESM-2.
- **Abramson, J. et al. (2024).** "Accurate structure prediction of biomolecular interactions with AlphaFold 3." *Nature* 630:493.
- **Lakshminarayanan, B. et al. (2017).** "Simple and Scalable Predictive Uncertainty Estimation Using Deep Ensembles." *NeurIPS*.

---

**Document prepared for review by Prof. J. Cesar Ignacio-Espinoza and Prof. Ran Libeskind-Hadas.**
**Available for discussion at the team's earliest convenience.**
