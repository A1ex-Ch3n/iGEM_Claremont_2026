# 05 — Structure Prediction

Predict 3D structures of (a) phage RBPs identified in Module 03 and (b) host receptors identified from bacterial genome annotation. Provides geometric features for variant design and a synthetic affinity prior for Module 06.

## Methods

### AlphaFold 3 (Abramson et al. 2024, *Nature* 630:493)
- MSA + diffusion-based structure generator
- Predicts RBP trimer + receptor structures
- Reports per-residue pLDDT + interface ipTM

### Boltz-2 (Passaro et al. 2025, bioRxiv)
- AlphaFold-style structure prediction with **explicit affinity head**
- Trained on PDBbind for binding ΔG estimation
- Used as zero-shot synthetic prior for cold-start (Cycle 0)

## Inputs

`inputs/` reads from:
- `03_rbp_identification/outputs/<phage>_rbp_sequences.faa`
- `02_annotation/outputs/host_proteins/<host>.faa` (TonB / ExbB / ExbD candidates)

## Outputs

`outputs/`
- `af3/<protein_id>.cif` — predicted 3D structure (AlphaFold 3)
- `af3/<protein_id>_metrics.json` — pLDDT / ipTM scores
- `boltz2/<rbp>__<receptor>.json` — predicted ΔG affinity per RBP-receptor pair
- `boltz2/zero_shot_priors.csv` — aggregate table fed into Module 06

## Status

⏳ Not started — runs on Laguna HPC. ESM-experienced + Laguna-trained team members will collaborate on submission scripts.
