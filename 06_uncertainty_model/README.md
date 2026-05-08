# 06 — Uncertainty-Aware Regression Model

Train a regression model that predicts ELISA binding score (EC50 / apparent Kd) from RBP variant + receptor sequence pair, **with calibrated uncertainty estimates**. Uncertainty drives the active learning acquisition function in Module 07.

## Method

**Deep ensemble** (Lakshminarayanan et al. 2017, *NeurIPS*) — 5 independently-trained MLPs over ESM-2 embeddings:
- Each member: 3-layer MLP, hidden dim 256, dropout 0.1
- Loss: Gaussian negative log-likelihood (each member outputs mean + sigma)
- Training: Adam, lr 1e-4, early stopping on validation NLL
- **Predictive mean** = average of 5 outputs
- **Epistemic uncertainty** = variance across 5 outputs

### Why deep ensemble vs Gaussian Process / MC Dropout
- Better calibration than MC Dropout (Beluch 2018; Ovadia 2019)
- Scales better than GP on 1280-D ESM-2 inputs
- Recently benchmarked best for protein engineering UQ (Greenman et al. 2025, *NAR Genom Bioinform*)

## Calibration check

After every cycle: produce a reliability diagram. If 80%-confidence predictions don't actually contain the truth ~80% of the time, apply temperature scaling.

## Inputs

`inputs/` reads from:
- `04_protein_embedding/outputs/*.npy` — ESM-2 embeddings of RBP variants
- `05_structure_prediction/outputs/boltz2/zero_shot_priors.csv` — synthetic prior for cold start
- `08_cycle_data/outputs/cycle_<N>_elisa.csv` — accumulating wet lab data

## Outputs

`outputs/`
- `cycle_<N>/ensemble_member_<i>.pt` — trained MLP weights
- `cycle_<N>/predictions.csv` — `(variant, receptor, predicted_score, uncertainty)` for all unmeasured pairs
- `cycle_<N>/reliability_diagram.png` — calibration check
- `cycle_<N>/training_log.json`

## Status

⏳ Not started — implementation begins after Module 04 (ESM-2 embeddings) outputs are stable.
