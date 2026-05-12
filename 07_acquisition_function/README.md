# 07 — Acquisition Function & Cycle Infrastructure

The active learning core: given Module 06's predictions + uncertainties, choose which RBP variants the wet lab should test next. Wraps the cycle infrastructure (data ingest → retrain trigger → recommendation output, 48 h SLA).

## Method

### BALD — Bayesian Active Learning by Disagreement (Houlsby et al. 2011)

Score each unmeasured (variant, receptor) pair:

```
BALD(x) = H[E_p[y|x]] - E_p[H[y|x]]
        = predictive entropy - expected entropy
        ≈ variance of ensemble member predictions
```

Equivalent to mutual information between observation and model parameters. Selects variants where ensemble members maximally disagree — measuring them gives the highest information gain.

### Why BALD over greedy / UCB
- In small-data regime, greedy gets stuck in local optima
- BALD maximizes information gain regardless of predicted value
- Active learning for protein engineering shown effective by Yang et al. 2025 (*Nat Commun*, ALDE) — note ALDE uses Thompson sampling, not BALD; we adopt BALD for its theoretical grounding in mutual information

## Control arm (mandatory)

Every cycle's recommended batch (4–5 BALD top picks) **also includes 1 randomly-selected variant**. Plus we maintain a parallel virtual log of "what would random selection have picked" for retrospective replay (Hie et al. 2022 standard practice).

## Cycle infrastructure

48-hour SLA: ELISA data arrives → ensemble retrains → BALD recomputed → next-cycle CSV out.

If SLA missed: wet lab proceeds from a pre-prepared "safe pick" backup list (PI + bio team intuition picks).

## Tools

- Python 3.11+ / PyTorch
- MLflow for experiment + model versioning
- Docker for environment reproducibility
- HPC submission scripts for Laguna

## Inputs

`inputs/` reads from `06_uncertainty_model/outputs/cycle_<N>/predictions.csv`.

## Outputs

`outputs/`
- `cycle_<N+1>/recommendations.csv` — top 4–5 BALD picks + 1 random + their scores
- `cycle_<N+1>/safe_pick_backup.csv` — handpicked fallback list
- `cycle_<N+1>/random_replay.csv` — virtual "what if random" picks for retrospective analysis
- `retrospective_replay/` — final analysis comparing AL vs random learning curves

## Status

✅ **Done** — implemented 2026-05-11.

Files:
- `processes/bald.py` — core BALD scoring + batch selection
- `processes/run_bald.py` — CLI cycle orchestrator (48 h SLA)
- `processes/tests/` — 18 tests (all pass)

First live run (synthetic data): `outputs/cycle_1/` — see `recommendations.csv`.

To run for Cycle N (replace N with actual cycle number):
```bash
python 07_acquisition_function/processes/run_bald.py \
  --cycle N \
  --measured_csv 08_cycle_data/outputs/cycle_<N-1>/elisa_processed.csv \
  --n_bald 4 --n_random 1
```

Module 06 schema change: `predictions.csv` now includes `epistemic_std` column
(Std of ensemble member means — true BALD score, not proxied by total std).
