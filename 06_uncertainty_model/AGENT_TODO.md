# AGENT TODO — Module 06 / `06_uncertainty_model/`

## Read first (mandatory)

1. `/INTERFACE.md` §Module 06 — your output spec.
2. `/CLAUDE.md` — notebook-first, bilingual.
3. `/docs/planning/iGEM_2026_项目大纲_中文版.md` §3.5 — the Chinese project plan's deep explanation of why deep ensembles vs alternatives.
4. Invoke `superpowers:test-driven-development` and `superpowers:verification-before-completion`.

## Your scope

Build the uncertainty-aware regression model that will eventually predict ELISA binding scores from RBP × receptor embedding pairs. Tonight you have **no real ELISA data** (first batch comes ~6/1 from wet lab). You build the architecture + training pipeline + tests using **synthetic labels derived from Module 05's Boltz-2 priors** — purely as a placeholder so the pipeline shape is correct and tomorrow we can plug in real labels.

You consume:
- `04_protein_embedding/outputs/embeddings_esm2_t6_8M_*.npz` — RBP and receptor embeddings.
- `05_structure_prediction/outputs/affinity_priors.csv` — synthetic labels (Boltz-2 predicted ΔG).

You produce a trained deep ensemble + predictions + calibration plot + tests.

## Goal (definition of done)

By morning, `06_uncertainty_model/` contains:

1. ✅ `processes/01_train_deep_ensemble_synthetic.ipynb` — trains a 5-member ensemble on synthetic data (or pure-noise data if upstream isn't ready).
2. ✅ `processes/02_calibration_check.ipynb` — produces a reliability diagram + calibration metrics.
3. ✅ `outputs/cycle_0/ensemble_member_<i>.pt` for `i ∈ {0..4}`.
4. ✅ `outputs/cycle_0/predictions.csv` per INTERFACE schema (synthetic predictions on held-out synthetic data).
5. ✅ `outputs/cycle_0/calibration.png`.
6. ✅ `outputs/cycle_0/training_log.json` and `model_meta.json`.
7. ✅ `outputs/MANIFEST.csv`.
8. ✅ ≥4 passing tests (one extra for ensemble diversity).
9. ✅ `AGENT_REPORT.md` — including the swap-in plan: when real ELISA labels arrive, exactly which lines of the notebook change.

## Setup

```bash
cd /path/to/agent-06-uncertainty-model
conda activate igem2026

# PyTorch (CPU build is fine for tonight's small ensemble)
pip install torch
pip install torchmetrics  # for calibration metrics

# That's it. scikit-learn, numpy, pandas, matplotlib, seaborn already in env.

# Verify
python -c "import torch; print(torch.__version__, 'CUDA:', torch.cuda.is_available())"
```

If you'd like to also benchmark against a Gaussian Process baseline:
```bash
pip install gpytorch
```
This is optional for tonight; if installed, mention in AGENT_REPORT for tomorrow's comparison.

## Step-by-step plan

### Step 1 — Synthetic data generator (45 min)

Real ELISA data isn't available. Build a synthetic generator that produces a small dataset matching the eventual real shape:

```python
def generate_synthetic_dataset(
    rbp_embeddings: np.ndarray,         # (N_rbp, embed_dim)
    receptor_embeddings: np.ndarray,    # (N_rec, embed_dim)
    boltz2_priors: pd.DataFrame,        # rows: rbp_id, receptor_id, predicted_pKd
    noise_std: float = 0.5,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """
    Returns:
      X: (n_pairs, 2*embed_dim) — concat(rbp_emb, receptor_emb)
      y: (n_pairs,) — synthetic binding score derived from Boltz-2 priors + noise
      meta: DataFrame with rbp_id, receptor_id, true_y, source
    """
```

Logic:
- For (rbp, receptor) pairs present in `boltz2_priors`: y = predicted_pKd + Gaussian noise.
- For (rbp, receptor) pairs NOT in `boltz2_priors`: y = sample from N(5.0, 1.0) — mimics typical pKd range.
- This is deliberately leaky so the ensemble can learn *something*; the goal tonight is pipeline plumbing, not biological accuracy.

If Module 04 / 05 outputs aren't there: generate fully-synthetic embeddings (`np.random.randn(20, 320)` for 20 RBPs and `np.random.randn(4, 320)` for 4 receptors) and synthetic Boltz-2 priors. Document the fallback.

### Step 2 — Build `01_train_deep_ensemble_synthetic.ipynb` (120 min)

Cells:

1. **Markdown** — Title + bilingual purpose. Cite Lakshminarayanan 2017.
2. **Code** — Imports, paths, version printouts, **set torch + numpy seeds = 42**.
3. **Markdown** — Method (bilingual): why deep ensemble vs MC Dropout vs GP. Cite Greenman 2025 NAR Genom Bioinform (best UQ benchmark for protein engineering).
4. **Code** — `class EnsembleMember(nn.Module)`: 3-layer MLP with hidden_dim=256, output mean + log_sigma (Gaussian NLL training).
5. **Code** — `train_member(X_train, y_train, X_val, y_val, seed: int) → (model, log)`. Adam optimizer, lr=1e-4, early stop on val NLL with patience=10, max 200 epochs.
6. **Code** — `class DeepEnsemble`: container that wraps 5 trained members. Methods:
   - `train(X_train, y_train, X_val, y_val, n_members=5)` — trains members with seeds 0..4.
   - `predict(X) → (mean, std)` — returns ensemble mean and std across members.
7. **Markdown** — Sample-then-batch.
8. **Code** — Sample run: load synthetic data from Step 1; do 80/20 train/val split (stratified by RBP if possible); train ensemble.
9. **Code** — Sanity assertions:
   - 5 members trained; each has different final weights.
   - Predictions on val set are finite.
   - std > 0 for at least 90% of predictions (means members actually disagree somewhere).
10. **Code** — Save members to `outputs/cycle_0/ensemble_member_<i>.pt`.
11. **Code** — Generate predictions CSV per INTERFACE §Module 06 schema. Include all (RBP × receptor) combinations seen in synthetic data, with `predicted_score`, `std`, `lower_95`, `upper_95`, `model_version`.
12. **Code** — Save `training_log.json` (per-member training curves) and `model_meta.json` (arch, hyperparams, train_size, repo commit sha, timestamp).
13. **Code** — Update MANIFEST.

### Step 3 — Build `02_calibration_check.ipynb` (75 min)

Cells:

1. **Markdown** — Title + bilingual purpose. Why calibration matters for active learning (BALD relies on well-calibrated uncertainty).
2. **Code** — Load saved ensemble + predictions from Step 2.
3. **Code** — Compute calibration metrics:
   - Expected Calibration Error (ECE) — bin predictions by predicted std, compare predicted std to actual residual.
   - Reliability diagram — scatter of predicted std vs absolute residual (binned).
4. **Code** — Plot reliability diagram, save to `outputs/cycle_0/calibration.png` (300 DPI, white background).
5. **Code** — If ECE > 0.1, suggest temperature scaling — implement and re-evaluate.
6. **Markdown** — Conclusion: is the ensemble well-calibrated? If yes, we're ready for BALD acquisition (Module 07). If no, document and add follow-up to AGENT_REPORT.

### Step 4 — Tests (75 min)

`processes/tests/`:
- `test_schema.py` — Open `predictions.csv`, assert columns. Open `model_meta.json`, assert required keys (`arch`, `hidden_dim`, `train_size`, `val_size`, `repo_commit_sha`, `timestamp_utc`).
- `test_smoke.py` — Generate 50 synthetic samples (random embeddings, random labels), train a 2-member ensemble for 10 epochs, assert it produces predictions with finite mean + std.
- `test_diversity.py` — Train 5 members with different seeds on the same data; assert their predictions on a held-out batch have std > 0 for ≥ 50% of points (pure agreement = ensemble useless).
- `test_calibration.py` — Generate well-calibrated synthetic data (predictions = labels + Gaussian noise with known std); train ensemble; assert ECE < 0.2.

Run: `pytest 06_uncertainty_model/processes/tests/ -v`

### Step 5 — Commit + report

Commits:
- `06_uncertainty_model: PyTorch ensemble architecture`
- `06_uncertainty_model: training notebook (synthetic data)`
- `06_uncertainty_model: calibration notebook`
- `06_uncertainty_model: cycle_0 outputs (synthetic placeholder)`
- `06_uncertainty_model: tests (incl. diversity + calibration checks)`
- `06_uncertainty_model: AGENT_REPORT (incl. swap-in plan for real ELISA)`

The AGENT_REPORT must contain a "When ELISA arrives" section listing exactly:
- Which file path the real labels go into.
- Which cell in the notebook reads labels.
- What format the labels need to be in (e.g., `08_cycle_data/outputs/cycle_0/elisa_processed.csv` with columns `rbp_id, receptor_id, EC50_nM, log_EC50`).

## References (cite in notebook markdown cells)

- **Deep ensembles (the method)**: Lakshminarayanan, B. et al. (2017) "Simple and Scalable Predictive Uncertainty Estimation Using Deep Ensembles." *NeurIPS 2017*. arXiv:1612.01474. **PRIMARY METHODOLOGY REFERENCE.**
- **Why deep ensembles outperform MC Dropout for UQ**:
  - Beluch, W.H. et al. (2018) "The Power of Ensembles for Active Learning in Image Classification." *CVPR 2018*.
  - Ovadia, Y. et al. (2019) "Can You Trust Your Model's Uncertainty? Evaluating Predictive Uncertainty Under Dataset Shift." *NeurIPS 2019*.
- **UQ benchmark for protein engineering** (key for our domain): Greenman, K.P. et al. (2025) "Benchmarking uncertainty quantification methods for protein engineering." *NAR Genom Bioinform*.
- **Calibration metrics**: Guo, C. et al. (2017) "On Calibration of Modern Neural Networks." *ICML 2017*.
- **Active learning context** (downstream Module 07 motivation): Houlsby, N. et al. (2011) "Bayesian Active Learning for Classification and Preference Learning." arXiv:1112.5745. (BALD acquisition function — relies on the ensemble uncertainty you produce.)
- **Active learning for protein engineering** (the closest published work):
  - Hie, B.L. et al. (2022) "Efficient evolution of human antibodies from general protein language models." *Cell* 185:1-15.
  - Yang, J. et al. (2025) "Active Learning-assisted Directed Evolution (ALDE)." *Nat Commun*. (Standardized benchmark of acquisition functions over deep ensembles + GP.)
- **MC Dropout (alternative we're not using, for reader context)**: Gal, Y. & Ghahramani, Z. (2016) "Dropout as a Bayesian Approximation." *ICML 2016*.

## Anti-goals

- ❌ Don't try to use real ELISA data — it doesn't exist yet (wet lab Cycle 0 starts ~6/1).
- ❌ Don't claim the synthetic-data calibration result generalizes to real data — it doesn't. State that explicitly in AGENT_REPORT.
- ❌ Don't train large models tonight — 3-layer MLP with hidden_dim=256 is enough.
- ❌ Don't write the BALD acquisition function — that's Module 07's job (which is excluded from tonight's parallel build).
- ❌ Don't push to remote.

## Time budget

~4-5 hours. Split: 1 hour data generator + 2 hours training notebook + 1 hour calibration notebook + 1 hour tests.

## If stuck

- Module 04 / 05 outputs not available → use fully-synthetic random embeddings (documented fallback in Step 1). The pipeline shape is what matters tonight.
- Ensemble members all converge to the same prediction (zero std) → check that random seeds actually differ, that initialization is random (not zero-init), and that you're using `torch.manual_seed` per member separately. Ensemble diversity is a hard requirement.
- Calibration is terrible (ECE > 0.3) → likely because synthetic data is too easy / too hard. Tweak `noise_std` in the generator; document.
- Out of memory → reduce `hidden_dim` from 256 to 64, this module needs to fit on a Mac CPU.
