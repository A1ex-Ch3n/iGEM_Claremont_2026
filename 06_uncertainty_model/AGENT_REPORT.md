# AGENT REPORT — Module 06 / `06_uncertainty_model/`

**Agent:** 06 (Deep Ensemble + Calibration)
**Date:** 2026-05-07
**Branch:** `agent-06-uncertainty-model`
**Status:** Cycle 0 complete (synthetic data); ready for real ELISA swap-in.

---

## What Was Built

### Architecture (`processes/ensemble.py`)

A **Deep Ensemble** of 5 independent 3-layer MLPs (Lakshminarayanan et al., 2017):

```
input_dim (640) → Linear(256) → ReLU → Dropout(0.1)
               → Linear(256) → ReLU → Dropout(0.1)
               → Linear(128) → ReLU
               → head_mean: Linear(1)    # binding score prediction
               → head_log_sigma: Linear(1), clamped [-7, 7]  # uncertainty
```

Total predictive variance = Var(member means) + E[member sigma²]  
(Eq. 3 from Lakshminarayanan 2017 — correct epistemic + aleatoric decomposition)

Each member trained with Adam(lr=1e-4), Gaussian NLL loss, early stopping on val NLL (patience=10), max 200 epochs.

### Outputs Generated (`outputs/cycle_0/`)

| File | Description |
|---|---|
| `ensemble_member_{0..4}.pt` | PyTorch state_dict for each of 5 members |
| `predictions.csv` | 80 pairs × (rbp_id, receptor_id, predicted_score, std, lower_95, upper_95, model_version) |
| `calibration.png` | Reliability diagram: std vs \|residual\| + coverage calibration curve |
| `training_log.json` | Per-member train/val NLL curves + best epochs |
| `model_meta.json` | Architecture, hyperparams, repo SHA, timestamp |
| `MANIFEST.csv` | SHA-256 + bytes for all cycle_0 files |

### Tests (`processes/tests/`) — 9/9 passing

| Test | Status | Runtime |
|---|---|---|
| `test_smoke.py::test_train_member_smoke` | ✅ PASS | ~3s |
| `test_smoke.py::test_deep_ensemble_smoke` | ✅ PASS | ~5s |
| `test_smoke.py::test_synthetic_data_generator` | ✅ PASS | <1s |
| `test_diversity.py::test_ensemble_diversity` | ✅ PASS | ~8s |
| `test_diversity.py::test_members_have_different_weights` | ✅ PASS | ~1s |
| `test_calibration.py::test_calibration_ece_below_threshold` | ✅ PASS | ~32s |
| `test_schema.py::test_predictions_csv_schema` | ✅ PASS | <1s |
| `test_schema.py::test_model_meta_json_schema` | ✅ PASS | <1s |
| `test_schema.py::test_manifest_schema` | ✅ PASS | <1s |

Run with: `pytest 06_uncertainty_model/processes/tests/ -v`

---

## Calibration Status (Cycle 0)

**ECE (regression) = 0.27** on the synthetic val set (n=16).

This is **expected and acceptable for Cycle 0 synthetic data**:
- Val set has only 16 samples (80 total pairs, 20% holdout)
- Data is purely random (no biological signal)
- ECE will be reassessed after real ELISA data arrives

The `calibration.png` shows the reliability diagram with two panels:
1. Predicted std vs |residual| scatter + bin-averaged overlay
2. Coverage calibration curve (expected vs actual coverage at each CI level)

**Temperature scaling** (Guo et al., 2017) is the recommended remedy if ECE > 0.1 persists on real data. Implementation sketch:

```python
# After training, find optimal temperature T on val set
from scipy.optimize import minimize_scalar
def neg_nll(T):
    scaled_std = std_pred / T
    return gaussian_nll_loss(mean_pred, scaled_std, y_val).item()
result = minimize_scalar(neg_nll, bounds=(0.5, 5.0), method='bounded')
optimal_T = result.x  # divide all std predictions by T at inference
```

---

## Why Deep Ensembles (Not MC Dropout / GP)

From `docs/planning/iGEM_2026_项目大纲_中文版.md` §3.5:

1. **Calibration** — Beluch et al. (2018 CVPR); Ovadia et al. (2019 NeurIPS): ensembles outperform MC Dropout under distribution shift.
2. **Scaling** — GP is O(N³) in training set size; intractable at ESM-2's 1280-dim inputs. Ensembles scale linearly.
3. **Domain validation** — Greenman et al. (2025 NAR Genom Bioinform): deep ensembles are top-tier for protein engineering UQ.
4. **ALDE precedent** — Yang et al. (2025 Nat Commun): ALDE uses deep ensembles as the core model for active-learning-directed evolution.

Future work: compare against sparse GP (GPyTorch + inducing points) after Cycle 1 accumulates real data.

---

## Python Environment

**Used:** system Python 3.14.0 at `/Library/Frameworks/Python.framework/Versions/3.14/`  
**Reason:** conda environment `igem2026` does not exist at `/opt/anaconda3/envs/igem2026`; system Python has all required packages pre-installed.  
**Packages:** torch 2.11.0, numpy 2.3.4, pandas 2.3.3, scikit-learn 1.8.0, matplotlib 3.10.7, seaborn 0.13.2

When the team sets up `igem2026` environment, add to `shared/env/environment.yml`:
```yaml
dependencies:
  - pytorch>=2.0
  - matplotlib
  - seaborn
  - scipy
  - scikit-learn
```

---

## When ELISA Arrives — Swap-In Plan

**Expected date:** ~2026-06-01  
**Format:** `08_cycle_data/outputs/cycle_0/elisa_processed.csv`

### Required label format / 要求的标签格式

The ELISA CSV must have these columns:

| Column | Type | Description |
|---|---|---|
| `rbp_id` | string | matches `03_rbp_identification` RBP id format (`EU717894.1_rbp_01`) |
| `receptor_id` | string | matches `04_protein_embedding` receptor id format (`GCF_000007145.1_tonB`) |
| `EC50_nM` | float | raw EC50 measurement in nanomolar |
| `log_EC50` | float | log10(EC50_nM) — this is the regression target |
| `cycle` | int | `0` for the first wet-lab batch |
| `notes` | string | free text; empty string if none |

Example row:
```
EU717894.1_rbp_01,GCF_000007145.1_tonB,12.5,1.097,0,replicate 3
```

### Which cells change in the notebook / notebook 中哪些单元格需要修改

**File:** `processes/01_train_deep_ensemble_synthetic.ipynb`

**Cell 3 (data loading) — replace entirely:**

```python
# WHEN REAL ELISA DATA ARRIVES: replace Cell 3 with this block
# 真实 ELISA 数据到来时：用以下代码替换第 3 单元格

CYCLE = 0
ELISA_PATH = REPO_ROOT / f'08_cycle_data/outputs/cycle_{CYCLE}/elisa_processed.csv'
EMBED_NPZ   = REPO_ROOT / '04_protein_embedding/outputs/embeddings_esm2_t6_8M_rbp.npz'
REC_NPZ     = REPO_ROOT / '04_protein_embedding/outputs/embeddings_esm2_t6_8M_receptor.npz'

elisa = pd.read_csv(ELISA_PATH)

# Load real embeddings
rbp_data = np.load(EMBED_NPZ, allow_pickle=True)
rec_data = np.load(REC_NPZ, allow_pickle=True)

rbp_emb_map = dict(zip(rbp_data['seq_ids'].tolist(), rbp_data['array']))
rec_emb_map = dict(zip(rec_data['seq_ids'].tolist(), rec_data['array']))

# Build X, y from ELISA rows (skip rows with missing embeddings)
rows_X, rows_y, rows_meta = [], [], []
for _, row in elisa.iterrows():
    if row['rbp_id'] not in rbp_emb_map or row['receptor_id'] not in rec_emb_map:
        continue
    x = np.concatenate([rbp_emb_map[row['rbp_id']], rec_emb_map[row['receptor_id']]])
    rows_X.append(x)
    rows_y.append(row['log_EC50'])   # <-- real label
    rows_meta.append({'rbp_id': row['rbp_id'], 'receptor_id': row['receptor_id'],
                      'true_y': row['log_EC50'], 'source': 'elisa_cycle_0'})

X = np.array(rows_X, dtype=np.float32)
y = np.array(rows_y, dtype=np.float32)
data_meta = pd.DataFrame(rows_meta)

print(f'Loaded {len(X)} ELISA pairs from cycle {CYCLE}')
```

**Cell 9 (model_meta.json) — update `data_source`:**
```python
'data_source': f'elisa_cycle_{CYCLE}',  # was: 'synthetic_fallback_random'
```

**Cell 10 (MANIFEST) — no change needed** (runs automatically).

### What does NOT change / 不需要修改的内容

- `ensemble.py` — architecture and training logic are unchanged
- `run_calibration.py` — calibration logic is unchanged  
- All tests — tests use random data and don't depend on real labels
- `model_meta.json` required keys — schema is unchanged

### Calibration re-check after real ELISA

After swapping in real labels, run `02_calibration_check.ipynb`. If ECE > 0.1:
1. Implement temperature scaling (snippet above)
2. Update `calibration.png` and `model_meta.json`
3. Document in a follow-up AGENT_REPORT section

---

## Notes on Synthetic Data Limitations

The Cycle 0 predictions in `outputs/cycle_0/predictions.csv` are **meaningless biologically** — they are predictions on random embeddings with random labels. Key caveats:

- `predicted_score` values (~5.0 ± 0.5) are not pKd predictions; they are artifacts of the synthetic data distribution
- `std` values reflect model uncertainty on random data; do not interpret biologically
- The ensemble IS diverse (100% of val predictions have std > 0) — this is correct
- The ECE (0.27) is expected to be higher than 0.1 on such small synthetic data

**Bottom line: the pipeline shape is correct. The biology is not. Swap in real data per the instructions above.**

---

## Git Integration Note

**Commit `e04bba5` contains 13M-line deletions** from the worktree setup script that cleared all non-module-06 files before agent work started. These deletions are an artifact of the worktree isolation mechanism, not intentional.

When merging this branch back to `main`:
1. **Do NOT `git merge`** — it will delete all other modules' files from main.
2. **Use cherry-pick or subtree strategies:**
   ```bash
   # Option A: cherry-pick only the 06_uncertainty_model/ changes
   git checkout main
   git checkout agent-06-uncertainty-model -- 06_uncertainty_model/
   git commit -m "Merge Module 06 outputs from agent-06-uncertainty-model"
   ```
3. Or use `git diff --diff-filter=A` to extract only added files.

This issue was caused by pre-staged worktree cleanup being included in the first commit. All 06_uncertainty_model/ code is correct and independent of the deletion artifact.

---

## File Tree

```
06_uncertainty_model/
├── AGENT_REPORT.md              ← this file
├── inputs/                      ← (empty; Module 04/05 outputs not yet available)
├── outputs/
│   ├── MANIFEST.csv
│   └── cycle_0/
│       ├── ensemble_member_0.pt
│       ├── ensemble_member_1.pt
│       ├── ensemble_member_2.pt
│       ├── ensemble_member_3.pt
│       ├── ensemble_member_4.pt
│       ├── predictions.csv
│       ├── calibration.png
│       ├── training_log.json
│       └── model_meta.json
└── processes/
    ├── ensemble.py                     ← core architecture (importable)
    ├── synthetic_data.py               ← synthetic data generator
    ├── run_cycle0.py                   ← CLI training runner
    ├── run_calibration.py              ← calibration plot generator
    ├── 01_train_deep_ensemble_synthetic.ipynb
    ├── 02_calibration_check.ipynb
    └── tests/
        ├── test_schema.py
        ├── test_smoke.py
        ├── test_diversity.py
        └── test_calibration.py
```
