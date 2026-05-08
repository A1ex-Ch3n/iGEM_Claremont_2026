# AGENT TODO — Module 05 / `05_structure_prediction/`

## Read first (mandatory)

1. `/INTERFACE.md` §Module 05 — your output spec.
2. `/LAGUNA.md` — production GPU runs go here; tonight is local CPU smoke + setup only.
3. `/CLAUDE.md` — notebook-first, bilingual.
4. Invoke `superpowers:test-driven-development` and `superpowers:verification-before-completion`.

## Your scope

3D structure prediction for phage RBPs and host receptors, plus zero-shot binding-affinity prediction for (RBP, receptor) pairs. Two tools, deliberately complementary:

- **Boltz-2** (Passaro et al. 2025) — primary tonight. Open-source, easy install, has explicit affinity head trained on PDBbind. Outputs binding ΔG estimates we feed to Module 06 as cold-start synthetic prior.
- **AlphaFold 3** (Abramson et al. 2024) — secondary. Higher-quality static structures but the official install requires model-weight access approval, which may not be granted by tonight. **Document the install path; do NOT block on actually running AF3 tonight.** Fallback: ColabFold for structures, accept that Boltz-2 covers affinity.

You read sequences from `03_rbp_identification/outputs/<phage_acc>_rbp_sequences.faa` and `04_protein_embedding/inputs/xcc_receptors.faa` (if Module 04 produced it; else extract directly). You write structures + affinity tables to `05_structure_prediction/outputs/`.

## Goal (definition of done)

By morning, `05_structure_prediction/` contains:

1. ✅ `processes/01_run_boltz2.ipynb` — runs Boltz-2 on 1-2 (RBP, receptor) pairs as proof of concept.
2. ✅ `processes/02_run_af3.ipynb` — **install + setup documentation only** if AF3 install non-trivial; with a clear fallback to ColabFold runnable.
3. ✅ `outputs/boltz2/<rbp_id>__<receptor_id>/{model.pdb, affinity.json}` for at least 1 pair (e.g., phiL7's top RBP × Xcc TonB).
4. ✅ `outputs/affinity_priors.csv` per INTERFACE schema (even if only 1-2 rows).
5. ✅ `outputs/MANIFEST.csv`.
6. ✅ ≥3 passing tests.
7. ✅ `AGENT_REPORT.md` — must include "Laguna runbook" with `sbatch` commands for batch Boltz-2 + AF3 runs (templates in `LAGUNA.md`).

## Setup

```bash
cd /path/to/agent-05-structure-prediction
conda activate igem2026

# Boltz-2 — open source from Passaro lab
pip install boltz -U

# Verify it works
boltz --help

# Pre-download Boltz-2 weights (~ 1-2 GB)
# First run will auto-download; force pre-download by running on a tiny dummy:
mkdir -p /tmp/boltz_test
echo -e ">test_seq\nMKLLILTCLVAVALARPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNELSKDIGSESTEDQAMEDIKQMEAESISSSEEIVPNSVEQKHIQKEDVPSERYLGYLEQLLRLKKYKVPQLEIVPNSAEERLHSMKEGIHAQQKEPMIGVNQELAYFYPELFRQFYQLDAYPSGAWYYVPLGTQYTDAPSFSDIPNPIGSENSEKTTMPLW" > /tmp/boltz_test/seq.fasta
# Note: full Boltz-2 invocation pattern documented in notebook

# AlphaFold 3 — try official install but DO NOT BLOCK ON IT
# Official: https://github.com/google-deepmind/alphafold3
# Requires model weight access (apply at https://docs.google.com/forms/d/e/1FAIpQLSfWZAgo1aYk0O4MuAXZj8sLCt4uemNsnflXLdU5DG_eJ45UFw/viewform)
# Fallback: ColabFold (bundled AF2 + AF3-ish multimer)
pip install colabfold-batch
# OR locally via Docker: docker pull deepmind/alphafold3
```

If Boltz-2 install fails, document and fall back to AlphaFold-Multimer's ipTM as an affinity proxy (poor man's version; document limitations).

## Step-by-step plan

### Step 1 — Build `01_run_boltz2.ipynb` (120 min)

Cells:

1. **Markdown** — Title + bilingual purpose. Cite Passaro 2025.
2. **Code** — Imports + paths + version printouts.
3. **Markdown** — Method explanation: Boltz-2 architecture (joint structure + affinity prediction). Why we use it as a synthetic prior, not ground truth.
4. **Code** — Helper: `prepare_pair_input(rbp_seq, receptor_seq, output_dir) → str`. Writes Boltz-2 input YAML/JSON in the format the CLI expects.
5. **Code** — Helper: `run_boltz2(input_path, output_dir) → dict`. Wraps the CLI (`subprocess.run(["boltz", "predict", ...])`). Returns parsed metadata: model PDB path, affinity JSON path, runtime.
6. **Code** — Helper: `parse_affinity_json(json_path) → dict`. Returns `{predicted_dG_kcal_mol, predicted_pKd, confidence}` per INTERFACE.
7. **Markdown** — Sample-then-batch.
8. **Code** — Sample run: 1 pair = phiL7 top RBP + Xcc TonB. Use ESM-2-extracted sequences from Module 04's `inputs/xcc_receptors.faa` if present; else hardcode Xcc TonB UniProt Q8P5W2.
9. **Code** — Sanity assertion: PDB output exists and is parseable; affinity JSON has all required keys; predicted_pKd is finite.
10. **Code** — Append row to `outputs/affinity_priors.csv`.
11. **Code** — (Optional) Run a second pair: phiL7 top RBP × Xcc ExbB to demonstrate batching pattern.

**CRITICAL**: Boltz-2 on CPU is slow (10-30 min per pair). Time-box the sample run to 30 min total. If it doesn't finish, save the input and document the Laguna sbatch path; don't burn the whole night on one prediction.

### Step 2 — Build `02_run_af3.ipynb` (60 min)

This is mostly documentation. Structure:

1. **Markdown** — Title + bilingual purpose. Cite Abramson 2024.
2. **Markdown** — Install instructions (3 paths):
   - Path A: Official AF3 from Google (requires weight access).
   - Path B: AlphaFold-Multimer (older, bundled with ColabFold; sufficient for ipTM-based affinity proxy).
   - Path C: ColabFold via API.
3. **Code** — Whichever path you got working: a callable function `run_af3_or_fallback(seq_a, seq_b, output_dir) → dict` returning structure path + ipTM/pLDDT.
4. **Code** — If got working: predict 1 sample pair. Else: write a stub that raises `NotImplementedError` with a clear message pointing to Laguna setup.
5. **Markdown** — "How to run for real" — `LAGUNA.md` Template B.

### Step 3 — Receptor sequence acquisition fallback (30 min)

If Module 04 hasn't produced `inputs/xcc_receptors.faa`, you need to get receptor sequences yourself for the sample Boltz-2 run. Quickest path:

```python
from Bio import Entrez
Entrez.email = "[email protected]"
# Xcc TonB UniProt Q8P5W2 → corresponding NCBI protein
# Or pull from your locally-fetched Xcc proteome by gene name regex
```

Save to `inputs/xcc_receptors_minimal.faa` for the 4 receptors (TonB, ExbB, ExbD1, ExbD2). This is the same data Module 04 needs, so document for cross-module sync.

### Step 4 — `affinity_priors.csv` (15 min)

Aggregate output table per INTERFACE §Module 05. Even with just 1-2 rows from tonight's smoke run, the columns must be exactly right because Module 06 reads this file.

### Step 5 — Tests (45 min)

`processes/tests/`:
- `test_schema.py` — Open `affinity_priors.csv`, assert columns. Open one Boltz-2 affinity JSON, assert keys.
- `test_smoke.py` — Mock the `run_boltz2` helper to return a fake JSON, assert the parsing/aggregation pipeline produces a valid CSV row.
- `test_sanity.py` — `outputs/boltz2/<rbp>__<recv>/model.pdb` exists, has `ATOM` records, parses with Bio.PDB.PDBParser.

Run: `pytest 05_structure_prediction/processes/tests/ -v`

### Step 6 — Commit + report

Commits:
- `05_structure_prediction: Boltz-2 install + smoke test`
- `05_structure_prediction: Boltz-2 notebook + 1 pair prediction`
- `05_structure_prediction: AF3 setup documentation`
- `05_structure_prediction: affinity_priors + manifest`
- `05_structure_prediction: tests`
- `05_structure_prediction: AGENT_REPORT (with Laguna runbook)`

## References (cite in notebook markdown cells)

- **Boltz-2**: Passaro, S. et al. (2025) "Boltz-2: All-atom protein folding with binding affinity prediction." *bioRxiv*. https://github.com/jwohlwend/boltz. **PRIMARY TONIGHT.**
- **AlphaFold 3**: Abramson, J. et al. (2024) "Accurate structure prediction of biomolecular interactions with AlphaFold 3." *Nature* 630:493-500. DOI: 10.1038/s41586-024-07487-w.
- **ColabFold (fallback)**: Mirdita, M. et al. (2022) "ColabFold: making protein folding accessible to all." *Nature Methods* 19:679-682.
- **AlphaFold-Multimer (poor man's affinity via ipTM)**: Evans, R. et al. (2021) "Protein complex prediction with AlphaFold-Multimer." *bioRxiv* 2021.10.04.463034.
- **Receptor identification context**: Wang, W.-T. et al. (2003) *Molecular Microbiology* 49:1097-1110. (Source for which Xcc receptors to model.)
- **PDB format spec**: https://www.wwpdb.org/documentation/file-format
- **mmCIF format spec**: https://mmcif.wwpdb.org/

## Anti-goals

- ❌ Don't waste GPU/CPU cycles trying to predict structures of all 5 phiL7 RBP candidates tonight — 1 pair is enough to validate the pipeline. Batch on Laguna tomorrow.
- ❌ Don't trust Boltz-2 affinities as ground truth — they're a prior to seed Module 06. Module 06's deep ensemble will be retrained on ELISA data once available.
- ❌ Don't block the whole night waiting for AF3 install. Boltz-2 + ColabFold covers our needs.
- ❌ Don't push to remote.

## Time budget

~4-5 hours. Boltz-2 install + 1 prediction is the main cost. Skip AF3 if it eats > 1 hour.

## If stuck

- Boltz-2 install fails on macOS → check `brew install pkg-config` then retry; or use Docker (`docker pull jwohlwend/boltz`).
- Boltz-2 prediction takes > 30 min on CPU → kill it, save the input config, document in AGENT_REPORT for Laguna.
- AF3 weights not granted → it's expected; fall back to ColabFold and note in AGENT_REPORT that AF3 requires manual onboarding (Alex / PI to apply).
- Module 04's receptor extraction not done → use the Bio.Entrez fallback (described in Step 3) and write the receptors yourself.
