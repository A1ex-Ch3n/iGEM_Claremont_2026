# Module 05 — Structure Prediction

**Status:** ✅ Sprint deliverable complete (2026-05-07 overnight build)

**Role:** 3D structure prediction for phage RBPs and host receptors, plus
binding-confidence estimation for (RBP, receptor) pairs. Feeds `affinity_priors.csv`
into Module 06 (deep ensemble, Cycle 0 cold start).

---

## Tools

| Tool | Version | Role | Status |
|---|---|---|---|
| **Boltz-2** | 2.0.3 | Structure + ipTM confidence (primary) | ✅ Installed, 1 pair running |
| **AlphaFold 3** | — | High-quality structures (requires weight access) | ⏳ Waiting for Google approval |
| **ColabFold** | — | AF2-Multimer fallback | Available via pip |

---

## Quick Start

```bash
# Install dependencies / 安装依赖
pip3 install -r 05_structure_prediction/requirements.txt

# Fix SSL on macOS Python 3.14 (one-time) / macOS SSL 修复（一次性）
bash "/Applications/Python 3.14/Install Certificates.command"

# Run Boltz-2 on a single pair (CPU, slow) / 在 CPU 上运行单对（较慢）
boltz predict 05_structure_prediction/inputs/boltz_input_*.fasta \
    --out_dir 05_structure_prediction/outputs/boltz2/ \
    --accelerator cpu --use_msa_server --model boltz2

# Run tests / 运行测试
python3 -m pytest 05_structure_prediction/processes/tests/ -v
```

---

## Directory Structure

```
05_structure_prediction/
├── AGENT_TODO.md          ← original sprint task spec
├── AGENT_REPORT.md        ← this sprint's decisions + Laguna runbook
├── README.md              ← this file
├── requirements.txt       ← Python dependencies
├── inputs/
│   ├── EU717894.1_rbp_sequences.faa                    ← phiL7 RBP candidates
│   ├── xcc_receptors_minimal.faa                       ← Xcc receptor sequences
│   ├── boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB.fasta  ← pair 1
│   └── boltz_input_EU717894.1_rbp_01__GCF_000007145.1_exbB.fasta  ← pair 2
├── outputs/
│   ├── MANIFEST.csv
│   ├── affinity_priors.csv                             ← fed to Module 06
│   ├── boltz2/<rbp_id>__<receptor_id>/
│   │   ├── model.pdb (after GPU run)
│   │   └── affinity.json
│   └── af3/<seq_id>/ (empty until AF3 weights approved)
└── processes/
    ├── 01_run_boltz2.ipynb  ← primary workflow
    ├── 02_run_af3.ipynb     ← AF3 setup documentation + ColabFold fallback
    └── tests/
        ├── test_schema.py   ← INTERFACE.md column validation
        ├── test_smoke.py    ← mock-based end-to-end pipeline test
        └── test_sanity.py   ← biology-level sanity checks
```

---

## Inputs

| Source | Path | Status |
|---|---|---|
| phiL7 RBP sequences | `inputs/EU717894.1_rbp_sequences.faa` | ✅ Manually fetched |
| Xcc receptor sequences | `inputs/xcc_receptors_minimal.faa` | ✅ Manually fetched |
| Module 03 output | `../03_rbp_identification/outputs/EU717894.1_rbp_sequences.faa` | ⏳ Pending |
| Module 04 output | `../04_protein_embedding/inputs/xcc_receptors.faa` | ⏳ Pending |

---

## Outputs

| File | Description | Schema |
|---|---|---|
| `outputs/affinity_priors.csv` | Aggregated binding confidence per pair | See INTERFACE.md §Module 05 |
| `outputs/boltz2/<pair>/model.pdb` | Predicted complex structure | PDB format |
| `outputs/boltz2/<pair>/affinity.json` | Confidence metrics (ipTM, pTM, pLDDT) | JSON |
| `outputs/MANIFEST.csv` | Checksums + metadata per output file | INTERFACE.md universal |

---

## Notes

- **Boltz-2 affinity head is for protein–small-molecule, not protein–protein.**
  For RBP–receptor (protein–protein) pairs, `predicted_dG_kcal_mol` and `predicted_pKd`
  are NaN. Use `confidence` (= ipTM) as the binding signal for Module 06.
- **GPU recommended for production.** CPU runs take 10-30 min per pair with reduced
  step counts. See `AGENT_REPORT.md §Laguna Runbook` for `sbatch` templates.
