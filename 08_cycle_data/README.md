# 08 — Cycle Data

Per-cycle wet lab data + model checkpoints + analysis artifacts. The single source of truth for what was measured when, by whom, and what the model knew at each cycle.

## Folder structure

Within `outputs/`:
```
outputs/
├── cycle_0/
│   ├── elisa_raw.csv              # raw ELISA OD450 readings
│   ├── elisa_processed.csv        # 4PL-fit EC50 / apparent Kd
│   ├── variants.csv               # cloned constructs (sequence, plate well, batch)
│   ├── controls.csv               # T7 gp17 positive control + BSA negative
│   ├── notes.md                   # batch-specific notes (failures, anomalies)
│   └── model_state/               # snapshot of Module 06 model at end of cycle 0
├── cycle_1/
│   └── ... (same structure)
└── cycle_2/
    └── ... (same structure)
```

## What goes in `inputs/`

Symbolic links or downstream-ready exports from Benchling (when implemented). For now, manual CSV uploads from wet lab team.

## Data quality requirements

Every cycle must have:
- ≥ 3 technical replicates per (variant, receptor) condition
- T7 gp17 positive control on every plate
- BSA-only negative control on every plate
- WT-RBP at fixed concentration for inter-plate normalization

## Versioning

All cycle data is version-controlled via git. Model checkpoints versioned via MLflow with run IDs cross-referenced in `model_state/mlflow_run_id.txt`.

## Status

⏳ Not started — populates as cycles execute (Cycle 0 starts ~6/1).
