# Step 6 — In-Silico Knockdown (fastISM)

## Purpose
Use fastISM (In-Silico Mutagenesis) to perform sequence perturbation on the highest-probability phage–host pairs identified in Step 5. Generates a binding sensitivity heatmap showing which nucleotide positions in key phage proteins most influence predicted infectivity.

## Inputs
- Top-scoring phage–host pairs from `05_predictive_modeling/outputs/`
- Corresponding phage protein sequences from `02_annotation/outputs/phage_proteins/`

## Processes
- `processes/fastism_knockdown.py` — to be written; wraps fastISM library
- Outputs sensitivity score per position per high-priority protein

## Outputs
- `outputs/binding_sensitivity_heatmap_<phage>_<host>.png`
- `outputs/sensitivity_scores.csv`

## Owner
TBD

## Status
Not started. Blocked on Step 5 (need final probability matrix to select top pairs).
