# Step 4 — Protein Embedding (ESM-2)

## Purpose
Embed high-priority phage and host protein sequences using the ESM-2 transformer model. Embeddings replace raw sequence as input to the Step 5 classifier, capturing structural and evolutionary context.

## Inputs
- Top-ranked phage proteins from `03_feature_weighting/outputs/proteins_weights.csv` (CRITICAL and HIGH priority rows)
- Corresponding `.faa` files from `02_annotation/outputs/phage_proteins/` and `host_proteins/`

## Processes
- ESM-2 wrapper script (to be written)
- Recommended model: `esm2_t33_650M_UR50D` or `esm2_t36_3B_UR50D` depending on GPU budget

## Outputs
- `outputs/<protein_id>.npy` — per-protein 1280-dim (or 2560-dim) embedding vectors (gitignored)
- `outputs/embedding_index.csv` — maps accession + protein_id to `.npy` filename

## Owner
TBD

## Status
Not started. Blocked on Step 3 (need ranked protein list) and Step 2 (need full `.faa` pool).
