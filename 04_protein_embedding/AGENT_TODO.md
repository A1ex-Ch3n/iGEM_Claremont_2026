# AGENT TODO — Module 04 / `04_protein_embedding/`

## Read first (mandatory)

1. `/INTERFACE.md` §Module 04 — your output spec. Pay close attention to the `.npz` key requirements.
2. `/LAGUNA.md` — for the production-scale ESM-2 3B run.
3. `/CLAUDE.md` — notebook-first, bilingual.
4. Invoke `superpowers:test-driven-development` and `superpowers:verification-before-completion`.

## Your scope

Convert protein sequences → numerical embeddings using ESM-2. Two regimes:
- **CPU smoke (tonight, default)**: ESM-2 8M (`esm2_t6_8M_UR50D`). Tiny model, runs on Mac CPU in seconds.
- **Laguna production (tomorrow, documented but not executed tonight)**: ESM-2 650M / 3B. Requires GPU. SLURM template in `LAGUNA.md`.

You read RBP / receptor sequences from `03_rbp_identification/outputs/<phage_acc>_rbp_sequences.faa`. You write `.npz` files to `04_protein_embedding/outputs/`.

You are also responsible for **extracting receptor sequences** from `02_annotation/outputs/host_proteins/<acc>.faa` (TonB / ExbB / ExbD1 / ExbD2 by gene name regex) and embedding them too — these will be the receptor side of the (RBP, receptor) prediction pairs in Module 06.

## Goal (definition of done)

By morning, `04_protein_embedding/` contains:

1. ✅ `processes/01_embed_esm2.ipynb` — main embedding notebook with model selector.
2. ✅ `processes/02_extract_receptors.ipynb` — pulls TonB / ExbB / ExbD from host proteome by name.
3. ✅ `outputs/embeddings_esm2_t6_8M_phiL7_rbps.npz` — phiL7 RBP candidates embedded.
4. ✅ `outputs/embeddings_esm2_t6_8M_xcc_receptors.npz` — Xcc receptor proteins embedded.
5. ✅ `outputs/embedding_index.csv` — manifest of all embedding files per INTERFACE.
6. ✅ `outputs/MANIFEST.csv`.
7. ✅ ≥3 passing tests.
8. ✅ `AGENT_REPORT.md` — including a "Laguna runbook" subsection with the exact `sbatch` command for the 650M run tomorrow.

## Setup

```bash
cd /path/to/agent-04-protein-embedding
conda activate igem2026

# fair-esm (Meta's official ESM library)
pip install fair-esm

# Test it loads
python -c "import esm; print(esm.__version__)"

# Pre-download ESM-2 8M weights (~ 30 MB)
python -c "
import esm
model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
print('ESM-2 8M loaded; embed_dim =', model.embed_dim)
"
```

For tomorrow's Laguna run, the env there will need:
```bash
pip install torch --index-url https://download.pytorch.org/whl/cu121  # CUDA 12.1 build
pip install fair-esm
```

## Step-by-step plan

### Step 1 — Build `01_embed_esm2.ipynb` (90 min)

Cells:

1. **Markdown** — Title + bilingual purpose. Cite Lin et al. 2023.
2. **Code** — Imports (`esm`, `torch`, `numpy`, `pandas`, `Bio.SeqIO`), paths, version printouts. **Set torch + numpy seeds (e.g. 42).**
3. **Markdown** — Model menu (bilingual): table of available ESM-2 sizes, embed dims, GPU memory needs, recommended use case.

   | Model | Layers | Embed dim | Params | GPU memory | Use for |
   |---|---|---|---|---|---|
   | `esm2_t6_8M_UR50D` | 6 | 320 | 8M | < 1 GB | CPU smoke / unit tests |
   | `esm2_t12_35M_UR50D` | 12 | 480 | 35M | 2 GB | Laptop GPU |
   | `esm2_t30_150M_UR50D` | 30 | 640 | 150M | 4 GB | Single small GPU |
   | `esm2_t33_650M_UR50D` | 33 | 1280 | 650M | 12 GB | Production default |
   | `esm2_t36_3B_UR50D` | 36 | 2560 | 3B | 24 GB | Laguna A100 |

4. **Code** — Helper: `load_model(model_name) → (model, alphabet)`. Use `esm.pretrained.<model_name>()`. Move to GPU if available.
5. **Code** — Helper: `embed_sequences(model, alphabet, seqs: List[Tuple[str, str]], pooling: str = "mean", batch_size: int = 8) → Dict`. Returns `{seq_ids, array, lengths, meta}` per INTERFACE §Module 04.
   - For pooling="mean": average over residues (excluding BOS/EOS tokens).
   - For pooling="per_residue": keep full residue tensor, pad to max length, return `lengths` for masking.
6. **Code** — Helper: `save_npz(out_path, embedding_dict)`. Per INTERFACE: `seq_ids`, `array`, `lengths`, `meta` (JSON-string).
7. **Markdown** — Sample-then-batch.
8. **Code** — Sample run: load `03_rbp_identification/outputs/EU717894.1_rbp_sequences.faa`, embed top 5 RBP candidates with ESM-2 8M, mean-pooled. Inspect output shape.
9. **Code** — Sanity assertion: `array.shape == (5, 320)`; `lengths` are reasonable; `meta` parses to JSON.
10. **Code** — Save to `outputs/embeddings_esm2_t6_8M_phiL7_rbps.npz`.
11. **Code** — Update `outputs/embedding_index.csv` and `outputs/MANIFEST.csv`.
12. **Markdown** — Laguna handoff: how to swap `esm2_t6_8M_UR50D` → `esm2_t33_650M_UR50D` for production. Reference `LAGUNA.md` Template A.

### Step 2 — Build `02_extract_receptors.ipynb` (60 min)

Module 06 needs receptor embeddings (TonB / ExbB / ExbD1 / ExbD2 from Xcc) to predict RBP-receptor binding. This notebook:

1. **Markdown** — Bilingual purpose. Cite Wang et al. 2003 for receptor identification in phiL7.
2. **Code** — Read `02_annotation/outputs/host_proteins/GCF_000007145.1.faa` (or whichever Xcc file Module 02 produced).
3. **Code** — For each receptor gene name, find the matching protein via:
   - First try: exact gene-name match in FASTA header (if header contains `gene=tonB`).
   - Fallback: BLAST-like local alignment of known receptor sequences (UniProt Q8P5W2 for Xcc TonB; you can fetch from UniProt as a one-off).
   - Fallback: HMM scan with TonB Pfam domain (PF03544 for TonB-dependent receptor plug domain).
4. **Code** — Write extracted receptor sequences to a small FASTA: `inputs/xcc_receptors.faa`.
5. **Code** — Embed these receptors with the same `embed_sequences` helper from notebook 01.
6. **Code** — Save to `outputs/embeddings_esm2_t6_8M_xcc_receptors.npz`.

If Module 02's Xcc annotation isn't ready, write the helper to read from `00_raw_data/bacteria/GCF_000007145.1/proteins.faa` directly (acceptable fallback because the file structure is identical even if header style differs).

### Step 3 — Sanity check determinism (15 min)

ESM-2 should be deterministic. Add a final cell:
- Embed the same 3 sequences twice.
- Assert `np.allclose(emb_run1, emb_run2)`.

If non-deterministic (e.g., due to torch dropout being on), set `model.eval()` and `torch.use_deterministic_algorithms(True)`.

### Step 4 — Tests (60 min)

`processes/tests/`:
- `test_schema.py` — Open the saved `.npz`, verify required keys (`seq_ids`, `array`, `lengths`, `meta`); verify `meta` is parseable JSON with `model`, `pooling`, `date_utc` keys.
- `test_smoke.py` — Embed 2 random 50-aa sequences with 8M model, assert output shape (2, 320), assert no NaN.
- `test_determinism.py` — Same input twice → same output.
- `test_sanity.py` — Embed phiL7 RBP candidate 1 (assume ~412 aa), assert mean-pooled output is finite, assert L2 norm ∈ [0.1, 100].

Run: `pytest 04_protein_embedding/processes/tests/ -v`

### Step 5 — Commit + report

Commits:
- `04_protein_embedding: ESM-2 install + smoke test`
- `04_protein_embedding: embed_esm2 notebook (CPU 8M default)`
- `04_protein_embedding: extract_receptors notebook`
- `04_protein_embedding: embeddings for phiL7 RBPs + Xcc receptors`
- `04_protein_embedding: tests`
- `04_protein_embedding: AGENT_REPORT (incl. Laguna runbook)`

## References (cite in notebook markdown cells)

- **ESM-2 (the model)**: Lin, Z. et al. (2023) "Evolutionary-scale prediction of atomic-level protein structure." *Science* 379(6637):1123-1130. DOI: 10.1126/science.ade2574. **PRIMARY METHODOLOGY REFERENCE.**
- **PLM-interact (transfer learning prior we will explore later)**: Liu, Y. et al. (2025) "PLM-interact: extending protein language models." *Nat Commun*.
- **PLM transfer to phage proteins**: Hie, B.L. et al. (2024) "Efficient evolution of human antibodies from general protein language models." *Nat Biotechnol*. (Demonstrates ESM transfer to non-canonical protein families.)
- **Receptor identification for phiL7**: Wang, W.-T. et al. (2003) "Involvement of *tonB-exbBD1D2* operon in infection of *Xanthomonas campestris* phage phiL7." *Molecular Microbiology* 49:1097-1110. DOI: 10.1046/j.1365-2958.2003.03565.x. (Source for which receptor genes to extract.)
- **TonB-dependent receptor Pfam**: Pfam PF03544 (TonB plug domain) — https://www.ebi.ac.uk/interpro/entry/pfam/PF03544/

## Anti-goals

- ❌ Don't try to run ESM-2 3B on Mac CPU; it will OOM or take hours. Document the Laguna path; don't run it tonight.
- ❌ Don't use HuggingFace's ESM-2 wrapper unless `fair-esm` install fails (slight numerical drift between the two implementations; we standardize on `fair-esm`).
- ❌ Don't ignore padding tokens in mean-pooling — that biases short sequences toward zero. Use `lengths` to mask.
- ❌ Don't push to remote.

## Time budget

~4 hours. Embedding the actual phiL7 RBPs is fast (< 30 sec); most time goes to writing reusable helpers and tests.

## If stuck

- `fair-esm` install fails on Apple Silicon → use `transformers` library: `pip install transformers && from transformers import AutoTokenizer, AutoModel; AutoModel.from_pretrained("facebook/esm2_t6_8M_UR50D")`. Document the drift in AGENT_REPORT.
- Module 02's annotation isn't done → mock with a tiny FASTA (3 random 100-aa sequences) so your pipeline is at least testable.
- Receptor extraction by gene-name fails (FASTA headers don't contain gene names) → BLAST-search a known TonB sequence (UniProt Q8P5W2) against the host proteome using `Bio.Blast` or `diamond blastp`.
