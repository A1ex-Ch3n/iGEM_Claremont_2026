# Step 2 — Standardized Annotation

## Purpose
Predict open reading frames and protein sequences for all phage genomes (PHANOTATE) and bacterial host genomes (Prodigal). Also runs pharokka for richer phage functional annotation.

## Inputs
Points to `01_data_ground_truth/outputs/phage_genomes/` and `01_data_ground_truth/outputs/host_genomes/`.

## Processes (`processes/`)
| Sub-folder | Tool | What it does |
|------------|------|-------------|
| `phage_phanotate/batch_phanotate.py` | PHANOTATE | ORF prediction on phage FASTAs → coordinates + `.faa` protein sequences |
| `host_prodigal/batch_prodigal.py` | Prodigal (pyrodigal) | ORF prediction on host FASTAs → `proteins.faa` + `genes.gff` |
| `pharokka/` | pharokka | Full annotation bundle for phages (PHROG hits, tRNA, CRISPR) |

**Note on tool split:** PHANOTATE is tuned for phage genomes (no stop-codon bias assumptions); Prodigal/pyrodigal is for bacterial hosts. Do not swap them.

## Outputs (`outputs/`)
| Path | Description |
|------|-------------|
| `phage_proteins/<acc>.faa` | Translated phage protein sequences (gitignored) |
| `host_proteins/<acc>.faa` | Translated host protein sequences (gitignored) |
| `pharokka_runs/<acc>/` | Full pharokka output bundle per phage (gitignored) |
| `phage_phanotate_coords/<acc>.txt` | Raw PHANOTATE coordinate files |

**5-sample pharokka results** for AB720063.2, AB720064.1, AP008979.1, EU717894.1, JN882298.1 are already in `pharokka_runs/` from the 2025-04-28 annotation sprint.

## Owners
Weitao (PHANOTATE/phage), Olivia (Prodigal/host)

## Gaps / Next steps
- `batch_phanotate.py` is incomplete — `.faa` translation loop was truncated. Needs finishing (parse PHANOTATE tab-delimited coords, translate via Biopython).
- Scale pharokka to the full 777-genome pool.
- `archive/legacy_master_pipeline.py` used Prodigal on phage genomes — do not use; it is kept for reference only.
