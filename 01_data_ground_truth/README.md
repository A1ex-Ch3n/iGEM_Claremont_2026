# Step 1 — Data & Ground Truth

## Purpose
Download phage and host genomes from NCBI and construct the phage–host interaction matrix (ground truth labels for the classifier).

## Inputs (`inputs/`)
- `complete_phage_data_Xanthomonas.csv` — manual seed list of *Xanthomonas*-infecting phages (from Olivia)
- `complete_phage_data_Pseudomonas.csv` — Pseudomonas phage seed list

## Processes (`processes/`)
| Script | What it does |
|--------|-------------|
| `fetch_positive_pairs.py` | NCBI Entrez search → positive affinity (=1) phage–host pairs for *Xanthomonas* phages |
| `fetch_negative_pairs.py` | Generates Affinity=0 negatives via 3 modules: cross-genus, GenBank comment mining, pathovar inference |
| `download_genomes.py` | Batch-downloads nucleotide FASTAs for phage accessions in the interaction matrix |

## Outputs (`outputs/`)
| Path | Description |
|------|-------------|
| `interaction_matrix/final_interaction_matrix.csv` | **Canonical** 2-D phage × host affinity matrix (1/0/NaN) |
| `interaction_matrix/phage_host_matrix_with_ids.csv` | Long-format pairs with host RefSeq accessions |
| `interaction_matrix/phage_host_matrix_raw.csv` | Raw long-format pairs (host names only) |
| `interaction_matrix/xanthomonas_phages_accession_list.csv` | Phage NCBI accession list |
| `negative_samples/negative_cross_genus.csv` | Module A cross-genus negatives |
| `negative_samples/negative_pv_inference.csv` | Module C same-species different-pathovar negatives |
| `negative_samples/negative_data_combined.csv` | All negatives combined |
| `phage_genomes/xanthomonas_pool/` | 777 phage genome FASTAs (gitignored; tracked via MANIFEST) |
| `phage_genomes/NC_*.fasta` | Sample genomes for testing annotation pipeline |
| `host_genomes/ncbi_dataset/` | *Xanthomonas* host genome FASTAs from NCBI (gitignored) |
| `download_summary.csv` | Download success/error log |

## Owners
Sarah (interaction matrix, positive & negative fetching), Weitao (genome download)

## Gaps / Next steps
- `fetch_negative_pairs.py` should write negatives directly into `interaction_matrix/final_interaction_matrix.csv` (currently separate CSVs)
- Add `MANIFEST.csv` listing all 777 accessions in `phage_genomes/xanthomonas_pool/` for git tracking
