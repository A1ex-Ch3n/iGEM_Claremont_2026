# Olivia — 04_dry_lab Contents

This directory contains Olivia's dry-lab work for iGEM Claremont 2026: an
integrated phage genomics pipeline that extracts proteomic/genomic features
from phage genomes and prepares them for an infectivity-prediction model.

## `pipeline/` — analysis code and documentation
- `master_pipeline.py` — top-level driver that orchestrates the full workflow.
- `batch_prodigal.py` — batch ORF/CDS prediction across genomes using Prodigal.
- `compute_pI_acidity.py` — computes per-genome isoelectric point (pI) and
  acidity features from predicted proteins.
- `plot_pI_acidity_density.py` — generates the pI/acidity density plot.
- `Integrated_Phage_Genomics_Pipeline_Guide.md` / `.pdf` — pipeline guide.
- `Prodigal_Manual.md` / `.pdf` — reference manual for Prodigal usage.

## `data/` — input datasets
- `complete_phage_data_Pseudomonas.csv` — labeled Pseudomonas phage dataset.
- `complete_phage_data_ Xanthomonas.csv` — labeled Xanthomonas phage dataset.
- `complete_phage_data_ Xanthomonas.csv.fasta` — concatenated FASTA for
  Xanthomonas phages.
- `ncbi_dataset/data/` — per-accession genome folders downloaded from NCBI
  (e.g. `AB720063`, `AP008979`, `KX181651`, …).

## `results/` — pipeline outputs
- `f02_pI_acidity_per_genome.csv` — per-genome pI and acidity feature table.
- `f02_pI_acidity_density.png` — density plot of pI/acidity across genomes.

## `step3_weighing/` — modeling step
- `README.md` — describes the linear-regression weighing model that consumes
  the per-genome features to predict phage infectivity.

## `workflow/` — diagrams
- `dry_lab_workflow_chart.jpeg` — visual chart of the dry-lab workflow.
