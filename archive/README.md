# Archive

Contains files that have been superseded or reorganized but are preserved for reference.

| File / Folder | Original location | Why archived |
|--------------|-------------------|-------------|
| `legacy_master_pipeline.py` | `olivia/04_dry_lab/pipeline/master_pipeline.py` | Mixed Step 1 (NCBI download) and Step 2 (Prodigal annotation on phage genomes). Phage annotation is now done by PHANOTATE (`02_annotation/processes/phage_phanotate/`); host annotation by Prodigal (`02_annotation/processes/host_prodigal/`). Download logic lives in `01_data_ground_truth/processes/download_genomes.py`. |
| `tasks/0428_linear_regression_cn_plan.md` | `TASK/0428_linear_regression/cn_plan.md` | Canonical 6-factor regression spec. Content has been translated and expanded into `03_feature_weighting/README.md`. Kept here in its original Mandarin form. |
