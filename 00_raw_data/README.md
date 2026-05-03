# 00_raw_data — Raw Genome Database

This folder is the single source of truth for all raw biological sequence data used in the iGEM 2026 phage–host infectivity pipeline. Every downstream step (annotation, feature engineering, ML) reads from here.

---

## Folder Structure

```
00_raw_data/
├── README.md                    ← this file
├── data_needs.md                ← detailed explanation of what files are needed and why
├── bacteria_list.csv            ← all bacteria with accessions, GCF mappings, and notes
├── phage_list.csv               ← all 777 phages with source tags (canonical / negative_pool)
│
├── bacteria/                    ← one subdirectory per bacterium
│   └── <seq_accession>/         ← named after the original sequence accession from the matrix
│       ├── genome.fna           ← whole-genome nucleotide sequence
│       ├── proteins.faa         ← NCBI-annotated protein sequences
│       ├── genes.gff            ← gene coordinates and annotation
│       └── cds.fna              ← CDS nucleotide sequences (for CAI, Factor 6)
│
├── phage/                       ← one subdirectory per phage
│   └── <accession>/             ← versioned accession (e.g., AB720063.2)
│       ├── genome.fna           ← whole-genome nucleotide sequence
│       ├── proteins.faa         ← NCBI-annotated protein sequences
│       └── cds.fna              ← CDS nucleotide sequences
│
└── processes/
    ├── fetch_bacteria.py        ← download script for bacteria
    ├── fetch_phages.py          ← download script for phages (777 total)
    └── accession_to_gcf.json    ← lookup table: seq accession → GCF assembly accession
```

---

## Organisms

### Bacteria — 38 resolved (+ 2 invalid flagged)

| Source | Count | Description |
|--------|-------|-------------|
| `confirmed_matrix` | 34 | Have accessions directly in `final_interaction_matrix.csv` |
| `newly_resolved` | 4 | Had `[No Complete Genome Found]` in matrix; accessions found via literature |
| `accession_invalid` | 2 | Matrix accessions are not genome assemblies (plasmid / patent sequence) — see below |
| `unresolved_TODO` | 11 groups | Cannot be resolved without team input — see `data_needs.md` |

**Newly resolved (add these to the interaction matrix):**

| Organism | Resolved Accession |
|----------|-------------------|
| *X. arboricola* pv. *juglandis* IVIA 1317-1a | CP076725 |
| *X. translucens* pv. *translucens* DSM 18974 | LT604072 |
| *X. arboricola* pv. *pruni* Xcp1 | CP090954 |
| *X. oryzae* pv. *oryzae* LN4 | CP045452 |

**Invalid accessions requiring correction in the interaction matrix:**

| Organism | Matrix Accession | Problem |
|----------|-----------------|---------|
| *Agrobacterium sp.* | KY000037 | Ti plasmid sequence, not a genome assembly |
| *Pseudomonas sp.* | PY746849 | Patent sequence (US 12571058 B2), not a genome assembly |

### Phages — 777 total

| Source | Count | Description |
|--------|-------|-------------|
| `canonical_xanthomonas` | 334 | Core Xanthomonas-infecting phages; appear in the interaction matrix |
| `negative_pool` | 443 | Non-Xanthomonas phages (*E. coli*, *Pseudomonas*, *Salmonella*, etc.); used as presumed-negative ML training samples |

> **Important:** The 443 negative-pool phages are NOT yet rows in `final_interaction_matrix.csv`. Before ML training, explicit `0` labels must be added for all (negative_pool phage × Xanthomonas host) pairs. Do not treat their absence as missing data.

---

## Setup — Install NCBI datasets CLI

The fetch scripts require the NCBI `datasets` CLI tool.

```bash
# Install into the project conda environment
conda install -n igem2026 -c conda-forge ncbi-datasets-cli

# Verify
datasets --version   # should show 18.x or later
```

If the `igem2026` environment does not yet exist:
```bash
conda env create -f shared/env/environment.yml
conda activate igem2026
conda install -c conda-forge ncbi-datasets-cli
```

---

## Running the Fetch Scripts

All scripts are run from the **repository root**:

```bash
# Download all bacteria (34 resolved genomes)
python 00_raw_data/processes/fetch_bacteria.py

# Download all 777 phages (can take 20–40 min depending on connection)
python 00_raw_data/processes/fetch_phages.py

# Download only canonical Xanthomonas phages (334)
python 00_raw_data/processes/fetch_phages.py --source canonical_xanthomonas

# Download only negative-pool phages (443)
python 00_raw_data/processes/fetch_phages.py --source negative_pool

# Dry run (shows what would be downloaded without downloading)
python 00_raw_data/processes/fetch_bacteria.py --dry-run
python 00_raw_data/processes/fetch_phages.py --dry-run

# Retry a single failed accession
python 00_raw_data/processes/fetch_bacteria.py --accession NZ_CP014028
python 00_raw_data/processes/fetch_phages.py --accession AB720063.2
```

Scripts are **resumable**: already-complete directories (all expected files present) are silently skipped.

---

## Known Issues and Required Follow-up

See `data_needs.md` for the full guide. Key action items:

1. **Fix duplicate accessions in the interaction matrix** — `NZ_JBWJFR000000000` is used for both *X. arboricola* and *X. arboricola* pv. *pruni*; `NZ_CP150073` is used for both *X. oryzae* and *X. oryzae* pv. *oryzae*. These pairs share one genome, making their matrix columns indistinguishable.

2. **Replace invalid accessions** — `KY000037` (Ti plasmid) and `PY746849` (patent sequence) are not genome assemblies. The matrix must be updated with proper genome accessions for *Agrobacterium sp.* and *Pseudomonas sp.* before these bacteria can be downloaded.

3. **Resolve 11 unresolved bacteria** — mostly *Xanthomonas* lab strains with no deposited genome. See `data_needs.md` Category A–C for resolution strategies.

4. **Add explicit 0-labels for negative-pool phages** — 443 × (number of Xanthomonas hosts) pairs need `0` in the interaction matrix before ML training.

---

## Data Flow to Downstream Steps

```
00_raw_data/bacteria/<acc>/genome.fna  ──→  02_annotation/  (Prodigal input)
00_raw_data/bacteria/<acc>/genome.fna  ──→  03_feature_weighting/  (GC content, CAI)
00_raw_data/bacteria/<acc>/proteins.faa ─→  03_feature_weighting/  (pI, GRAVY, length, similarity)
00_raw_data/bacteria/<acc>/cds.fna     ──→  03_feature_weighting/  (Factor 6: CAI codon table)

00_raw_data/phage/<acc>/genome.fna     ──→  02_annotation/  (PHANOTATE input)
00_raw_data/phage/<acc>/genome.fna     ──→  03_feature_weighting/  (GC content)
00_raw_data/phage/<acc>/cds.fna        ──→  03_feature_weighting/  (Factor 6: CAI computation)
```

Step 02 annotation (`02_annotation/`) re-derives `proteins.faa` and `.gff` for phages using PHANOTATE/pharokka. The NCBI copies here are a reference baseline, not the primary annotation source for phages.
