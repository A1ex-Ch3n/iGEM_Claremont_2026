# INTERFACE.md — Inter-Module Data Contract

**Status:** Locked for the 2026-05-08 overnight parallel build. Any deviation requires a `BREAKING_INTERFACE.md` note in the offending agent's commit.

**Audience:** All 7 parallel agents (00-06). Each agent must read this file in full before writing any code. Downstream agents use this as the spec for what their inputs look like; upstream agents use it as the spec for what their outputs must conform to.

---

## Universal conventions

### Path anchoring

Inside notebooks: `REPO_ROOT = Path.cwd().resolve().parents[1]` (notebooks live in `<module>/processes/`).
Inside `.py`: `REPO_ROOT = Path(__file__).resolve().parents[2]`.

**Never hard-code absolute paths.** Never use `~`. Never assume CWD.

### Identifier conventions

| Concept | Format | Example |
|---|---|---|
| Phage NCBI accession | `<base>.<version>` | `EU717894.1` |
| Bacterial assembly | `GCF_*` or `GCA_*` | `GCF_000007145.1` |
| ORF id (per genome) | `<acc>_orf_<5-digit>` | `EU717894.1_orf_00031` |
| RBP candidate id | `<acc>_rbp_<2-digit>` | `EU717894.1_rbp_01` |
| Receptor id (host protein) | `<host_acc>_<gene_name>` | `GCF_000007145.1_tonB` |
| Variant id (designed RBP) | `<parent_rbp>_<change>_<idx>` | `EU717894.1_rbp_01_trunc_03` |

### Filename conventions

- Snake_case, no spaces, no `·`, no Chinese characters in filenames.
- ISO dates in filenames: `YYYY-MM-DD`.
- Cycle tags: `cycle_<N>` where N is integer ≥ 0.
- Model versions in filenames: `<tool>_<version>` (e.g., `esm2_t6_8M_UR50D`).

### MANIFEST.csv (every `outputs/` directory must have one)

Required columns: `filename, sha256, bytes, n_records, created_utc, source_acc, source_module, notes`

`n_records` = number of FASTA records / CSV rows / structures in that file. NaN if N/A.
`created_utc` = ISO 8601 with `Z` suffix.
`source_module` = the module that produced it (e.g., `02_annotation`).

### CSV conventions

- UTF-8, Unix line endings, RFC 4180.
- First row is header.
- No trailing whitespace in cells.
- Empty cells = empty string `""` (not `NaN` or `null`).
- Floats: 6 decimal places max.

### FASTA conventions

- One record per protein/gene.
- Header format: `><id> | source=<acc> | length=<aa> | <free_text>`
- Example: `>EU717894.1_orf_00031 | source=EU717894.1 | length=412 | tail spike candidate (PhageRBPdetect rank 1)`
- Wrap at 60 chars/line.
- Use `*` for stop codon if present in source; strip leading `M` only if explicitly noted.

### NumPy / npz conventions

For embeddings and arrays:
- Save as `.npz` (compressed if > 50 MB).
- Required keys: `seq_ids` (array of strings), `array` (the data), `meta` (JSON-string with model name, version, date, pooling strategy).
- Embedding arrays: `dtype=float32`, shape `(N_sequences, embed_dim)` for pooled, `(N_sequences, max_len, embed_dim)` for per-residue (with separate `lengths` key).

---

## Per-module contracts

Below: each module's **inputs (what it reads)** and **outputs (what it writes)**. If your module's actual outputs deviate, fix your module — don't fix this file.

---

### Module 00 — `00_raw_data/`

**Role:** Raw genome dump. Source of truth for nucleotide sequences. Read-only after this module finalizes.

**Inputs:** None (downloads from NCBI). Existing files preserved per user decision.

**Outputs:**

```
00_raw_data/
├── README.md
├── MANIFEST.csv                                    ← global checksum index
├── data_needs.md                                   ← what we plan to fetch (intent doc)
├── phage_list.csv                                  ← columns: accession, description, length_bp, source, host_genus, in_interaction_matrix
├── bacteria_list.csv                               ← columns: accession, gcf_accession, organism_name, source, matrix_interactions, notes
├── phage/<accession>/                              ← one subdirectory per phage
│   ├── genome.fna                                  ← single FASTA, nucleotide
│   ├── cds.fna                                     ← coding sequences (nucleotide)
│   └── proteins.faa                                ← translated proteins
├── bacteria/<accession>/                           ← same structure as phage/
│   ├── genome.fna
│   ├── cds.fna
│   └── proteins.faa
└── processes/
    ├── 01_verify_dataset.ipynb                     ← notebook checking integrity
    └── tests/
```

**Schema validators (all must pass):**

- Each `phage/*/genome.fna` is parseable with `Bio.SeqIO`, contains exactly 1 sequence record (some phages have 2 if circular).
- Each `bacteria/*/genome.fna` may contain >1 record (chromosomes + plasmids).
- `proteins.faa` records have headers matching `>(\S+) ` regex.
- `MANIFEST.csv` row count == number of files in `phage/` + `bacteria/` subtrees.

---

### Module 01 — `01_data_ground_truth/`

**Role:** NCBI download orchestration + curated reference interaction matrix.

**Inputs:** External NCBI Entrez/Datasets API; existing `inputs/` seed lists.

**Outputs:**

```
01_data_ground_truth/
├── inputs/
│   └── reference_targets.csv                       ← which references we want fetched
│                                                     columns: category, assembly_acc, nucleotide_acc, label
├── outputs/
│   ├── MANIFEST.csv
│   ├── interaction_matrix.csv                      ← see schema below
│   ├── reference_genomes_index.csv                 ← columns: label, accession, target_dir, status, sha256_genome
│   └── download_log_<YYYY-MM-DD>.csv               ← columns: accession, status, attempts, error_msg, fetched_utc
└── processes/
    ├── 01_fetch_reference_genomes.ipynb            ← (already exists, refactor to match this contract)
    ├── 02_build_interaction_matrix.ipynb           ← canonicalize the historical matrix
    └── tests/
```

**`interaction_matrix.csv` schema:**

| Column | Type | Description |
|---|---|---|
| `phage_acc` | string | Phage NCBI accession with version (e.g., `EU717894.1`) |
| `host_acc` | string | Host NCBI accession or assembly id |
| `host_organism` | string | Genus + species + pathovar |
| `label` | int | 1 = lyses, 0 = does not lyse, -1 = unknown/untested |
| `source` | string | `literature_curated`, `experimental_in_house`, `inferred_taxonomy` |
| `confidence` | float | 0.0-1.0 |
| `notes` | string | Free text |

**Note:** Reference genomes themselves are written into `00_raw_data/{phage,bacteria}/<acc>/`, NOT into `01_data_ground_truth/outputs/`. Module 01 only writes the *index* of what was fetched.

---

### Module 02 — `02_annotation/`

**Role:** Predict ORFs and translate to proteins. Phage uses PHANOTATE; bacteria uses Prodigal.

**Inputs:** `00_raw_data/{phage,bacteria}/<acc>/genome.fna`

**Outputs:**

```
02_annotation/
├── outputs/
│   ├── MANIFEST.csv
│   ├── phage_proteins/<acc>.faa                    ← protein FASTA from PHANOTATE
│   ├── phage_orfs/<acc>.gff3                       ← ORF coordinates (GFF3)
│   ├── host_proteins/<acc>.faa                     ← protein FASTA from Prodigal
│   ├── host_orfs/<acc>.gff3                        ← ORF coordinates (GFF3)
│   └── annotation_summary.csv                      ← columns: acc, type, n_orfs, mean_orf_len, tool, tool_version, runtime_s
└── processes/
    ├── 01_run_phanotate.ipynb                      ← phage annotation
    ├── 02_run_prodigal.ipynb                       ← host annotation
    └── tests/
```

**FASTA header format (REQUIRED):**

```
>EU717894.1_orf_00031 | source=EU717894.1 | length=412 | start=12345 | end=13580 | strand=+ | tool=PHANOTATE_1.5.0
```

ORF numbering: 5-digit zero-padded, in genome order, starting at 00001.

**Sanity-check expectations** (use as test assertions):
- phiL7 (EU717894.1): 50-65 phage ORFs (Lee et al. 2009 reports 59).
- Xcc ATCC 33913 (GCF_000007145.1): 4000-4500 host ORFs (da Silva et al. 2002 reports 4181).

---

### Module 03 — `03_rbp_identification/`

**Role:** Identify receptor-binding protein candidates from phage proteomes.

**Inputs:** `02_annotation/outputs/phage_proteins/<acc>.faa`

**Outputs:**

```
03_rbp_identification/
├── outputs/
│   ├── MANIFEST.csv
│   ├── <phage_acc>_rbp_candidates.csv              ← all candidates with scores
│   ├── <phage_acc>_rbp_sequences.faa               ← top-K candidates as FASTA (default K=5)
│   └── all_rbp_candidates.csv                      ← unified across all phages
└── processes/
    ├── 01_run_phagerbpdetect.ipynb
    └── tests/
```

**`<phage_acc>_rbp_candidates.csv` schema:**

| Column | Type | Description |
|---|---|---|
| `orf_id` | string | matches `02_annotation` ORF id |
| `phage_acc` | string | source phage accession |
| `length_aa` | int | protein length in amino acids |
| `hmm_score` | float | Pfam HMM track score (NaN if not detected) |
| `hmm_match_pfam` | string | Matched Pfam domain id (NaN if none) |
| `ml_score` | float | XGBoost-on-ESM track probability [0,1] |
| `combined_score` | float | aggregated score (HMM rule: if HMM hit, use 1.0; else use ml_score) |
| `evidence_track` | string | `hmm`, `ml`, `both` |
| `rank` | int | 1 = top candidate |
| `passes_threshold` | bool | combined_score ≥ 0.5 |

**`<phage_acc>_rbp_sequences.faa` header format:**

```
>EU717894.1_rbp_01 | source_orf=EU717894.1_orf_00031 | length=412 | combined_score=0.94 | evidence=both
```

**Sanity check:** phiL7 should yield 1-3 RBP candidates with combined_score > 0.7 (the tail spike gp25 is well-characterized).

---

### Module 04 — `04_protein_embedding/`

**Role:** Generate ESM-2 embeddings for RBP variants and host receptors.

**Inputs:**
- Primary: `03_rbp_identification/outputs/<phage_acc>_rbp_sequences.faa`
- Secondary: receptor sequences extracted from `02_annotation/outputs/host_proteins/`

**Outputs:**

```
04_protein_embedding/
├── outputs/
│   ├── MANIFEST.csv
│   ├── embeddings_esm2_t6_8M_<set>.npz             ← CPU-friendly default
│   ├── embeddings_esm2_t33_650M_<set>.npz          ← Laguna run
│   └── embedding_index.csv                         ← columns: seq_id, source_file, model, n_residues, embed_dim, npz_file
└── processes/
    ├── 01_embed_esm2.ipynb
    └── tests/
```

**`.npz` keys (REQUIRED):**

- `seq_ids` — array of bytes/strings, shape `(N,)`
- `array` — float32, shape `(N, embed_dim)` for mean-pooled or `(N, max_len, embed_dim)` for per-residue
- `lengths` — int32, shape `(N,)` — original sequence lengths (used to mask padding in per-residue)
- `meta` — JSON string with: `{model, model_revision, pooling, date_utc, repo_commit_sha}`

**Default model selection:** `esm2_t6_8M_UR50D` for CPU testing (fast, low memory), `esm2_t33_650M_UR50D` for production.

---

### Module 05 — `05_structure_prediction/`

**Role:** 3D structure prediction (AlphaFold 3) and zero-shot binding affinity (Boltz-2).

**Inputs:**
- `03_rbp_identification/outputs/<phage_acc>_rbp_sequences.faa`
- Receptor sequences (TonB / ExbB / ExbD1 / ExbD2 from Xcc)

**Outputs:**

```
05_structure_prediction/
├── outputs/
│   ├── MANIFEST.csv
│   ├── af3/<seq_id>/
│   │   ├── model.cif                               ← top-ranked model
│   │   └── metrics.json                            ← {plddt, ipTM, ptm, ranking_score}
│   ├── boltz2/<rbp_id>__<receptor_id>/
│   │   ├── model.pdb
│   │   └── affinity.json                           ← {predicted_dG_kcal_mol, predicted_pKd, confidence, model_version}
│   └── affinity_priors.csv                         ← aggregated table fed into Module 06
└── processes/
    ├── 01_run_boltz2.ipynb                         ← primary (open source, runs locally)
    ├── 02_run_af3.ipynb                            ← optional (AF3 install hard; ColabFold fallback documented)
    └── tests/
```

**`affinity_priors.csv` schema:**

| Column | Type | Description |
|---|---|---|
| `rbp_id` | string | matches Module 03 RBP id |
| `receptor_id` | string | matches Module 04 receptor id format |
| `model` | string | `boltz2_1.0`, `af3_multimer_v3`, etc. |
| `predicted_dG_kcal_mol` | float | binding free energy estimate |
| `predicted_pKd` | float | -log10 Kd estimate |
| `confidence` | float | model's self-reported confidence [0,1] |
| `interface_ipTM` | float | interface predicted TM score (AF3 only, NaN otherwise) |

---

### Module 06 — `06_uncertainty_model/`

**Role:** Train a deep ensemble regressor that predicts RBP-receptor binding score with calibrated uncertainty.

**Inputs:**
- `04_protein_embedding/outputs/embeddings_*.npz`
- `05_structure_prediction/outputs/affinity_priors.csv` (used as cold-start synthetic labels for Cycle 0 demo)
- `08_cycle_data/outputs/cycle_<N>/elisa_processed.csv` (real ELISA labels — not available tonight)

**Outputs:**

```
06_uncertainty_model/
├── outputs/
│   ├── MANIFEST.csv
│   ├── cycle_<N>/
│   │   ├── ensemble_member_<i>.pt                  ← i ∈ {0, 1, 2, 3, 4}
│   │   ├── predictions.csv                         ← schema below
│   │   ├── calibration.png                         ← reliability diagram
│   │   ├── training_log.json
│   │   └── model_meta.json                         ← {arch, hidden_dim, train_size, val_size, ...}
└── processes/
    ├── 01_train_deep_ensemble_synthetic.ipynb     ← demo on Boltz-2 priors as fake labels
    ├── 02_calibration_check.ipynb
    └── tests/
```

**`predictions.csv` schema:**

| Column | Type | Description |
|---|---|---|
| `rbp_id` | string |  |
| `receptor_id` | string |  |
| `predicted_score` | float | ensemble mean |
| `std` | float | ensemble std (epistemic uncertainty) |
| `lower_95` | float | mean - 1.96*std |
| `upper_95` | float | mean + 1.96*std |
| `model_version` | string | git commit sha + cycle tag |

---

## Cross-cutting requirements

### Tests (every module)

Each module's `processes/tests/` must contain at minimum:

1. **`test_schema.py`** — load each output file, assert columns/keys match this INTERFACE.md.
2. **`test_smoke.py`** — run the notebook end-to-end on the smallest possible input, assert it completes.
3. **`test_sanity.py`** — module-specific sanity check (e.g., phiL7 → ≈59 ORFs).

Tests run via `pytest <module>/processes/tests/`. Each test must complete in < 60 seconds on CPU.

### Reproducibility

Every notebook must:
- Print versions of all key libraries in cell 2.
- Set random seeds (numpy, torch, python random) explicitly.
- Save `repo_commit_sha` into output metadata.

### Logging

Use Python `logging` module (not `print`) for any code in `.py` files. Notebooks may use `print` for human-readable output.

### Dependencies

Each module owns a `<module>/requirements.txt` if it needs Python packages beyond `shared/env/environment.yml`. Don't pollute the shared env with module-specific tools.
