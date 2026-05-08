# AGENT TODO — Module 02 / `02_annotation/`

## Read first (mandatory)

1. `/INTERFACE.md` §Module 02 — your output spec.
2. `/CLAUDE.md` — notebook-first, bilingual.
3. The legacy `processes/phage_phanotate/` and `processes/host_prodigal/` if they exist — read for context, then replace with notebooks.
4. Invoke `superpowers:test-driven-development` and `superpowers:verification-before-completion`.

## Your scope

Genome → ORF → protein FASTA. Two distinct tools for two distinct problems:
- **Phage genomes** → PHANOTATE (handles overlapping ORFs, ~10% of phage genes overlap)
- **Bacterial genomes** → Prodigal (assumes non-overlapping, optimized for prokaryotes)

DO NOT swap tools. PHANOTATE on bacteria gives garbage; Prodigal on phages misses ~15% of genes (McNair 2019).

You read from `00_raw_data/{phage,bacteria}/<acc>/genome.fna`. You write to `02_annotation/outputs/`.

## Goal (definition of done)

By morning, `02_annotation/` contains:

1. ✅ `processes/01_run_phanotate.ipynb` — annotates phage genomes, produces `outputs/phage_proteins/<acc>.faa` + `outputs/phage_orfs/<acc>.gff3`.
2. ✅ `processes/02_run_prodigal.ipynb` — annotates host genomes, produces `outputs/host_proteins/<acc>.faa` + `outputs/host_orfs/<acc>.gff3`.
3. ✅ At least 1 phage annotated (phiL7 / EU717894.1).
4. ✅ At least 1 bacterium annotated (Xcc / GCF_000007145.1 if Module 01 fetched it; otherwise pick the smallest bacterium in `00_raw_data/bacteria/`).
5. ✅ `outputs/annotation_summary.csv` with one row per annotated genome.
6. ✅ `outputs/MANIFEST.csv`.
7. ✅ ≥3 passing tests.
8. ✅ `AGENT_REPORT.md`.

## Setup

```bash
cd /path/to/agent-02-annotation
conda activate igem2026

# PHANOTATE — pip from author's repo
pip install phanotate
# OR (fallback if pip fails):
git clone https://github.com/deprekate/PHANOTATE.git /tmp/PHANOTATE
cd /tmp/PHANOTATE && pip install . && cd -

# Prodigal — already in conda env via pyrodigal (Python wrapper)
# Verify:
python -c "import pyrodigal; print(pyrodigal.__version__)"

# pharokka (optional, for richer functional annotation)
# Heavy install; skip if pip fails — Module 03 (PhageRBPdetect) is the main consumer
# and PhageRBPdetect uses Pfam HMMs not pharokka categories
```

If PHANOTATE install fails on macOS (some compilation issues), document in `SETUP.md` and proceed with `pyrodigal` + a manual ORF-finding fallback for the smoke test only.

## Step-by-step plan

### Step 1 — Build the PHANOTATE notebook (`01_run_phanotate.ipynb`) (90 min)

Cells:

1. **Markdown** — Title + bilingual purpose. Cite McNair et al. 2019.
2. **Code** — Imports, paths, version printouts.
3. **Markdown** — Method explanation (bilingual): why PHANOTATE for phages.
4. **Code** — Helper: `run_phanotate(genome_fna_path, output_dir) → dict`. Wraps the CLI call (`subprocess`) or Python API. Returns metadata: n_orfs, runtime_s, version.
5. **Code** — Helper: `parse_phanotate_output_to_faa(phanotate_output) → faa_path`. Translates DNA ORFs → protein FASTA per INTERFACE format.
6. **Code** — Helper: `parse_phanotate_output_to_gff3(phanotate_output) → gff3_path` per Sequence Ontology GFF3 spec.
7. **Markdown** — Sample-then-batch pattern explanation.
8. **Code** — Sample run: 1 phage (phiL7 = EU717894.1).
9. **Code** — Validation: assert FASTA header format matches INTERFACE; assert n_orfs ∈ [50, 65] (phiL7 expected ~59 per Lee 2009).
10. **Code** — Batch run: optionally annotate all phages in `00_raw_data/phage/` (parallelize with `concurrent.futures.ProcessPoolExecutor`, save partial progress every 50 genomes).
11. **Code** — Append to MANIFEST.csv + annotation_summary.csv.

FASTA header MUST follow INTERFACE format:
```
>EU717894.1_orf_00031 | source=EU717894.1 | length=412 | start=12345 | end=13580 | strand=+ | tool=PHANOTATE_<version>
```

### Step 2 — Build the Prodigal notebook (`02_run_prodigal.ipynb`) (60 min)

Same structure as PHANOTATE notebook but for bacteria:

- Use `pyrodigal` (Python bindings, faster than CLI).
- Run in `single` mode (not `meta` — meta is for metagenomes).
- Use `pyrodigal.OrfFinder()` with default training on each genome (each bacterium gets its own model).
- Same FASTA header format.

Sample run: 1 bacterium (Xcc if available; else the smallest in `bacteria/`).
Sanity assertion: Xcc should have 4000-4500 ORFs (da Silva 2002 reports 4181).

### Step 3 — `annotation_summary.csv` (15 min)

Schema (per INTERFACE §Module 02):
`acc, type, n_orfs, mean_orf_len, tool, tool_version, runtime_s`

One row per annotated genome.

### Step 4 — Tests (45 min)

`processes/tests/`:
- `test_schema.py` — Open the smallest annotated `<acc>.faa`, assert every record's header matches the regex `^>(\S+) \| source=(\S+) \| length=(\d+) \|`. Open `annotation_summary.csv`, assert columns.
- `test_smoke.py` — Run `01_run_phanotate.ipynb`'s helper on a 5-kb dummy phage genome, assert returns ≥1 ORF.
- `test_sanity.py` — `outputs/phage_proteins/EU717894.1.faa` exists, has ≥50 records and ≤65 records, mean protein length ∈ [50, 800] aa. `outputs/host_proteins/GCF_000007145.1.faa` (or whichever bacterium ran) has ≥3000 records.

Run: `pytest 02_annotation/processes/tests/ -v`

### Step 5 — Commit + report

Commits:
- `02_annotation: PHANOTATE notebook + phiL7 annotation`
- `02_annotation: Prodigal notebook + Xcc annotation`
- `02_annotation: annotation_summary + manifest`
- `02_annotation: tests`
- `02_annotation: AGENT_REPORT`

## References (cite in notebook markdown cells)

- **PHANOTATE**: McNair, K. et al. (2019) "PHANOTATE: a novel approach to gene identification in phage genomes." *Bioinformatics* 35(22):4537-4542. DOI: 10.1093/bioinformatics/btz265
- **Prodigal**: Hyatt, D. et al. (2010) "Prodigal: prokaryotic gene recognition and translation initiation site identification." *BMC Bioinformatics* 11:119. DOI: 10.1186/1471-2105-11-119
- **pyrodigal** (Python binding used for performance): Larralde, M. (2022) "Pyrodigal: Python bindings and interface to Prodigal..." *Journal of Open Source Software* 7(72):4296. DOI: 10.21105/joss.04296
- **GFF3 spec**: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
- **phiL7 expected ORF count**: Lee, C.N. et al. (2009) *Appl Environ Microbiol* 75:7828.
- **Xcc expected ORF count**: da Silva, A.C.R. et al. (2002) *Nature* 417:459.

## Anti-goals

- ❌ Don't run PHANOTATE on bacteria.
- ❌ Don't run Prodigal on phages.
- ❌ Don't modify `00_raw_data/`.
- ❌ Don't try to annotate all 777 phages tonight if PHANOTATE per-genome is > 30s on this hardware — save the batch for Laguna. Just do phiL7 + maybe 5 more for the smoke test.
- ❌ Don't push to remote.

## Time budget

~4 hours. PHANOTATE notebook is the long pole. Skip optional pharokka enrichment if behind schedule.

## If stuck

- PHANOTATE install failure → document in SETUP.md, fall back to `pyrodigal` in "single" mode on phages with the explicit caveat that ~15% of overlapping ORFs will be missed (logged in AGENT_REPORT under "DOWNSTREAM IMPACT").
- Out of memory on Xcc Prodigal → use `meta` mode instead (faster, slightly less accurate; document choice).
