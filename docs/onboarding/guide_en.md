# Onboarding Guide — iGEM Claremont 2026

**Active-Learning Phage Engineering for *Xanthomonas* Biocontrol**

> Companion to `slides_en.pdf`. Read this if you missed the talk, or to find the actual file paths, commands, and references behind a slide.

| | |
|---|---|
| **Core engineer** | Alex Chen |
| **Wet lab leads** | Sarah, Olivia, Weitao, Carol |
| **Contributors** | Ryan, Leah |
| **PI** | Prof. J. Cesar Ignacio-Espinoza |
| **Faculty advisor** | Prof. Ran Libeskind-Hadas |
| **Branch** | `active-learning-pipeline` (main is empty — do not work from main) |
| **Today's status** | Modules 00–07 complete and tested. Module 08 opens with the first ELISA delivery (~June 1). |

---

## 1. Project at a glance

We are building a **closed-loop active-learning pipeline** that pairs a machine-learning model of phage RBP (receptor-binding protein) × bacterial-receptor interaction with iterative wet-lab validation. After each ELISA round, the model retrains, the BALD acquisition function ranks every untested variant by epistemic uncertainty, and the wet lab measures the next 4–5. The whole loop is designed to extract maximum information from each expensive ELISA measurement — addressing the central pain point of phage-host prediction: **data scarcity**.

Reference dry-lab scaffold:
- Host: *Xanthomonas campestris* pv. *campestris* (Xcc) ATCC 33913 — NCBI `GCF_000007145.1`.
- Phage: phiL7 — NCBI `EU717894`.
- Receptor system: TonB-ExbB-ExbD1 essential, ExbD2 **not** required (Hung 2003, *BBRC* 302:878–884, PMID 12646254).

Wet lab self-isolates *Xanthomonas* + phage from California brassicas (per PI consultation 2026-05-07), bypassing the multi-month USDA APHIS PPQ-526 permit.

iGEM tracks targeted: **Best Agriculture Project · Best Model · Best Composite Part**.

---

## 2. The science — why we're doing this

### *Xanthomonas* and the host-range problem

The genus *Xanthomonas* contains plant pathogens of >400 species (Ryan 2011, *Nat Rev Microbiol*). California-relevant pathovars include Xcc (black rot of brassicas, present in commercial vegetable production), *X. citri* subsp. *citri* (citrus canker), and *X. perforans* / *X. euvesicatoria* (bacterial spot of tomato). Current control relies on copper bactericides, against which resistance is now widespread (Aiello 2019, *Plant Dis*).

Bacteriophage biocontrol offers an attractive alternative — phages are host-specific, self-amplifying, environmentally degradable. Field trials match copper efficacy for bacterial spot of tomato (Iriarte 2018) and suppress Xcc in cabbage (Holtappels 2022). But the same specificity that makes phages safe also makes them hard to deploy: each strain requires the right phage.

Predicting which phage infects which bacterium is unreliable. The state-of-the-art model, PhageHostLearn (Boeckaerts 2024, *Nat Commun*), reaches **AUC ≈ 0.82 within one genus** and degrades to **0.67–0.70 across genera** (PAML benchmark, Mutalik 2025). The fundamental bottleneck is data scarcity — quantitative (phage, host) binding labels are slow and expensive to generate.

### phiL7 + Xcc — what we already know

phiL7 is a 44,080 bp lytic Siphoviridae characterized by Lee et al. 2009 (*AEM* 75:7828). It infects *X. campestris* through the **TonB-ExbB-ExbD1** outer-membrane import complex. Hung et al. 2003 demonstrated this experimentally via Tn5 mutagenesis: TonB, ExbB, and ExbD1 are essential for phage penetration; ExbD2 — though co-transcribed in the same operon — is **not required** (strain CH620, ΔexbD2, retains full sensitivity).

This last point is more useful than it sounds. Generating ΔexbD2 alongside ΔtonB / ΔexbB / ΔexbD1 gives us a **free negative control**: if our knockout system is working, ΔexbD2 should still allow infection. That's a built-in validation of the entire knockout pipeline.

### rbp_01 — Lee 2009 and our HMM rediscovery

Lee 2009 explicitly searched for a tail-fiber homolog of OP1 ORF25 in phiL7 and **was unable to find one** by sequence similarity:

> "We were unable to identify a homolog of the OP1 tail fiber protein (ORF25) thought to be involved in host range determination…"

p20 (1105 aa, tail protein III) was suggested as host-range related, but no protein was named a "tail spike". Our **rbp_01 (712 aa)** comes from PhageRBPdetect's Tail_spike_N HMM (Boeckaerts 2022, *Viruses*) — a structural-profile method that finds proteins too diverged for BLAST to detect. This is not a contradiction with Lee 2009; HMMs and BLAST have different sensitivity profiles, and rbp_01 is a complementary discovery using a more sensitive tool.

### Why active learning

Active learning (AL) is a mathematical response to data scarcity (Lindley 1956; Settles 2009). Rather than passively training on whatever data exists, an AL system uses its current uncertainty to pick the experiment that, if performed, would most reduce that uncertainty. Acquisition functions like BALD (Houlsby 2011) formalize this as information-gain maximization.

Recent demonstrations:
- **Hie 2024** (*Nat Biotechnol* 42:275) — human antibody affinity maturation guided by ESM-1b/1v language-model likelihood + AL, ~20 variants per antibody, up to 160-fold affinity improvement.
- **Yang 2025 ALDE** (*Nat Commun*) — DNN ensemble + Thompson sampling (note: *not* BALD) on enzyme directed evolution; yield 12% → 93% in 2 wet-lab rounds.

**Neither has been applied to phage RBP × bacterial receptor.** That's our methodological contribution.

---

## 3. Pipeline architecture

### The three layers

```
LAYER 0 — Priors          Boltz-2 ipTM, ESM-2 embeddings, PLM-interact (optional)
                          (informative prior before any wet-lab data)
                                       │
                                       ▼
LAYER 1 — AL loop         Ensemble → BALD → recommendations.csv
                          → wet-lab ELISA → retrain → next recommendations
                                       │
                                       ▼
LAYER 2 — Causality       ΔtonB / ΔexbB / ΔexbD1 knockouts + ΔexbD2 negative control
                          (decouples receptor-specific binding from defense)
```

Binding (Layer 1) is necessary but **not sufficient** for productive infection (Farquharson 2021, T4 × *E. coli*: RBP binds 85 % of strains but plaques only on 11 %). Layer 2 quantifies how much of the model's binding signal translates to infection.

### Module map (00 → 08)

![Pipeline flow](figures/pipeline_flow.png)

| Module | Tool | Status | Key output |
|--------|------|--------|------------|
| 00 Raw Data | NCBI fetch | ✅ | 777 phage + 34 bacteria; `MANIFEST.csv` (binaries gitignored) |
| 01 Ground Truth | curated | ✅ | `interaction_matrix.csv`: 2,236 phage–host pairs |
| 02 Annotation | PHANOTATE / pyrodigal | ✅ | phiL7: 80 ORFs; Xcc: 4,344 ORFs |
| 03 RBP ID | PhageRBPdetect HMM + XGBoost | ✅ | `EU717894.1_rbp_candidates.csv`; rbp_01 = 712 aa |
| 04 Embedding | ESM-2 (8M local; 650M / 3B on Laguna) | ✅ | `embeddings_esm2_*.npz` |
| 05 Structure | Boltz-2 (Laguna); AF3 pending | ✅ | `affinity.json` (`ipTM = 0.365`), PDB, PAE matrix |
| 06 Uncertainty | 5-MLP deep ensemble | ✅ | `predictions.csv` with `std` and **`epistemic_std`** |
| 07 BALD | Var_k[μ_k] acquisition | ✅ | `recommendations.csv` + `safe_pick_backup.csv` |
| 08 Cycle Data | wet-lab ELISA | ⏳ ~June 1 | `elisa_processed.csv` |

### The data-contract convention

![Data contract](figures/data_contract.png)

Every module has the same three subdirectories:

- **`inputs/`** — read-only pointers to upstream `outputs/` (or external seeds like NCBI accession lists). **Never write generated data here.**
- **`processes/`** — the only place code lives. Scripts and notebooks read from `inputs/` and write to their own `outputs/`.
- **`outputs/`** — canonical products consumed by the next module. Large artifact trees are gitignored — a `MANIFEST.csv` with SHA-256 + size per file makes the set reproducible.

Full spec — including identifier formats, FASTA header conventions, MANIFEST schema — is in `INTERFACE.md` at the repo root.

### Notebook-first workflow + bilingual comments

From `CLAUDE.md`:

- New code is authored as **Jupyter notebooks first** (`<NN>_<short_name>.ipynb`) for exploration.
- Once a notebook (a) runs end-to-end and (b) matches a verification cell, **freeze it as `.py`** in the same folder; rename the notebook to `<NN>_<short_name>__frozen.ipynb` with a pointer to the `.py`.
- **Bilingual comments** mandatory — short inline: `# English / 中文` on the same line; longer: separate paragraphs in markdown cells.
- Module 07 is the exception — production orchestration code, written as `.py` from day 1.
- `nbstripout` is enforced repo-wide so notebook outputs never pollute git diffs. Run `pip install nbstripout && nbstripout --install` once per clone.

---

## 4. The ML core in depth

### Module 06 — Deep ensemble (Lakshminarayanan 2017, *NeurIPS*)

Five independently-initialized MLPs (3 hidden layers, ReLU + dropout, two output heads for `mean` and `log_sigma`) are trained on ESM-2 embeddings to predict the ELISA binding score. Each member produces a Gaussian $(\mu_k, \sigma_k^2)$; training loss is the Gaussian negative log-likelihood.

```python
# From 06_uncertainty_model/processes/ensemble.py:57–68
def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Returns (mean, sigma) — both shape (batch,).
    返回 (均值, 标准差)，各形状为 (batch,)。
    """
    h = self.net(x)
    mean = self.head_mean(h).squeeze(-1)
    # Clamp log_sigma to [-7, 7] for numerical stability
    # 将 log_sigma 截断到 [-7, 7] 防止数值爆炸
    log_sigma = self.head_log_sigma(h).squeeze(-1).clamp(-7.0, 7.0)
    sigma = torch.exp(log_sigma)
    return mean, sigma
```

The total predictive variance decomposes:

$$\sigma^2_{\text{total}}(x) \;=\; \underbrace{\mathrm{Var}_k[\mu_k(x)]}_{\text{epistemic — reducible}} + \underbrace{\mathbb{E}_k[\sigma_k^2(x)]}_{\text{aleatoric — noise}}$$

`predictions.csv` exports both: `std` (total) and `epistemic_std` (BALD input).

**Why deep ensembles?** Better calibration than MC Dropout (Ovadia 2019); scales to ESM-2's 1280-D inputs (Gaussian Processes don't). The Greenman 2025 audit (*PLoS Comput Biol* 21:e1012639) is important here:

> *"There is no single best UQ method"* across protein-engineering benchmarks. Ensembles tend to be most accurate but worst-calibrated.

We pick ensembles for scalability + ALDE precedent, and patch ECE with temperature scaling if it drifts past 0.1.

#### Cycle 0 output preview

`06_uncertainty_model/outputs/cycle_0/predictions.csv` (first rows):

| rbp_id | receptor_id | predicted_score | std | epistemic_std |
|--------|-------------|-----------------|-----|---------------|
| rbp_00 | rec_00      | 4.872           | 0.736 | 0.125 |
| rbp_01 | rec_02      | **5.128**       | 0.726 | **0.190** |
| rbp_07 | rec_02      | 5.177           | 0.721 | **0.218** |

80 pairs scored; `epistemic_std` ranges 0.04–0.22 across the pool. Synthetic data; will be replaced by real ELISA after Cycle 0.

The calibration plot (`outputs/cycle_0/calibration.png`) is near-diagonal on synthetic data — a meaningful version appears after Cycle 0 ELISA.

### Module 07 — BALD acquisition (Houlsby 2011 extension)

Original Houlsby 2011 formulation (Gaussian Process Classifier):

$$\mathrm{BALD}(x) \;=\; I(y;\theta \mid x, D) \;=\; H[y \mid x, D] - \mathbb{E}_{\theta \sim p(\theta \mid D)}\![H[y \mid x, \theta]]$$

For a Gaussian deep ensemble regression target:

$$\mathrm{BALD}(x) \;\approx\; \mathrm{Var}_k[\mu_k(x)] \;\;\Longrightarrow\;\; \text{score} = \text{epistemic\_std}$$

**Extension note.** Houlsby 2011 was originally classification with a GP. Extending the information-theoretic objective to deep ensemble regression is straightforward but is our extension — *not* a direct citation. When presenting to academic reviewers, mention this explicitly.

```python
# From 07_acquisition_function/processes/bald.py:38–135 (abridged)
def bald_score(epistemic_std: np.ndarray) -> np.ndarray:
    """BALD score = epistemic_std (rank-equivalent to epistemic variance)."""
    if np.any(epistemic_std < 0):
        raise ValueError("epistemic_std must be non-negative; got negatives.")
    return epistemic_std.copy()

# Inside select_batch():
ranked = pool.sort_values(bald_col, ascending=False).reset_index(drop=True)
ranked["bald_rank"] = ranked.index + 1
bald_picks = ranked.iloc[:n_bald].copy()
bald_picks["selection_reason"] = [f"bald_top_{i+1}" for i in range(len(bald_picks))]

# Random control: sample from remaining (NOT in BALD top-K)
remaining = ranked[~ranked.index.isin(set(bald_picks.index))]
random_picks = ranked.loc[_random.Random(seed).sample(list(remaining.index), k=n_random)]
random_picks["selection_reason"] = "random_control"
```

#### Cycle 1 output preview

`07_acquisition_function/outputs/cycle_1/recommendations.csv`:

| rbp_id | receptor_id | predicted_score | epistemic_std | bald_rank | selection_reason |
|--------|-------------|-----------------|---------------|-----------|------------------|
| rbp_07 | rec_02      | 5.177           | **0.218**     | 1         | bald_top_1       |
| rbp_03 | rec_01      | 4.810           | 0.197         | 2         | bald_top_2       |
| rbp_05 | rec_02      | 4.906           | 0.196         | 3         | bald_top_3       |
| rbp_01 | rec_02      | 5.128           | 0.190         | 4         | bald_top_4       |
| rbp_03 | rec_03      | 5.128           | 0.127         | 19        | random_control   |

4 BALD picks + 1 random control (Hie 2022 pattern). These are synthetic placeholders — real Cycle 1 = rbp_01 variants × TonB after Cycle 0 ELISA.

### ALDE caveat — Yang 2025 ≠ BALD validation

Yang 2025 ALDE is the closest published analog to our project but **does not validate BALD**. It uses Thompson sampling and one-hot encoding (not ESM-2). So ALDE supports the broader claim *"AL + UQ + DNN ensemble works in protein engineering"*; for BALD specifically, the citation chain is Houlsby 2011 → our extension to ensembles.

### A worked example of the data-contract pattern

`07_acquisition_function/processes/run_bald.py` reads from `inputs/` pointers and writes to `outputs/`:

```python
# Inputs (read-only pointers)
predictions_path = (
    REPO_ROOT
    / "06_uncertainty_model" / "outputs" / f"cycle_{cycle-1}"
    / "predictions.csv"
)
measured_pairs = load_measured_pairs(args.measured_csv) if args.measured_csv else set()

df = pd.read_csv(predictions_path)
df["bald_score"] = bald_score(df["epistemic_std"].values)

recommendations, random_replay, safe_pick_backup = select_batch(
    df, n_bald=args.n_bald, n_random=args.n_random,
    measured_pairs=measured_pairs, seed=args.seed,
)

# Outputs (canonical artifacts + provenance)
out_dir = REPO_ROOT / "07_acquisition_function" / "outputs" / f"cycle_{cycle}"
out_dir.mkdir(parents=True, exist_ok=True)
recommendations.to_csv(out_dir / "recommendations.csv", index=False)
safe_pick_backup.to_csv(out_dir / "safe_pick_backup.csv", index=False)
random_replay.to_csv(out_dir / "random_replay.csv", index=False)
(out_dir / "run_meta.json").write_text(json.dumps(run_meta, indent=2))
```

No hard-coded absolute paths. `REPO_ROOT` is anchored via `Path(__file__).resolve().parents[2]`.

---

## 5. The current Boltz-2 result

We ran Boltz-2 (Passaro 2025) on Laguna (job 59986, NVIDIA L40S) for the pair:

- **Chain A:** rbp_01, 712 aa, from `03_rbp_identification/outputs/EU717894.1_rbp_sequences.faa`
- **Chain B:** TonB, 604 aa, from `04_protein_embedding/outputs/...`

The three numbers from `affinity.json`:

| Metric | Value | Interpretation |
|--------|-------|---------------|
| `interface_ipTM` | **0.365** | Low. Model is uncertain how rbp_01 and TonB dock. Expected for a novel system with no PDB template; ELISA resolves this. |
| `chain_A_ptm` | **0.808** | High. rbp_01 monomer is structurally well-constrained — reliable for variant design. |
| `confidence_score` | 0.683 | Overall complex quality, moderate. |
| `predicted_dG` | `null` | **Affinity head is small-molecule only.** |

The low ipTM is not a model failure — it defines the experiment. That structural uncertainty IS the question the ELISA + active learning loop is designed to answer.

**Boltz-2 affinity head caveat (Passaro 2025).** The affinity head was trained on small-molecule × protein binding data (PubChem, ChEMBL, BindingDB). For protein-protein pairs (RBP × TonB), it outputs `NaN`. Always use **ipTM** as the structural confidence proxy — never describe Boltz-2 as providing a "zero-shot affinity prior" for protein-protein.

### PAE heatmap

![rbp_01 × TonB PAE](../../pae_heatmap.png)

Generated from `pae_*.npz` (1316×1316 float32). Residues 0–711 = Chain A (rbp_01); 712–1315 = Chain B (TonB). The **off-diagonal block** (rows 712–1315 × cols 0–711) is the interface region — low values (dark blue) indicate residue pairs where the model is confident about their relative position. Light bands in this block correspond to interface uncertainty — the parts ELISA will most informatively constrain.

To regenerate the heatmap from raw `.npz` data, see the code snippet at the bottom of `docs/planning/PI_briefing_2026-05-11.md`.

### File locations

All Boltz-2 outputs live under:

```
05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/
├── affinity.json                                          ← summary scores
└── boltz_results_boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/
    └── predictions/boltz_input_EU717894.1_rbp_01__GCF_000007145.1_tonB/
        ├── boltz_input_..._model_0.pdb                    ← 3D atomic structure
        ├── pae_..._model_0.npz                            ← PAE matrix
        ├── plddt_..._model_0.npz                          ← per-residue pLDDT
        └── confidence_..._model_0.json                    ← per-chain confidence
```

Open the PDB in PyMOL (`pymol.org`) or ChimeraX (`rbvi.ucsf.edu/chimerax/`) to inspect the predicted complex. PI_briefing.md has full untruncated paths.

---

## 6. The 48-hour cycle — wet/dry interface

![48-hour cycle](figures/cycle_48h.png)

### Wet lab → dry lab (end of cycle)

Wet lab delivers to `08_cycle_data/outputs/cycle_<N>/`:

| File | Required columns |
|------|-----------------|
| `elisa_processed.csv` | `variant_id, receptor_id, ec50_nM, hill_slope, r2, plate_id, date` |
| `plaque_results.csv`  | `variant_id, strain_id, pfu_per_ml, plaque_morphology, date` (WT and ΔReceptor) |
| `qc_report.md`        | SDS-PAGE image path, Bradford concentration, expression notes |

**Minimum for retraining:** ≥3 valid EC50 per variant with R² > 0.9. Failed variants are marked `ec50_nM = NaN` with a `failed_reason` — the model handles missing data.

### Dry lab → wet lab (48-hour SLA)

Within 48 hours of receiving ELISA data, dry lab produces, in `07_acquisition_function/outputs/cycle_<N+1>/`:

| File | Purpose |
|------|---------|
| `recommendations.csv` | 4 BALD picks + 1 random control; primary task list |
| `primer_sequences.txt` | NEB Q5-compatible SDM primers (auto-generated) |
| `uncertainty_bands.png` | Calibration: predicted vs measured from previous cycle |
| `safe_pick_backup.csv` | Top-10 BALD; used only if the 48-h SLA misses |
| `run_meta.json` | Provenance: git SHA, timestamp, pool size, score stats |

**Blinded control:** wet lab does NOT know which row is the random control — the recommendations CSV does not label `random_control` separately when sent to the wet lab. This preserves the retrospective AL-vs-random comparison.

### Cloning execution

- **Cycle 0** — gene synthesis (IDT/Twist), 4–6 variants chosen by structure-based design + expert picks. ~2 weeks lead time, ~$150/variant.
- **Cycles 1+** — site-directed mutagenesis (NEB Q5 kit) on existing constructs. ~4 days, ~$50/variant — 3× cheaper and 3.5× faster than gene synthesis.

Vector: pET-28a (Addgene 69864-3). Host: BL21(DE3). Induction: 0.5 mM IPTG, 18 °C overnight (favors soluble trimer assembly per Studier & Moffatt 1986).

### ELISA protocol (cell-based, Boeckaerts 2024 + Latka 2021)

1. Coat 96-well plate with heat-inactivated *Xanthomonas* (10⁸ CFU/well).
2. Block 3% BSA, 1 h.
3. Apply serial dilutions of purified His6-RBP (1 nM – 1 µM).
4. Detect with HRP-anti-His6, TMB substrate, OD₄₅₀.
5. Fit 4PL → **EC50** is the active-learning target variable.

Controls per plate: BSA-only (non-specific binding), WT-RBP at fixed concentration (inter-plate normalizer), heat-denatured RBP (folding-specific binding), 3 technical replicates per concentration.

### Receptor knockouts — pK18mobsacB (Schäfer 1994)

Markerless deletion via suicide vector + sucrose counter-selection:

1. Build deletion plasmid: ~500 bp upstream + ~500 bp downstream of the target gene flanking the multiple cloning site of pK18mobsacB (Addgene 87097).
2. Electroporate into Xcc isolate; select on kanamycin (single crossover).
3. Counter-select on sucrose; sacB lethality kills non-resolvers. PCR + sequencing confirm deletion.

Targets and expected phenotypes (from Hung 2003):

| Strain | Expected phiL7 outcome |
|--------|------------------------|
| WT | Sensitive |
| ΔtonB | Resistant |
| ΔexbB | Resistant |
| ΔexbD1 | Resistant |
| **ΔexbD2** | **Still sensitive** (negative control) |

Timeline: 4–6 weeks per gene, parallelizable.

### Validation tiers

| Tier | What you measure | Story |
|------|------------------|-------|
| 1 | ELISA only (WT host) | "We found variants that bind better." |
| 2 | + Plaque assay (WT) | "Binding → infection confirmed." |
| 3 | + Δreceptor panel | **"Receptor-specific causality."** (paper-grade) |

Recommendation (PI briefing 2026-05-11): commit to Tier 3 if knockouts start May 17.

### Failure modes and quality gates

| Failure | Action |
|---------|--------|
| Dry lab misses 48-h SLA | Wet lab uses `safe_pick_backup.csv` (pre-approved by PI) |
| Variant insoluble | Mark NaN; attempt backup truncation variant |
| ELISA R² < 0.9 | Down-weight in retraining |
| All 5 BALD picks fail QC | 1 expert pick + 2 random; retrain on partial data |
| Calibration ECE > 0.1 | Temperature scaling before next BALD run |

---

## 7. How to reproduce — environment, commands, tests

### Quick-start (~10 min total)

```bash
# 1. Clone + branch (active-learning-pipeline, NOT main)
git clone https://github.com/A1ex-Ch3n/iGEM_Claremont_2026.git
cd iGEM_Claremont_2026
git checkout active-learning-pipeline

# 2. Conda environment (one-time, ~5 min)
conda env create -f shared/env/environment.yml
conda activate igem2026
pip install nbstripout && nbstripout --install  # once per clone

# 3. Minimum genomes for local dev (~5 MB)
python 00_raw_data/processes/fetch_phages.py   --accession EU717894.1
python 00_raw_data/processes/fetch_phages.py   --accession NC_001604.1     # T7, used by tests
python 00_raw_data/processes/fetch_bacteria.py --accession GCF_000007145.1

# 4. One-time PhageRBPdetect data + HMM press (Module 03)
bash 03_rbp_identification/inputs/setup_inputs.sh

# 5. Launch JupyterLab
jupyter lab
```

The full 777-phage + 34-bacteria dataset (~630 MB) is gitignored. Fetch it only on Laguna before batch jobs.

### Per-module entry points

| Module | Open this first |
|--------|-----------------|
| 00 | `00_raw_data/processes/01_verify_dataset.ipynb` |
| 01 | `01_data_ground_truth/processes/01_fetch_reference_genomes.ipynb` |
| 02 | `02_annotation/processes/01_run_phanotate.ipynb` (phage) + `02_run_prodigal.ipynb` (host) |
| 03 | `03_rbp_identification/processes/01_run_phagerbpdetect.ipynb` |
| 04 | `04_protein_embedding/processes/01_embed_esm2.ipynb` |
| 05 | `05_structure_prediction/processes/01_run_boltz2.ipynb` (Laguna for production) |
| 06 | `06_uncertainty_model/processes/run_cycle0.py` (3 s on CPU; verified 2026-05-17) |
| 07 | `07_acquisition_function/processes/run_bald.py` (<1 s) |

### Running tests

```bash
pytest 00_raw_data/processes/tests/             -v    # 15+ pass (3 expected fails)
pytest 01_data_ground_truth/processes/tests/    -v    # 22/22
pytest 02_annotation/processes/tests/           -v    # 26/26
pytest 03_rbp_identification/processes/tests/   -v    # 25+ pass (2 expected — hmmpress)
pytest 04_protein_embedding/processes/tests/    -v    # 17/17
pytest 05_structure_prediction/processes/tests/ -v    # 28/28 (1 expected skip — GPU run)
pytest 06_uncertainty_model/processes/tests/    -v    # 9/9
pytest 07_acquisition_function/processes/tests/ -v    # 18/18
```

### Laguna HPC — when

| Task | Local? | Why Laguna? |
|------|--------|-------------|
| ESM-2 8M (3 RBPs) | Yes | seconds on CPU |
| ESM-2 650M / 3B (777 phages) | No | GPU memory + time |
| Boltz-2 protein-protein | No | ~15 min/pair, A100/L40S only |
| AF3 batch | No | GPU + weights approval (Google form) |
| Ensemble retrain + BALD | Yes | seconds on CPU |

Full Laguna recipe (CUDA 12.4, torch 2.5.1+cu121, conda `boltz2` env, NVIDIA L40S, OnDemand portal) lives in `LAGUNA.md`. The boltz install upgrades torch to cu130 — there's a one-line fix in `LAGUNA.md` to pin back to cu121.

### Rebuilding this onboarding deck

```bash
# Figures
python docs/onboarding/figures/make_figures.py

# Slides (English)
pandoc docs/onboarding/slides_en.md \
  -t beamer --pdf-engine=xelatex --slide-level=2 \
  -o docs/onboarding/slides_en.pdf

# Slides (Chinese)
pandoc docs/onboarding/slides_zh.md \
  -t beamer --pdf-engine=xelatex --slide-level=2 \
  -o docs/onboarding/slides_zh.pdf
```

One-time TeX dependencies: `tlmgr install beamer booktabs mdframed xecjk ctex pgfplots fontaxes needspace zref`. Mac CJK font: `Hiragino Sans GB` (built-in).

---

## 8. Conventions reference

### Identifier formats (from `INTERFACE.md`)

| Concept | Format | Example |
|---------|--------|---------|
| Phage NCBI accession | `<base>.<version>` | `EU717894.1` |
| Bacterial assembly | `GCF_*` / `GCA_*` | `GCF_000007145.1` |
| ORF id (per genome) | `<acc>_orf_<5-digit>` | `EU717894.1_orf_00031` |
| RBP candidate id | `<acc>_rbp_<2-digit>` | `EU717894.1_rbp_01` |
| Receptor id | `<host_acc>_<gene_name>` | `GCF_000007145.1_tonB` |
| Variant id | `<parent_rbp>_<change>_<idx>` | `EU717894.1_rbp_01_trunc_03` |

### Path anchoring

- Inside `.ipynb`: `REPO_ROOT = Path.cwd().resolve().parents[1]` (notebook lives in `<module>/processes/`).
- Inside `.py`: `REPO_ROOT = Path(__file__).resolve().parents[2]`.
- **Never** hard-code absolute paths. **Never** use `~`. **Never** assume CWD.

### Cycle versioning

Every model checkpoint, predictions CSV, and ELISA dataset is tagged `cycle_<N>` and lives in a per-cycle directory. Provenance (`run_meta.json`) includes the git SHA at run time. MLflow tracks runs (planned for Cycle 0 onward).

### Tool split (Module 02)

PHANOTATE for **phage** ORF calling; Prodigal / pyrodigal for **bacterial** hosts. **Never swap.** Prodigal assumes non-overlapping ORFs and loses ~10–15% of phage genes; PHANOTATE handles overlapping ORFs with dynamic programming.

---

## 9. References — compact per-module

| Module | Paper | Role |
|--------|-------|------|
| 00 / 01 | da Silva 2002 *Nature*; Lee 2009 *AEM*; Hung 2003 *BBRC* | reference genomes + receptor |
| 02 | McNair 2019 *Bioinformatics* (PHANOTATE); Hyatt 2010 (Prodigal); Bouras 2023 (pharokka) | gene calling |
| 03 | Boeckaerts 2022 *Viruses* (PhageRBPdetect) | RBP HMM + XGBoost |
| 04 | Lin 2023 *Science* (ESM-2); Liu 2025 *Nat Commun* (PLM-interact) | embeddings + PPI prior |
| 05 | Abramson 2024 *Nature* (AF3); Passaro 2025 (Boltz-2) | structure + ipTM |
| 06 | Lakshminarayanan 2017 *NeurIPS*; Greenman 2025 *PLoS Comput Biol* | deep ensembles + UQ audit |
| 07 | Houlsby 2011 (BALD); Yang 2025 *Nat Commun* (ALDE); Hie 2024 *Nat Biotechnol* | acquisition |
| 08 | Schäfer 1994 (pK18mobsacB); Gibson 2009; Latka 2021 *mBio* | wet lab methods |

Full annotated reading guide with 🔴/🟡/⚪ priority tags: `docs/reference/papers.md` (19 papers).

### Five literature corrections (May 2026 audit)

These were found and fixed during a full read of the 19 core papers (Alex, 2026-05-11). Don't reintroduce them.

1. **ExbD2 is NOT essential** — Hung 2003 (the real source) shows ΔexbD2 retains phiL7 sensitivity. The previous "Wang 2003 *Mol Microbiol*" citation was a hallucinated reference.
2. **Boltz-2 affinity head is small-molecule only** — outputs NaN for protein-protein pairs. Use ipTM as a structural confidence proxy, not as an affinity prior.
3. **Greenman 2025 journal is *PLoS Comput Biol*** (not *NAR Genomics*); the conclusion is *"no single best UQ method"*, not "deep ensembles win".
4. **Hie 2024 used ESM-1b/1v** (not ESM-2); each antibody saw **~20** variants (not ~50); mechanism is language-model likelihood filtering, not BALD-style closed loop.
5. **Lee 2009 explicitly searched and could not find** an OP1 ORF25 homolog in phiL7. Our rbp_01 is an HMM-based complementary discovery — *not* a contradiction with Lee 2009. HMMs detect structurally similar but sequence-diverged proteins that BLAST misses.

Full audit: `docs/reference/paper_reading_notes.md`.

---

## 10. Glossary + further reading

### Fast lookup

| Term | Meaning |
|------|---------|
| RBP | Receptor-binding protein — the phage's "key" to the host cell |
| HMM | Hidden Markov Model — sequence-profile method (Boeckaerts 2022) |
| ESM-2 | Evolutionary-Scale Modeling v2 — protein language model (Lin 2023) |
| ipTM | Interface predicted TM-score (0–1, structural confidence) |
| ELISA | Enzyme-linked immunosorbent assay (binding readout) |
| EC50 | Half-maximal effective concentration (from 4PL fit) |
| SDM | Site-directed mutagenesis (NEB Q5) |
| BALD | Bayesian Active Learning by Disagreement (Houlsby 2011) |
| epistemic | Reducible model uncertainty (= BALD target) |
| aleatoric | Irreducible measurement noise |
| pK18mobsacB | Suicide vector for markerless deletion (Schäfer 1994) |

### Further reading

- `docs/reference/glossary.md` — full glossary.
- `docs/planning/iGEM_2026_Project_Plan.md` — the PI-facing English project plan.
- `docs/planning/iGEM_2026_项目大纲_中文版.md` — Chinese version with deeper dry-lab module mechanics + the 6-layer data-scarcity strategy.
- `docs/planning/PI_briefing_2026-05-11.md` — bilingual status briefing with every output path and the work log (May 7–12).
- `docs/protocols/` — wet-lab Benchling SOPs (cultivation, transformation, plaque, infection curves, lysate amp).
- `INTERFACE.md` — locked data contract between modules.
- `LAGUNA.md` — HPC setup, SLURM templates, CUDA gotchas.

---

## Appendix — Reproduction one-liner

```bash
# From a fresh checkout — runs Modules 06 and 07 end-to-end (~5 seconds)
conda activate igem2026
python 06_uncertainty_model/processes/run_cycle0.py
python 07_acquisition_function/processes/run_bald.py --cycle 1 --n_bald 4 --n_random 1
# Inspect:
head 06_uncertainty_model/outputs/cycle_0/predictions.csv
head 07_acquisition_function/outputs/cycle_1/recommendations.csv
```

That's the active-learning loop in two commands. The hard part starts on June 1.
