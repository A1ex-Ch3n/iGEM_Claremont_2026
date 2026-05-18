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

> [!info] Key vocabulary for this section / 本节关键术语
> - **Phage (bacteriophage)** — a virus that infects bacteria. We work with *lytic* phages — they kill the host after replicating. (Contrast: *lysogenic* phages integrate into the genome and stay dormant.)
> - **RBP (receptor-binding protein)** — the protein on a phage's tail that latches onto a specific bacterial surface receptor. It's the lock-and-key step that decides whether the phage can even attach.
> - **Active learning** — a way of training a machine-learning model where, after each round of training, *the model* picks which experiment to run next, instead of a human picking arbitrarily. We use it because each ELISA measurement is expensive — we want every one of them to be maximally informative.
> - **Closed-loop** — the model recommends → wet lab measures → result feeds back → model retrains. The loop closes itself.
> - **ELISA** — Enzyme-Linked Immunosorbent Assay. A standard 96-well plate experiment that yields a binding curve. The numerical target our model predicts.
> - **BALD** — the algorithm (an *acquisition function*) that ranks every unmeasured variant by how much measuring it would reduce the model's uncertainty. Detailed in §4.
> - **Epistemic uncertainty** — the kind of uncertainty that *more data could reduce*. Different from random measurement noise. BALD targets this.

We are building a **closed-loop active-learning pipeline** that pairs a machine-learning model of phage RBP (receptor-binding protein) × bacterial-receptor interaction with iterative wet-lab validation. After each ELISA round, the model retrains, the BALD acquisition function ranks every untested variant by epistemic uncertainty, and the wet lab measures the next 4–5. The whole loop is designed to extract maximum information from each expensive ELISA measurement — addressing the central pain point of phage-host prediction: **data scarcity**.

Reference dry-lab scaffold:
- Host: *Xanthomonas campestris* pv. *campestris* (Xcc) ATCC 33913 — NCBI `GCF_000007145.1`.
- Phage: phiL7 — NCBI `EU717894`.
- Receptor system: TonB-ExbB-ExbD1 essential, ExbD2 **not** required (Hung 2003, *BBRC* 302:878–884, PMID 12646254).

Wet lab self-isolates *Xanthomonas* + phage from California brassicas (per PI consultation 2026-05-07), bypassing the multi-month USDA APHIS PPQ-526 permit.

iGEM tracks targeted: **Best Agriculture Project · Best Model · Best Composite Part**.

---

## 2. The science — why we're doing this

> [!info] Key vocabulary for this section / 本节关键术语
> - **Pathovar (pv.)** — a strain of a plant-pathogenic bacterium defined by which plant it infects. *X. campestris* pv. *campestris* (Xcc) infects brassicas; *X. campestris* pv. *vesicatoria* infects pepper. Same species, different hosts.
> - **Lytic phage** — a phage that bursts (lyses) the host cell after replicating. The kind we want for biocontrol.
> - **Siphoviridae** — a phage family with long, flexible, non-contractile tails. phiL7 belongs to this family.
> - **AUC** — short for Area Under the ROC Curve. A single number from 0.5 (no better than random) to 1.0 (perfect) that summarises how well a model can tell positive pairs apart from negatives. 0.82 ≈ "right about 82 % of the time."
> - **BLAST** — a sequence-search method that asks "find me proteins whose amino-acid letters are similar to my query." Fast and standard, but it misses proteins that have diverged so much that the letters no longer match — even when the proteins still have the same shape and function.
> - **HMM (Hidden Markov Model)** — a more powerful sequence search that learns *which positions* in a protein family matter and which are free to change. Catches proteins BLAST misses. Module 03 uses HMMs to find rbp_01.
> - **TonB / ExbB / ExbD1** — three bacterial proteins that together form an energy-coupled import channel in the outer membrane. phiL7 hijacks this channel as its entry point. Identified by Hung 2003.

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

> [!info] Key vocabulary for this section / 本节关键术语
> - **Module** — a numbered subdirectory in this repo (`00_…` through `08_…`). Each module does one job; outputs from module *N* become inputs to module *N+1*.
> - **ORF (open reading frame)** — a stretch of DNA that codes for a single protein, marked by a start codon (ATG) and a stop codon. The first job of any genome pipeline is to find all the ORFs.
> - **Gitignored** — listed in the project's `.gitignore` file so git doesn't track it. We exclude large binary files (raw genomes, model weights, structure predictions) from git but record their identity in a `MANIFEST.csv`, so anyone can re-fetch the exact same set.
> - **MANIFEST.csv** — a per-output-folder ledger (filename, SHA-256 checksum, size, record count) that makes gitignored artifacts reproducible — you can verify you're looking at the same file Alex was.
> - **Conda environment** — a project-specific Python install with every library pinned to an exact version. `conda activate igem2026` switches into ours.
> - **Jupyter notebook** — an interactive document with code cells, output, and prose mixed together. Good for exploration. Once code is stable we "freeze" it into a plain `.py` module.
> - **Bilingual comments** — every code comment in this repo is written in both English and Chinese, so the whole team can read it.

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

> [!info] Key vocabulary for this section / 本节关键术语
> - **Neural network** — a chain of mathematical functions stacked into "layers" that turns one vector of numbers (the input) into another (the prediction). Trained by adjusting the numbers inside each layer to minimise a *loss* — how wrong the predictions are on training examples.
> - **MLP (multilayer perceptron)** — the simplest neural-network architecture. It's just a stack of fully-connected layers separated by a *nonlinearity*. What our Module 06 uses.
> - **Layer / hidden layer** — one step in the chain. A "layer" applies a matrix multiplication plus a constant offset to its input. A *hidden* layer is any layer whose output isn't the final prediction — it's an intermediate representation the network learned. Our MLP has 3 hidden layers of widths 256, 256, 128.
> - **ReLU** — the standard nonlinearity. It keeps positive numbers unchanged and replaces negative numbers with zero. Without a nonlinearity between layers, stacking layers would mathematically collapse into a single linear function — no learning capacity beyond a straight line.
> - **Dropout** — during training, randomly turn off a fraction of the connections in a layer each step. This forces the network to spread information across many connections instead of memorising the training set; it generalises better to unseen variants.
> - **Gaussian distribution** — the bell curve. Fully described by two numbers: a **mean** (where it's centred) and a **standard deviation** (σ, how spread out it is). When we say "the model predicts a Gaussian," we mean it outputs *two* numbers per input — its best guess and how confident it is.
> - **NLL (negative log-likelihood) loss** — the loss function for a probabilistic model. It rewards the model for putting high probability on the *true* outcome. Being wrong-but-uncertain is penalised less than being wrong-and-overconfident.
> - **Deep ensemble** — train the *same* network architecture several times with different random initialisations. Where the data is informative, the members converge to similar predictions; where the data is sparse, they disagree. The disagreement is *epistemic uncertainty*.
> - **Calibration / ECE** — does a "90 % credible interval" actually contain the true value 90 % of the time? Expected Calibration Error (ECE) measures the gap. Ensembles tend to be accurate but not always well-calibrated; we monitor ECE and patch with temperature scaling if it drifts.
> - **Mutual information** — a number measuring "how much does knowing $X$ tell me about $Y$?" Zero if they're independent. BALD picks the experiment whose result would give the most mutual information about the model's parameters.
> - **Epistemic vs aleatoric** — two flavours of uncertainty. *Epistemic* is the model's ignorance — it *shrinks with more training data*. *Aleatoric* is irreducible measurement noise (e.g., ELISA pipetting variance from one well to another). Only epistemic uncertainty is fixable by collecting more data, which is why BALD targets epistemic only.

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

4 BALD picks + 1 random control (Hie 2024 control-arm pattern). These are synthetic placeholders — real Cycle 1 = rbp_01 variants × TonB after Cycle 0 ELISA.

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

> [!info] Key vocabulary for this section / 本节关键术语
> - **Structure prediction** — given an amino-acid sequence (or two), predict the protein's 3D shape: where each atom sits in space. Replaces months of X-ray crystallography with a "good enough" starting point.
> - **Boltz-2 / AlphaFold 3** — two state-of-the-art structure-prediction tools. Boltz-2 is open-weight; AF3 requires applying to Google.
> - **PDB file** — the standard plain-text format for atomic 3D structures (one line per atom). Open in **PyMOL** (free) or **ChimeraX** (free) to visualise.
> - **ipTM** — *interface predicted TM-score.* A single number from 0 to 1 that tells you how confident the model is about *the interface between two chains*. Rough scale: ≥ 0.6 = trustworthy interface; 0.4–0.6 = ambiguous; < 0.4 = the model is essentially guessing how the two proteins dock together.
> - **chain pTM** — the same TM-score idea but for *a single chain on its own*. Tells you whether the monomer is well-predicted, independent of how the chains stick together. Our `chain_A_ptm = 0.808` means rbp_01's monomer structure is reliable.
> - **PAE (predicted aligned error)** — a 2D matrix. Entry [i, j] answers: "if I align my prediction at residue *i*, how far off (in Ångström) am I likely to be at residue *j*?" Low PAE = high confidence in the *relative* position of those two residues. The off-diagonal block between two chains tells you about the interface.
> - **pLDDT** — per-residue confidence in the local geometry of the predicted backbone (range 0–100 or 0–1 depending on the tool). High pLDDT = the model is sure about that residue's local position.

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

> [!info] Key vocabulary for this section / 本节关键术语
> - **Gibson assembly** — a one-tube reaction (5′ exonuclease + polymerase + ligase) that joins overlapping DNA fragments into a single circular plasmid. How we build new RBP-variant constructs in Cycle 0.
> - **SDM (site-directed mutagenesis)** — making a single specific change to an existing plasmid (e.g., changing residue K450 to alanine). Cheaper and faster than ordering a whole new gene from a vendor — the natural fit for the small, targeted mutations BALD picks.
> - **BL21(DE3) / pET-28a / IPTG** — the standard *E. coli* recipe for producing recombinant protein. BL21(DE3) is the bacterial strain; pET-28a is the plasmid backbone; **IPTG** is the chemical that flips on the gene's expression.
> - **His6-tag + Ni-NTA** — we tack six histidine residues onto the protein's N-terminus. Those histidines bind tightly to a nickel-charged resin, so we can pull our protein out of a cell lysate in a single column step (called **immobilised metal affinity chromatography**, IMAC).
> - **SDS-PAGE** — gel electrophoresis that separates proteins by size in a polyacrylamide gel. Used to verify purity and identify your protein's band.
> - **EC50 / 4PL fit** — fit the ELISA binding curve (signal vs. RBP concentration) to a 4-parameter logistic. **EC50** is the concentration at which the binding signal is half-maximal — our quantitative "binding strength" and the target variable the model predicts.
> - **Plaque assay / PFU** — drop diluted phage onto a lawn of bacteria; each infectious phage particle kills the cells around it, leaving a clear spot called a **plaque**. Counting plaques gives **PFU/mL** (plaque-forming units per millilitre) — how infectious your phage stock is.
> - **pK18mobsacB / sacB selection** — a "suicide plasmid" for making clean gene deletions in *Xanthomonas*. The trick: the plasmid carries `sacB`, which makes cells die in sucrose. So we can select for cells that have *lost* the plasmid through a second recombination event — leaving only the desired in-frame deletion behind.
> - **Electroporation** — applying a brief electric pulse to make bacterial cell walls temporarily permeable so plasmid DNA can enter.

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

> [!info] Key vocabulary for this section / 本节关键术语
> - **Branch** — a parallel line of git history. All pipeline code lives on the `active-learning-pipeline` branch; `main` is intentionally empty. If you check out `main` and see nothing, that's expected — switch branches.
> - **Conda environment** — a Python install with every library pinned to an exact version. `conda activate igem2026` switches you into ours; until you do, you're using the system Python and will see import errors.
> - **pytest** — Python's standard test runner. Each module ships a small test suite that locks in its behaviour. A passing run is the contract — if you change code and the tests still pass, you didn't break anything downstream.
> - **JupyterLab** — the browser-based interface for opening and running `.ipynb` notebooks.
> - **Laguna (CARC)** — the USC Center for Advanced Research Computing's GPU cluster. We use it for jobs too big or slow for a laptop: Boltz-2 structure prediction, ESM-2 650M / 3B inference, AlphaFold 3.
> - **SLURM** — the job scheduler that Laguna runs. You submit a small shell script with `sbatch`; SLURM finds a free GPU node and runs it there. Templates live in `scripts/`.
> - **CUDA** — Nvidia's GPU compute platform. Its version has to match between the cluster's GPU driver and PyTorch's build — currently we pin to CUDA 12.1 / `torch==2.5.1+cu121`. `LAGUNA.md` has the recipe.

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

> [!info] Key vocabulary for this section / 本节关键术语
> - **NCBI accession** — a unique string identifying a sequence record in the NCBI public database. `EU717894.1` is phiL7's genome; `GCF_000007145.1` is Xcc's. The `.1` is a version suffix.
> - **REPO_ROOT** — a Python variable we set in every script so file paths are anchored to the project root (wherever you cloned the repo), not to your shell's current directory. That makes scripts work no matter where you run them from.
> - **MANIFEST.csv** — a per-output-folder ledger (filename, SHA-256 checksum, byte size, record count, UTC timestamp) that makes large gitignored outputs reproducible. If your SHA-256 matches the manifest, you have the same file.
> - **INTERFACE.md** — the locked data contract between modules. Defines identifier formats, FASTA header conventions, CSV column schemas, and what counts as a breaking change.
> - **MLflow** — a tracking system for ML experiments: every model checkpoint, hyperparameter, and metric is logged with a unique run ID. Planned for Cycle 0 onward.

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

---

## Citation list

Full bibliographic details for every paper mentioned in this guide. Project shorthand (e.g., "Hung 2003") shown in parentheses; the canonical annotated reading guide with priority tags lives at `docs/reference/papers.md`. Where the project's literature audit (May 2026, `docs/reference/paper_reading_notes.md`) found a factual issue with how a paper had previously been cited, the correction is noted inline.

### Receptor biology & reference genomes

1. **da Silva, A.C.R. et al.** (2002). "Comparison of the genomes of two *Xanthomonas* pathogens with differing host specificities." *Nature* **417**, 459–463. (`da Silva 2002`) — source of the Xcc ATCC 33913 reference genome (NCBI `GCF_000007145.1` / GenBank `AE008922`).
2. **Lee, C.-N. et al.** (2009). "Genomic characterization of the intron-containing T7-like phage phiL7 of *Xanthomonas campestris* pv. *campestris*." *Applied and Environmental Microbiology* **75**(24), 7828–7838. NCBI accession EU717894. (`Lee 2009`) — phiL7 genome paper. *Audit note:* Lee 2009 explicitly searched for and was unable to find an OP1 ORF25 tail-fiber homolog; our rbp_01 is an HMM-based complementary discovery, not a contradiction.
3. **Hung, C.-H. et al.** (2003). "Involvement of *tonB-exbBD1D2* operon in infection of *Xanthomonas campestris* phage ϕL7." *Biochemical and Biophysical Research Communications* **302**(4), 878–884. PMID 12646254. DOI: 10.1016/S0006-291X(03)00255-9. (`Hung 2003`) — TonB, ExbB, and ExbD1 are essential for phiL7 entry; ExbD2 is **not** required (strain CH620). *Audit note:* an earlier project document mis-cited "Wang 2003" as the source; that paper cannot be verified to exist.

### *Xanthomonas* pathology & phage biocontrol

4. **Ryan, R.P. et al.** (2011). "*Xanthomonas* genomics and molecular plant–microbe interactions." *Nature Reviews Microbiology* **9**, 344–355. (`Ryan 2011`) — >400-species host range; economic-impact context.
5. **Aiello, D. et al.** (2019). Reports of copper resistance in plant-pathogenic *Xanthomonas* spp. *Plant Disease*. (`Aiello 2019`) — supports the "copper resistance widespread" claim.
6. **Iriarte, F.B. et al.** (2018). "Combination of plant defense elicitors and bacteriophage for biocontrol of bacterial spot of tomato." *Frontiers in Plant Science* **9**, 1–12. (`Iriarte 2018`) — field-trial efficacy for phage biocontrol on tomato.
7. **Holtappels, D. et al.** (2022). "The future of phage biocontrol in integrated plant protection for sustainable crop production." *Microbial Biotechnology* **15**(3), 597–610. (`Holtappels 2022`) — Xcc biocontrol context.
8. **Farquharson, E.L. et al.** (2021). "Phage resistance is driven by reduced infection efficiency of receptor mutants." (`Farquharson 2021`) — T4 × *E. coli*: RBP binds 85 % of strains, plaques on only 11 %; the "binding ≠ infection" confounder our Layer 2 knockouts address.

### Phage host-range prediction & RBP identification

9. **Boeckaerts, D. et al.** (2022). "Identification of phage receptor-binding protein sequences with hidden Markov models and an extreme gradient boosting classifier." *Viruses* **14**(6), 1329. (`Boeckaerts 2022`) — PhageRBPdetect, the tool we use in Module 03. Precision-recall AUC 93.8 %, F1 84.0 %.
10. **Boeckaerts, D. et al.** (2024). "Prediction of *Klebsiella* phage-host specificity at the strain level (PhageHostLearn)." *Nature Communications* **15**, art. 4768 (48675). DOI: 10.1038/s41467-024-48675-6. (`Boeckaerts 2024`) — strain-level SOTA. *Audit note:* the headline AUC of 0.82 is "up to 81.8 % at 100 % identity threshold"; cross-strain it drops to ~0.70.
11. **Mutalik group** (2025). "Phage Anti-Microbial Landscape (PAML) benchmark." *bioRxiv* (preprint). (`Mutalik 2025`) — independent confirmation that PhageHostLearn's cross-genus AUC drops to 0.67–0.70.

### Genome annotation tools

12. **McNair, K. et al.** (2019). "PHANOTATE: a novel approach to gene identification in phage genomes." *Bioinformatics* **35**(22), 4537–4542. (`McNair 2019`) — dynamic-programming ORF caller for phages.
13. **Hyatt, D. et al.** (2010). "Prodigal: prokaryotic gene recognition and translation initiation site identification." *BMC Bioinformatics* **11**, 119. (`Hyatt 2010`) — bacterial ORF caller (used via pyrodigal binding).
14. **Bouras, G. et al.** (2023). "Pharokka: a fast scalable bacteriophage annotation tool." *Bioinformatics* **39**(1), btac776. (`Bouras 2023`) — PHROG / VFDB / CARD functional annotation for phages.

### Protein language models and structure prediction

15. **Lin, Z. et al.** (2023). "Evolutionary-scale prediction of atomic-level protein structure with a language model." *Science* **379**(6637), 1123–1130. (`Lin 2023`) — ESM-2. Trained on ~65 M unique sequences; 650M params → 1280-D embedding; 3B params → 2560-D.
16. **Liu, D. et al.** (2025). "PLM-interact: extending protein language models to predict protein–protein interactions." *Nature Communications* **16**, art. 64512. DOI: 10.1038/s41467-025-64512-w. (`Liu 2025`) — AUPR +16–28 % when transferred to mouse / fly / worm / yeast / *E. coli* PPI. *Audit note:* not tested on phage–bacteria; our project may be the first such transfer.
17. **Abramson, J. et al.** (2024). "Accurate structure prediction of biomolecular interactions with AlphaFold 3." *Nature* **630**, 493–500. DOI: 10.1038/s41586-024-07487-w. (`Abramson 2024`) — AF3; requires Google form for weights.
18. **Passaro, J.M. et al.** (2025). "Boltz-2: towards accurate and efficient binding affinity prediction." *bioRxiv* (preprint). (`Passaro 2025`) — what we ran on Laguna. *Audit note:* affinity head is trained on small-molecule × protein data (PubChem / ChEMBL / BindingDB) and outputs `NaN` for protein–protein pairs; we use **ipTM** as a structural-confidence proxy, NOT a quantitative affinity.

### Uncertainty quantification

19. **Lakshminarayanan, B. et al.** (2017). "Simple and Scalable Predictive Uncertainty Estimation Using Deep Ensembles." *NeurIPS* **30**. arXiv:1612.01474. (`Lakshminarayanan 2017`) — the 5-MLP architecture and Gaussian NLL loss in Module 06.
20. **Ovadia, Y. et al.** (2019). "Can You Trust Your Model's Uncertainty? Evaluating Predictive Uncertainty Under Dataset Shift." *NeurIPS* **32**. (`Ovadia 2019`) — deep ensembles vs MC Dropout vs temperature scaling under distribution shift.
21. **Greenman, K.P. et al.** (2025). "Benchmarking uncertainty quantification for protein engineering." *PLOS Computational Biology* **21**(1), e1012639. DOI: 10.1371/journal.pcbi.1012639. (`Greenman 2025`) — *Audit note:* journal is *PLoS Comput Biol* (NOT *NAR Genomics & Bioinformatics*); conclusion is "no single best UQ method," NOT "deep ensembles are best."

### Active learning and directed evolution

22. **Lindley, D.V.** (1956). "On a measure of the information provided by an experiment." *Annals of Mathematical Statistics* **27**, 986–1005. (`Lindley 1956`) — Bayesian optimal experimental design.
23. **Settles, B.** (2009). *Active Learning Literature Survey.* Computer Sciences Technical Report 1648, University of Wisconsin–Madison. (`Settles 2009`) — canonical AL reference.
24. **Houlsby, N. et al.** (2011). "Bayesian Active Learning for Classification and Preference Learning." arXiv:1112.5745. (`Houlsby 2011`) — BALD. *Audit note:* original paper applies BALD to a Gaussian Process Classifier (classification); extending the information-theoretic objective to deep-ensemble *regression* is our extension, not a direct citation.
25. **Hie, B.L. et al.** (2024). "Efficient evolution of human antibodies from general protein language models." *Nature Biotechnology* **42**, 275–283. Preprinted April 2022 on bioRxiv; published January 2024. (`Hie 2024`) — *Audit note:* uses **ESM-1b / ESM-1v** (not ESM-2); ~20 variants per antibody (not ~50); mechanism is language-model-likelihood filtering, not a BALD closed loop. Earlier project documents that cite "Hie 2022 *Cell*" refer to this same paper.
26. **Yang, J. et al.** (2025). "Active learning-assisted directed evolution (ALDE)." *Nature Communications* **16**, art. 55987. DOI: 10.1038/s41467-025-55987-8. (`Yang 2025`) — *Audit note:* ALDE uses **Thompson sampling** acquisition and **one-hot** encoding (not BALD + ESM-2). It validates "AL + UQ + DNN ensembles work in protein engineering" but does NOT validate BALD specifically.
27. **Wittmann, B.J. et al.** (2021). "Informed training set design enables efficient machine learning-assisted directed evolution." *Cell Systems* **12**(11), 1026–1045. (`Wittmann 2021`) — exploration vs exploitation tradeoff analysis (referenced in the project plan, supporting the "4 BALD + 1 random" batch design).

### Wet lab methods

28. **Schäfer, A. et al.** (1994). "Small mobilizable multi-purpose cloning vectors derived from the *Escherichia coli* plasmids pK18 and pK19." *Gene* **145**(1), 69–73. (`Schäfer 1994`) — pK18mobsacB; Addgene #87097.
29. **Gibson, D.G. et al.** (2009). "Enzymatic assembly of DNA molecules up to several hundred kilobases." *Nature Methods* **6**, 343–345. (`Gibson 2009`) — Gibson assembly for variant cloning.
30. **Studier, F.W. & Moffatt, B.A.** (1986). "Use of bacteriophage T7 RNA polymerase to direct selective high-level expression of cloned genes." *Journal of Molecular Biology* **189**(1), 113–130. (`Studier & Moffatt 1986`) — original pET / BL21(DE3) expression system.
31. **Latka, A. et al.** (2021). "Engineering the modular receptor-binding proteins of *Klebsiella* phages switches their *in vitro* host range." *mBio* **12**(6), e02329-21. (`Latka 2021`) — RBP modular truncation strategy; ELISA assay format reference. *Audit note:* the *Klebsiella* RBP is a CPS depolymerase — the structural N-anchor / C-specificity-head principle transfers, but the biochemical mechanism does not directly apply to phiL7 × TonB.

---

## Document log

A short record of how this file evolved. Larger context (the pre-build plan, the AGENT_PROMPT, the speaker prep) is also in `docs/onboarding/`.

| Date (UTC) | Author | Change |
|------------|--------|--------|
| 2026-05-17 | Claude (onboarding agent), reviewed by Alex Chen | Initial draft. 10 sections + reproduction appendix; figures generated via `figures/make_figures.py`; rendered companion deck `slides_en.pdf`; demo runbook `DEMO.md`. |
| 2026-05-17 | Claude (onboarding agent) | Path-verification pass: corrected `00_raw_data/phage_genomes/` → `00_raw_data/phage/<acc>/`; replaced raw `hmmpress` call with `setup_inputs.sh`. Committed as `5304f69` on `active-learning-pipeline`. |
| 2026-05-18 | Claude (onboarding agent) | Added Obsidian-style `> [!info]` callouts at the start of each of sections 1–8, scaffolding key vocabulary for first-time readers (no CS background assumed). Section 4 (the ML core) gets the deepest vocabulary box — neural network / MLP / hidden layer / ReLU / dropout / Gaussian / NLL / deep ensemble / calibration / ECE / mutual information / epistemic vs aleatoric. |
| 2026-05-18 | Claude (onboarding agent) | Drive-by fix: unified "Hie 2022" → "Hie 2024" (same paper; preprint 2022 → published 2024, per audit). Added this Citation list (31 entries) and Document log. |

