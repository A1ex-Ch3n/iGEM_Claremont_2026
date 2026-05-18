# CONCEPT_DEEPDIVE.md — Seven concepts Alex must own

For each concept: 30-second pitch (PI/wet lab); 2-minute (dry lab); full mechanism; analogies; one verbatim sentence Alex should memorize; known failure modes.

---

## 1. Why active learning solves data scarcity (vs random / supervised)

### 30-second pitch
Every ELISA costs ~$50 and 4 days. If we pick variants at random, most measurements land in regions the model already understands. Active learning instead picks variants where the model is **most uncertain** — those measurements teach the model the most. In adjacent fields (Hie 2024 antibodies, Yang 2025 enzymes), this gives 2–5× more progress per experiment.

### 2-minute explanation
A supervised model trained on a fixed dataset is stuck with whatever bias that dataset has. Active learning closes the loop: the model uses its **own uncertainty** as a signal for what to measure next. Formally, you have a posterior over model parameters `p(θ | D)` after seeing data `D`. The next measurement at input `x` will be informative if `y(x)` reduces the entropy of that posterior. That's the optimal experimental design objective Lindley (1956) wrote down. BALD operationalizes it: pick `x` to maximize mutual information `I(y; θ | x, D)`. In practice, that means "pick where ensemble members disagree most." Compared to random: random hits informative variants with low probability when the variant space is huge (20^712 for a 712-aa protein). Active learning hits them by design — Yang 2025 showed yield going from 12% to 93% in two cycles where random would have taken many more.

### Full mechanism
1. Train model on `D`. Get a posterior `p(θ | D)`.
2. For each candidate `x` in the pool, compute `acq(x) = I(y; θ | x, D)`. Closed form for Gaussians.
3. Pick top-K by `acq(x)`.
4. Measure them in the lab — get `y_obs`.
5. Update: `D ← D ∪ {(x, y_obs)}`. Retrain. Repeat.
The "closed loop" is the iteration. The "active" is that the model — not the human — picks `x`.

### Analogies
- **Doctor diagnosing**: don't order every test in the catalogue; order the one whose result would most narrow the diagnosis.
- **20 Questions**: each question should split the remaining possibilities in half, not test arbitrary hypotheses.
- **Sonar in fog**: ping where you're least sure — not where you already see clearly.

### Memorize verbatim
> "Active learning picks the next experiment that maximally reduces the model's uncertainty — so every ELISA buys us as much information as possible."

### Where it could fail
- **Calibration is wrong** — model is "confident" where it shouldn't be, so BALD picks the wrong place. Mitigation: temperature scaling, monitor ECE per cycle.
- **Model bias is fundamental** — wrong feature space (e.g., ESM-2 misses something receptor-specific). BALD then asks the wrong questions efficiently. Mitigation: random control arm + Tier-3 causal validation.
- **Aleatoric noise dominates** — if ELISA noise is too high, no epistemic signal stands out. Mitigation: 3 technical replicates per plate, ECE-gated retrain.

---

## 2. ESM-2 — what a protein language model embedding actually IS

### 30-second pitch
ESM-2 is a neural network trained to "fill in the blank" on ~65 million protein sequences — like GPT, but for proteins. After training, you feed it any new protein and pull out the model's internal numbers — a 1280-dimensional **embedding** that captures what the protein "looks like" to a model that has seen 65M others. Similar-function proteins land near each other in that 1280-D space, even when their amino-acid letters don't match.

### 2-minute explanation
ESM-2 (Lin 2023, *Science*) is a transformer trained on UniRef50/90 with **masked language modelling**: randomly hide 15% of amino acids, ask the model to predict them from the context. To do that well, the model has to internally represent biochemistry — hydrophobic clusters, secondary structure motifs, evolutionary conservation. After training, you don't actually use the predicted amino acids. You use the **per-residue hidden states** — these are 1280-D vectors (for the 650M model) that have learned to encode position-aware structural and functional context. We **mean-pool** across residues to get one vector per protein. For our pipeline, that vector replaces hand-crafted features (charge, hydrophobicity, etc.) and is the input to the deep ensemble.

### Full mechanism
Architecture: 33-layer transformer (650M params), 1280-D hidden state, learned positional encoding. Training objective: cross-entropy on masked amino acid prediction. Training data: UniRef50/90, ~65M unique sequences seen during training. Inference: forward pass through encoder, extract `hidden_state[layer=33]` of shape `(L, 1280)`. Mean-pool over L: `embedding = hidden_state.mean(dim=0) → (1280,)`. We then concatenate `embedding_RBP ⊕ embedding_receptor → (2560,)` as the input to the MLP.

### Analogies
- **Spotify recommendation embeddings**: songs you've never heard land near similar ones based on patterns the model learned from millions of other songs.
- **Word embeddings (word2vec)**: "king − man + woman ≈ queen" — the geometry of the space encodes semantics.
- **Spectrogram of a song**: the raw waveform is uninformative; the frequency-time representation makes patterns visible.

### Memorize verbatim
> "ESM-2 is a transformer trained to fill in masked amino acids on 65 million proteins — what we use is the internal representation, not the predictions, because that representation captures structural and functional context."

### Where it could fail
- **Domain shift**: ESM-2 trained on UniRef which is mostly cellular proteins. Phage tail spikes are unusual structures (β-helices, trimer interlocks) — embedding may underweight what matters.
- **Mean-pooling discards positional info**: a point mutation 200 aa from the binding site has the same global mean as a mutation in the binding site. Per-residue or attention-weighted alternatives exist.
- **Scale matters but is sublinear**: 8M → 650M gives a big jump; 650M → 3B is smaller. PLM-interact (Liu 2025) suggests fine-tuning on PPI data is more impactful than raw scale.

---

## 3. Boltz-2 ipTM — the BIGGEST misconception risk

### 30-second pitch
Boltz-2 predicts 3D structures of protein complexes and gives a score called **ipTM** (interface predicted TM-score) between 0 and 1. **ipTM is a confidence measure** — how sure the model is about how the two chains dock. It is **NOT a binding affinity**. Boltz-2 does have an affinity head, but it was trained on **small molecules only** (PubChem, ChEMBL, BindingDB) — for protein-protein pairs the affinity output is `NaN`. Our ipTM for rbp_01 × TonB is 0.365: low confidence in the interface geometry. That doesn't mean "they don't bind" — it means the model doesn't know how.

### 2-minute explanation
TM-score (template modeling score) measures structural similarity, 0–1, where >0.5 typically indicates the same fold. **pTM** is TM-score predicted from the structure model's internal confidence — a self-assessment of "if I had the true structure, how similar would my prediction be?" **ipTM** restricts that to the **interface region** between two chains in a complex — how well-predicted is the docking geometry. Boltz-2 (Passaro 2025) reports ipTM as part of `affinity.json`. The model was trained largely on PDB structures plus distillation; for protein-protein it produces a structural prediction with ipTM. Separately, Boltz-2 has an **affinity head** that predicts binding free energy, but this head was trained on small-molecule binding databases. When we feed a protein-protein pair, the affinity head outputs `NaN` because there's no learned mapping. Our `predicted_dG = null`. This is documented in the Boltz-2 paper and confirmed in our build report.

### Full mechanism
- Boltz-2 input: two amino acid sequences (or one + a ligand).
- Model: diffusion-based atomic structure prediction with confidence head.
- Outputs:
  - `*.pdb` — atomic coordinates of the predicted complex
  - `chain_X_ptm` — per-chain monomer pTM (rbp_01: 0.808 = trustworthy fold)
  - `interface_ipTM` — interface confidence (0.365 for us = low)
  - `confidence_score` — weighted overall (0.683 for us)
  - `pae_*.npz` — pairwise alignment error matrix (per-residue confidence in relative positions)
  - `plddt_*.npz` — per-residue local confidence (0.30–0.98, mean 0.76)
  - `predicted_dG = null` for protein-protein

For our rbp_01 × TonB: ipTM 0.365 means the model is genuinely uncertain where on TonB rbp_01 docks. The monomers are well-predicted (chain pTMs high). This is consistent with the biology — phiL7 is known to enter through TonB (Hung 2003), but the specific binding interface has never been crystallized.

### Analogies
- **GPS confidence circle**: ipTM is the radius of "I'm 95% sure you're somewhere in here" — small circle = high ipTM. The circle says nothing about whether you'll actually find a parking spot.
- **Weather forecast probability**: "70% chance of rain" is a confidence; it's not "how much rain". ipTM is confidence in the prediction; it's not how strong the binding will be.
- **Identity parade**: ipTM 0.365 is "I think it's one of these three faces" — the model has narrowed it but not committed.

### Memorize verbatim
> "Boltz-2's affinity head is small-molecule only — it outputs NaN for protein-protein pairs. We use ipTM as a structural confidence proxy, not as binding affinity."

### Where it could fail
- **Mis-citing ipTM as affinity** — biggest credibility risk; rehearse the verbatim sentence.
- **Low ipTM ≠ no binding** — Hung 2003 confirmed phiL7 enters via TonB. Low ipTM means we can't predict geometry, not that there's no interaction.
- **High ipTM ≠ strong binding** — even a 0.9 ipTM tells us the docking is confident, not the affinity. Confidence and strength are orthogonal.

---

## 4. Deep ensemble — why 5 networks give epistemic uncertainty

### 30-second pitch
We train the same neural network architecture 5 times with different random starts. The 5 trained models agree where data tells them what to do, and **disagree** where data doesn't. The spread of their predictions = epistemic uncertainty = "what the model doesn't know." It's the simplest and most reliable way to quantify model ignorance in a neural net.

### 2-minute explanation
A single neural network gives a point prediction — no honest uncertainty signal. **Bayesian neural nets** put a distribution over weights but are computationally hard. **MC Dropout** approximates that but is poorly calibrated (Ovadia 2019). **Deep ensembles** (Lakshminarayanan 2017) are dead simple: train K models with different random seeds; the variance across their predictions approximates the posterior predictive variance. In our setup, each member outputs a **Gaussian** — both `μ_k(x)` and `σ_k(x)`. The **epistemic** variance is `Var_k[μ_k(x)]` — the spread of the member means. The **aleatoric** variance is `E_k[σ_k²(x)]` — the average within-member variance. Total = epistemic + aleatoric (Lakshminarayanan eq. 3, implemented in `ensemble.py:286–294`). Choice of K=5: diminishing returns past 5 in the original paper; ALDE (Yang 2025) also uses 5. Greenman 2025 audit: ensembles tend toward best accuracy but worst calibration — patch with temperature scaling if ECE > 0.1.

### Full mechanism
For each member k ∈ {0, …, 4}:
1. Set random seed = k (controls weight init AND DataLoader shuffle).
2. Initialize 3-layer MLP (256 / 256 / 128 hidden, ReLU, 0.1 dropout).
3. Two output heads: `head_mean → μ_k`, `head_log_sigma → log σ_k`, clamped log_sigma to [-7, 7].
4. Train with Gaussian NLL loss for max 200 epochs, early-stop on val NLL (patience 10).
5. Save state_dict.

Prediction (`DeepEnsemble.predict()`):
- Forward all 5 members → `means[5, N]`, `sigmas[5, N]`.
- `epistemic_var = means.var(axis=0)`
- `aleatoric_var = (sigmas**2).mean(axis=0)`
- `total_std = sqrt(epistemic_var + aleatoric_var)`
- `epistemic_std = sqrt(epistemic_var)` ← THIS is what BALD uses.

### Analogies
- **Five doctors, second opinions**: when they all agree, the diagnosis is settled; when they disagree, the disagreement itself is the signal that more tests are needed.
- **Five GPS systems on different satellites**: when they agree on your position, low uncertainty; when they diverge, you're in a region with poor coverage.
- **Bootstrap resampling**: classical statistics' way of getting a confidence interval from a single dataset.

### Memorize verbatim
> "Five MLPs, same architecture, different random seeds — they agree where data tells them what to do, and disagree where it doesn't. That disagreement is epistemic uncertainty."

### Where it could fail
- **Members collapse to same solution** — checked in `run_cycle0.py:194–199` via `frac_diverse > 0.5`. If they collapse, BALD scores all zero.
- **Calibration drift** — ensembles tend to be over-confident on out-of-distribution inputs (Greenman 2025).
- **Aleatoric vs epistemic confusion** — if you accidentally use `total_std` for BALD, you'll pick variants where the **measurement is noisy**, not where the **model is uncertain**.

---

## 5. BALD math — Std_k[μ_k] for regression (derivation)

### 30-second pitch
BALD = "Bayesian Active Learning by Disagreement." It picks the experiment that gives you the most information about the model's parameters. For a deep ensemble doing regression with Gaussian outputs, this reduces to: pick the input where the ensemble member means disagree the most. Mathematically: `BALD(x) ≈ Var_k[μ_k(x)] → score = epistemic_std`.

### 2-minute explanation
Houlsby et al. 2011 introduced BALD for Gaussian Process **Classification**. Their objective: `argmax_x I(y ; θ | x, D)`, where `I` is mutual information between the label and the model parameters. Equation 2 of the paper rewrites this as `H[y | x, D] − E_{θ ~ p(θ|D)}[H[y | x, θ]]`. **First term** = predictive entropy (total uncertainty). **Second term** = expected per-parameter-sample entropy (aleatoric only — fixing θ removes parameter uncertainty). Difference = epistemic uncertainty.

For our setup — regression with Gaussian deep ensemble — substitute Gaussian entropies. The per-member Gaussian has entropy `½ log(2π e σ_k²)`. The mixture-of-Gaussians (averaged over K members) has predictive variance `Var_k[μ_k] + E_k[σ_k²]`. Plug in:
- `H[y | x, D] ≈ ½ log(2π e (Var_k[μ_k] + E_k[σ_k²]))`
- `E_θ[H[y | x, θ]] = ½ log(2π e) + ½ E_k[log σ_k²]`
- Difference (ignoring `log` Jensen's inequality slack) ≈ `½ log(1 + Var_k[μ_k] / E_k[σ_k²])`, which is monotone in `Var_k[μ_k]` when aleatoric is roughly constant across `x`.

In our regime, `E_k[σ_k²]` is roughly stable across the pool (homoscedastic per-pair noise), so ranking by `Var_k[μ_k]` is rank-equivalent to ranking by the full BALD score. We take the **standard deviation** (`epistemic_std`) for the same reason — monotone, same units as the prediction.

**This extension from GPC classification to deep ensemble regression is OUR work** — we cite Houlsby for the original objective and treat the reduction as our contribution.

### Full derivation (with explicit Gaussian)
Per-member predictive: `y | x, θ_k ~ N(μ_k(x), σ_k²(x))`.
Mixture predictive: `y | x, D ~ (1/K) Σ_k N(μ_k(x), σ_k²(x))`.
- `Var[y | x, D] = E_k[σ_k²] + Var_k[μ_k]` (law of total variance).
- Predictive entropy of a mixture is **upper-bounded** by Gaussian-with-the-same-variance entropy: `H[y | x, D] ≤ ½ log(2π e (E_k[σ_k²] + Var_k[μ_k]))`.
- `E_θ[H[y | x, θ]] = E_k[½ log(2π e σ_k²)]`.
- `BALD(x) = H[y | x, D] − E_θ[H[y | x, θ]] ≈ ½ log( (E_k[σ_k²] + Var_k[μ_k]) / exp(E_k[log σ_k²]) )`.

When `σ_k(x) ≈ σ` (homoscedastic per-x):
- Numerator → `σ² + Var_k[μ_k]`.
- Denominator → `σ²`.
- BALD → `½ log(1 + Var_k[μ_k] / σ²)`, monotone in `Var_k[μ_k]`. So rank by `Var_k[μ_k]` = rank by `epistemic_std`. QED.

### Analogies
- **Crowd-sourcing**: ask the question where the crowd is most split — measuring the truth on that one shifts everyone's beliefs the most.
- **Investing**: rebalance toward the position whose value you're least certain about — the marginal price discovery matters most there.
- **Triangulation**: surveyors don't measure where two lines of sight already intersect; they measure where the intersection is fuzziest.

### Memorize verbatim
> "BALD for our deep-ensemble regression is `epistemic_std = Std_k[μ_k]` — the standard deviation of the ensemble member means. Houlsby 2011 derived it for GPC classification; the deep-ensemble regression form is our extension."

### Where it could fail
- **Heteroscedastic noise** — if `σ_k(x)` varies a lot, the homoscedasticity assumption breaks and ranking by `Var_k[μ_k]` is no longer rank-equivalent to the full BALD objective.
- **Posterior approximation** — deep ensembles are a crude approximation to a true Bayesian posterior. Calibration audits (Ovadia 2019, Greenman 2025) flag this.
- **Batch correlation** — when picking the top-4 BALD, the picks may be redundant if epistemic uncertainty correlates across similar inputs. Batch-BALD (Kirsch 2019) addresses this — we don't use it yet.

---

## 6. The 48-hour closed-loop cycle — why the SLA matters

### 30-second pitch
The wet lab cycle is 10–14 days end-to-end (clone → express → ELISA). The dry lab cycle is 48 hours (retrain model on new ELISA → BALD → next recommendations). If dry lab slips, the wet lab either idles (wasting time) or runs from the old recommendations (wasting information). The 48-h SLA is what makes the loop **closed** — not just sequential.

### 2-minute explanation
Naively, you'd say "the model can take however long it needs." But the wet lab works in fixed cycles — SDM primers ordered, plasmids transformed, induced, purified, ELISA, regression. That pipeline runs ~10–14 days. If new model recommendations arrive **before** the wet lab plans the next round, they go into Cycle N+1. If they arrive **after**, you have two bad options: idle the wet lab (lose a week of throughput) or use stale recommendations (lose the information from the just-finished cycle). 48 hours is the buffer that keeps wet lab moving without ever using stale data. Operationally, the SLA is enforced by: pre-computed `safe_pick_backup.csv` (the top-10 BALD picks from the previous cycle — used only if SLA misses), MLflow tracking of training time, and a documented escalation path. If we slip Cycle 1, we slip downstream; if we slip three cycles, we miss the iGEM wiki freeze.

### Full mechanism
Cycle N → Cycle N+1 timeline:
- **Day 0**: ELISA data finalized; CSV uploaded to `08_cycle_data/outputs/cycle_<N>/elisa_processed.csv`.
- **Day 0+0h**: Dry lab pulls data, regenerates ESM-2 embeddings if new variants, retrains ensemble (Module 06).
- **Day 0+24h**: Retraining done, calibration check (ECE ≤ 0.1 or apply temperature scaling). Predictions on full pool exported.
- **Day 0+30h**: BALD scoring + selection (`run_bald.py`).
- **Day 0+36h**: Cross-check with PI (sanity review of top picks).
- **Day 0+48h** (SLA): `recommendations.csv` + `primer_sequences.txt` delivered to wet lab.
- **Day 0+48h+**: Wet lab orders SDM primers (NEB Q5 design tool wired into `primer_sequences.txt`), begins Cycle N+1.

Quality gates: ECE > 0.1 → temp scaling; >50% variant SDM failure rate previous cycle → backup picks; BALD score variance shrinking faster than expected → flag for PI review (model may be overfitting).

### Analogies
- **Just-in-time manufacturing**: Toyota's kanban system — downstream pulls only when ready, no inventory of stale info.
- **Recipe testing on a TV show**: chef prepares two versions while waiting for feedback on the first; if feedback is late, the second is cooked blind.
- **CI/CD pipeline**: code commit → tests → deploy. If tests take longer than the next commit, you either queue or you skip.

### Memorize verbatim
> "The 48-hour SLA keeps the loop closed — wet lab never waits on us, and we never push stale recommendations. Safe-pick backup is the parachute, not the default."

### Where it could fail
- **Retraining doesn't converge** — old models loaded, fresh hyperparameters retried, escalation to PI before pushing.
- **ECE blow-up** — temperature scaling applied; if still bad, fall back to safe-pick.
- **ELISA data arrives malformed** — schema validation in pipeline rejects malformed CSVs upfront with a clear error.
- **Cumulative slippage** — three slipped cycles = wiki freeze risk. We treat slippage as a critical-path PI escalation.

---

## 7. Lee 2009 HMM-vs-BLAST — why rbp_01 is complementary, not contradictory

### 30-second pitch
Lee 2009 — the phiL7 genome paper — explicitly searched for a tail-fiber homolog of OP1 ORF25 and **said in print they couldn't find one** using sequence similarity. Our pipeline found rbp_01 using a Hidden Markov Model. HMMs catch what BLAST misses — proteins so diverged the letters no longer match but the structural pattern does. Lee 2009 used the right tools for 2009; we have a more sensitive tool. The findings are complementary, not contradictory.

### 2-minute explanation
BLAST (and the search tools Lee 2009 used) finds protein homologs by **sequence similarity** — alignment of amino acid letters, scored by substitution matrices like BLOSUM62. Sensitive when sequences are >25–30% identical; rapidly loses signal below. HMMs (Hidden Markov Models, profile HMM via tools like HMMER) instead build a probabilistic model of a **protein family** from a multiple sequence alignment — capturing which positions are conserved, which tolerate substitution, which insert/delete gaps. Sequence-divergent but structurally conserved proteins **stay detectable** by HMM long after BLAST loses them. PhageRBPdetect (Boeckaerts 2022) ships HMM profiles for tail-spike domains — including `Tail_spike_N`. rbp_01 hits this profile with score 342 (very high). It's plausibly the same functional family of tail spikes Lee was looking for, just sequence-diverged enough that 2009 BLAST didn't find it. **Saying "Lee 2009 missed it" is wrong** — they used the right tool for the time. Saying "rbp_01 contradicts Lee 2009" is also wrong — Lee searched for an OP1-ORF25 *homolog* and couldn't find one; we're not claiming rbp_01 is that homolog, just that it's a candidate tail spike found by a more sensitive method.

### Full mechanism
BLAST: scores `S = Σ_i s(a_i, b_i)` over aligned residues, then computes E-value vs random database. Threshold: E < 0.01 typically. Sensitivity drops sharply when sequence identity falls below ~30%.

Profile HMM (HMMER's `hmmsearch`): given a profile built from N homologs, each position has emission probabilities (which amino acids are allowed there) and transition probabilities (insertion/deletion). Score is log-likelihood against the profile, calibrated to bit score. Domain-level hits: even short conserved motifs can register if the profile captures them.

Crucially: profile HMMs use **information from N homologs**, not just the query against one target. Boeckaerts 2022 built `Tail_spike_N` from manually curated tail-spike domain alignments across hundreds of phages. Our rbp_01 hits at score 342 — orders of magnitude beyond the trust threshold.

Lee 2009 ran BLAST against OP1 ORF25 — a single template. rbp_01 vs OP1 ORF25: not a BLAST hit. rbp_01 vs Tail_spike_N HMM (built from many homologs): hit. These are not in conflict.

### Analogies
- **Police sketch artist vs face recognition**: sketch is one image (BLAST), face-recognition uses thousands of training faces (HMM) — recognizes someone in disguise that the sketch doesn't.
- **Tasting a wine you've never tasted**: knowing one specific Bordeaux (BLAST template) vs having tasted 1000 wines (HMM profile) — the latter catches family resemblance even with novel grape ratios.
- **Translation**: knowing one bilingual sentence vs knowing the grammar — the grammar generalizes to sentences you've never seen.

### Memorize verbatim
> "Lee 2009 used BLAST and explicitly said they couldn't find an OP1 ORF25 homolog. We used a Hidden Markov Model — a more sensitive tool — and found rbp_01. The two results are complementary, not contradictory."

### Where it could fail
- **HMM false positive**: rbp_01 might hit the profile incidentally without being a true tail spike. Mitigation: Boltz-2 monomer pTM 0.808 supports a legitimate fold; ELISA will measure direct binding to TonB.
- **rbp_01 not the binding RBP**: phiL7 might have multiple host-range proteins; rbp_01 could be a structural tail protein that doesn't actually bind TonB. Plaque assay on ΔtonB tests this.
- **Sequence divergence is real but irrelevant**: rbp_01 is the right family but doesn't bind Xanthomonas TonB specifically. The whole pipeline is designed to test this.
