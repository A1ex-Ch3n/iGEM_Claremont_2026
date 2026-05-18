---
title: "Active-Learning Phage Engineering for *Xanthomonas* Biocontrol"
subtitle: "iGEM Claremont 2026 — Team Onboarding"
author: "Core engineer: Alex Chen · PI: J. Cesar Ignacio-Espinoza · Faculty advisor: Ran Libeskind-Hadas"
date: "May 2026"
mainfont: "Hiragino Sans GB"
CJKmainfont: "Hiragino Sans GB"
monofont: "Menlo"
fontsize: 10pt
aspectratio: 169
classoption:
  - "aspectratio=169"
theme: "default"
colortheme: "seahorse"
innertheme: "rectangles"
outertheme: "infolines"
slide_level: 2
header-includes:
  - \usepackage{xeCJK}
  - \setCJKmainfont{Hiragino Sans GB}
  - \usepackage{fancyvrb}
  - \usepackage{booktabs}
  - \usepackage{xcolor}
  - \usepackage{mdframed}
  - \definecolor{cmcblue}{HTML}{1A4480}
  - \definecolor{wetlab}{HTML}{A01A4D}
  - \setbeamercolor{frametitle}{bg=cmcblue,fg=white}
  - \setbeamercolor{title}{fg=cmcblue}
  - \setbeamertemplate{navigation symbols}{}
  - \newenvironment{noteblock}{\begin{mdframed}[backgroundcolor=cmcblue!8,linecolor=cmcblue,linewidth=0.5pt,roundcorner=2pt,innertopmargin=4pt,innerbottommargin=4pt]}{\end{mdframed}}
  - \newenvironment{warnblock}{\begin{mdframed}[backgroundcolor=red!5,linecolor=red!60,linewidth=0.5pt,roundcorner=2pt,innertopmargin=4pt,innerbottommargin=4pt]}{\end{mdframed}}
---

# Part 0 — Roadmap

## TL;DR

- **Closed-loop active learning** for phage RBP engineering against *Xanthomonas campestris* pv. *campestris*.
- **Modules 00–07 done** (computational); Module 08 (wet-lab data) opens ~June 1.
- First Boltz-2 result on `rbp_01 (712 aa) × TonB`: **ipTM = 0.365, chain_A_ptm = 0.808**.
- BALD recommends the next 4 variants + 1 random control in **< 1 second**.

\vspace{0.5em}
\small *Targeting: Best Agriculture Project · Best Model · Best Composite Part.*

## Roadmap

1. **The science** — *Xanthomonas*, phage host-range, the data-scarcity problem.   `[WET LAB] [PI]`
2. **Pipeline architecture** — Modules 00–08, data-contract convention.            `[DRY LAB] [PI]`
3. **The ML core** — deep ensemble, BALD, code + sample outputs.                   `[DRY LAB] [PI]`
4. **Current Boltz-2 result** — what the 0.365 ipTM means.                          all
5. **The 48-hour cycle** — wet/dry handoffs, ELISA, knockouts.                     all
6. **Reproducing & demoing** — quick-start, tests, Laguna.                         `[DRY LAB] [WET LAB]`
7. **Risks, decisions, asks**                                                       `[PI]`
8. **References & appendix**                                                        all

## How to read this deck

Tags on each section header tell you who needs it:

- **`[WET LAB]`** — biological motivation, lab protocols, handoffs you'll execute.
- **`[DRY LAB]`** — code, algorithms, file conventions, reproducibility.
- **`[PI]`** — overall framing, decisions needed, risks.

Tagless sections are for everyone.

## Suggested cuts

- **30-min wet-lab onboarding:** slides 1–4, 7–12, 17–22 (modules at a glance),
  37–44 (the 48-hour cycle and protocols).
- **45-min PI / advisor briefing:** slides 1–4, 5–12 (science), 13–14, 23, 28–32 (ML core),
  33–36 (Boltz-2), 51–55 (risks + decisions), 56 (deliverables).
- **90-min full team deep dive:** all slides.

Demo (≈15 min) is independent; run `DEMO.md` after slide 48 or at the end.

# Part 1 — The science · `[WET LAB] [PI]`

## *Xanthomonas* — what & why we care

- Genus of plant-pathogenic gamma-proteobacteria; >400 host plant species (Ryan 2011, *Nat Rev Microbiol*).
- California-relevant pathovars:
  - **Xcc** — black rot of brassicas (cabbage, broccoli, kale).
  - *X. citri* subsp. *citri* — citrus canker.
  - *X. perforans* / *X. euvesicatoria* — bacterial spot of tomato.
- Current control: **copper bactericides**; resistance now widespread (Aiello 2019, *Plant Dis*).
- Field trials show phage biocontrol matches copper efficacy (Iriarte 2018; Holtappels 2022).

## Phage biocontrol — the host-range problem

- Phages are **host-specific** — that's the feature (no off-target microbes) AND the bug (each strain needs the right phage).
- Predicting which phage infects which bacterium is unreliable:
  - **PhageHostLearn** (Boeckaerts 2024, *Nat Commun*) reaches **AUC 0.82 within genus**.
  - **Cross-genus AUC drops to 0.67–0.70** (PAML benchmark, Mutalik 2025).
- The bottleneck: labeled (phage, host) interaction data is scarce, expensive, species-bound.

## Our reference scaffold

\begin{center}
\begin{tabular}{ll}
\toprule
Host & \textbf{Xcc ATCC 33913} \quad NCBI: GCF\_000007145.1 \\
Phage & \textbf{phiL7} \quad NCBI: EU717894 \\
Source & da Silva 2002 (\emph{Nature}); Lee 2009 (\emph{AEM}) \\
\bottomrule
\end{tabular}
\end{center}

\vspace{0.4em}

- Dry-lab uses these public references.
- **Wet lab self-isolates** from California brassicas — bypasses 4-month USDA APHIS PPQ-526 permit (per PI consultation 2026-05-07).

## phiL7 receptor system — Hung 2003

**Hung, C.-H. et al. (2003) *BBRC* 302:878–884, PMID 12646254.**

- phiL7 enters Xcc through the **TonB-ExbB-ExbD1** complex (energy-coupled outer-membrane import).
- Tn5 mutagenesis + complementation: **TonB, ExbB, ExbD1 essential**.
- ΔexbD2 (strain CH620) retains full sensitivity → **ExbD2 is NOT required**.

\begin{noteblock}
\textbf{Free negative control.} The Hung 2003 audit (May 2026) corrected an
earlier mis-citation that listed ExbD2 as essential. Building $\Delta$exbD2
alongside $\Delta$tonB / $\Delta$exbB / $\Delta$exbD1 validates the entire
knockout system: $\Delta$exbD2 should \emph{still} allow infection.
\end{noteblock}

## phiL7's RBP — Lee 2009 and our HMM rediscovery

**Lee, C.N. et al. (2009) *AEM* 75:7828** — the phiL7 genome paper.

- Lee 2009 **explicitly searched** for an OP1 tail-fiber (ORF25) homolog and **found none** by sequence similarity:

  > *"We were unable to identify a homolog of the OP1 tail fiber protein (ORF25)…"*

- p20 (1105 aa) was suggested as host-range related; no protein was named a "tail spike".
- **Our rbp_01 (712 aa)** comes from PhageRBPdetect's Tail_spike_N HMM — a structural-profile method that finds proteins too diverged for BLAST.
- **Not a contradiction with Lee 2009 — a complementary discovery using a more sensitive tool.**

## The data-scarcity bottleneck

- We have **one** experimentally confirmed (phage, host) interaction for our system: phiL7 × Xcc.
- The interaction matrix (Module 01) has 2,236 phage–host pairs from literature curation — none are quantitative binding affinities.
- For a 712 aa protein, the variant space is $20^{712}$ — astronomically larger than any wet-lab budget.
- **We need a learning strategy that gets the most out of each expensive ELISA measurement.**

## Why active learning is the right framing

- **Active learning (AL):** the model picks its next training point by maximizing expected information gain (Lindley 1956; Settles 2009).
- Demonstrated 2–5× speedup over random selection in adjacent domains:
  - **Hie 2024** (*Nat Biotechnol*) — human antibody affinity maturation with ESM-1b/1v + AL (~20 variants per antibody).
  - **Yang 2025 ALDE** (*Nat Commun*) — DNN ensemble + **Thompson sampling** (not BALD) on enzyme directed evolution, 12% → 93% yield in 2 rounds.
- **None applied to phage RBP × bacterial receptor.** That's our methodological contribution.

## What this maps to in iGEM

\begin{tabular}{ll}
\toprule
\textbf{Track} & \textbf{Our deliverable} \\
\midrule
Best Agriculture Project & Phage biocontrol pipeline against Xcc \\
Best Model               & Closed-loop AL + UQ + causal validation \\
Best Composite Part      & Registered RBP-His6 expression library \\
\bottomrule
\end{tabular}

\vspace{0.6em}

Plus: in-house *Xanthomonas* + phage isolates, sequenced + deposited.

# Part 2 — Pipeline architecture · `[DRY LAB] [PI]`

## The three layers

\begin{center}
\begin{tabular}{ll}
\toprule
\textbf{Layer 0 — Priors}      & Boltz-2 ipTM, ESM-2, PLM-interact (optional) \\
\textbf{Layer 1 — AL loop}     & Ensemble $\to$ BALD $\to$ ELISA $\to$ retrain \\
\textbf{Layer 2 — Causality}   & $\Delta$tonB / $\Delta$exbB / $\Delta$exbD1 knockouts \\
\bottomrule
\end{tabular}
\end{center}

\vspace{0.6em}

Binding (Layer 1) is necessary but not sufficient for infection (Farquharson 2021). Layer 2 decouples receptor-specific binding from defense-system contributions.

## Pipeline at a glance — Modules 00 → 08

\begin{center}
\includegraphics[width=\textwidth]{docs/onboarding/figures/pipeline_flow.png}
\end{center}

\small Color coding: data (blue) · features (yellow) · ML (green) · wet lab (pink).

## The data-contract convention

\begin{center}
\includegraphics[width=0.85\textwidth]{docs/onboarding/figures/data_contract.png}
\end{center}

\small
- `inputs/` = pointers to upstream outputs (read-only).
- `processes/` = the **only** place code lives.
- `outputs/` = canonical artifacts + `MANIFEST.csv` (large files gitignored).

Full spec: `INTERFACE.md` (paths, identifiers, FASTA headers, MANIFEST schema).

## Notebook-first workflow + bilingual comments (CLAUDE.md)

- New code is authored as **Jupyter notebooks first** (`<NN>_<short_name>.ipynb`).
- Once a notebook (a) runs end-to-end and (b) matches a verification cell → **freeze to `.py`**.
- **Bilingual comments** mandatory — English / 中文 on the same line for short comments.
- Module 07 is the exception — production orchestration, written as `.py` from day 1.
- `nbstripout` enforced repo-wide so notebook outputs never pollute git diffs.

## Module 00 — Raw Data

- **777 phage + 34 bacteria** complete genomes from NCBI.
- Genome binaries **gitignored** (~630 MB); re-fetchable via:
  - `python 00_raw_data/processes/fetch_phages.py`
  - `python 00_raw_data/processes/fetch_bacteria.py`
- `MANIFEST.csv` carries SHA-256 + size per genome for reproducibility.
- Canonical reference pair: `EU717894.1.fna` (phiL7) + `GCF_000007145.1.fna` (Xcc).

## Module 01 — Ground Truth interaction matrix

`01_data_ground_truth/outputs/interaction_matrix.csv` — first 4 rows:

\scriptsize
\begin{tabular}{lllll}
\toprule
phage\_acc & host\_acc & label & source & notes \\
\midrule
EU717894.1 & GCF\_000007145.1 & 1 & literature\_curated & Hung 2003 (receptor); Lee 2009 (genome) \\
NC\_054459.1 & NZ\_CP150073    & 1 & literature\_curated & X. oryzae pv. oryzae \\
ON758385.1 & — & 1 & literature\_curated & \emph{Xanthomonas} sp. \\
ON711490.1 & — & 1 & literature\_curated & Xcc XC114 \\
\bottomrule
\end{tabular}

\normalsize
- 2,236 pairs (315 positive + 1,921 negative + 1 ground truth).
- 22/22 tests pass.

## Module 02 — Annotation

- **Phage:** PHANOTATE (McNair 2019) — dynamic-programming over overlapping ORFs.
- **Bacteria:** pyrodigal binding of Prodigal (Hyatt 2010).
- **Never swap them.** Prodigal assumes non-overlapping ORFs → loses ~10–15 % of phage genes.
- Optional second pass: pharokka (PHROG / VFDB / CARD) for functional categories.

\begin{tabular}{lll}
\toprule
Genome & Tool & ORFs found \\
\midrule
phiL7 (EU717894.1) & PHANOTATE & 80 \\
Xcc (GCF\_000007145.1) & pyrodigal & 4,344 \\
\bottomrule
\end{tabular}

## Module 03 — RBP identification (PhageRBPdetect, Boeckaerts 2022)

`03_rbp_identification/outputs/EU717894.1_rbp_candidates.csv` — top 3 of 80:

\scriptsize
\begin{tabular}{llllll}
\toprule
orf\_id & length\_aa & hmm\_score & hmm\_match & combined\_score & rank \\
\midrule
EU717894.1\_orf\_00001 & \textbf{712} & 342.0 & unknown\_C54  & 1.000 & 1 \\
EU717894.1\_orf\_00021 & 918 & 235.1 & unknown\_C112 & 1.000 & 2 \\
EU717894.1\_orf\_00003 & 224 &  56.7 & unknown\_C294 & 1.000 & 3 \\
\bottomrule
\end{tabular}

\normalsize
- **rbp_01** (= orf_00001, 712 aa, HMM score 342) — primary Cycle 0 target.
- Hits the Tail_spike_N HMM — exactly what Lee 2009 couldn't find via BLAST.
- 25+ tests pass (2 expected fails until you run `hmmpress` locally).

## Module 04 — Protein embedding (ESM-2, Lin 2023)

- ESM-2 = masked-language-model trained on ~65 M protein sequences (UniRef50/90).
- Mean-pooling over residues → sequence-level vector.

\begin{tabular}{lll}
\toprule
Variant & Dim & Where it runs \\
\midrule
\texttt{esm2\_t6\_8M\_UR50D}    & 320  & Local CPU (proof of concept; current outputs) \\
\texttt{esm2\_t33\_650M\_UR50D} & 1280 & Laguna A100 / L40S (production target) \\
\texttt{esm2\_t36\_3B\_UR50D}   & 2560 & Laguna only (final benchmark) \\
\bottomrule
\end{tabular}

\vspace{0.5em}
\small Optional layer: PLM-interact (Liu 2025) — ESM-2 fine-tuned on human PPIs. AUPR +16–28 % when transferred to mouse / fly / worm / yeast / *E. coli* PPI. **Untested on phage-bacteria** — our project may be the first.

## Module 05 — Structure prediction (Boltz-2 + AF3)

- **Boltz-2** (Passaro 2025) — predicts complex 3D structure + ipTM interface confidence.
- **AF3** (Abramson 2024) — higher-quality static structures; trimer support.

\begin{warnblock}
\textbf{Critical caveat — Boltz-2 affinity head is small-molecule only.} It was trained on PubChem / ChEMBL / BindingDB. For \textbf{protein–protein pairs} (e.g., RBP × TonB) the affinity head outputs \texttt{NaN}. Use \textbf{ipTM} as a \emph{structural confidence proxy}, NOT a quantitative binding affinity.
\end{warnblock}

# Part 3 — The ML core · `[DRY LAB] [PI]`

## Module 06 — Deep ensemble for predictive uncertainty

- **Lakshminarayanan 2017** (*NeurIPS*) — 5 MLPs, independently trained, predictive mean + uncertainty.
- Each member predicts a Gaussian $(\mu_k, \sigma_k^2)$ via Gaussian NLL loss.

**Why deep ensembles?**
- Better calibration than MC Dropout (Ovadia 2019).
- Scales to ESM-2's 1280-D inputs (Gaussian Processes don't).

\begin{noteblock}
\textbf{Greenman 2025 audit:} \emph{"There is no single best UQ method"} across protein-engineering benchmarks. Ensembles tend to be most accurate but worst-calibrated. We pick ensembles because (a) scalability, (b) ALDE's precedent, (c) ECE can be patched with temperature scaling.
\end{noteblock}

## ensemble.py — per-member prediction (snippet)

\scriptsize
```python
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
\normalsize
\small Source: `06_uncertainty_model/processes/ensemble.py:57–68`.

## Epistemic vs aleatoric — the key decomposition

For a Gaussian ensemble of $K$ members, total predictive variance decomposes as:

$$\sigma_{\text{total}}^2(x) \;=\; \underbrace{\mathrm{Var}_k[\mu_k(x)]}_{\text{epistemic (reducible)}} \;+\; \underbrace{\mathbb{E}_k[\sigma_k^2(x)]}_{\text{aleatoric (noise)}}$$

- **Epistemic** = what *the model* doesn't know → shrinks with more data → **BALD's target**.
- **Aleatoric** = measurement / experimental noise → fundamentally irreducible.

\vspace{0.4em}
`predictions.csv` exports both as `std` (total) and `epistemic_std` (BALD input).

## Module 06 output — Cycle 0 (synthetic)

`06_uncertainty_model/outputs/cycle_0/predictions.csv` — head:

\scriptsize
\begin{tabular}{llrrr}
\toprule
rbp\_id & receptor\_id & predicted\_score & std & epistemic\_std \\
\midrule
rbp\_00 & rec\_00 & 4.872 & 0.736 & 0.125 \\
rbp\_00 & rec\_01 & 4.959 & 0.755 & 0.114 \\
rbp\_01 & rec\_02 & \textbf{5.128} & 0.726 & \textbf{0.190} \\
rbp\_01 & rec\_03 & 5.094 & 0.715 & 0.049 \\
rbp\_07 & rec\_02 & 5.177 & 0.721 & \textbf{0.218} \\
\bottomrule
\end{tabular}

\normalsize
- 80 (RBP × receptor) pairs scored.
- `epistemic_std` ranges 0.04–0.22 across the pool.
- `model_version = aa99d51_cycle_0` (git SHA + cycle).

## Calibration plot — Cycle 0 (synthetic)

\begin{center}
\includegraphics[width=0.62\textwidth]{06_uncertainty_model/outputs/cycle_0/calibration.png}
\end{center}

\small Reliability diagram: predicted credible-interval coverage vs observed. On synthetic data, calibration is near-diagonal — a meaningful version will be plotted after Cycle 0 ELISA arrives.

## Module 07 — BALD intuition

\begin{center}
\includegraphics[width=0.85\textwidth]{docs/onboarding/figures/bald_intuition.png}
\end{center}

\small Pick where ensemble members disagree the most — those measurements reduce model uncertainty fastest.

## BALD math (regression extension)

Original (Houlsby 2011, GPC classification):
$$\mathrm{BALD}(x) \;=\; I(y;\theta \mid x, D) \;=\; H[y \mid x, D] - \mathbb{E}_{\theta \sim p(\theta \mid D)}\!\bigl[H[y \mid x, \theta]\bigr]$$

For a Gaussian deep ensemble + regression target:
$$\mathrm{BALD}(x) \;\approx\; \mathrm{Var}_k[\mu_k(x)] \;\;\Longrightarrow\;\; \text{score} = \text{epistemic\_std}$$

\begin{noteblock}
\textbf{Extension note.} Houlsby 2011 originally applied BALD to a Gaussian Process \emph{Classifier}. The information-theoretic objective extends naturally to deep ensemble \emph{regression}, but this is our extension — not a direct citation. Mention this when presenting to academic audiences.
\end{noteblock}

## bald.py — `bald_score()` + `select_batch()`

\scriptsize
```python
def bald_score(epistemic_std: np.ndarray) -> np.ndarray:
    """BALD score = epistemic_std (rank-equivalent to epistemic variance)."""
    if np.any(epistemic_std < 0):
        raise ValueError("epistemic_std must be non-negative; got negatives.")
    return epistemic_std.copy()

# Inside select_batch():  rank by BALD score (descending), then take top-K
ranked = pool.sort_values(bald_col, ascending=False).reset_index(drop=True)
ranked["bald_rank"] = ranked.index + 1
bald_picks = ranked.iloc[:n_bald].copy()
bald_picks["selection_reason"] = [f"bald_top_{i+1}" for i in range(len(bald_picks))]

# Random control: sample from remaining (NOT in BALD top-K)
remaining = ranked[~ranked.index.isin(set(bald_picks.index))]
random_picks = ranked.loc[_random.Random(seed).sample(list(remaining.index), k=n_random)]
random_picks["selection_reason"] = "random_control"
```
\normalsize
\small Source: `07_acquisition_function/processes/bald.py:38–135`.

## Module 07 output — Cycle 1 recommendations

`07_acquisition_function/outputs/cycle_1/recommendations.csv`:

\scriptsize
\begin{tabular}{llrrrl}
\toprule
rbp\_id & receptor\_id & predicted\_score & epistemic\_std & bald\_rank & selection\_reason \\
\midrule
rbp\_07 & rec\_02 & 5.177 & \textbf{0.218} & 1 & bald\_top\_1 \\
rbp\_03 & rec\_01 & 4.810 & 0.197 & 2 & bald\_top\_2 \\
rbp\_05 & rec\_02 & 4.906 & 0.196 & 3 & bald\_top\_3 \\
rbp\_01 & rec\_02 & 5.128 & 0.190 & 4 & bald\_top\_4 \\
rbp\_03 & rec\_03 & 5.128 & 0.127 & 19 & random\_control \\
\bottomrule
\end{tabular}

\normalsize
- 4 BALD picks + 1 random control (Hie 2022 control-arm pattern).
- These are **synthetic placeholders**. Real Cycle 1 = rbp_01 variants × TonB after Cycle 0 ELISA.

## ALDE caveat — Yang 2025 ≠ BALD validation

- **Yang 2025 ALDE** (*Nat Commun*) is the closest published work — active learning + DNN ensemble + protein engineering.
- But: ALDE uses **Thompson sampling** as the acquisition function and **one-hot encoding** as features. Not BALD, not ESM-2.
- So ALDE validates: *"AL + UQ + DNN ensemble works in protein engineering"*.
- BALD specifically still needs its own citation chain (Houlsby 2011 + our extension).

# Part 4 — The current Boltz-2 result

## What we ran (job 59986)

\begin{tabular}{ll}
\toprule
Chain A & rbp\_01 (712 aa tail spike candidate, EU717894.1) \\
Chain B & TonB (604 aa, GCF\_000007145.1) \\
GPU     & NVIDIA L40S on CARC Laguna \\
Runtime & $\sim$3 minutes \\
\bottomrule
\end{tabular}

\vspace{0.6em}

History note: an earlier run (job 59949) used the wrong protein — a 85 aa P25 ORF mis-labeled as `rbp_01`. **Job 59986 uses the correct 712 aa rbp_01** from `03_rbp_identification/outputs/EU717894.1_rbp_sequences.faa`.

## The three numbers that matter

\begin{tabular}{lrl}
\toprule
\textbf{Metric} & \textbf{Value} & \textbf{Interpretation} \\
\midrule
\texttt{interface\_ipTM} & \textbf{0.365} & Low — model uncertain how they dock. \\
\texttt{chain\_A\_ptm}   & \textbf{0.808} & High — rbp\_01 monomer well-predicted. \\
\texttt{confidence\_score} & 0.683 & Moderate overall complex quality. \\
\bottomrule
\end{tabular}

\vspace{0.6em}

The low ipTM **is not a failure** — it defines the experiment. That uncertainty IS what the ELISA + active learning loop is designed to answer. High chain A ptm means rbp_01 is structurally well-constrained → reliable for variant design.

## PAE heatmap — interface block

\begin{center}
\includegraphics[width=0.55\textwidth]{pae_heatmap.png}
\end{center}

\scriptsize PAE matrix 1316×1316 (rbp_01 residues 0–711, TonB 712–1315). Low PAE (dark) = high relative confidence. The **off-diagonal block** (rows 712–1315 × cols 0–711) is the rbp_01 × TonB interface — light = low confidence.

Load via `np.load(.../pae_*.npz)['pae']` (full path on next slide).

## Where the Boltz-2 outputs live

\scriptsize
\texttt{05\_structure\_prediction/outputs/boltz2/\\
EU717894.1\_rbp\_01\_\_GCF\_000007145.1\_tonB/}

Inside that directory:

- `affinity.json` — interface_ipTM, chain_A_ptm, confidence_score.
- `boltz_results_.../predictions/.../*.pdb` — full atomic coordinates (open in PyMOL / ChimeraX).
- `pae_*.npz` — PAE matrix.
- `plddt_*.npz` — per-residue pLDDT (range 0.30–0.98, mean 0.76).

Full paths in `docs/planning/PI_briefing_2026-05-11.md` (no abbreviations).

\normalsize
**Boltz-2 reminder:** `predicted_dG = null` because the affinity head is small-molecule only.

# Part 5 — The 48-hour cycle · `[WET LAB] [DRY LAB]`

## Cycle structure

\begin{center}
\includegraphics[width=0.95\textwidth]{docs/onboarding/figures/cycle_48h.png}
\end{center}

\small Dry lab turns ELISA → recommendations in **48 hours**; wet lab cycle (SDM → expression → ELISA) is **10–14 days**.

## Wet lab → dry lab handoff

End of each cycle, wet lab delivers to `08_cycle_data/outputs/cycle_<N>/`:

\begin{tabular}{lll}
\toprule
File & Required columns & Notes \\
\midrule
\texttt{elisa\_processed.csv} & variant\_id, receptor\_id, ec50\_nM, hill\_slope, r2, plate\_id, date & \\
\texttt{plaque\_results.csv}  & variant\_id, strain\_id, pfu\_per\_ml, date & WT and $\Delta$Receptor \\
\texttt{qc\_report.md}        & SDS-PAGE image path, concentration, expression notes & one per cycle \\
\bottomrule
\end{tabular}

\vspace{0.4em}
**Minimum to retrain:** ≥3 valid EC50 per variant with R² > 0.9. Failed variants marked `ec50_nM = NaN` + `failed_reason` — the model handles missing data.

## Dry lab → wet lab handoff (48-hour SLA)

After receiving ELISA data:

\begin{tabular}{ll}
\toprule
File & Purpose \\
\midrule
\texttt{recommendations.csv}    & 4 BALD + 1 random; primary task list \\
\texttt{primer\_sequences.txt}  & NEB Q5-compatible SDM primers \\
\texttt{uncertainty\_bands.png} & Calibration: predicted vs measured (previous cycle) \\
\texttt{safe\_pick\_backup.csv} & Top-10 BALD, used only if 48-h SLA misses \\
\texttt{run\_meta.json}         & Git SHA + timestamp + pool stats (provenance) \\
\bottomrule
\end{tabular}

\vspace{0.4em}
Wet lab **does not know** which row is the random control — preserves blinded baseline for retrospective AL vs random comparison at project end.

## Cloning execution

- **Cycle 0** — gene synthesis (IDT/Twist), 4–6 variants from structure-based design. ~2 weeks lead time, ~$150/variant.
- **Cycle 1+** — site-directed mutagenesis (NEB Q5), point mutations on existing constructs. ~4 days, ~$50/variant (3× cheaper, 3.5× faster).
- BALD on a small data regime selects point mutations → SDM is the natural execution method.

Vector: pET-28a (Addgene 69864-3). Host: BL21(DE3). Induction: 0.5 mM IPTG, 18 °C overnight (favors soluble trimer assembly).

## ELISA protocol — cell-based binding

Adapted from Boeckaerts 2024 + Latka 2021:

1. Coat 96-well plate with heat-inactivated *Xanthomonas* (10⁸ CFU/well).
2. Block 3 % BSA, 1 h.
3. Serial dilution His6-RBP (1 nM – 1 µM).
4. HRP-anti-His6 detection, TMB substrate, read OD~450~.
5. **4PL fit → EC50** is the active-learning target variable.

Controls per plate: BSA-only (background); WT-RBP (inter-plate normalizer); heat-denatured RBP (folding-specific binding); 3 technical replicates.

## Receptor knockouts — pK18mobsacB (Schäfer 1994)

Markerless deletion via suicide vector + sucrose counter-selection:

1. Build deletion plasmid: ~500 bp upstream + ~500 bp downstream + pK18mobsacB.
2. Electroporate Xcc isolate; select on kanamycin (single crossover).
3. Counter-select on sucrose (sacB kills non-resolvers); PCR-confirm.

**Targets, from Hung 2003:**

- ΔtonB, ΔexbB, ΔexbD1 → expected to **block** phiL7 infection.
- ΔexbD2 → expected to **retain** infection (built-in negative control).

Timeline: 4–6 weeks per gene (parallelizable).

## Validation tiers

\begin{tabular}{lll}
\toprule
\textbf{Tier} & \textbf{What you measure} & \textbf{Story} \\
\midrule
1 & ELISA only (WT host) & "We found variants that bind better." \\
2 & + Plaque assay (WT) & "Binding $\to$ infection confirmed." \\
3 & + $\Delta$tonB / $\Delta$exbB / $\Delta$exbD1 & "Receptor-specific causality." (paper-grade) \\
\bottomrule
\end{tabular}

\vspace{0.6em}

**Recommendation (PI briefing):** commit to Tier 3 if knockouts start May 17. ΔexbD2 negative control validates the entire knockout system for free.

## Quality gates and failure modes

\scriptsize
\begin{tabular}{ll}
\toprule
\textbf{Failure} & \textbf{Action} \\
\midrule
Dry lab misses 48-h SLA   & Wet lab uses \texttt{safe\_pick\_backup.csv} \\
Variant insoluble         & Mark NaN; attempt backup truncation \\
ELISA R² $< 0.9$          & Down-weight in retraining \\
All 5 picks fail QC       & 1 expert pick + 2 random; retrain on partial data \\
Calibration ECE $> 0.1$   & Apply temperature scaling before next BALD \\
\bottomrule
\end{tabular}

# Part 6 — Reproducing & demoing · `[DRY LAB] [WET LAB]`

## Quick-start

\scriptsize
```bash
# Clone + branch / 克隆并切换分支
git clone https://github.com/A1ex-Ch3n/iGEM_Claremont_2026.git
cd iGEM_Claremont_2026
git checkout active-learning-pipeline

# Environment (one-time)
conda env create -f shared/env/environment.yml
conda activate igem2026
pip install nbstripout && nbstripout --install

# Minimum genomes for local dev (~5 MB total)
python 00_raw_data/processes/fetch_phages.py   --accession EU717894.1
python 00_raw_data/processes/fetch_phages.py   --accession NC_001604.1
python 00_raw_data/processes/fetch_bacteria.py --accession GCF_000007145.1

# Launch JupyterLab
jupyter lab
```

\normalsize
Full per-module checklists: `GETTING_STARTED.md`.

## Per-module entry points

\scriptsize
\begin{tabular}{ll}
\toprule
Module & Open this first \\
\midrule
00 & \texttt{00\_raw\_data/processes/01\_verify\_dataset.ipynb} \\
01 & \texttt{01\_data\_ground\_truth/processes/01\_fetch\_reference\_genomes.ipynb} \\
02 & \texttt{02\_annotation/processes/01\_run\_phanotate.ipynb} \\
03 & \texttt{03\_rbp\_identification/processes/01\_run\_phagerbpdetect.ipynb} \\
04 & \texttt{04\_protein\_embedding/processes/01\_embed\_esm2.ipynb} \\
05 & \texttt{05\_structure\_prediction/processes/01\_run\_boltz2.ipynb} (Laguna for production) \\
06 & \texttt{06\_uncertainty\_model/processes/run\_cycle0.py} \\
07 & \texttt{07\_acquisition\_function/processes/run\_bald.py} \\
\bottomrule
\end{tabular}

## Where outputs live

\scriptsize
\begin{tabular}{ll}
\toprule
File & Path (from repo root) \\
\midrule
Interaction matrix    & \texttt{01\_data\_ground\_truth/outputs/interaction\_matrix.csv} \\
RBP candidates        & \texttt{03\_rbp\_identification/outputs/EU717894.1\_rbp\_candidates.csv} \\
ESM-2 embeddings      & \texttt{04\_protein\_embedding/outputs/embeddings\_esm2\_t6\_8M\_phiL7\_rbps.npz} \\
Boltz-2 PDB / ipTM    & \texttt{05\_structure\_prediction/outputs/boltz2/EU717894.1\_rbp\_01\_\_GCF\_000007145.1\_tonB/} \\
Ensemble predictions  & \texttt{06\_uncertainty\_model/outputs/cycle\_0/predictions.csv} \\
Calibration plot      & \texttt{06\_uncertainty\_model/outputs/cycle\_0/calibration.png} \\
BALD recommendations  & \texttt{07\_acquisition\_function/outputs/cycle\_1/recommendations.csv} \\
PAE heatmap           & \texttt{pae\_heatmap.png} (repo root) \\
\bottomrule
\end{tabular}

## Live demo plan (full runbook → `DEMO.md`)

1. **Module 03** — `pytest 03_rbp_identification/processes/tests/ -v` (~15 s).
2. **Module 04** — inspect cached `.npz` (shape `(3, 320)`; seq_ids match rbp_01/02/03) (<1 s).
3. **Module 06** — `python 06_uncertainty_model/processes/run_cycle0.py` (~3 s, verified).
4. **Module 07** — `python 07_acquisition_function/processes/run_bald.py --cycle 1 --n_bald 4 --n_random 1` (<1 s).
5. **Module 05 (read-only)** — open `pae_heatmap.png` and the PDB in PyMOL.

\vspace{0.4em}
\small Each step has "what to say" + "if it fails" in `DEMO.md`.

## Laguna HPC — when to push to GPU

\scriptsize
\begin{tabular}{lll}
\toprule
Task & Local? & Why Laguna? \\
\midrule
ESM-2 8M embedding (3 RBPs)    & Yes & seconds on CPU \\
ESM-2 650M / 3B (777 phages)   & No  & GPU memory + time \\
Boltz-2 protein-protein        & No  & ~15 min/pair, A100/L40S only \\
AF3 batch                      & No  & GPU + weights approval \\
Ensemble train + BALD          & Yes & seconds on CPU \\
\bottomrule
\end{tabular}

\normalsize
Setup recipe + SLURM templates: `LAGUNA.md` (CUDA 12.4, torch 2.5.1+cu121).

## Tests — current pass rates

\begin{tabular}{ll}
\toprule
Module & Tests \\
\midrule
00 Raw Data            & 15+ pass (3 expected fails — GCF/T7 not in seed list) \\
01 Ground Truth        & 22/22 \\
02 Annotation          & 26/26 \\
03 RBP ID              & 25+ pass (2 expected — HMM needs local \texttt{hmmpress}) \\
04 Embedding           & 17/17 \\
05 Structure           & 28/28 (1 expected skip — needs completed GPU run) \\
06 Uncertainty Model   & 9/9 \\
07 BALD                & 18/18 \\
\bottomrule
\end{tabular}

# Part 7 — Risks, decisions, asks · `[PI]`

## Pending PI decisions

\scriptsize
\begin{tabular}{lll}
\toprule
\textbf{Decision} & \textbf{Deadline} & \textbf{Default} \\
\midrule
pK18mobsacB vs CRISPRi for knockouts  & May 17 & \textbf{pK18mobsacB} (Hung 2003 used it in Xcc) \\
AlphaFold 3 weights — who applies?     & This week & Either Alex or PI; institutional email preferred \\
Validation tier commitment             & Before June 1 & \textbf{Tier 3} if knockouts start May 17 \\
Phage enrichment source (LA County)    & Before June 1 & TBD — needs PI input on lab access \\
Manuscript ambition                    & Before Cycle 0 & Concurrent submission \emph{Bioinformatics} / \emph{NAR Genomics} \\
\bottomrule
\end{tabular}

## Critical risks + mitigations

\scriptsize
\begin{tabular}{lll}
\toprule
\textbf{Risk} & \textbf{Likelihood / Impact} & \textbf{Mitigation} \\
\midrule
Strain isolation fails                & L / H & Dual source (brassica + citrus); 2 sampling rounds \\
Phage isolation yields nothing lytic  & M / H & Sewage / runoff enrichment; Phage Directory backup \\
RBP expression insoluble              & M / M & T7 gp17 positive control first; GCN4 trimer tag backup \\
ELISA dynamic range insufficient      & M / H & 2 weeks of Cycle 0 reserved for ELISA optimization \\
AL underperforms random               & L / H (reportable) & Standard control arm; honest reporting \\
Dry lab misses 48-h SLA               & M / M & \texttt{safe\_pick\_backup.csv} pre-approved by PI \\
Receptor knockout fails               & M / M & Multiple targets in parallel; CRISPRi fallback \\
\bottomrule
\end{tabular}

## What "AL underperforms random" looks like

- After Cycle 2: BALD trajectory's test R² **statistically equivalent** to random.
- We commit to honest reporting:
  - Document the negative result with the same rigor as a positive.
  - Useful contribution to the field — *Hie 2024 worked; ALDE worked; phage-RBP didn't (yet)*.
- Concrete metrics archived per cycle:
  - Test-set R² on held-out variants.
  - Calibration ECE.
  - Information gain per experiment (KL divergence posterior vs prior).

## Timeline

\scriptsize
\begin{tabular}{ll}
\toprule
Date & Milestone \\
\midrule
2026-05-17 & Wet lab launches: brassica sampling, pK18mobsacB construction, Cycle 0 synthesis order \\
2026-06-01 & Cycle 0 starts: ELISA optimization + first binding measurements \\
2026-06-14 & Cycle 1: ensemble retrained, first BALD picks delivered, SDM \\
2026-06-28 & Cycle 2: round-2 recommendations, receptor knockouts complete \\
2026-07-12 & Analysis + extension: trained model on soil-derived phage genomes \\
2026-09-12 & Wiki freeze prep \\
2026-10-21 & Wiki freeze \\
2026-11-13 & Grand Jamboree \\
\bottomrule
\end{tabular}

## iGEM deliverables checklist

- **Wiki:** Engineering DBTL · Modeling · Hardware/Software · Human Practices.
- **Composite Part:** ≥4–6 RBP-His6 constructs from in-house phage isolates.
- **Software repo:** open-source closed-loop AL pipeline (this repo).
- **Promotion video:** demo of one closed-loop cycle.
- **Community contribution:** in-house *Xanthomonas* + phage isolates → iGEM Registry; sequenced genomes → NCBI.

## Open questions for the team

1. Wet lab — what's the LA County sewage / agricultural runoff access?
2. Wet lab — do we want a positive expression control (T7 gp17) before committing the variant library?
3. Dry lab — who owns the ESM-2 3B Laguna job (1280→2560-D, final benchmark)?
4. Everyone — manuscript submission timing relative to iGEM wiki freeze?
5. PI — APHIS permit posture if our isolate ID reveals a quarantine-status pathovar.

# Part 8 — References & appendix

## Key papers per module — compact

\scriptsize
\begin{tabular}{lll}
\toprule
Module & Paper & Role \\
\midrule
00 / 01 & da Silva 2002 \emph{Nature}; Lee 2009 \emph{AEM}; Hung 2003 \emph{BBRC} & reference genomes + receptor \\
02 & McNair 2019 \emph{Bioinformatics} (PHANOTATE); Hyatt 2010 (Prodigal); Bouras 2023 (pharokka) & gene calling \\
03 & Boeckaerts 2022 \emph{Viruses} (PhageRBPdetect) & RBP HMM + XGBoost \\
04 & Lin 2023 \emph{Science} (ESM-2); Liu 2025 \emph{Nat Commun} (PLM-interact) & embeddings + PPI prior \\
05 & Abramson 2024 \emph{Nature} (AF3); Passaro 2025 (Boltz-2) & structure + ipTM \\
06 & Lakshminarayanan 2017 \emph{NeurIPS}; Greenman 2025 \emph{PLoS Comput Biol} & deep ensembles + UQ audit \\
07 & Houlsby 2011 (BALD); Yang 2025 \emph{Nat Commun} (ALDE); Hie 2024 \emph{Nat Biotechnol} & acquisition \\
08 & Schäfer 1994 (pK18mobsacB); Gibson 2009; Latka 2021 \emph{mBio} & wet lab methods \\
\bottomrule
\end{tabular}

\vspace{0.4em}
\small Full annotated reading guide: `docs/reference/papers.md` (19 papers).

## Five literature corrections (May 2026 audit)

\scriptsize
\begin{tabular}{ll}
\toprule
\textbf{Was} & \textbf{Corrected to} \\
\midrule
"TonB-ExbB-ExbD1-ExbD2 all essential" (Wang 2003)            & TonB-ExbB-ExbD1 essential; ExbD2 NOT (Hung 2003 BBRC) \\
"Boltz-2 affinity head for RBP × receptor"                    & affinity head is small-molecule only $\to$ \texttt{NaN}; use ipTM \\
"Greenman 2025 \emph{NAR Genomics}; deep ensemble = best UQ"  & Greenman 2025 \emph{PLoS Comput Biol}; "no single best UQ" \\
"Hie 2024 used ESM-2, $\sim$50 antibody variants"             & Hie 2024 used ESM-1b/1v, $\sim$20 variants per antibody \\
"Lee 2009 named the tail spike"                               & Lee 2009 \emph{searched and found no OP1 ORF25 homolog} \\
\bottomrule
\end{tabular}

\vspace{0.5em}
\small Full audit: `docs/reference/paper_reading_notes.md`.

## Navigation reference

\small
\begin{tabular}{ll}
\toprule
Want to know… & Read this \\
\midrule
What each module produces  & \texttt{<module>/README.md} \\
End-to-end status snapshot & \texttt{docs/planning/PI\_briefing\_2026-05-11.md} \\
Project plan (EN, full)    & \texttt{docs/planning/iGEM\_2026\_Project\_Plan.md} \\
Project plan (ZH, full)    & \texttt{docs/planning/iGEM\_2026\_项目大纲\_中文版.md} \\
Paper list (annotated)     & \texttt{docs/reference/papers.md} \\
Literature corrections     & \texttt{docs/reference/paper\_reading\_notes.md} \\
Data contracts             & \texttt{INTERFACE.md} \\
HPC setup                  & \texttt{LAGUNA.md} \\
Wet-lab SOPs               & \texttt{docs/protocols/} (5 Benchling PDFs) \\
\bottomrule
\end{tabular}

## Glossary — fast lookup

\scriptsize
\begin{tabular}{ll}
\toprule
Term & Meaning \\
\midrule
RBP        & Receptor-binding protein — the phage's "key" to the host cell \\
HMM        & Hidden Markov Model — sequence-profile method (Boeckaerts 2022) \\
ESM-2      & Evolutionary-scale protein language model (Lin 2023) \\
ipTM       & Interface predicted TM-score (0–1, structural confidence) \\
ELISA      & Enzyme-linked immunosorbent assay (binding readout) \\
EC50       & Half-maximal effective concentration (from 4PL fit) \\
SDM        & Site-directed mutagenesis (NEB Q5 kit) \\
BALD       & Bayesian Active Learning by Disagreement (Houlsby 2011) \\
Epistemic  & Reducible model uncertainty (= BALD target) \\
Aleatoric  & Irreducible measurement noise \\
\bottomrule
\end{tabular}

\vspace{0.4em}
\small Full glossary: `docs/reference/glossary.md`.

## Thank you / 谢谢

\vspace{0.6em}
\Large
**Questions, comments, push-back welcome.**

\normalsize
\vspace{0.5em}

- Dry lab: Alex Chen — `CChen29@cmc.edu`
- PI: Prof. J. Cesar Ignacio-Espinoza
- Faculty advisor: Prof. Ran Libeskind-Hadas
- Wet lab leads: Sarah, Olivia, Weitao, Carol
- Contributors: Ryan, Leah

\vspace{0.6em}
\small Live demo (15 min) follows. See `docs/onboarding/DEMO.md`.
