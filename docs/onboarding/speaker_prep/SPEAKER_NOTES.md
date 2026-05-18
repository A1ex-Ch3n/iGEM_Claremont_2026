# SPEAKER_NOTES.md — Per-slide notes for Alex

**Companion to** `docs/onboarding/slides_en.md` (60+ slides, 8 Parts).
**Target total:** 60–90 min talk + 15 min Q&A.
**Convention:** Each Part below has a time budget, then per-slide notes. Vocabulary slides ("Key vocabulary — for this part") share a single short note block.

---

## Part 0 — Roadmap (slides: TL;DR, Roadmap, How to read, Suggested cuts)

### Time: **4 min**

### Slide: TL;DR
- **Say this:**
  - "One sentence — we built a closed-loop active-learning pipeline so each ELISA measurement makes the model smarter, not just bigger."
  - "Modules 00–07 are done and tested in code. Module 08 — the wet-lab side — opens around June 1."
  - "First real Boltz-2 structure prediction is in: rbp_01 (712 aa) docked against TonB. ipTM 0.365 — moderate. I'll explain what that number means later."
  - "BALD picks the next 4 variants + 1 random control in under a second. The first cycle has already run on synthetic data."
- **Don't say:** "We have a working model that predicts binding affinity." (We have predictions on **synthetic** data only; real binding labels arrive Cycle 0.)
- **If asked "Is this paper-quality?":** "The pipeline architecture is. The data hasn't arrived yet — that's exactly what Cycle 0 produces."
- **Transition:** "Here's how I'll walk you through it…"

### Slide: Roadmap
- **Say this:**
  - "Eight parts. Each has an audience tag — `[WET LAB]`, `[DRY LAB]`, `[PI]`, or all."
  - "Skim sections that aren't tagged for you, but you're welcome to interrupt anywhere."
- **Transition:** "Read the tags. Let's start with the science."

### Slide: How to read this deck / Suggested cuts
- **Say this:** brief — "If we run short, I'll skip the optional sections marked in `TIMING.md`. Demo at the end is 15 min, can be cut if Q&A runs long."
- **Transition:** "OK — into the science."

---

## Part 1 — The science (slides: vocab + 7 content slides) · `[WET LAB] [PI]`

### Time: **10 min** (this is the part that hooks the wet lab + PI)

### Slide: Key vocabulary — Part 1
- **Say this:**
  - "Quick vocabulary block. I'll define things again as I go, but these are the seven words you'll hear most."
  - Spend ~30 sec; don't read every line — point at RBP, BLAST vs HMM, TonB. Move on.

### Slide: *Xanthomonas* — what & why
- **Say this:**
  - "*Xanthomonas* is a genus of plant pathogens — over 400 host species. Locally relevant: black rot of cabbage, citrus canker, bacterial spot of tomato."
  - "The current control is **copper sprays**, and resistance is widespread in California (Aiello 2019). That's the agricultural pain point."
  - "Field trials show **phage biocontrol matches copper efficacy** (Iriarte 2018, Holtappels 2022). So phages aren't speculative — they work. The problem is that you need the **right** phage for the right strain."
- **Don't say:** "Phages can replace copper." (They're a complement, not a wholesale replacement; PI may push back on overstatement.)
- **If asked "How big is the market?":** "Brassica crops in California alone are ~$300M annually; black rot losses run 10–20%. I can pull the USDA numbers if you want the exact figure."
- **Transition:** "And the reason you need the right phage is the host-range problem…"

### Slide: Phage biocontrol — the host-range problem
- **Say this:**
  - "Phages are **host-specific** — that's both the feature and the bug. No off-target microbes, but you have to match strain to strain."
  - "Best published predictor — PhageHostLearn, Boeckaerts 2024 *Nat Commun* — gets AUC up to **0.82 within one genus**. Cross-genus it drops to **0.67–0.70** (Mutalik 2025)."
  - "Bottleneck: labeled (phage, host) interaction data is scarce."
- **Don't say:** "PhageHostLearn is 82% accurate." (AUC is not accuracy. And 0.82 is best-case at 100% identity threshold; it drops at lower thresholds.)
- **If asked "Why not just sequence more phages?":** "Sequence is cheap; **quantitative binding labels** are the bottleneck. That's what active learning addresses."
- **Transition:** "Our reference system to start dry lab is…"

### Slide: Our reference scaffold
- **Say this:**
  - "Dry lab uses public references — Xcc ATCC 33913 plus phiL7 from NCBI. These are well-characterized."
  - "Wet lab will **self-isolate** from California brassicas — bypasses the USDA APHIS PPQ-526 permit, which is a 4-month process per PI consultation."
- **If asked "Why self-isolate when references exist?":** "(a) permit timeline, (b) novel community-contribution to iGEM Registry, (c) genuinely local strains."
- **Transition:** "And the receptor on Xcc that phiL7 actually binds…"

### Slide: phiL7 receptor system — Hung 2003
- **Say this:**
  - "Critical paper: Hung 2003, *BBRC*. They did Tn5 mutagenesis on Xcc and asked what mutations block phiL7."
  - "Three genes essential: **TonB, ExbB, ExbD1**. These form an energy-coupled outer-membrane import complex."
  - "And — this matters — **ΔexbD2 still allows infection.** ExbD2 is in the same operon but it is **not required**. So if we build ΔexbD2 alongside the other three knockouts, we get a built-in negative control for free: ΔexbD2 should still let phiL7 in."
- **Don't say:** "TonB-ExbB-ExbD1-ExbD2 are all essential." (This is the old mis-citation we corrected in May. Do NOT regress to it.)
- **If asked "Why does the receptor matter?":** "Because Tier 3 validation — the gold-standard story for the paper and for iGEM Best Model — requires showing the binding gain we engineer is **receptor-specific**, not generic stickiness."
- **Transition:** "And the RBP on phiL7's side — this is where we made our first discovery."

### Slide: phiL7's RBP — Lee 2009 and our HMM rediscovery
- **Say this (this is the most important slide of Part 1):**
  - "Lee 2009 — the phiL7 genome paper — **explicitly searched** for a tail-fiber homolog of OP1's ORF25 and **said in print they couldn't find one**. I'll quote the exact line." [point at the quote]
  - "Our pipeline found a 712-aa protein — rbp_01 — that hits PhageRBPdetect's **Tail_spike_N HMM profile**."
  - "**This is not a contradiction with Lee 2009.** They used sequence-similarity search — BLAST-style. We used a Hidden Markov Model profile, which catches structural homology that BLAST misses. HMMs and BLAST have different sensitivity profiles. rbp_01 is **complementary** to Lee 2009, not a correction of it."
- **Don't say:** "Lee 2009 missed it" or "Lee 2009 was wrong." (They did the right experiment with the tools of the time. Phrase it as we have a more sensitive tool.)
- **If asked "How sure are you rbp_01 is a real RBP?":** "Two signals: HMM hit score 342 against Tail_spike_N, and Boltz-2 monomer pTM 0.808 — the monomer is well-predicted. The binding to TonB is what we're going to measure to find out."
- **Transition:** "And now you see why this matters — we have basically one labeled interaction…"

### Slide: The data-scarcity bottleneck
- **Say this:**
  - "One experimentally confirmed interaction for our system. 2,236 pairs from literature curation but none are quantitative binding affinities — they're presence/absence."
  - "Variant space for a 712-aa protein is 20^712 — astronomically larger than any wet-lab budget. We can't brute-force this."
  - "So we need a strategy that gets the most information per ELISA measurement."
- **Transition:** "And that strategy is active learning."

### Slide: Why active learning is the right framing
- **Say this:**
  - "Active learning is **not new** — Lindley wrote about Bayesian optimal experimental design in 1956."
  - "Recent precedent in protein engineering: Hie 2024 — antibody affinity maturation with ESM-1b/1v, about 20 variants per antibody. Yang 2025 (ALDE) — DNN ensemble + Thompson sampling on enzyme directed evolution, yield 12% → 93% in two rounds."
  - "**Nothing published** applies this to phage RBP × bacterial receptor. That's our methodological contribution."
- **Don't say:** "Hie 2024 used ESM-2" (ESM-1b/1v) or "ALDE used BALD" (Thompson sampling). PI will catch either error.
- **If asked "Why this combination if nobody's done it?":** "Each piece is validated in adjacent domains. Our contribution is the integration plus the wet-lab closed loop on a phage system."
- **Transition:** "And in iGEM terms…"

### Slide: What this maps to in iGEM
- **Say this:**
  - "Three medal tracks: Best Agriculture Project (phage biocontrol), Best Model (closed-loop AL + UQ), Best Composite Part (RBP-His6 library from in-house isolates)."
- **Transition:** "OK, that's the why. Now let me show you the how."

---

## Part 2 — Pipeline architecture (slides: vocab + 9 content slides) · `[DRY LAB] [PI]`

### Time: **12 min**

### Slide: Key vocabulary — Part 2
- **Say this:** "If you only catch one word: **MANIFEST.csv**. Every output folder has one — it's how we make gitignored large files reproducible."

### Slide: The three layers
- **Say this:**
  - "Layer 0 — **priors** from pre-trained models. ESM-2 sees 65M proteins, Boltz-2 gives interface confidence. Free knowledge."
  - "Layer 1 — the **AL loop**: ensemble → BALD → ELISA → retrain."
  - "Layer 2 — **causality**: knock out the receptor, show the binding gain disappears. This is what makes it paper-grade not just blog-grade."
- **If asked "Why three layers?":** "Binding is necessary but not sufficient for infection — Farquharson 2021 showed T4 binds 85% of *E. coli* strains but plaques on 11%. Layer 2 separates 'phage docks' from 'phage works'."
- **Transition:** "Here's the module flow…"

### Slide: Pipeline at a glance — Modules 00 → 08
- **Say this:**
  - "Eight modules. Blue = data, yellow = features, green = ML, pink = wet lab. Read left to right."
  - "Each module is a self-contained folder with `inputs/`, `processes/`, `outputs/`."

### Slide: The data-contract convention
- **Say this:**
  - "`inputs/` is read-only pointers to upstream outputs. `processes/` is the **only** place code lives. `outputs/` has canonical artifacts plus a `MANIFEST.csv`."
  - "This is how parallel agents — and Sarah and Olivia working in parallel — don't clobber each other."
- **If asked "Is this enforced?":** "Convention, not enforcement. Reviewing PRs we check that nothing writes to `inputs/`."
- **Transition:** "Code style…"

### Slide: Notebook-first workflow + bilingual comments
- **Say this:**
  - "Author as Jupyter notebooks first; freeze to `.py` once stable. Module 07 is the exception — it's production orchestration."
  - "Every comment is **bilingual** — English / 中文 on the same line for short comments. This is in CLAUDE.md."
- **If asked "Why bilingual?":** "Team mix. It's friction at first but it's a hard constraint — I'd rather slow notebooks down than have anyone on the team blocked by language."
- **Transition:** "Now per module…"

### Slide: Module 00 — Raw Data
- **Say this:** "777 phage + 34 bacteria genomes from NCBI. Binaries are gitignored — ~630 MB. `MANIFEST.csv` tracks SHA-256 + size so it's still reproducible."
- **If asked "Why 777 phages?":** "Curated from PhagesDB + NCBI — broad-enough sample to train a generalist model later, but every one is a complete genome."

### Slide: Module 01 — Ground Truth interaction matrix
- **Say this:**
  - "2,236 phage-host pairs: 315 positive + 1,921 negative + 1 ground-truth pair (phiL7 × Xcc)."
  - "All literature-curated — Hung 2003 for the receptor, plus species-level matches. **No quantitative affinities** — these are presence/absence labels."
  - "22 tests pass."
- **Don't say:** "We have 2,236 binding affinities." (Presence/absence, not affinity.)

### Slide: Module 02 — Annotation
- **Say this:**
  - "Phage genes: **PHANOTATE** (McNair 2019) — dynamic programming over overlapping ORFs."
  - "Bacterial genes: **pyrodigal** — Python binding of Prodigal (Hyatt 2010)."
  - "**Never swap them.** Prodigal assumes non-overlapping ORFs and loses 10–15% of phage genes."
  - "phiL7: 80 ORFs (Lee 2009 reported 59 with an older caller). Xcc: 4,344 ORFs."

### Slide: Module 03 — RBP identification
- **Say this:**
  - "PhageRBPdetect (Boeckaerts 2022) — HMM track + XGBoost track. We use the HMM track because the ML track was blocked."
  - "Three candidates for phiL7. **rbp_01 is 712 aa, HMM score 342** against Tail_spike_N — this is our Cycle 0 target."
- **If asked "Why HMM and not ML track?":** "ML track requires features we don't have wired up yet — HMM is the published recommendation when in doubt."

### Slide: Module 04 — Protein embedding (ESM-2)
- **Say this:**
  - "ESM-2 is a language model — masked-language-modelling on 65M proteins. The output **embedding** is a vector that captures structural and functional context."
  - "Three sizes: **8M (320-D)** local CPU — what's in the outputs now. **650M (1280-D)** A100/L40S — the production target. **3B (2560-D)** Laguna only — final benchmark."
  - "Optional layer: **PLM-interact** (Liu 2025) — ESM-2 fine-tuned on PPIs, +16–28% AUPR on cross-species PPI. Never tested on phage-bacteria; our project may be the first."
- **Don't say:** "ESM-2 predicts structure." (It predicts amino-acid distributions. Structure-like signal emerges; it isn't a direct output.)

### Slide: Module 05 — Structure prediction (Boltz-2 + AF3)
- **Say this (RED FLAG SLIDE):**
  - "Boltz-2 (Passaro 2025) predicts complex 3D structure and gives an ipTM score — interface confidence."
  - "**The Boltz-2 affinity head is small-molecule only.** Trained on PubChem, ChEMBL, BindingDB. For protein-protein pairs it outputs **NaN**."
  - "So we use **ipTM** as a **structural confidence proxy** — it tells us how confident the model is in the interface geometry. It is **not** a quantitative binding affinity. I will keep repeating this."
- **Don't say:** "Boltz-2 gives us zero-shot binding affinity for RBP × receptor." (Wrong. Will be caught.)
- **If asked "Why use it at all if no affinity?":** "Two reasons: (1) structural confidence is a real signal for variant prioritization, (2) the predicted complex coordinates feed PyMOL/ChimeraX and let us design site-directed mutations rationally."
- **Transition:** "OK, into the ML core."

---

## Part 3 — The ML core (slides: 2 vocab + 7 content slides) · `[DRY LAB] [PI]`

### Time: **15 min** (this is the part the PI and dry lab will judge you on)

### Slide: Key vocabulary — neural-net basics + UQ/AL
- **Say this:** "Two vocabulary blocks here — neural net basics, then uncertainty / active learning. I'll move fast; flag anything you want me to slow on."

### Slide: Module 06 — Deep ensemble for predictive uncertainty
- **Say this:**
  - "Lakshminarayanan 2017: train 5 MLPs independently — same architecture, **different random seeds**. Each predicts a Gaussian — mean and sigma."
  - "Why deep ensemble: better calibration than MC Dropout (Ovadia 2019), scales to ESM-2's 1280-D inputs unlike Gaussian Processes."
  - "**Caveat — Greenman 2025**, *PLoS Comp Biol*: 'no single best UQ method.' Ensembles tend to be most accurate but worst-calibrated. We pick ensembles for scalability + ALDE precedent + ECE can be patched with temperature scaling."
- **Don't say:** "Greenman 2025 says deep ensemble is best." (They say no method is universally best.)
- **If asked "Why 5 and not 10?":** "Compute budget. Lakshminarayanan showed diminishing returns past ~5. We can scale up post-Cycle-0 if calibration is poor."

### Slide: ensemble.py — per-member prediction (snippet)
- **Say this:**
  - "Each member is a 3-layer MLP — 256 / 256 / 128 — with two output heads: mean and log-sigma. Log-sigma is clamped to [-7, 7] for numerical stability — directly from Lakshminarayanan §3.1."
  - "Training loss is Gaussian NLL — penalises wrong-and-confident more than wrong-and-uncertain. That's what teaches the model to **know what it doesn't know**."
- **If asked "Why log-sigma instead of sigma directly?":** "Forces positivity through exp, avoids gradient pathologies near zero."

### Slide: Epistemic vs aleatoric — the key decomposition
- **Say this (KEY SLIDE):**
  - "Total predictive variance splits into two: **epistemic** — `Var_k[μ_k]`, the variance of the member means, which **shrinks as we collect more data** — and **aleatoric** — `E_k[σ_k²]`, the expected per-member variance, which is **irreducible measurement noise**."
  - "BALD targets **epistemic only** — because that's what experiments can reduce. There is no point doing an experiment to learn about the pipette."
  - "`predictions.csv` exports both — `std` (total) and `epistemic_std` (the BALD input)."
- **If asked "How do you know which is which?":** "It's a mathematical decomposition — Lakshminarayanan equation 3 — directly visible in `ensemble.py` line 286–294."

### Slide: Module 06 output — Cycle 0
- **Say this:**
  - "80 (RBP × receptor) pairs scored on synthetic data."
  - "Epistemic_std ranges 0.04 to 0.22. Compare rbp_01 × rec_02 — epistemic 0.190 — to rbp_01 × rec_03 — epistemic 0.049. Same RBP, very different information content."
  - "Model version is `aa99d51_cycle_0` — git SHA plus cycle number. Every prediction is traceable to a commit."
- **Don't say:** "These are real binding predictions." (Synthetic until Cycle 0 ELISA.)

### Slide: Calibration plot — Cycle 0
- **Say this:** "Calibration looks near-diagonal — but this is **synthetic data**. The real plot lands after Cycle 0 ELISA. If ECE > 0.1, temperature scaling kicks in before the next BALD run."

### Slide: Module 07 — BALD intuition
- **Say this:** "Pick where ensemble members disagree most. That's where a new measurement reduces uncertainty fastest. **One sentence: BALD picks the experiments the model is most confused about.**"

### Slide: BALD math (regression extension)
- **Say this (BIG SLIDE):**
  - "Houlsby 2011: BALD is mutual information between the prediction and the model parameters — `H[y|x,D] - E_θ[H[y|x,θ]]`. The acquisition picks the input that maximises this."
  - "Houlsby derived this for Gaussian Process **Classification**. For a Gaussian deep ensemble on a **regression** target, the same objective reduces to `Var_k[μ_k(x)]` — the variance of member means."
  - "We use the square root — `epistemic_std` — because it has the same units as the prediction and is rank-equivalent."
  - "**This regression extension is ours, not in Houlsby's paper.** I want to be honest about that — academic audiences will ask."
- **Don't say:** "Houlsby 2011 proved this for regression." (They proved it for GPC classification.)
- **If asked "Why not Thompson sampling like ALDE?":** "Two reasons. (1) BALD is **information-theoretic** — picks for max uncertainty reduction, which is what we want when each experiment costs ~$50 and 4 days. (2) Thompson is exploitation-biased; we're in early cycles, exploration matters."

### Slide: bald.py — `bald_score()` + `select_batch()`
- **Say this:**
  - "30 lines. Rank by BALD descending, take top-K, then sample n_random from **remaining unmeasured** for the control arm."
  - "**Control arm idea is from Hie 2024**: blind the wet lab to which one is the random control — preserves the retrospective AL-vs-random comparison at project end."

### Slide: Module 07 output — Cycle 1 recommendations
- **Say this:** "Top 4 BALD picks plus 1 random. Wet lab doesn't know which is which — they're shuffled. These are still synthetic; real Cycle 1 picks land after Cycle 0 ELISA."

### Slide: ALDE caveat — Yang 2025 ≠ BALD validation
- **Say this (CRUCIAL SLIDE — read it carefully):**
  - "Yang 2025 (ALDE) is the closest published work: active learning + DNN ensemble + protein engineering."
  - "**But** — ALDE uses **Thompson sampling**, not BALD. And **one-hot encoding**, not ESM-2."
  - "So ALDE validates the **general framework** — AL + UQ + DNN ensemble works in protein engineering. **It does not validate our specific method.** BALD-on-deep-ensemble-regression-with-ESM-2 still needs its own citation chain — Houlsby + Lakshminarayanan + our extension."
- **Don't say:** "ALDE validated BALD." (They validated Thompson sampling. Different acquisition function.)
- **Transition:** "OK — now the Boltz-2 result everybody's been waiting for."

---

## Part 4 — Current Boltz-2 result (slides: vocab + 4 content slides)

### Time: **8 min**

### Slide: Key vocabulary — structure prediction
- **Say this:** "Three numbers matter: **ipTM** (interface — how good the docking is), **chain pTM** (monomer alone), **PAE** (per-residue error matrix). I'll explain as we go."

### Slide: What we ran (job 59986)
- **Say this:**
  - "Chain A: rbp_01, 712 aa tail spike. Chain B: TonB, 604 aa. NVIDIA L40S on Laguna, runtime ~3 minutes."
  - "**History note**: an earlier job (59949) ran on the wrong protein — a 85-aa P25 ORF I mis-labeled as rbp_01. Job 59986 is the correct one."

### Slide: The three numbers that matter
- **Say this (CRITICAL):**
  - "**ipTM = 0.365** — that's **low**. Model is uncertain how rbp_01 and TonB dock."
  - "**chain_A_ptm = 0.808** — that's **high**. The rbp_01 monomer fold itself is well-predicted."
  - "**Confidence_score = 0.683** — moderate overall."
  - "**The low ipTM is not a failure. It defines the experiment.** That uncertainty is exactly what the ELISA + active learning loop is designed to resolve."
  - "High chain pTM means rbp_01's fold is reliable — so when we design site-directed mutations, the structural scaffold is trustworthy even if the interface isn't."
- **Don't say:** "ipTM 0.365 means the binding is weak." (ipTM is about confidence in the structure, NOT the affinity. Mixing them up will lose credibility.)
- **If asked "Could ipTM 0.365 mean they don't really bind?":** "Possible — but Hung 2003 already showed phiL7 enters Xcc through TonB. So the interaction exists; we just can't predict the geometry. That's a model limitation, not a biological negative."

### Slide: PAE heatmap — interface block
- **Say this:** "PAE is a 1316×1316 matrix — per-residue alignment error. The off-diagonal block — rows from TonB by columns from rbp_01 — is the interface region. **Light** means low confidence. You can see we're light in the interface, dark on the diagonals — model knows each monomer well, doesn't know the docking."

### Slide: Where the Boltz-2 outputs live
- **Say this:** "All paths are in `PI_briefing_2026-05-11.md`. `affinity.json` has the three numbers, `*.pdb` opens in PyMOL or ChimeraX. **Reminder: `predicted_dG = null` because affinity head is small-molecule only.**"

---

## Part 5 — The 48-hour cycle (slides: vocab + 7 content slides) · `[WET LAB] [DRY LAB]`

### Time: **12 min** (this is the part wet lab will judge you on)

### Slide: Key vocabulary — Part 5
- **Say this:** "These are the wet-lab terms. Sarah, Olivia, Weitao, Carol — flag anything that doesn't match how you'd say it."

### Slide: Cycle structure
- **Say this:**
  - "Two clocks: **dry lab** turns ELISA into recommendations in **48 hours**. Wet lab cycle — SDM → expression → ELISA — is **10–14 days**."
  - "The two clocks don't have to match — they hand off at well-defined points."

### Slide: Wet lab → dry lab handoff
- **Say this:** "Three files at the end of each cycle: `elisa_processed.csv` with EC50s and R², `plaque_results.csv` for WT and ΔReceptor strains, `qc_report.md` with the SDS-PAGE image. **Minimum to retrain: 3 valid EC50s with R² > 0.9.** Failures: mark `NaN` + `failed_reason` — model handles missing data."
- **If asked "What if a whole batch fails?":** "Then we run the safe-pick backup from the previous cycle, retrain on partial data, and the model just doesn't update much. The pipeline is robust to one bad cycle."

### Slide: Dry lab → wet lab handoff (48-h SLA)
- **Say this:**
  - "Five files. **`recommendations.csv`** is the task list — 4 BALD + 1 random, **shuffled so wet lab doesn't know which is which**. That blind preserves the AL-vs-random comparison."
  - "`safe_pick_backup.csv` — top-10 BALD picks — is used **only if** the 48-hour SLA misses."

### Slide: Cloning execution
- **Say this:**
  - "Cycle 0 uses **gene synthesis** — 4–6 variants from structure-based design, ~$150 each, 2-week lead time."
  - "Cycle 1+ uses **site-directed mutagenesis** (NEB Q5) — single point mutations on existing constructs, ~$50, 4 days. **3× cheaper and 3.5× faster.**"
  - "BALD on small data selects point mutations — SDM is the natural execution method."

### Slide: ELISA protocol — cell-based binding
- **Say this:** "Adapted from Boeckaerts 2024 and Latka 2021. Coat plate with heat-inactivated *Xanthomonas*, block, serial-dilute His6-RBP, HRP-anti-His6 detection, TMB substrate, OD~450~. **4PL fit → EC50** is the active-learning target variable."

### Slide: Receptor knockouts — pK18mobsacB
- **Say this:**
  - "**Markerless deletion** via suicide vector + sucrose counter-selection — clean genome edits, no antibiotic scar."
  - "Four targets from Hung 2003. **ΔtonB / ΔexbB / ΔexbD1** should block infection. **ΔexbD2 should retain infection** — that's our built-in negative control for the entire knockout system."

### Slide: Validation tiers
- **Say this:**
  - "Tier 1: ELISA only — 'we found variants that bind better.' Tier 2: + plaque on WT — 'binding leads to infection.' Tier 3: + receptor knockouts — '**receptor-specific** causality.' Tier 3 is paper-grade."
  - "**PI recommendation**: commit to Tier 3 if knockouts start May 17."

### Slide: Quality gates and failure modes
- **Say this:** brief — point at table, "Every failure mode has a documented action. None of these are panics — they're rehearsed."

---

## Part 6 — Reproducing & demoing (slides: vocab + 6 content slides) · `[DRY LAB] [WET LAB]`

### Time: **6 min** if no live demo; **15 min** including live demo

### Slide: Key vocabulary — Part 6
- **Say this:** "Hand-wave through unless someone asks. The terms matter when you actually clone the repo."

### Slide: Quick-start
- **Say this:** "Six commands and you have the env. **`conda activate igem2026` before anything else** — otherwise imports fail."

### Slide: Per-module entry points
- **Say this:** "Each module has a numbered notebook to open first. `01_…` is the canonical example notebook — it's the bilingual-comment template."

### Slide: Where outputs live
- **Say this:** "Paths are absolute from repo root. Module 07's BALD output is in `outputs/cycle_1/recommendations.csv` — that's the file wet lab will pull each cycle."

### Slide: Live demo plan (full runbook → `DEMO.md`)
- **Say this:** "I'll do a 5-minute live demo: Module 03 tests, Module 04 embedding shape check, Module 06 cycle 0 (~3 sec), Module 07 BALD (<1 sec). If something errors, `DEMO.md` has the recovery."
- **If demo fails:** Pull up `recommendations.csv` directly from disk — "Here's the file the demo would have produced." Don't apologize for more than one sentence; pivot to the file.

### Slide: Laguna HPC — when to push to GPU
- **Say this:** brief — "ESM-2 8M and BALD run on a laptop. ESM-2 650M, Boltz-2, AF3 are Laguna. Recipe is in LAGUNA.md."

### Slide: Tests — current pass rates
- **Say this:** "Every module has tests. Module 07 — 18/18. Module 06 — 9/9. Total ~140+ across the pipeline. The 'expected fails' are external dependencies — like hmmpress needing a local install."

---

## Part 7 — Risks, decisions, asks (slides: 5 content slides) · `[PI]`

### Time: **6 min** (target your PI directly here)

### Slide: Pending PI decisions
- **Say this:**
  - "Five decisions, in priority order. The big one is **pK18mobsacB vs CRISPRi for knockouts** — May 17. Default is pK18mobsacB because Hung 2003 used it successfully in Xcc."
- **If asked "Why not CRISPRi by default?":** "CRISPRi is knockdown, not knockout — leaves residual expression. For Tier 3 causality we want a clean null."

### Slide: Critical risks + mitigations
- **Say this:** "Every risk has a mitigation. The two I'd flag: **strain isolation failing** — dual-source plan, brassica + citrus. And **AL underperforming random** — that's actually a publishable negative result if we report it honestly."

### Slide: What "AL underperforms random" looks like
- **Say this:**
  - "After Cycle 2, BALD's test R² statistically equal to random's. We commit to honest reporting — document the negative with the same rigor as a positive. Hie 2024 worked; ALDE worked; phage-RBP didn't (yet) is itself a useful contribution."
  - "Per-cycle archive: test R², calibration ECE, information gain per experiment."

### Slide: Timeline
- **Say this:** "Key dates: Cycle 0 starts ~June 1. Cycle 1 BALD picks ~June 14. Knockouts done by June 28. Wiki freeze October 21. Jamboree November 13."

### Slide: iGEM deliverables checklist
- **Say this:** brief — "Five deliverables for iGEM: wiki, composite part, software repo, video, community contribution (in-house isolates → Registry)."

### Slide: Open questions for the team
- **Say this:** "Five open questions — flag any of these now if you have input, or after the talk."

---

## Part 8 — References & appendix (4 citation slides + glossary + log + thank-you)

### Time: **1 min** — don't read citations. "Slides are in the deck — full annotated list in `papers.md`."

### Slide: Thank you
- **Say this:** "Questions, push-back, comments — bring them on. Demo to follow if we have time."
- **If the room goes silent:** Ask a specific person — "Sarah, anything from the wet-lab handoff you'd push back on?" — to break the ice.

---

# 中文版 / SPEAKER_NOTES (ZH)

> 与上面英文逐 slide 对应。Alex 可视听众灵活切换中英文。

## Part 0 — Roadmap

### 时间: **4 分钟**

### Slide: TL;DR
- **要讲:**
  - "一句话:我们做了一个闭环主动学习管道,目的是让每一次 ELISA 测量都让模型变聪明,而不只是变大。"
  - "模块 00–07 已完成且通过测试,模块 08(湿实验数据)预计 6 月 1 日左右启动。"
  - "第一个真实 Boltz-2 结构预测出来了:rbp_01(712 aa)对接 TonB,ipTM 0.365——中等水平,这个数字稍后我会展开。"
  - "BALD 在一秒内挑出下一轮 4 个变体加 1 个随机对照。第一个 cycle 已经在合成数据上跑过了。"
- **不要说:** "我们已经有能预测结合亲和力的模型。"(目前预测都在合成数据上;真实结合标签要等 Cycle 0。)
- **如果被问 "这够发论文吗?":** "管道架构够了,数据还没到——那正是 Cycle 0 要产出的。"
- **过渡:** "下面我会带你们走一遍……"

### Slide: Roadmap / How to read / Suggested cuts
- **要讲:** "八个部分,每个有标签——`[WET LAB]`、`[DRY LAB]`、`[PI]`,或者全员。跳读没标的部分都可以,但欢迎随时打断。"

---

## Part 1 — 科学背景

### 时间: **10 分钟**(这部分是钓 PI 和湿实验的关键)

### Slide: 关键术语 — Part 1
- **要讲:** "快速过一下七个会反复出现的词。"30 秒;指着 RBP、BLAST vs HMM、TonB 几个关键点,不要逐行念。

### Slide: *Xanthomonas* — 是什么 / 为什么
- **要讲:**
  - "*Xanthomonas* 是植物病原菌属——超过 400 种宿主。本地相关:卷心菜黑腐病、柑橘溃疡病、番茄细菌性斑点。"
  - "当前控制手段是**铜剂喷洒**,加州抗药性已经普遍(Aiello 2019)。这是农业上的痛点。"
  - "田间试验显示**噬菌体生物防治可以匹敌铜剂效果**(Iriarte 2018, Holtappels 2022)。所以噬菌体不是假说——它能工作。问题是要找到针对特定菌株的合适噬菌体。"
- **不要说:** "噬菌体可以替代铜剂。"(它是补充,不是完全替代。)
- **如果被问 "市场多大?":** "光是加州十字花科蔬菜每年约 3 亿美元;黑腐病损失 10–20%。要精确数字我可以查 USDA。"

### Slide: 噬菌体生物防治 — 宿主范围问题
- **要讲:**
  - "噬菌体是**宿主特异**的——这既是优势也是缺点。没有脱靶微生物,但必须菌株对菌株匹配。"
  - "最好的已发表预测器——PhageHostLearn,Boeckaerts 2024 *Nat Commun*——属内 AUC 最高 **0.82**。跨属下降到 **0.67–0.70**(Mutalik 2025)。"
  - "瓶颈:已标注的 (噬菌体, 宿主) 相互作用数据稀缺。"
- **不要说:** "PhageHostLearn 准确率 82%。"(AUC 不是准确率,且 0.82 是最佳条件——100% identity threshold——下的值。)
- **如果被问 "为什么不多测一些噬菌体?":** "测序便宜;**定量结合标签**才是瓶颈,这正是主动学习的用武之地。"

### Slide: 参考骨架
- **要讲:** "干实验用公开参考——Xcc ATCC 33913 + phiL7。湿实验**自分离**——绕过 USDA APHIS PPQ-526 许可(4 个月,PI 已沟通)。"
- **如果被问 "为什么不直接用参考菌株?":** "三个原因:(a) 许可时间;(b) iGEM Registry 的原创社区贡献;(c) 本地真实菌株。"

### Slide: phiL7 受体系统 — Hung 2003
- **要讲(关键 slide):**
  - "关键论文:Hung 2003 *BBRC*。Xcc 的 Tn5 突变库,找出阻断 phiL7 的突变。"
  - "三个基因 essential:**TonB, ExbB, ExbD1**。这三个组成跨外膜的能量耦合输入通道。"
  - "**注意——ΔexbD2 仍然允许感染**。ExbD2 在同一操纵子里但**不是必需的**。所以如果我们在敲除 TonB/ExbB/ExbD1 的同时也做 ΔexbD2,就免费得到一个负对照:ΔexbD2 应该仍然让 phiL7 进入。"
- **不要说:** "TonB-ExbB-ExbD1-ExbD2 都是必需的。"(这是 5 月修正掉的旧错误,不要回退。)
- **如果被问 "受体重要在哪?":** "因为 Tier 3 验证——也就是论文和 iGEM Best Model 的金标准——必须证明工程化得到的结合增益是**受体特异**的,不是泛粘性。"

### Slide: phiL7 的 RBP — Lee 2009 与我们的 HMM 重发现
- **要讲(Part 1 最重要 slide):**
  - "Lee 2009——phiL7 基因组论文——**明确搜索**了 OP1 ORF25 的 tail fiber 同源物,并**在文中写明找不到**。"[指出原文引语]
  - "我们的管道用 PhageRBPdetect 的 **Tail_spike_N HMM** 发现了 712 aa 的 rbp_01。"
  - "**这与 Lee 2009 不矛盾。**他们用的是序列同源搜索——BLAST 类。我们用的是 HMM——能捕捉 BLAST 看不到的结构同源。两个方法的敏感度不同。rbp_01 是对 Lee 2009 的**补充**,不是更正。"
- **不要说:** "Lee 2009 漏掉了 rbp_01" 或 "Lee 2009 错了"。
- **如果被问 "你怎么确定 rbp_01 是真正的 RBP?":** "两个信号:HMM hit score 342 against Tail_spike_N,以及 Boltz-2 monomer pTM 0.808——单体折叠预测可靠。它和 TonB 的实际结合,正是我们要测的。"

### Slide: 数据稀缺瓶颈
- **要讲:** "我们这个系统只有 **一对** 实验确证的相互作用。文献整理出 2236 对噬菌体-宿主,但都是**有/无**标签,不是定量亲和力。712 aa 蛋白的变体空间是 20^712——远超任何湿实验预算。所以我们需要从每次 ELISA 测量中榨出最多信息的策略。"

### Slide: 为什么主动学习是对的框架
- **要讲:** "主动学习不是新东西——Lindley 1956 已经写过贝叶斯实验设计。蛋白质工程近期先例:Hie 2024 用 ESM-1b/1v 做抗体亲和成熟,每个抗体筛 20 个左右;Yang 2025 (ALDE) 用 DNN ensemble + Thompson sampling 做酶定向进化,两轮把 yield 从 12% 提到 93%。**没有任何已发表工作把它用到噬菌体 RBP × 细菌受体上——这是我们的方法学贡献。**"
- **不要说:** "Hie 2024 用了 ESM-2"(他们用 ESM-1b/1v)或 "ALDE 用了 BALD"(用 Thompson sampling)。

### Slide: 对应 iGEM 三个赛道
- **要讲:** "Best Agriculture Project(噬菌体生物防治)· Best Model(闭环 AL + UQ)· Best Composite Part(in-house 噬菌体的 RBP-His6 库)。"

---

## Part 2 — 管道架构

### 时间: **12 分钟**

### Slide: 关键术语 — Part 2
- **要讲:** "如果你只记住一个词——**MANIFEST.csv**。每个 output 文件夹都有,这是让 gitignore 的大文件可复现的方法。"

### Slide: 三层架构
- **要讲:**
  - "Layer 0 — **先验**:预训练模型给的免费知识。ESM-2 见过 65M 蛋白,Boltz-2 给界面置信度。"
  - "Layer 1 — **AL 闭环**:ensemble → BALD → ELISA → 重训练。"
  - "Layer 2 — **因果性**:敲除受体,展示结合增益消失。这一层把 blog 级别提升到 paper 级别。"
- **如果被问 "为什么需要三层?":** "结合是感染的必要条件但不充分——Farquharson 2021 显示 T4 噬菌体能结合 85% *E. coli* 菌株但只在 11% 上成噬菌斑。Layer 2 把'对接上'和'真的感染'分开。"

### Slide: Pipeline 概览 / 数据契约 / 笔记本优先
- **要讲:** "8 个模块。`inputs/` 只读、`processes/` 唯一放代码、`outputs/` 是产物 + MANIFEST.csv。代码用 Jupyter notebook 开发,稳定后冻结成 `.py`。所有注释**双语**——英文 / 中文同行。"

### Slides: Module 00 / 01 / 02 / 03 / 04 / 05
- 见英文版要点,逐 slide 简短讲解。
- **关键防线(Module 03):** "**rbp_01 是 712 aa,HMM 分 342,Tail_spike_N profile** ——Cycle 0 的目标蛋白。"
- **关键防线(Module 05):** "**Boltz-2 的 affinity head 是针对小分子的——蛋白-蛋白对输出 NaN**。我们用 **ipTM** 作为结构置信度 proxy,**不是定量亲和力**。这一点我会反复强调。"

---

## Part 3 — ML 核心

### 时间: **15 分钟**(PI 和干实验最看重的部分)

### Slide: Module 06 — Deep ensemble
- **要讲:**
  - "Lakshminarayanan 2017:5 个 MLP 独立训练——同样架构,**不同随机种子**。每个输出高斯——均值和 sigma。"
  - "为什么 deep ensemble:比 MC Dropout 校准更好(Ovadia 2019),能扩展到 ESM-2 的 1280 维输入(GP 不行)。"
  - "**Greenman 2025 警告**:'没有单一最好的 UQ 方法。'Ensemble 通常 accuracy 最高但 calibration 最差。我们选 ensemble 是因为可扩展性 + ALDE 先例 + ECE 可以用温度缩放补救。"
- **不要说:** "Greenman 2025 说 deep ensemble 最好。"

### Slide: ensemble.py 片段
- **要讲:** "每个成员是 3 层 MLP——256/256/128——两个输出头:均值和 log_sigma。log_sigma 截断到 [-7, 7] 防止数值爆炸。损失是 Gaussian NLL——'错且自信'比'错但不确定'惩罚更重。这就是让模型**知道自己不知道什么**的训练方式。"

### Slide: Epistemic vs aleatoric — 关键分解
- **要讲(关键 slide):**
  - "总预测方差分两部分:**epistemic** ——`Var_k[μ_k]`,各成员均值的方差,**随数据增加而缩小**——以及 **aleatoric** ——`E_k[σ_k²]`,各成员方差的期望,**不可减少的测量噪声**。"
  - "BALD 只针对 **epistemic**——因为实验能减少的只有这部分。没必要做实验去了解移液器。"
  - "`predictions.csv` 同时导出两者——`std`(总)和 `epistemic_std`(BALD 输入)。"

### Slide: Module 06 输出 / 校准图
- **要讲:** "80 个 (RBP × receptor) 对。epistemic_std 范围 0.04–0.22。同样 rbp_01:rec_02 epistemic 0.190 vs rec_03 epistemic 0.049——信息含量差很多。校准图看起来贴对角线但是**合成数据**,真的图等 Cycle 0 ELISA 之后。"

### Slide: BALD 直觉 / BALD 数学
- **要讲(关键 slide):**
  - "BALD = 模型最困惑的实验。"
  - "Houlsby 2011:BALD 是预测和模型参数之间的互信息——`H[y|x,D] - E_θ[H[y|x,θ]]`。"
  - "Houlsby 在论文中是**高斯过程分类**。我们把这个 objective 套用到 **deep ensemble 回归**上,化简为 `Var_k[μ_k(x)]`——各成员均值的方差。"
  - "我们用平方根 `epistemic_std`,单位和预测一致,排名等价。"
  - "**这个回归扩展是我们的,不在 Houlsby 原文里。**面对学术听众我会主动说清楚。"
- **如果被问 "为什么不用 Thompson sampling 像 ALDE?":** "两个原因。(1) BALD 是**信息论**的——挑最大不确定性减少,这正是每次实验 50 美元 + 4 天时我们想要的。(2) Thompson 偏 exploitation;我们在早期 cycle,exploration 很重要。"

### Slide: bald.py 片段
- **要讲:** "30 行。按 BALD 降序排,取 top-K,从**剩余未测**的里采样 n_random 做对照。**对照盲法来自 Hie 2024**:湿实验不知道哪个是随机对照,保留项目结束时 AL vs random 的回顾性对比。"

### Slide: Module 07 输出
- **要讲:** "top 4 BALD + 1 random,顺序打乱。仍然是合成的;真实 Cycle 1 picks 等 Cycle 0 ELISA 之后。"

### Slide: ALDE caveat
- **要讲(关键 slide):**
  - "Yang 2025 (ALDE) 是最接近的已发表工作:AL + DNN ensemble + 蛋白质工程。"
  - "**但是**——ALDE 用 **Thompson sampling**,不是 BALD。用 **one-hot encoding**,不是 ESM-2。"
  - "所以 ALDE 验证的是**通用框架**——AL + UQ + DNN ensemble 在蛋白质工程有效。**不是验证我们的具体方法。**BALD + deep ensemble regression + ESM-2 还需要自己的引用链——Houlsby + Lakshminarayanan + 我们的扩展。"

---

## Part 4 — 当前 Boltz-2 结果

### 时间: **8 分钟**

### Slide: 关键术语
- **要讲:** "三个数字重要:**ipTM**(界面——对接质量)、**chain pTM**(单体)、**PAE**(每残基误差矩阵)。"

### Slide: 我们跑的 (job 59986)
- **要讲:** "Chain A:rbp_01,712 aa tail spike。Chain B:TonB,604 aa。NVIDIA L40S on Laguna,~3 分钟。**历史注**:更早一个 job (59949) 跑错了蛋白——85 aa P25 ORF 被误标为 rbp_01。59986 是正确的。"

### Slide: 三个关键数字
- **要讲(关键 slide):**
  - "**ipTM = 0.365** ——**低**。模型对 rbp_01 和 TonB 怎么对接不确定。"
  - "**chain_A_ptm = 0.808** ——**高**。rbp_01 单体折叠预测可靠。"
  - "**confidence_score = 0.683** ——中等总体质量。"
  - "**低 ipTM 不是失败——它定义了实验。** 这个不确定性正是 ELISA + 主动学习闭环要回答的。chain pTM 高意味着 rbp_01 折叠是可信的,所以我们设计点突变时结构骨架可信,即便界面不可信。"
- **不要说:** "ipTM 0.365 意味着结合弱。"(ipTM 是结构置信度,不是亲和力。)
- **如果被问 "ipTM 0.365 会不会意味着它们其实不结合?":** "可能——但 Hung 2003 已经证明 phiL7 通过 TonB 进入 Xcc。所以相互作用是存在的;我们只是预测不出几何形状。这是模型局限,不是生物学否定。"

### Slide: PAE 热图 / 输出位置
- **要讲:** "PAE 是 1316×1316 矩阵。off-diagonal block(TonB 行 × rbp_01 列)是界面区域。**亮 = 低置信**。可以看到界面亮,对角暗——模型知道单体,不知道对接。路径在 PI_briefing 里。`predicted_dG = null`——再次说明 affinity head 是小分子的。"

---

## Part 5 — 48 小时闭环

### 时间: **12 分钟**(湿实验最看重的部分)

### Slide: 关键术语
- **要讲:** "这些是湿实验术语。Sarah / Olivia / Weitao / Carol,有不一致就指出来。"

### Slide: Cycle 结构 / 双向 handoff
- **要讲:** "两个时钟:**干实验** 48 小时把 ELISA 转成推荐;**湿实验** 10–14 天一个 cycle(SDM → 表达 → ELISA)。两个时钟不需要匹配,在明确的点上交接。"
- **湿实验交付**:`elisa_processed.csv` / `plaque_results.csv` / `qc_report.md`。**最小重训练数据:3 个有效 EC50,R² > 0.9**。失败 → `NaN` + `failed_reason`,模型能处理缺失。
- **干实验交付**:5 个文件。**`recommendations.csv` 是任务清单——4 BALD + 1 random,打乱顺序,湿实验不知道哪个是哪个**——保留 AL vs random 比较的盲法。

### Slide: 克隆 / ELISA / 敲除 / Tier
- **要讲:** "Cycle 0:**基因合成**——4-6 个变体,$150 / 2 周。Cycle 1+:**位点定向突变**(NEB Q5)——$50 / 4 天,3 倍便宜 3.5 倍快。ELISA 用 4PL 拟合得 EC50——这是 AL 目标变量。敲除用 pK18mobsacB(Hung 2003)——**ΔtonB / ΔexbB / ΔexbD1 应该阻断感染,ΔexbD2 应该保留感染**(免费负对照)。**PI 推荐:如果敲除 5 月 17 开始,承诺 Tier 3。**"

### Slide: 质量门与失败模式
- **要讲:** "每个失败模式都有文档化的应对。这不是恐慌——是预演。"

---

## Part 6 — 复现与演示

### 时间: **6 分钟** 不演示;**15 分钟** 含演示

- **要讲:** "六行命令把环境装好。**`conda activate igem2026` 之前任何 import 都会失败。** 每个模块有编号的入口 notebook。Module 07 的 `recommendations.csv` 是湿实验每个 cycle 要拉的文件。"
- **演示计划:** 5 分钟——Module 03 测试 / Module 04 embedding 形状检查 / Module 06 cycle 0(~3 秒)/ Module 07 BALD(<1 秒)。
- **演示失败:** 直接打开 `recommendations.csv` 说"这是演示本来要产出的文件",不超过一句道歉,转回到文件。

---

## Part 7 — 风险 / 决策 / 请求

### 时间: **6 分钟**(直接对准 PI)

- **要讲:** "五个待决定的事,按优先级。最大的:**pK18mobsacB vs CRISPRi 做敲除**——5 月 17。默认 pK18mobsacB,因为 Hung 2003 在 Xcc 里用成功过。"
- "两个我会标的风险:**菌株分离失败**——双源计划(brassica + citrus)。**AL 不比 random 好**——如果如实报告这本身就是一个可发表的负面结果。"
- "时间线:Cycle 0 ~6/1。Cycle 1 BALD picks ~6/14。敲除 6/28 之前。Wiki freeze 10/21。Jamboree 11/13。"

---

## Part 8 — 引用 / 附录

### 时间: **1 分钟** — 不要逐条念,"slide 在那里,完整带注释列表在 `papers.md`。"

### Thank-you
- **要讲:** "提问、反对、评论都欢迎。还有时间的话演示在后面。"
- **如果冷场:** 点名问一个具体人——"Sarah,湿实验交接有什么你想推翻的吗?"——破冰。
