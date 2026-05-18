# QA_PREP.md — Anticipated questions with model answers + fallbacks

Three audience sections. For each Q: a **confident 3–5 sentence answer** Alex can deliver; then a **fallback** if Alex doesn't actually know.

---

## [WET LAB Q&A] — Sarah, Olivia, Weitao, Carol

### W1. "What do I actually do when `recommendations.csv` arrives?"
**Confident answer:** Pull the CSV from `07_acquisition_function/outputs/cycle_<N>/`. It has 5 rows — 4 BALD picks + 1 random control, **shuffled so you don't know which is which**. Each row tells you `rbp_id`, `receptor_id`, and the predicted score. Pair with `primer_sequences.txt` (NEB Q5-compatible) and execute the SDM workflow as documented in the Benchling protocol. Hand back `elisa_processed.csv` when ELISA is done. If you finish early, do the extra picks from `safe_pick_backup.csv` — they're vetted alternates.
**Fallback:** "Good question — let me pull up the protocol for the exact handoff format" → open `docs/protocols/` or `08_cycle_data/README.md`.

### W2. "What if a variant doesn't express in BL21?"
**Confident answer:** Mark it `ec50_nM = NaN` and add `failed_reason = "insoluble"` (or "low_expression", whichever fits). The model handles missing data — it doesn't break the cycle. We have a planned fallback: GCN4 trimer-tag fusion if more than ~30% of variants are insoluble, since tail spikes naturally want to trimerize. PI briefing flagged this as Medium-likelihood/Medium-impact risk with a clear path.
**Fallback:** "We have a backup plan in PI briefing 2026-05-11 risk table — let me check the exact escalation."

### W3. "Why do we need a random control if BALD is the smart pick?"
**Confident answer:** Three reasons. **First**, it's the only way to honestly answer "did BALD beat random?" at project end — Hie 2024 established this as the standard practice. **Second**, it's a calibration sanity check — if random consistently outperforms BALD, the model is broken. **Third**, you stay blinded to which row is which, which preserves the comparison's validity. It's not wasted effort — it's the whole reason we can claim AL works (if it does).
**Fallback:** "It's the standard control arm pattern — let me show you Hie 2024's setup."

### W4. "What if ELISA fails for the whole batch?"
**Confident answer:** That triggers Quality Gate row 5 — model retrains on partial data (just the previous cycle's accumulated set), wet lab uses `safe_pick_backup.csv`, and we treat the failed cycle as a "lost week" in the timeline. We've reserved two weeks of Cycle 0 specifically for ELISA optimization, so a single batch failure isn't catastrophic. The escalation path is in `docs/onboarding/guide_en.md` §quality-gates.
**Fallback:** "Let me look that up — we documented an escalation per failure mode."

### W5. "Cycle 0 takes 2 weeks for gene synthesis but Cycle 1+ only 4 days for SDM. Why?"
**Confident answer:** Cycle 0 starts from scratch — full-length DNA from IDT or Twist. That's $150 and 2 weeks. Cycle 1+ takes the Cycle 0 plasmids and makes single point mutations via NEB Q5 site-directed mutagenesis — $50 and 4 days. **3× cheaper, 3.5× faster.** BALD on small data tends to pick point mutations, which matches SDM's economics. If BALD ever recommends a huge multi-site rewrite, we go back to synthesis for that variant.
**Fallback:** "Cost and speed — the numbers are in the cloning slide. Let me pull them up."

### W6. "Why heat-inactivated *Xanthomonas* for ELISA instead of live cells?"
**Confident answer:** Two reasons. **Safety/handling** — once we move beyond BSL-1 strains, live cells in ELISA wells add containment overhead. **Reproducibility** — heat-inactivation gives stable, storable plates with consistent surface display. The trade-off is that we're measuring binding to a denatured-receptor-presenting surface, which is fine for screening but Tier 3 plaque assays use live cells to confirm receptor-specific infection. This is the Boeckaerts 2024 + Latka 2021 protocol; we're following the published method.
**Fallback:** "The protocol is adapted from Boeckaerts 2024 — let me pull the exact rationale from the paper."

### W7. "ΔexbD2 is the negative control — what if it doesn't infect either?"
**Confident answer:** Then we know something's wrong with the knockout system itself — maybe a polar effect from the deletion disrupted neighboring genes, or the antibiotic selection introduced an off-target mutation. The whole point of ΔexbD2 is **system validation**: Hung 2003 showed it shouldn't matter, so if it does matter in our hands, our other knockouts are also suspect. We'd verify via PCR + complementation rescue before trusting the ΔtonB / ΔexbB / ΔexbD1 data.
**Fallback:** "Good — that would be a flag I'd want to debug together. Let me think about the polar effect possibility."

### W8. "How do I know my Cycle 0 ELISA plate gave good data?"
**Confident answer:** Three checks built into the protocol. **R² > 0.9** on the 4PL fit (anything lower marks the EC50 as low-confidence). **WT-RBP positive control** within 2× of the historical EC50 baseline (inter-plate normalizer). **BSA-only background** below threshold (~0.1 OD~450~). All three are in `qc_report.md` template. If any fails, the model down-weights that plate's data in retraining.
**Fallback:** "There's a QC template in `08_cycle_data/templates/` — let me bring it up."

### W9. "Can I add an extra variant of my own choosing?"
**Confident answer:** Yes — anything beyond the 5 picks goes in as `selection_reason = "expert_pick"`. The model treats it as additional training data. We do ask that you note the rationale in `qc_report.md` so we can analyze whether expert picks systematically agreed with or diverged from BALD's recommendations — that's also publishable.
**Fallback:** "Yes — we should standardize the tag. Let me follow up on the exact convention."

### W10. "What's the absolute minimum experimental output you need from one cycle?"
**Confident answer:** **3 valid EC50s with R² > 0.9.** Below that, the model doesn't get enough new signal to meaningfully update. If we're at 2 valid measurements, we'd retrain anyway but flag the cycle as low-confidence in PI briefing and skip BALD selection — just use safe-pick backup for the next cycle.
**Fallback:** "Let me check the spec — minimum data thresholds are in the cycle infrastructure doc."

---

## [DRY LAB Q&A] — Ryan, Leah, ESM-experienced teammate

### D1. "Why deep ensemble instead of MC Dropout?"
**Confident answer:** Calibration. Ovadia 2019 (NeurIPS) showed MC Dropout is poorly calibrated under dataset shift, while deep ensembles maintain better-calibrated predictive intervals. For active learning, calibration matters more than raw accuracy — a badly-calibrated UQ means BALD picks the wrong variants. The cost is 5× training compute, which is trivial here (3 seconds end-to-end for Cycle 0). MC Dropout's only real advantage — single model fits in memory — doesn't matter at our scale.
**Fallback:** "Calibration is the headline difference — let me pull up Ovadia 2019's calibration plots if you want details."

### D2. "Why ESM-2 8M locally and 650M on Laguna — what does scale buy us?"
**Confident answer:** Lin 2023 showed embeddings of structural quality improve roughly with parameters but sublinearly past 650M. 8M is **proof-of-concept** — verifies the embedding pipeline end-to-end without GPU. 650M (1280-D vs 320-D) gives genuinely better structural features and is the production target. 3B is for the final Cycle 4+ benchmark. The PLM-interact paper (Liu 2025) is even more interesting — fine-tuning 650M on PPI data beats raw 3B for interaction prediction, suggesting scale-vs-fine-tune is an interesting axis we'll explore.
**Fallback:** "Empirically — let me pull the ESM-2 scaling figure from Lin 2023."

### D3. "How do you handle the cold-start with no real ELISA labels?"
**Confident answer:** Two answers. **Today** — Cycle 0 trains on synthetic data with documented structure: low-rank target with noise, used only to verify the pipeline runs and the BALD scoring sorts as expected. The synthetic-data plot is in calibration.png. **June 1 onward** — the first real ELISA measurements (4 variants from gene synthesis, structure-based picks) become the seed training set. Cycle 1 retrains on that. We're effectively bootstrapping from one Cycle 0 batch — small but sufficient because each member is regularized by Gaussian NLL.
**Fallback:** "It's bootstrap from Cycle 0 — let me walk through the exact data flow."

### D4. "What's the data contract between modules?"
**Confident answer:** Each module folder has `inputs/`, `processes/`, `outputs/`. `inputs/` is read-only pointers to upstream `outputs/`. Code lives **only** in `processes/`. Outputs are canonical and accompanied by a `MANIFEST.csv` with SHA-256 + size + record count per file. Schemas are documented in `INTERFACE.md`. The most important contract: predictions.csv must have `rbp_id, receptor_id, predicted_score, std, epistemic_std` — BALD reads `epistemic_std` directly.
**Fallback:** "Documented in `INTERFACE.md` — let me show you."

### D5. "Why BALD's `epistemic_std` and not total `std`?"
**Confident answer:** Total `std` includes **aleatoric** uncertainty — irreducible measurement noise (ELISA pipetting variance). If we sorted by total `std`, BALD would prefer pairs where the measurement is noisy. That's the opposite of what we want — we want pairs where **the model is uncertain**, not where the experiment is uncertain. `epistemic_std` is the variance-of-member-means component; it's what actually shrinks with more data. This is in `bald.py:38–52` and the verbatim Lakshminarayanan 2017 eq. 3 decomposition.
**Fallback:** "Aleatoric vs epistemic — the decomposition is in `ensemble.py:286–294`. Let me pull it up."

### D6. "What's the test coverage story? How do you know the pipeline isn't silently broken?"
**Confident answer:** Each module has pytest suites — 140+ tests total across modules 00–07. Module 07 (BALD) has 18 tests covering edge cases like empty pool, all-zero scores, measured-pair exclusion, seed reproducibility. Module 06 has 9 tests including the ensemble diversity assertion (`frac_diverse > 0.5` in `run_cycle0.py:194`). Module 05 has 28 tests + 1 expected skip pending a completed GPU run. Pass rates are in the "Tests — current pass rates" slide.
**Fallback:** "Test counts per module are on the tests slide; happy to walk through specific suites."

### D7. "Why is Module 07 written as `.py` while the rest are notebooks?"
**Confident answer:** Production orchestration, not exploration. Module 07 runs on a 48-hour SLA — it needs to be a deterministic, idempotent CLI we can invoke from a cron or shell script. Notebooks are great for `inputs → tinker → outputs` exploration but bad for `cron job → CI/CD → deploy`. CLAUDE.md flags this exception explicitly. All other modules will graduate from notebook to `.py` once stable; Module 07 was born production.
**Fallback:** "It's a deliberate exception in CLAUDE.md — production code path needs CI-friendliness."

### D8. "What happens to old model checkpoints when you retrain?"
**Confident answer:** Each cycle gets its own subfolder: `06_uncertainty_model/outputs/cycle_<N>/`. The model_version string in predictions.csv is `<git_sha>_cycle_<N>`, so every prediction is traceable to a specific git commit + cycle. Old `ensemble_member_*.pt` checkpoints are kept indefinitely (small files; <10 MB total per cycle). Calibration plots per cycle make it easy to diff cycle-over-cycle. MLflow tracking will be wired in for Cycle 1.
**Fallback:** "Provenance is via git SHA + cycle number — let me show the model_version field."

### D9. "How would you extend this to a generalist phage-host model?"
**Confident answer:** The Module 04 embedding step is already generalist — any phage RBP and any bacterial receptor get the same ESM-2 treatment. To extend, we'd train Module 06's ensemble on the full 2,236-pair Module 01 interaction matrix (currently we use synthetic in Cycle 0). The Mutalik 2025 PAML benchmark gives us 0.67–0.70 cross-genus AUC as a baseline — anything above is a contribution. The post-Jamboree research roadmap explicitly includes "train on soil-derived phage genomes" as Phase 2.
**Fallback:** "It's the project's natural extension — let me show the post-Jamboree plan."

### D10. "Have you actually run this end-to-end?"
**Confident answer:** Yes — modules 00–07 run start-to-finish in under 5 minutes on a laptop, plus the 3-minute Boltz-2 run on Laguna. The demo runbook in `DEMO.md` is the exact sequence. The only step that requires external data is Module 08 (ELISA), which is what June 1 starts producing. All output CSVs are committed (synthetic for Module 06 onward), so the contract is locked.
**Fallback:** "Live demo coming up — happy to walk through any module."

---

## [PI Q&A] — Prof. Ignacio-Espinoza, Prof. Libeskind-Hadas

### P1. "What's the iGEM-medal-relevant differentiator over existing tools?"
**Confident answer:** Three. **First**, no published system pairs ESM-2 embeddings + deep ensemble UQ + BALD-driven closed-loop ELISA for phage-RBP × bacterial receptor — Hie 2024 did antibodies, Yang 2025 did enzymes, Boeckaerts 2024 did host-range prediction without active learning. We're the first to combine them in this domain. **Second**, the Tier 3 receptor-knockout validation gives causality, not just correlation — that's what makes it paper-grade beyond iGEM. **Third**, the in-house *Xanthomonas* and phage isolates are a community contribution to the iGEM Registry + NCBI.
**Fallback:** "The combination is novel — let me detail the closest prior work."

### P2. "How do you know the model is actually learning vs random?"
**Confident answer:** The random control arm. Every cycle, 1 of the 5 picks is sampled randomly from the unmeasured pool, blinded to wet lab. At project end, we compare BALD's trajectory of test R² and information gain vs the random arm's trajectory. If they're indistinguishable, BALD didn't help. We also track per-cycle ECE and KL-divergence of the posterior — those are model-internal sanity checks. The whole comparison is pre-registered (in PLAN.md + risks slide).
**Fallback:** "Control arm + pre-registered metrics — happy to walk through the experimental design."

### P3. "What's the failure mode that would invalidate the whole approach?"
**Confident answer:** Two scenarios. **One** — ELISA's quantitative signal-to-noise is too low to distinguish meaningful binding gains from pipetting variance (aleatoric dominates epistemic). We've reserved two weeks of Cycle 0 for ELISA optimization specifically to de-risk this. **Two** — ESM-2 embeddings turn out to lack receptor-relevant signal (e.g., specific surface loops matter and mean-pooling washes them out). The fallback is per-residue or attention-pooled embeddings, which we'd switch to between cycles. Either failure is itself a publishable negative result if we report honestly.
**Fallback:** "Two main risks — ELISA noise and ESM-2 signal. Let me detail the mitigations."

### P4. "Why these papers and not others?"
**Confident answer:** Each citation maps to a specific design decision. Lakshminarayanan 2017 for ensemble; Houlsby 2011 for BALD objective; Lin 2023 for ESM-2; Boeckaerts 2022 for PhageRBPdetect HMM; Hung 2003 for the receptor system that defines our Tier 3 controls. We do NOT lean on ALDE (Yang 2025) as a BALD validation — they used Thompson sampling, different acquisition. We're explicit about that distinction in slides and in `paper_reading_notes.md`. The May 2026 paper audit corrected five mis-citations that earlier drafts contained — those are now fixed.
**Fallback:** "Each paper backs a specific decision — happy to walk citation by citation."

### P5. "What's your iGEM gold-medal critical path?"
**Confident answer:** Three deliverables drive the medal. **Best Composite Part** — at least 4–6 registered RBP-His6 constructs by wiki freeze. **Best Model** — the closed-loop AL pipeline + at least 2 wet-lab cycles' worth of data showing BALD's behavior (positive or honestly negative). **Best Agriculture Project** — the in-house *Xanthomonas* + lytic phage from California brassicas, sequenced and deposited. Tier 3 receptor-knockout validation is the gold-medal-grade differentiator on the Model track. Timeline is tight but achievable if we start cloning May 17.
**Fallback:** "Three deliverables — let me walk the timeline."

### P6. "What happens to this code after iGEM?"
**Confident answer:** Open source on the existing GitHub repo (active-learning-pipeline branch is the canonical one — main is intentionally empty). The May 2026 audit + bilingual comments + per-module READMEs make it actually onboardable. We've discussed concurrent manuscript submission to *Bioinformatics* or *NAR Genomics* (pending PI decision before Cycle 0). The repo, the in-house isolate sequences, and the trained model checkpoints would all be public deliverables.
**Fallback:** "Manuscript timing is a pending PI decision — let me flag it."

### P7. "Have you talked to anyone outside the team about the methodology?"
**Confident answer:** Not formally yet. The methodology survey is in `docs/planning/iGEM_2026_Project_Plan.md` and was reviewed at the May 7 PI consultation. We have informal connections to the ALDE group's published code (Frances Arnold lab at Caltech) which we could leverage if needed. Before submitting the manuscript I'd want at least one external methodology review — happy to ask you for suggestions of reviewers.
**Fallback:** "Not yet — would value your suggestions for external reviewers."

### P8. "What's the single biggest risk you lose sleep over?"
**Confident answer:** Honest answer: **strain isolation failing in May/June.** Everything downstream — Tier 3 validation, composite part, agriculture project — chains on having an in-house *Xanthomonas* isolate plus a lytic phage. We have dual-source plans (brassica + citrus) and Phage Directory as a backup phage source, but if the May sampling rounds yield nothing, we'd have to pivot the wet-lab story heavily. That's the only "no good fallback" risk on the table. Dry-lab can keep moving regardless, but the medal story needs wet-lab.
**Fallback:** "Strain isolation is the critical-path risk — let me detail the dual-source plan."

---

## Bonus: Hostile-listener questions (cycle ~10 minutes of paranoia)

### H1. "You said ipTM 0.365 is moderate — but doesn't that mean Boltz-2 is essentially failing to predict the interface?"
**Answer:** Yes, that's right — Boltz-2's confidence is low on this interface. That's the **expected** result for a novel phage-receptor pair that's never been crystallized. The model has no homologous complex in training data. The value of the run isn't a confident structure — it's the **monomer fold** (chain pTM 0.808) which gives us a trustworthy scaffold for site-directed mutation design, and the **acknowledgement of interface uncertainty** which is exactly what active learning is meant to resolve.

### H2. "If your synthetic Cycle 0 calibration looks perfect, isn't that just self-consistency theater?"
**Answer:** Correct — synthetic-data calibration is a sanity check that the pipeline doesn't break, not a validation of model quality. It tells us `epistemic_std` exports correctly, `select_batch` excludes measured pairs, and the BALD ranking matches the spec. The real calibration test is Cycle 1+ with actual ELISA data. We flag this explicitly on the calibration slide and in the predictions.csv `data_source = synthetic_fallback_random` metadata field.

### H3. "Why should the PI trust your interpretation of Lee 2009 over the published consensus?"
**Answer:** The "published consensus" on phiL7's tail spike is precisely what Lee 2009 says — which is that they couldn't find one by BLAST homology to OP1 ORF25. We're not contradicting them; we're saying HMM search (a more sensitive method, published 2022 in Boeckaerts) found a candidate they couldn't have found with 2009 tools. The framing "rbp_01 = phiL7 tail spike" is a hypothesis to test — ELISA + plaque assay on ΔtonB tests it directly. We don't claim it's confirmed.

### H4. "Active learning has been around since 1956 and is taught in undergrad ML. What's actually new here?"
**Answer:** The combination + the application. BALD on deep ensemble regression with ESM-2 embeddings, closed-loop with quantitative ELISA, on phage-RBP × bacterial-receptor, with receptor-knockout causal validation — that exact stack has no published precedent we could find. Each component is known; the integration and the wet-lab grounding in a phage system are the contribution. If a reviewer points to prior work we missed, we'd want to update the lit review immediately.

### H5. "How do you handle epistemic-aleatoric decomposition when the ELISA EC50 itself comes from a noisy 4PL fit with its own confidence interval?"
**Answer:** We currently treat the EC50 point estimate as the label, with R² > 0.9 as a quality filter. The 4PL fit's parameter standard error could be incorporated as observation noise in the Gaussian NLL — that's a Cycle 2+ refinement. Right now we down-weight noisy plates (R² 0.7–0.9 enter with halved weight; <0.7 excluded). Per-measurement aleatoric is on the post-Cycle-1 roadmap.

### H6. "If Hie 2024 only used 20 variants per antibody, why are you committing to 4–5 per cycle × many cycles?"
**Answer:** Different problem. Hie 2024 was language-model likelihood filtering — a single batch of candidate variants ranked once, not a closed loop. Our ALDE-analog (Yang 2025) used larger batches (~90 variants per round, 2–3 rounds). We're between those — small batches per cycle because ELISA optimization in Cycle 0 limits throughput, larger total over many cycles. If wet-lab throughput grows, we'd scale batch size up.

### H7. "Your model_version is git SHA + cycle. What if I want to reproduce Cycle 1 six months from now?"
**Answer:** `git checkout <sha>` recreates the code state. `06_uncertainty_model/outputs/cycle_1/ensemble_member_*.pt` has the trained weights. The synthetic data uses a seeded RNG with seed=42 logged in model_meta.json. Real-data reproduction needs the ELISA CSV from `08_cycle_data/outputs/cycle_1/elisa_processed.csv`, which is gitignored but checksummed in MANIFEST.csv. Six months from now you'd reproduce from those three artifacts.

### H8. "You said wet lab doesn't know which pick is random. How is that enforced in practice?"
**Answer:** The `recommendations.csv` rows are shuffled and the `selection_reason` column is **not** sent to wet lab — only `rbp_id, receptor_id, predicted_score, primer_sequence` are. The full CSV with `selection_reason` is archived in dry-lab outputs for the retrospective comparison. It's an honor system — Sarah, Olivia, Weitao, and Carol have to not peek — but the file separation makes it operationally easy.

### H9. "Lakshminarayanan 2017 also uses adversarial training. Are you using it?"
**Answer:** No — we use the Gaussian NLL component but not the adversarial training augmentation. Their ablations showed NLL alone is most of the win for in-distribution prediction; adversarial training mainly helped with out-of-distribution robustness. For Cycle 0 our pool is in-distribution by construction, so we skip it. This is a deliberate choice, not an oversight — happy to add adversarial training if calibration is poor after Cycle 1.

### H10. "What does 'Module 03 ML track was blocked' actually mean?"
**Answer:** The XGBoost track of PhageRBPdetect requires precomputed feature vectors and a model file that we couldn't reproducibly install via conda/pip — it depends on an XGBoost version mismatch with our pinned env. Boeckaerts 2022 recommends HMM as the primary track anyway; the ML track is a secondary precision booster. We left it as a "Cycle 2+ enhancement" rather than block the pipeline on it. The 25+ tests on Module 03 confirm the HMM track works end-to-end.

### H11. "Your TonB is 604 aa — is that the full protein including transmembrane, or just the periplasmic domain?"
**Answer:** Full protein — sequence pulled from GCF_000007145.1 as annotated. The transmembrane anchor is in the N-terminal ~30 aa; the periplasmic domain that actually contacts RBPs is the bulk. Boltz-2 predicts the whole thing; the periplasmic side dominates the chain pTM and is where any meaningful interface with rbp_01 would form. We can crop to periplasmic-only as a sensitivity check if ipTM stays low after Cycle 1 — that's on the roadmap.
