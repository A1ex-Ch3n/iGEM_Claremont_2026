# Onboarding Package — Plan

**Author:** Claude (onboarding agent, session 2026-05-17)
**Target audience:** Wet lab (Sarah, Olivia, Weitao, Carol), dry lab contributors (Ryan, Leah), PI (J. Cesar Ignacio-Espinoza), faculty advisor (Ran Libeskind-Hadas).

---

## Format choice

**Slides:** Pandoc + Beamer + xelatex + xeCJK (verified working with `Hiragino Sans GB`).

Rationale (post-pivot):
- I tried Marp first; `npx @marp-team/marp-cli` hangs on the initial Chromium download (this machine has no headless browser cached). Killed.
- Pandoc + Beamer + xelatex was already on the machine; `tlmgr install beamer xecjk ctex fontspec` resolved the missing classes. A 3-slide bilingual hello-world rendered cleanly (21 KB PDF, CJK glyphs intact).
- `Hiragino Sans GB` is the working macOS CJK font for both matplotlib (Simplified Chinese) and xeCJK. Confirmed via `fc-list`.
- Markdown source is what the team edits; the `.tex` and `.pdf` are build artifacts.

Commit the `.md` source AND the rendered PDF. Build command goes in the guide.

**Chinese variant:** Simplified Chinese throughout (matches `docs/planning/iGEM_2026_项目大纲_中文版.md`). The Traditional Chinese in `paper_reading_notes.md` stays where it is.

**Written guides:** Plain Markdown (`guide_en.md` / `guide_zh.md`) — directly readable in any editor or rendered on the iGEM wiki.

**Demo runbook:** Plain Markdown (`DEMO.md`).

---

## Visuals to generate (matplotlib, into `docs/onboarding/figures/`)

1. `pipeline_flow.png` — modules 00→08 left-to-right with data-flow arrows and per-module 1-line outputs.
2. `data_contract.png` — the `inputs/ → processes/ → outputs/` block per module, with arrows showing how outputs feed downstream inputs.
3. `cycle_48h.png` — 48-hour dry-lab cycle: ELISA arrives → retrain ensemble → BALD → recommendations.csv → wet lab clones/expresses → ELISA → loop.
4. `bald_intuition.png` — toy plot showing predicted-vs-uncertainty across candidates, with BALD top picks highlighted.

Existing figures embedded as-is:
- `pae_heatmap.png` (repo root) — PAE matrix for rbp_01 × TonB.
- `06_uncertainty_model/outputs/cycle_0/calibration.png` — calibration plot.
- `archive/2026-05-pivot/05_predictive_modeling/outputs/baseline_taxonomy_knn/heatmap_combined.png` — historical baseline (only on a "where we came from" slide).

---

## Slide outline (EN — ~60 slides, 60–90 min deep dive)

Each section header is tagged so the presenter can skip per audience.

### Part 0 — Title & roadmap (5 slides)
1. Title + team + date.
2. Roadmap (5 parts).
3. TL;DR (3 bullets: closed loop, Modules 00–07 done, Cycle 0 starts ~June 1).
4. How to read this deck (the [WET LAB] / [DRY LAB] / [PI] tags).
5. **Suggested cuts:** 30-min wet-lab cut, 45-min PI cut, 90-min full team — list slide numbers per cut.

### Part 1 — The science [WET LAB] [PI] (8 slides)
5. *Xanthomonas* — what & why we care (Ryan 2011, Iriarte 2018, copper-resistance gap).
6. Phage biocontrol & the host-range problem.
7. Our reference scaffold: phiL7 + Xcc ATCC 33913.
8. The receptor system — Hung 2003 (TonB-ExbB-ExbD1 essential, ExbD2 NOT — free negative control).
9. RBP = the "key" — Lee 2009 explicitly searched and found NO OP1 ORF25 homolog; rbp_01 is a HMM complementary discovery.
10. The data-scarcity bottleneck (PhageHostLearn AUC 0.82 within genus → 0.67 across).
11. Why active learning is the right framing (Settles 2009, Hie 2024 antibody, Yang 2025 ALDE).
12. iGEM tracks targeted: Agriculture / Model / Composite Part.

### Part 2 — Pipeline architecture [DRY LAB] [PI] (10 slides)
13. The 3-layer architecture (Layer 0 priors → Layer 1 AL loop → Layer 2 causal validation).
14. Pipeline diagram (00→08) — `pipeline_flow.png`.
15. Data-contract convention (`inputs/ / processes/ / outputs/`) — `data_contract.png`.
16. Bilingual notebook-first workflow (CLAUDE.md convention).
17. Module 00 — Raw Data (777+34 genomes; MANIFEST).
18. Module 01 — Ground Truth (2,236 pairs; phiL7×Xcc confirmed) + interaction_matrix.csv head preview.
19. Module 02 — Annotation (PHANOTATE phage / pyrodigal bacteria — NEVER swap).
20. Module 03 — RBP ID (PhageRBPdetect; rbp_01 = 712 aa primary).
21. Module 04 — Embedding (ESM-2; 1280-dim production).
22. Module 05 — Structure (Boltz-2 ipTM, NOT affinity for protein-protein).

### Part 3 — The ML core [DRY LAB] [PI] (10 slides)
23. Module 06 — Deep Ensemble overview (Lakshminarayanan 2017) + Greenman 2025 "no single best UQ method" — why ensemble is still our choice.
24. (mean, sigma) per member — `ensemble.py` snippet.
25. Epistemic vs aleatoric uncertainty (Var_k[μ_k] vs E_k[σ_k²]).
26. Calibration plot (`calibration.png`).
27. Predictions table preview (head of predictions.csv).
28. Module 07 — BALD intuition (`bald_intuition.png`).
29. BALD math: I(y;θ|x,D) ≈ Var_k[μ_k] — with the Houlsby 2011 → deep-ensemble-regression extension note.
30. `bald.py` snippet: `bald_score()` + `select_batch()` with selection_reason audit.
31. Recommendations table preview (head of recommendations.csv).
32. ALDE caveat — Yang 2025 uses Thompson, not BALD; one-hot, not ESM-2.

### Part 4 — The current Boltz-2 result [WET LAB] [DRY LAB] [PI] (4 slides)
33. rbp_01 × TonB — what we ran (job 59986, 712aa correct after the 85aa P25 fix).
34. ipTM 0.365 + chain_A_ptm 0.808 — what those numbers mean.
35. PAE heatmap (`pae_heatmap.png`) — interface block reading guide.
36. Boltz-2 affinity head warning (NaN for protein-protein — Passaro 2025).

### Part 5 — The 48-hour cycle [WET LAB] [DRY LAB] (8 slides)
37. Cycle diagram (`cycle_48h.png`).
38. Wet lab → dry lab handoff (elisa_processed.csv schema).
39. Dry lab → wet lab handoff (recommendations.csv + primer_sequences.txt + safe_pick_backup.csv).
40. Cloning execution: Gibson Cycle 0, SDM Cycles 1+.
41. ELISA protocol summary (Boeckaerts 2024-style cell-based ELISA).
42. Receptor knockout system (pK18mobsacB; Hung 2003 confirms targets).
43. Validation tiers (ELISA only → +plaque → +ΔtonB/ΔexbB/ΔexbD).
44. Quality gates + failure modes.

### Part 6 — Reproducing & demoing [DRY LAB] [WET LAB] (6 slides)
45. Quick-start (`conda env create`, `git checkout active-learning-pipeline`).
46. Per-module entry points (the GETTING_STARTED checklists).
47. Where outputs live (the table from PI briefing).
48. Live demo plan (see DEMO.md): Module 03 → 04 → 06 → 07.
49. Laguna HPC primer (when to push to GPU).
50. Tests: `pytest <module>/processes/tests/` — pass rates per module.

### Part 7 — Risks, decisions, asks [PI] (6 slides)
51. Pending PI decisions (pK18mobsacB vs CRISPRi; AF3 weights; validation tier; phage source; manuscript ambition).
52. Critical risks + mitigations (strain isolation fails, RBP insolubility, AL underperforms random).
53. What "AL underperforms random" looks like and how we'd report it honestly.
54. Timeline (May 17 wet lab launch → Cycle 0 June 1 → Cycle 2 mid-July).
55. iGEM deliverables checklist.
56. Open questions for the team.

### Part 8 — References & appendix (4 slides)
57. The 19 papers per module (compact table).
58. Five literature corrections (the audit) — ExbD2, Lee 2009, Boltz-2 affinity, Greenman journal, Hie ESM-1b.
59. Glossary pointer + docs/ navigation.
60. Q&A / thank-you.

---

## Guide TOC (`guide_en.md`)

Mirrors slide structure but expanded. Each section ≈ 1–3 paragraphs.

1. Project at a glance (TL;DR from PI briefing).
2. The science (Xanthomonas, phage biocontrol, host-range problem, phiL7+Xcc background).
3. Pipeline architecture (diagram + data contract + each module's purpose & key output, with sample table previews).
4. The ML core in depth (ensemble, calibration, BALD math + code snippets).
5. Module 05 result (rbp_01 × TonB; how to read ipTM/PAE/pLDDT).
6. The 48-hour cycle protocol (file schemas + Gibson/SDM/ELISA/knockout summaries).
7. How to reproduce — environment, genome fetch, per-module commands, tests, Laguna pointer.
8. Conventions reference (notebook-first, bilingual comments, MANIFEST, INTERFACE).
9. References (compact paper table with 1-line role per paper) + 5 corrections.
10. Glossary + further reading.

`guide_zh.md`: same structure, full Chinese translation, all figures/paths/code identical.

---

## Demo runbook (`DEMO.md`)

Pre-meeting checklist (5 items):
- `git checkout active-learning-pipeline` + `git pull`
- `conda activate igem2026`
- Genome present: `EU717894.1.fna` and `GCF_000007145.1.fna`
- HMM pressed: `hmmpress 03_rbp_identification/inputs/phagerbpdetect_data/RBPdetect_phageRBPs.hmm`
- Browser tab open to PI_briefing.

Live steps (each: command, runtime, expected output, failure recovery):
1. **Module 03** — `pytest 03_rbp_identification/processes/tests/ -v` (~15 s, 25+ pass). Open `EU717894.1_rbp_candidates.csv` → show rbp_01 at 712 aa, HMM score 342.
2. **Module 04** — Run `01_embed_esm2.ipynb` cell 1 (load) + cell 2 (embed phiL7 RBPs, ~1 min CPU). Show `.npz` shape (3, 320).
3. **Module 06** — `python 06_uncertainty_model/processes/run_cycle0.py` (~25 s synthetic). Show `predictions.csv` head + `calibration.png`.
4. **Module 07** — `python 07_acquisition_function/processes/run_bald.py --cycle 1 --n_bald 4 --n_random 1` (<1 s). Show `recommendations.csv` — rbp_07×rec_02 BALD top.
5. **Module 05 (read-only)** — Open `pae_heatmap.png` and the PDB in PyMOL (if available); explain ipTM/PAE interpretation.

Each step has a "what to say" line (≤2 sentences) and a "if it fails" rescue (point at the AGENT_REPORT.md or the cached output already committed).

---

## Conventions I will respect

- Paper shorthand exactly as `docs/reference/papers.md` and `paper_reading_notes.md` document them — including the 5 corrections.
- Bilingual code comments shown verbatim in snippets.
- All file paths typed in full (no `...` abbreviations).
- Module 02: never call PHANOTATE/pyrodigal "interchangeable".
- Module 05: ipTM = structural confidence proxy; affinity = NaN for protein-protein.
- Module 07: BALD as extension of Houlsby 2011 (originally GPC); ALDE used Thompson, not BALD.
- Cite team names where appropriate (Alex / Sarah / Olivia / Weitao / Carol / Ryan / Leah / PI / advisor).
