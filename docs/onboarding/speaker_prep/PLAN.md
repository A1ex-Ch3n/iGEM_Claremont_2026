# Speaker Prep — Plan

**Audience:** mixed — wet lab (Sarah, Olivia, Weitao, Carol), dry lab (Ryan, Leah), PI (Ignacio-Espinoza), faculty advisor (Libeskind-Hadas). 60–90 min talk + Q&A tomorrow morning.

**Deliverables** (all under `docs/onboarding/speaker_prep/`, EN + ZH each):
- `SPEAKER_NOTES` — per-Part notes (time / say-this / don't-say / if-asked / transition), grouped by the 8 Parts of `slides_en.md`. Each Part covers all slides in it; vocabulary slides treated as one block.
- `CONCEPT_DEEPDIVE` — 7 hardest concepts (active learning, ESM-2 embedding, Boltz-2 ipTM, deep ensemble, BALD math, 48h cycle SLA, Lee 2009 HMM-vs-BLAST). Each: 30-sec / 2-min / full / analogies / verbatim sentence / failure modes.
- `QA_PREP` — 3 audience sections (wet/dry/PI), ~10/10/8 anticipated questions with confident model answer + honest fallback.
- `REHEARSAL` + `REHEARSAL_ANSWERS` — 3-tier self-quiz (cold facts / 1-min concept / design defense).
- `CHEAT_SHEET` — one-page printable safety net: numbers, citations, 3 verbatim sentences, freeze-recovery lines.
- `TIMING` — minute-by-minute target with 15/30/45/60 checkpoints, [CORE]/[OPTIONAL] tags, audience-priority annotations.

**Accuracy guardrails** (cross-checked against `paper_reading_notes.md` & code):
- Boltz-2 ipTM = structural confidence proxy, NOT binding affinity. Affinity head is NaN for protein-protein.
- BALD score = `epistemic_std = Std_k[μ_k]`, not total_std, not variance — Std. (`bald.py:38–52`).
- ALDE (Yang 2025) uses Thompson sampling + one-hot, NOT BALD + ESM-2 — only validates AL+UQ generally.
- Lee 2009 searched for OP1 ORF25 homolog and explicitly found none; rbp_01 = complementary HMM discovery (Tail_spike_N), not a contradiction.
- Hung 2003: TonB / ExbB / ExbD1 essential; ExbD2 NOT required → ΔexbD2 is the free negative control.
- Hie 2024 used ESM-1b/1v (not ESM-2), ~20 variants per antibody, two rounds, language-model likelihood filtering — NOT a BALD closed loop.
- Greenman 2025 is *PLoS Comput Biol* (not *NAR Genomics*); conclusion: "no single best UQ method."

**Workflow:** Files A→F in order; verify every numeric claim against source code or `paper_reading_notes.md`; final hostile-listener pass adds embarrassing-Q backstop; commit (do not push) on the worktree branch — Alex merges into `active-learning-pipeline`.
