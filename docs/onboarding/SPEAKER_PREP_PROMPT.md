# Speaker Prep Agent Prompt

Paste the block below into a fresh Claude Code session (auto-accept ON) **after** the
onboarding-package agent (see `AGENT_PROMPT.md`) has finished generating
`docs/onboarding/slides_en.*`, `guide_en.md`, `DEMO.md`, etc.

This agent does NOT produce audience-facing material. It produces **presenter-facing**
material: speaker notes, deep-dive on the 5–7 hardest concepts, anticipated Q&A by audience,
a self-quiz to verify understanding, and a one-page cheat sheet.

---

```
You are preparing Alex Chen (core engineer) to deliver a 60–90 minute onboarding talk
TOMORROW MORNING to a mixed audience: wet lab teammates (Sarah, Olivia, Weitao, Carol),
dry lab contributors (Ryan, Leah), and the PI/faculty advisor (Prof. Ignacio-Espinoza,
Prof. Libeskind-Hadas). Alex built the pipeline himself but needs to (a) lock in his own
understanding of every step's principle so he can teach in his own words, (b) anticipate
questions per audience, and (c) make sure listeners can themselves run and extend the codebase
after the talk.

You have auto-accept permissions. Read freely, generate figures if helpful, install tools
if needed. Do NOT modify the slides/guide/DEMO files — those are the audience artifacts.

REPOSITORY
- Working dir: /Users/alexy/Desktop/Claude Workspace/iGEM_Claremont_2026
- Branch: active-learning-pipeline (check out first)

INPUTS — read these before producing anything
1. The audience-facing onboarding package (already generated):
   - docs/onboarding/slides_en.*  (and slides_zh.*)
   - docs/onboarding/guide_en.md  (and guide_zh.md)
   - docs/onboarding/DEMO.md
   - docs/onboarding/PLAN.md
2. The source-of-truth project docs:
   - README.md, CLAUDE.md, GETTING_STARTED.md, LAGUNA.md
   - docs/planning/PI_briefing_2026-05-11.md
   - docs/planning/iGEM_2026_Project_Plan.md
   - docs/reference/papers.md
   - docs/reference/paper_reading_notes.md (CRITICAL — has the Lee 2009, Boltz-2 affinity,
     ALDE-vs-BALD corrections that Alex MUST get right)
3. The actual code for the trickiest modules:
   - 06_uncertainty_model/processes/ensemble.py + run_cycle0.py
   - 07_acquisition_function/processes/bald.py + run_bald.py + tests/
   - 03_rbp_identification/processes/01_run_phagerbpdetect.ipynb

DELIVERABLES — write all of these to docs/onboarding/speaker_prep/

A. SPEAKER_NOTES.md  (EN + ZH side-by-side per slide)
   For each slide in slides_en, write:
   - **Time**: how long to spend (target 60–90 min total)
   - **Say this**: 2–4 bullet points in Alex's natural voice (not robotic narration)
   - **Don't say**: common over-claims to avoid (e.g., don't call Boltz-2 ipTM "binding affinity";
     don't say BALD was validated by ALDE — ALDE used Thompson sampling)
   - **If asked X, say Y**: 1–3 likely interruptions and short answers
   - **Transition**: one-line bridge to the next slide

B. CONCEPT_DEEPDIVE.md  (the 5–7 hardest concepts, explained from first principles)
   For each concept, provide:
   - **The 30-second elevator pitch** (for the PI / wet lab who just want the gist)
   - **The 2-minute explanation** (for dry lab who wants the mechanism)
   - **The full derivation / mechanism** (for Alex to internalize)
   - **2–3 analogies** in different domains (so Alex can re-explain on the fly)
   - **The single sentence Alex should memorize verbatim**
   - **What's still unknown / where the method could fail**
   Concepts to cover (at minimum — add more if you find genuinely hard ones in the code):
   1. Why active learning solves the data-scarcity problem (vs random screening / supervised)
   2. ESM-2 — what a protein language model embedding actually IS, and why it works for RBPs
   3. Boltz-2 ipTM — what it measures, why it's a structural confidence proxy, why the
      affinity head is NaN for protein-protein (Alex MUST nail this — biggest misconception risk)
   4. Deep ensemble — why 5 independently-trained networks give epistemic uncertainty
      (vs MC-dropout, vs Bayesian NNs); the variance-decomposition: total² = epistemic² + aleatoric²
   5. BALD as Std_k[μ_k] for regression — derive from Houlsby 2011 mutual-information form,
      explain why it reduces to ensemble variance under Gaussian assumptions, and why this is
      "epistemic only" (not total)
   6. The 48-hour closed-loop cycle — why the SLA matters, what breaks if it slips
   7. Why Lee 2009 finding "no OP1 ORF25 homolog" + our HMM-found rbp_01 is complementary,
      not contradictory (BLAST vs HMM detection thresholds)

C. QA_PREP.md  (three audience sections — anticipated questions with model answers)
   - **[WET LAB Q&A]** (~10 questions) e.g.:
     "What do I actually do when recommendations.csv arrives?"
     "What if a variant doesn't express in BL21?"
     "Why do we need a random control if BALD is the smart pick?"
     "What if ELISA fails for the whole batch?"
   - **[DRY LAB Q&A]** (~10 questions) e.g.:
     "Why deep ensemble instead of MC-dropout?"
     "Why ESM-2 8M locally and 650M on Laguna — what does scale buy us?"
     "How do you handle the cold-start with no real ELISA labels?"
     "What's the data contract between modules?"
   - **[PI Q&A]** (~8 questions) e.g.:
     "What's the iGEM-medal-relevant differentiator over existing tools?"
     "How do you know the model is actually learning vs random?"
     "What's the failure mode that would invalidate the whole approach?"
     "Why these papers and not others?"
   For each question: give a CONFIDENT 3–5 sentence model answer, plus a "fallback if you don't
   know" honest response (e.g., "Good question — we're tracking that in cycle 2; let me follow up").

D. REHEARSAL.md  (self-quiz Alex does TONIGHT before bed and AGAIN over coffee tomorrow)
   - **Tier 1 — must answer in <10 seconds, cold** (~20 questions):
     numerical facts (777 phages, 712 aa rbp_01, ipTM=0.365, 5-member ensemble, 18 tests pass,
     2236 interaction pairs, etc.), tool-to-module mapping, paper-to-claim mapping.
   - **Tier 2 — must explain in 1 minute** (~10 questions): each concept in CONCEPT_DEEPDIVE.
   - **Tier 3 — must defend the design choice** (~5 questions):
     "Why BALD and not Thompson sampling?", "Why ESM-2 and not ProtBERT?",
     "Why Boltz-2 instead of just AlphaFold?", "Why ipTM instead of pLDDT?"
   Answers go in a separate file (REHEARSAL_ANSWERS.md) so Alex can self-test without peeking.

E. CHEAT_SHEET.md  (ONE page, printable — Alex's safety net during the talk)
   - All key numbers in one table (genome counts, pair counts, dimensions, ipTM, test counts)
   - All paper citations in shortform (Author Year — one-line claim)
   - The 3 sentences Alex MUST get exactly right (Boltz-2 ipTM ≠ affinity; BALD = epistemic_std;
     Lee 2009 complementary not contradictory)
   - Three "if I freeze" recovery lines ("Let me pull up the actual code…",
     "That's an open question — let me note it down…", "Good question, the short answer is X,
     the long answer is in the guide section Y")
   - One-line URL/path to slides, guide, DEMO.md for the audience to follow along

F. TIMING.md  (minute-by-minute target with checkpoints)
   - Total target: 75 minutes (with 15 min buffer for Q&A)
   - Checkpoint at minute 15, 30, 45, 60 — if behind, which sections to compress
   - Mark each slide as [CORE — cannot skip] or [OPTIONAL — skip if running long]
   - Mark each section by audience priority — e.g., "BALD math: critical for dry lab + PI;
     skim if wet-lab-heavy crowd"

ZH VERSIONS
For each deliverable above, also produce a `_zh.md` companion in the same folder. Maintain
the same structure; translate idiomatically (not literally). Alex's audience is bilingual;
he'll switch languages on the fly, so the prep material must be in both.

WORKFLOW
1. Check out branch. Read inputs. Walk the codebase to verify everything you're going to
   teach is accurate (don't rely only on docs — open the code).
2. Write a 1-paragraph plan to docs/onboarding/speaker_prep/PLAN.md and call advisor() to
   sanity-check coverage of the 7 concepts and 3 audiences.
3. Produce A → F in order. Each file: verify every claim against source code or papers.md.
4. Final pass: simulate being a hostile listener for 10 minutes — what's the question that
   would embarrass Alex? Add it to QA_PREP.md.
5. Commit to active-learning-pipeline with message "onboarding: add speaker prep materials".
   Do NOT push — leave that to Alex.

CRITICAL ACCURACY GUARDRAILS — if you violate any of these, Alex will lose credibility
in front of the PI. Verify each one against docs/reference/paper_reading_notes.md:
- Boltz-2 affinity head returns NaN for protein-protein. We use ipTM (interface confidence,
  0–1). Do NOT call it "binding affinity" or "zero-shot affinity prior."
- BALD score for our regression case = epistemic_std = Std_k[μ_k(x)] across ensemble members.
  Not total_std. Not Var. Std.
- ALDE (Yang 2025) uses Thompson sampling and one-hot encoding, NOT BALD with ESM-2. It
  validates "active learning + UQ works in protein engineering" generally — it does NOT
  validate our exact method.
- Lee 2009 explicitly searched for an OP1 ORF25 tail-fiber homolog in phiL7 and found none.
  Our rbp_01 (712 aa) was found via the Tail_spike_N HMM profile — HMM captures structural
  homology beyond BLAST's detection threshold. Frame this as COMPLEMENTARY discovery, not
  as "Lee 2009 missed it" or "Lee 2009 was wrong."
- Hung 2003 (BBRC 302:878-884) showed phiL7 receptor system = TonB + ExbB + ExbD1.
  ExbD2 is NOT required. (Older project docs had this wrong — do not regress.)

WHEN DONE
Print a one-paragraph summary: which concepts you flagged as highest-risk for Alex to nail,
which Q&A questions you expect to come up most often, and any spots in the slides/guide
where you noticed an inaccuracy that the slides agent should fix.
```
