# Onboarding Agent Prompt

Paste the block below into a fresh Claude Code session (auto-accept permissions ON) to
regenerate the onboarding package — slides + guide + live-demo runbook — for the iGEM
Claremont 2026 active-learning phage engineering pipeline.

The prompt is self-contained: the new session has no memory of prior conversations.

---

```
You are onboarding a new dry-lab + wet-lab team to the iGEM Claremont 2026 active-learning
phage engineering pipeline. Your job: read the entire codebase and produce a polished onboarding
package — slides + written guide + live-demo runbook — for three audiences (wet lab, dry lab,
PI/faculty). You have auto-accept permissions: install tools, run notebooks, generate figures,
write files freely.

REPOSITORY
- Working dir: /Users/alexy/Desktop/Claude Workspace/iGEM_Claremont_2026
- Active branch: active-learning-pipeline  (CHECK OUT FIRST — main is empty)
- Project: closed-loop ML pipeline that uses Bayesian Active Learning by Disagreement (BALD)
  to engineer phage receptor-binding proteins (RBPs) for Xanthomonas biocontrol. Pipeline has
  9 modules (00–08); Modules 00–07 are complete, Module 08 starts ~June 1 with wet-lab data.

READ THESE FIRST (in this order — they explain everything else)
1.  README.md
2.  CLAUDE.md                            ← conventions, team roles, data contracts
3.  GETTING_STARTED.md                    ← per-module guide with checklists
4.  docs/planning/PI_briefing_2026-05-11.md  ← current status, work log, all output paths
5.  docs/planning/iGEM_2026_Project_Plan.md
6.  docs/reference/papers.md              ← 19 key papers, per-module reading guide
7.  docs/reference/paper_reading_notes.md ← detailed notes + corrections (Lee 2009, Boltz-2, etc.)
8.  LAGUNA.md                             ← HPC setup

Then read each module's README.md and the key script(s) in processes/:
- 00_raw_data/                        (genome library + MANIFEST)
- 01_data_ground_truth/               (interaction matrix)
- 02_annotation/                      (PHANOTATE for phage / pyrodigal for bacteria — NEVER swap)
- 03_rbp_identification/              (PhageRBPdetect HMM — rbp_01 = 712aa tail spike)
- 04_protein_embedding/               (ESM-2)
- 05_structure_prediction/            (Boltz-2 on Laguna; ipTM, not affinity, for protein-protein)
- 06_uncertainty_model/processes/ensemble.py + run_cycle0.py  (5-member deep ensemble; epistemic_std)
- 07_acquisition_function/processes/bald.py + run_bald.py    (BALD = Std_k[μ_k]; Cycle 0 done)
- 08_cycle_data/                      (will fill in starting ~June 1)

DELIVERABLES — write all of these to docs/onboarding/

A. SLIDES (one unified deck, two language versions)
   - docs/onboarding/slides_en.{pdf,pptx,or beamer .tex — you choose the best format and install
     what's needed; commit the source AND the rendered PDF}
   - docs/onboarding/slides_zh.{same format, Chinese version}
   - Single deck per language, with section headers tagged [WET LAB] / [DRY LAB] / [PI]
     so the presenter can skip or expand sections per audience.
   - Depth: deep dive (60–90 min). Cover:
       * The scientific problem (Xanthomonas biocontrol, phage host-range engineering, why data scarcity)
       * Pipeline architecture (modules 00–08, data flow contract, the 48-hr cycle)
       * Each module's purpose, key tool, and one concrete output
       * Key algorithms: ESM-2 protein language model, Boltz-2 structure (ipTM not affinity!),
         5-member deep ensemble for epistemic UQ, BALD acquisition (Std_k[μ_k] for regression)
       * Wet-lab integration: ELISA → Module 08 → retrain → next recommendations.csv
       * Mention the key papers inline on the relevant slides — at minimum:
         - PhageRBPdetect (Boeckaerts 2022) — Module 03
         - ESM-2 (Lin 2023) — Module 04
         - Boltz-2 (Passaro 2025) — Module 05
         - Deep ensembles (Lakshminarayanan 2017) — Module 06
         - BALD (Houlsby 2011) + ALDE (Yang 2025, NOTE: Thompson sampling not BALD) — Module 07
         - Hung 2003 (TonB-ExbB-ExbD1 receptor system for phiL7)
         - Lee 2009 (phiL7 genome; explicitly searched and found NO OP1 ORF25 tail-spike
           homolog — rbp_01 is a complementary HMM-based discovery, not a contradiction)

B. WRITTEN GUIDES (separate EN + ZH, markdown)
   - docs/onboarding/guide_en.md
   - docs/onboarding/guide_zh.md
   - Self-contained: someone who missed the talk can read this and understand the project,
     run any module, and know where to find things. Mirror the slide structure but with
     more prose. Include a "How to reproduce" appendix.

C. LIVE DEMO RUNBOOK
   - docs/onboarding/DEMO.md
   - Copy-paste terminal commands to demo Modules 03 → 04 → 06 → 07 live during the talk.
     Each step: command, expected runtime, expected output (file path + 1-line description),
     and how to recover if it fails. Include a "before the meeting" prep checklist
     (conda env, branch checked out, genomes downloaded).

VISUALS — include all of these
1. Embed EXISTING figures from the repo. Find them by walking the tree; at minimum:
   - pae_heatmap.png (in repo root)
   - 06_uncertainty_model/outputs/cycle_0/calibration.png
   - any Boltz-2 PDB renders or ipTM figures from 05_structure_prediction/outputs/
2. GENERATE new diagrams:
   - Pipeline architecture flowchart (00 → 08 with data flow arrows)
   - Data-contract diagram (inputs/processes/outputs convention)
   - 48-hour cycle diagram (wet lab ↔ dry lab handoff)
   Use mermaid or matplotlib; render to PNG/SVG and embed in slides + guide.
3. CODE SNIPPETS (5–15 lines each) for:
   - bald.py: bald_score() and the top-k selection in select_batch()
   - ensemble.py: the (mean, sigma) per-member prediction and the epistemic_std computation
   - One example of the data-contract pattern (a process script reading inputs/ and writing outputs/)
4. SAMPLE OUTPUT TABLES (preview first 5–10 rows):
   - 07_acquisition_function/outputs/cycle_1/recommendations.csv
   - 06_uncertainty_model/outputs/cycle_0/predictions.csv
   - 01_data_ground_truth/outputs/interaction_matrix.csv (head)
   - 03_rbp_identification/outputs/EU717894.1_rbp_candidates.csv

CONVENTIONS TO RESPECT
- CLAUDE.md mandates bilingual comments and notebook-first workflow — explain these in the guide.
- Cite paper sources with the Lee 2009 / Hung 2003 / etc. shorthand the team already uses
  (see docs/reference/papers.md for the canonical list and corrections — particularly the
  Lee 2009, ALDE-vs-BALD, and Boltz-2 affinity-vs-ipTM nuances; do not re-introduce those errors).
- Team: Core engineer Alex Chen; Wet lab Sarah, Olivia, Weitao, Carol; Dry lab contributors
  Ryan, Leah; PI Prof. J. Cesar Ignacio-Espinoza; Faculty advisor Prof. Ran Libeskind-Hadas.

WORKFLOW
1. Check out branch, read the priority files above, walk the module tree.
2. Before producing anything, write a short plan (docs/onboarding/PLAN.md) listing the slide
   outline + guide TOC + demo steps. Call advisor() to sanity-check it.
3. Pick a slide format. Install whatever is needed (python-pptx, marp-cli, texlive, mermaid-cli).
   Pin choice in PLAN.md with a one-line rationale.
4. Build the EN slides + guide + DEMO.md first. Render slides to PDF and verify they open.
5. Translate to ZH (preserve all figures, code, paths, and paper citations).
6. Final pass: every file path mentioned in the slides/guide must actually exist in the repo
   (verify with ls/grep). No `...` abbreviations in paths — write them out in full.
7. Commit everything to active-learning-pipeline with a clear commit message. Do NOT push;
   leave the push to Alex.

WHEN DONE
Print a summary table: file path, purpose, size, language. Then call advisor() for a final review.
```
