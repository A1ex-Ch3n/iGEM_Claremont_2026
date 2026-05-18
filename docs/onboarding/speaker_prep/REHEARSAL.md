# REHEARSAL.md — Self-quiz (no peeking at answers)

**Use:** Tonight before bed (cold pass). Tomorrow with coffee (warm pass). Mark each Q as ✅ / 🤔 / ❌. Anything ❌ → re-read the relevant slide/concept. Answers in `REHEARSAL_ANSWERS.md`.

---

## Tier 1 — Must answer in <10 seconds, cold (20 questions)

### Numerical facts
1. How many phage genomes in Module 00?
2. How many bacterial genomes in Module 00?
3. How many pairs in the Module 01 interaction matrix? (total / positive / negative)
4. Length of rbp_01 in amino acids?
5. HMM bit score for rbp_01 against Tail_spike_N?
6. How many RBP candidates did Module 03 find for phiL7?
7. ESM-2 8M embedding dimension?
8. ESM-2 650M embedding dimension?
9. ESM-2 3B embedding dimension?
10. ipTM for rbp_01 × TonB Boltz-2 run (job 59986)?
11. chain_A_ptm for rbp_01 in same run?
12. Number of ensemble members in Module 06?
13. Number of tests passing in Module 07?
14. Default n_bald and n_random per cycle?
15. Dry-lab SLA for cycle turnaround (hours)?

### Tool-to-module mapping
16. Which tool annotates phage ORFs? Which tool annotates bacterial ORFs?
17. Which tool finds RBP candidates? Two sub-tracks?
18. Which two structure-prediction tools are in Module 05?

### Paper-to-claim mapping
19. Which paper established phiL7 receptor system? What did it find essential / not essential?
20. Which paper introduced BALD? On what type of model?

---

## Tier 2 — Must explain in 1 minute (10 questions)

These correspond to the 7 concepts in CONCEPT_DEEPDIVE.md plus 3 stretch.

21. Why active learning solves data scarcity (vs random / vs supervised).
22. What an ESM-2 embedding actually IS — explain to a wet-lab teammate.
23. What ipTM measures and why it is NOT a binding affinity.
24. Why a deep ensemble gives epistemic uncertainty.
25. Derive (or sketch) why BALD ≈ Var_k[μ_k] for our setup.
26. Why the 48-h SLA matters and what breaks if it slips.
27. Why Lee 2009 + our HMM rbp_01 are complementary, not contradictory.
28. What ΔexbD2 is doing in the knockout panel.
29. Explain the difference between epistemic and aleatoric uncertainty in one breath.
30. What the random control arm is for, in one sentence.

---

## Tier 3 — Must defend the design choice (5 questions)

31. Why BALD and not Thompson sampling (like ALDE)?
32. Why ESM-2 and not ProtBERT (or ProtTrans, or a hand-crafted feature set)?
33. Why Boltz-2 instead of just AlphaFold 3?
34. Why ipTM instead of pLDDT as the structural-confidence proxy?
35. Why deep ensemble instead of MC Dropout or a Bayesian neural net?

---

## Stretch (bonus — only if Tier 1–3 all green)

36. What ECE is and how temperature scaling fixes it.
37. What `safe_pick_backup.csv` is for and when wet lab uses it.
38. Why the wet lab doesn't know which row in recommendations.csv is the random control.
39. What "model_version = aa99d51_cycle_0" tells you and why it matters for reproducibility.
40. If someone asked you to extend the pipeline to predict T7 × E. coli K-12 RBP × receptor binding, what's the minimum change?

---

## Scoring (be honest)

| Tier | Q correct | Comment |
|------|-----------|---------|
| Tier 1 (20) | __/20 | <16 → re-skim slides; <12 → cancel the talk |
| Tier 2 (10) | __/10 | <7 → re-read CONCEPT_DEEPDIVE |
| Tier 3 (5)  | __/5  | <3 → expect to lose credibility in Q&A |

**Hard rule:** if you can't get the Lee 2009 / BALD math / Boltz-2 ipTM / ExbD2 questions right, **do not give the talk** — those are the four credibility-killers.
