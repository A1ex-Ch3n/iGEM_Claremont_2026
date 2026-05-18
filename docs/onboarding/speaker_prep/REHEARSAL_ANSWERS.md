# REHEARSAL_ANSWERS.md — Answer key

Open AFTER attempting REHEARSAL.md. Self-grade.

---

## Tier 1 — Numerical / mapping

1. **777** phage genomes.
2. **34** bacterial genomes.
3. **2,236** pairs total = **315 positive + 1,921 negative + 1 ground-truth**.
4. **712** amino acids (rbp_01 = orf_00001).
5. HMM bit score **342** against `Tail_spike_N` (`unknown_C54`).
6. **3** RBP candidates for phiL7 (rbp_01: 712 aa @ 342; rbp_02: 918 aa @ 235.1; rbp_03: 224 aa @ 56.7).
7. ESM-2 8M → **320-D**.
8. ESM-2 650M → **1280-D**.
9. ESM-2 3B → **2560-D**.
10. ipTM = **0.365**.
11. chain_A_ptm = **0.808**. (Confidence score = 0.683.)
12. **5** members (`n_members = 5` default, Lakshminarayanan 2017 recipe).
13. **18/18** tests passing in Module 07. (Module 06: 9/9.)
14. **n_bald=4, n_random=1** → 5 picks per cycle.
15. **48-hour** SLA (dry lab cycle).
16. **PHANOTATE** for phage ORFs (McNair 2019); **pyrodigal** (Python binding of Prodigal, Hyatt 2010) for bacteria. Never swap.
17. **PhageRBPdetect** (Boeckaerts 2022). HMM track + XGBoost ML track. We use HMM.
18. **AlphaFold 3** (Abramson 2024) + **Boltz-2** (Passaro 2025).
19. **Hung 2003 *BBRC* 302:878–884 PMID 12646254.** Essential: **TonB, ExbB, ExbD1**. **ExbD2 NOT required** (ΔexbD2 retains infection).
20. **Houlsby 2011 arXiv:1112.5745**. Originally applied to **Gaussian Process Classifier (GPC)**. We extend to deep-ensemble regression.

---

## Tier 2 — 1-minute explanations

(Full answers in `CONCEPT_DEEPDIVE.md`; quick check below.)

21. **AL solves data scarcity** — passive ML wastes measurements where the model is already confident. AL picks where the model is most uncertain → each ELISA reduces uncertainty maximally. 2–5× speedup in adjacent fields (Hie 2024 antibodies, Yang 2025 enzymes).

22. **ESM-2 embedding** — vector summary of a protein from a transformer trained to fill in masked amino acids on 65M proteins. The internal representation (1280-D for 650M model) captures structural + functional context. For wet lab: "It's the model's compressed understanding of what each protein looks like."

23. **ipTM = interface predicted TM-score**, 0–1, **confidence** in the docking geometry. NOT affinity. Boltz-2 affinity head is small-molecule only → NaN for protein-protein. ipTM 0.365 = "model uncertain how they dock" not "they don't bind".

24. **Deep ensemble** — 5 networks, same architecture, different random seeds → trained models agree where data tells them what to do, disagree elsewhere. Spread of predictions = epistemic uncertainty (what model doesn't know, shrinks with data).

25. **BALD ≈ Var_k[μ_k]** — Houlsby's mutual-information objective `H[y|x,D] - E_θ[H[y|x,θ]]`, under Gaussian-ensemble homoscedasticity, reduces to `½ log(1 + Var_k[μ_k]/σ²)` which is monotone in `Var_k[μ_k]`. We use the std for unit consistency.

26. **48-h SLA** — wet lab cycles 10–14 days, dry lab 48h. Sliding either causes (a) wet lab idle (wasted week) or (b) stale recommendations (wasted information). Safe-pick backup is the parachute. Three slipped cycles = wiki freeze risk.

27. **Lee 2009 + HMM** — Lee 2009 BLAST-searched OP1 ORF25 homolog in phiL7, found none. We HMM-searched Tail_spike_N family profile, found rbp_01. HMM is more sensitive than BLAST at low sequence identity. Not contradictory: different tools, different sensitivities.

28. **ΔexbD2** — Hung 2003 showed ExbD2 not essential; if we build ΔexbD2 alongside ΔtonB/ΔexbB/ΔexbD1, it should **retain** infection. That validates our knockout system for free — if ΔexbD2 doesn't infect, something's broken upstream.

29. **Epistemic** = reducible model uncertainty (shrinks with more data); **aleatoric** = irreducible measurement noise (ELISA pipetting). Total variance = epistemic + aleatoric. BALD targets epistemic only.

30. **Random control arm** — 1 in 5 picks sampled randomly, wet-lab blinded; lets us prove at project end that BALD actually beat random selection (Hie 2024 standard).

---

## Tier 3 — Defending design choices

31. **BALD over Thompson sampling**: BALD is information-theoretic (maximize uncertainty reduction); Thompson is exploitation-biased (sample from posterior, act). Early cycles benefit from exploration; BALD is the canonical exploration acquisition. ALDE (Yang 2025) used Thompson because they were further along in directed evolution where exploitation matters more. Different problem, different choice.

32. **ESM-2 over ProtBERT/etc.**: ESM-2 trained on 65M proteins (UniRef50/90, masked LM), 33-layer transformer at 650M; outperforms ProtBERT on most benchmarks (Lin 2023, *Science*). Plus PLM-interact (Liu 2025) provides a clear path to fine-tune on PPIs. ESM family has the strongest community + tooling. Hand-crafted features lose to learned embeddings on essentially all protein-ML benchmarks since 2020.

33. **Boltz-2 over AF3 only**: AF3 also predicts complexes but with weights-access gating and slower runs. Boltz-2 is open-weights (Github jwohlwend/boltz), gives the same ipTM + PAE family of outputs, runs in ~3 min on L40S. We use BOTH — AF3 for higher-quality static structures, Boltz-2 for iterating on RBP variants where speed + open weights matter.

34. **ipTM over pLDDT**: pLDDT is **per-residue local** confidence (backbone position). ipTM is **interface** confidence (how chains dock). For an RBP × receptor pair, we care about the interface, not the residue-local geometry of each monomer. The monomer pLDDT (mean 0.76) tells us each chain is well-predicted; ipTM tells us about the docking we actually care about.

35. **Deep ensemble over MC Dropout / BNN**: Calibration. Ovadia 2019 showed MC Dropout is poorly calibrated under shift; full BNNs are computationally hard at our scale. Deep ensembles are simple, scale to ESM-2's 1280-D inputs, and have published precedent in protein ML (ALDE). Greenman 2025 warns no single method is universally best but ensembles are a defensible default with temperature-scaling backstop.

---

## Stretch answers

36. **ECE (Expected Calibration Error)** — bin predicted probabilities, check actual frequency in each bin, average the gap. ECE > 0.1 = overconfident. **Temperature scaling** = divide logits by a learned scalar `T`; preserves ranking, recalibrates confidence. Single-parameter fix, applied post-hoc.

37. **`safe_pick_backup.csv`** = top-10 BALD picks (sorted by BALD score) used **only if** dry lab misses 48h SLA. Pre-vetted by PI/dry lab so wet lab can immediately proceed without waiting for a re-run.

38. **Blinding** — preserves the AL-vs-random comparison's statistical validity at project end. If wet lab knew which row was random, they might subconsciously do better work on BALD picks (Hawthorne effect), or worse on the random (skip controls).

39. **model_version = git_sha_cycle_N** — every prediction is reproducible from that string: `git checkout <sha>` + load `ensemble_member_*.pt` + same `seed=42` → identical predictions. Critical for debugging "why did Cycle 3 pick this variant?".

40. **Extend to T7 × E. coli K-12**: minimum change is the inputs to Module 01 (add the T7 × *E. coli* K-12 ground-truth pair) and Module 04 (embed T7 RBP + *E. coli* receptors). Modules 02 (PHANOTATE), 03 (PhageRBPdetect), 06 (ensemble), 07 (BALD) work unchanged — that's the entire point of the data-contract design. Cycle 0 retrains, picks new variants. No code changes.

---

## ZH version inline

For each Q, the EN answer above is the source of truth. ZH 翻译可参考 SPEAKER_NOTES_zh 和 CONCEPT_DEEPDIVE_zh.md。Alex 可在脑中双语对照,但**数字、引用、人名一律以英文版为准** ——这些不翻译,直接原文用。
