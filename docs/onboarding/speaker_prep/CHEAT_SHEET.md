# CHEAT_SHEET.md — One-page safety net

**Print this. Have it on the lectern. Don't read from it — but knowing it's there will keep you calm.**

---

## Key numbers — memorize cold

| Item | Value |
|------|-------|
| Phage genomes (Module 00) | **777** |
| Bacteria genomes (Module 00) | **34** |
| Interaction matrix pairs | **2,236** (315 + / 1,921 −) |
| rbp_01 length | **712 aa** |
| rbp_01 HMM bit score (Tail_spike_N) | **342** |
| RBP candidates for phiL7 | **3** |
| ESM-2 8M / 650M / 3B embed dim | **320 / 1280 / 2560** |
| Boltz-2 ipTM (rbp_01 × TonB, job 59986) | **0.365** |
| Boltz-2 chain_A_ptm | **0.808** |
| Boltz-2 confidence_score | **0.683** |
| Ensemble members | **5** |
| Tests pass: M01 / M02 / M03 / M04 / M05 / M06 / M07 | 22 / 26 / 25+ / 17 / 28 / 9 / **18** |
| Cycle picks | **4 BALD + 1 random** |
| Dry-lab SLA | **48 h** |
| Wet-lab cycle | **10–14 days** |
| SDM cost / time | **$50 / 4 days** |
| Gene synthesis cost / time | **$150 / 2 weeks** |

---

## Citations — shortform

| Tag | Claim |
|-----|-------|
| **Lee 2009** *AEM* 75:7828 | phiL7 genome; explicitly searched for OP1 ORF25 homolog, found none |
| **Hung 2003** *BBRC* 302:878 PMID 12646254 | TonB/ExbB/ExbD1 essential; **ExbD2 NOT** required |
| **Boeckaerts 2022** *Viruses* 14:1329 | PhageRBPdetect HMM + XGBoost; PR AUC 93.8% |
| **Boeckaerts 2024** *Nat Commun* 15:4768 | PhageHostLearn AUC up to 0.82 (within genus) |
| **Mutalik 2025** *bioRxiv* | PAML benchmark; cross-genus AUC 0.67–0.70 |
| **Lin 2023** *Science* 379:1123 | ESM-2; 65M proteins, masked LM |
| **Liu 2025** *Nat Commun* 16:64512 | PLM-interact; +16–28% AUPR cross-species PPI |
| **Abramson 2024** *Nature* 630:493 | AlphaFold 3 |
| **Passaro 2025** *bioRxiv* | Boltz-2; affinity head SMALL-MOLECULE ONLY |
| **Lakshminarayanan 2017** *NeurIPS* arXiv:1612.01474 | Deep ensembles + Gaussian NLL |
| **Ovadia 2019** *NeurIPS* | UQ under shift; ensemble > MC Dropout calibration |
| **Greenman 2025** *PLoS Comp Biol* 21:e1012639 | "No single best UQ method"; NOT *NAR Genomics* |
| **Houlsby 2011** arXiv:1112.5745 | BALD; originally for GPC classification |
| **Hie 2024** *Nat Biotechnol* 42:275 | Antibody affinity maturation; **ESM-1b/1v** (not ESM-2); ~20 variants |
| **Yang 2025** *Nat Commun* 16:55987 | ALDE; **Thompson sampling** + one-hot (NOT BALD/ESM-2) |
| **McNair 2019** *Bioinformatics* 35:4537 | PHANOTATE (phage ORFs, dynamic programming) |
| **Hyatt 2010** *BMC Bioinf* 11:119 | Prodigal (bacterial ORFs) — we use pyrodigal binding |
| **Schäfer 1994** *Gene* 145:69 | pK18mobsacB (suicide vector + sacB counter-selection) |
| **Farquharson 2021** | T4 × E. coli: binding ≠ infection (85% bind, 11% plaque) |

---

## ⚠️ Three sentences you must say verbatim

> **1.** "Boltz-2's affinity head is small-molecule only — it outputs NaN for protein-protein pairs. We use ipTM as a structural confidence proxy, not as binding affinity."

> **2.** "BALD for our deep-ensemble regression is `epistemic_std = Std_k[μ_k]` — the standard deviation of the ensemble member means. Houlsby 2011 derived it for GPC classification; the deep-ensemble regression form is our extension."

> **3.** "Lee 2009 used BLAST and explicitly said they couldn't find an OP1 ORF25 homolog. We used a Hidden Markov Model — a more sensitive tool — and found rbp_01. The two results are complementary, not contradictory."

---

## ⚠️ Three things you must NEVER say

> ❌ "Boltz-2 gives us zero-shot binding affinity."
> ❌ "ALDE (Yang 2025) validated BALD." (They used Thompson sampling.)
> ❌ "Lee 2009 missed it" or "Hie 2024 used ESM-2." (ESM-1b/1v; and Lee searched and reported negative result.)

---

## If you freeze — three recovery lines

> 🎯 **"Let me pull up the actual code…"** — then open `bald.py` or `ensemble.py`. The code itself answers most ML questions. Reading from your own repo looks competent, not panicked.

> 🎯 **"That's an open question — let me note it down."** — pull out a pen, write it on paper visibly. Don't bluff. Honesty + a note is better than a wrong answer.

> 🎯 **"Short answer is X. Long answer is in `docs/onboarding/guide_en.md` section Y — I can walk through it after."** — keeps momentum, defers without ducking.

---

## Audience artifacts — paths to share

- Slides (EN): `docs/onboarding/slides_en.pdf`
- Slides (ZH): `docs/onboarding/slides_zh.pdf`
- Companion guide (EN): `docs/onboarding/guide_en.md`
- Companion guide (ZH): `docs/onboarding/guide_zh.md`
- Demo runbook: `docs/onboarding/DEMO.md`
- Paper audit: `docs/reference/paper_reading_notes.md`
- Project plan (EN): `docs/planning/iGEM_2026_Project_Plan.md`
- PI briefing snapshot: `docs/planning/PI_briefing_2026-05-11.md`

---

## Live-demo emergency: minimum viable demo

If everything breaks:
1. `cat 07_acquisition_function/outputs/cycle_1/recommendations.csv` — show the actual recommendations file.
2. `cat 06_uncertainty_model/outputs/cycle_0/predictions.csv | head` — show predictions with epistemic_std.
3. Open `05_structure_prediction/outputs/boltz2/EU717894.1_rbp_01__GCF_000007145.1_tonB/affinity.json` — show ipTM 0.365.

These three files prove the pipeline works. Nothing has to run live.

---

## Bilingual switch — short cues for ZH transition

- "下面用中文复述一下…" / "Let me say that in Chinese…"
- "回到英文 / Back to English."
- "用 Sarah 听得最舒服的语言…" / "In the language Sarah hears most easily…"

You're allowed to mid-sentence switch — the audience is bilingual.
