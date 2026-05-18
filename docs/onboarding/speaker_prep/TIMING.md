# TIMING.md — Minute-by-minute target

**Target total:** 75 minutes of talk + 15 minutes Q&A = 90 minutes door-to-door.
**Hard floor:** 60 minutes if you're rushing (skip OPTIONAL sections).
**Demo (15 min):** independent — run after Part 6 OR at the end during Q&A buffer.

---

## Per-part allocation

| Part | Min | Cum | Tag | Audience priority |
|------|-----|-----|-----|---------|
| 0. Roadmap | 4 | 4 | [CORE] | All |
| 1. Science | 10 | 14 | [CORE] | WET LAB + PI ★★★ |
| 2. Pipeline architecture | 12 | 26 | [CORE except Modules 04 ALT path] | DRY LAB + PI |
| 3. ML core | 15 | 41 | [CORE] | DRY LAB + PI ★★★ |
| 4. Boltz-2 result | 8 | 49 | [CORE] | All |
| 5. 48-h cycle | 12 | 61 | [CORE] | WET LAB + DRY LAB ★★★ |
| 6. Reproducing | 6 | 67 | [OPTIONAL — skip if running long] | DRY LAB + WET LAB |
| 7. Risks & decisions | 6 | 73 | [CORE for PI] | PI ★★★ |
| 8. References & thank-you | 2 | 75 | [CORE thank-you only] | All |
| Q&A buffer | 15 | 90 | — | All |
| Demo (optional) | 15 | +15 | [OPTIONAL] | DRY LAB + WET LAB |

---

## Checkpoints — adjust if behind

### Minute 15 — should be at end of Part 1
- **On track:** ending Part 1, transitioning to architecture.
- **5 min behind:** skip "Suggested cuts" + "How to read"; compress vocab block; rejoin at Part 1.
- **>10 min behind:** cut "data-scarcity bottleneck" slide (it's redundant with intro); rejoin at "Why active learning".

### Minute 30 — should be mid-Part 2 (Module 03 / 04)
- **On track:** ending Module 03 RBP slide, starting embeddings.
- **5 min behind:** skip "Module 02 — Annotation" and "Module 04 PLM-interact paragraph"; rejoin at Module 05 Boltz-2.
- **>10 min behind:** compress all of Module 00–04 to one minute total; jump to Module 05.

### Minute 45 — should be mid-Part 3 (BALD math)
- **On track:** at BALD math slide.
- **5 min behind:** skip "ensemble.py snippet" + "calibration plot"; rejoin at BALD math.
- **>10 min behind:** lean entirely on the "BALD intuition" picture; skip math derivation; rejoin at bald.py snippet.

### Minute 60 — should be mid-Part 5 (handoff)
- **On track:** at "Wet lab → dry lab handoff" or "Dry lab → wet lab handoff".
- **5 min behind:** compress receptor knockout slide; rejoin at validation tiers.
- **>10 min behind:** skip cloning execution + ELISA protocol slides (they're in the guide); jump to validation tiers + quality gates.

### Minute 75 — should be at thank-you
- **On track:** wrapping up Part 7 / 8, opening Q&A.
- **Over:** start Q&A immediately. Demo runs only if asked.

---

## Audience-priority annotations (which to compress for which audience)

**Heavy wet-lab attendance** (more wet than dry/PI):
- Compress: Part 3 BALD math + ensemble code snippets (cut to "5 networks vote, picks where they disagree").
- Expand: Part 5 handoff details, ELISA protocol, knockout strategy.
- Time: shift 5 min from Part 3 → Part 5.

**Heavy dry-lab/PI attendance**:
- Compress: Part 5 cloning execution + ELISA protocol (one slide overview).
- Expand: Part 3 BALD math + ALDE caveat + Greenman 2025 caveat.
- Time: shift 5 min from Part 5 → Part 3.

**PI-only briefing (45 min)**:
- Cut Part 6 entirely.
- Compress Part 5 to "two handoffs in 48h SLA".
- Lead with Part 7 risk + decisions.
- Order: 1, 5 (Hung 2003 receptor), 3 (briefly), 4 (Boltz-2 ipTM), 5 (48h cycle + Tier 3), 7 (decisions).

---

## Cuts you can make WITHOUT losing the talk

These slides are skippable in any priority mode:
- "How to read this deck" (Part 0)
- "Suggested cuts" (Part 0)
- "Notebook-first workflow + bilingual comments" (Part 2)
- "Calibration plot — Cycle 0" (Part 3) — caveat: needs to be mentioned verbally
- "ensemble.py snippet" (Part 3) — visual aid only
- "Where the Boltz-2 outputs live" (Part 4) — point at PI briefing instead
- "Quality gates and failure modes" (Part 5) — wave at it
- "Per-module entry points" (Part 6) — in `GETTING_STARTED.md`
- All 4 citation appendix slides (Part 8) — point at `papers.md`

If you cut all of the above, you save ~10 minutes.

---

## DO NOT skip — even if you're running late

- Module 03 (rbp_01 discovery) — your key dry-lab result
- Module 05 Boltz-2 (ipTM caveat) — biggest credibility risk slide
- BALD math + ALDE caveat — your methodological contribution + accuracy guardrail
- Boltz-2 three numbers slide (ipTM 0.365) — your live result
- Receptor knockouts (Tier 3 + ΔexbD2 negative control) — PI-relevant
- Risks + decisions slides — what the PI is here to decide

---

## If Q&A goes long

- Demo can be deferred to after the talk (small group).
- Stay on questions until 95-min mark, then announce: "I'll stay for as long as anyone wants — let's officially close so people can leave."
- Don't end abruptly mid-question. End on the natural beat.

---

## ZH version — TIMING_zh.md mirror

The same plan in Chinese is in `TIMING_zh.md`. Stick to one language for delivery; switch on the fly for ZH speakers in Q&A.
