# Onboarding Package

Author: Claude (onboarding agent, 2026-05-17). Target audience: wet lab + dry lab + PI / faculty advisor.

## Contents

| File | Purpose | Language |
|------|---------|----------|
| `slides_en.md` + `slides_en.pdf` | Deep-dive deck (~60 min), tagged `[WET LAB] / [DRY LAB] / [PI]` per section | English |
| `slides_zh.md` + `slides_zh.pdf` | Same deck, Simplified Chinese | 简体中文 |
| `guide_en.md`  | Self-contained written guide with reproduction appendix | English |
| `guide_zh.md`  | Same guide, Simplified Chinese | 简体中文 |
| `DEMO.md`      | Live-demo runbook (Modules 03 → 04 → 06 → 07) | English |
| `PLAN.md`      | The agent's pre-build plan + format rationale | English |
| `figures/make_figures.py`           | Regenerates the four diagrams below | — |
| `figures/pipeline_flow.png`         | Modules 00 → 08 with data flow | — |
| `figures/data_contract.png`         | inputs / processes / outputs convention | — |
| `figures/cycle_48h.png`             | 48-hour wet ↔ dry handoff | — |
| `figures/bald_intuition.png`        | BALD selection intuition (toy plot) | — |

## Rebuild

```bash
# From repo root
python docs/onboarding/figures/make_figures.py

pandoc docs/onboarding/slides_en.md -t beamer --pdf-engine=xelatex --slide-level=2 \
  -o docs/onboarding/slides_en.pdf
pandoc docs/onboarding/slides_zh.md -t beamer --pdf-engine=xelatex --slide-level=2 \
  -o docs/onboarding/slides_zh.pdf
```

One-time TeX dependencies (TinyTeX / TeX Live):

```bash
tlmgr install beamer booktabs mdframed xecjk ctex pgfplots fontaxes needspace zref
```

Mac CJK font: `Hiragino Sans GB` (built-in).

## Suggested presenter cuts

- **30-min wet-lab onboarding** — slides 1–4, 7–12, 17–22, 37–44.
- **45-min PI / advisor briefing** — slides 1–4, 5–12, 13–14, 23, 28–32, 33–36, 51–55, 56.
- **90-min full team deep dive** — all slides.

The demo (≈15 min) is independent — see `DEMO.md`.
