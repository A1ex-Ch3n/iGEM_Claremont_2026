# `docs/` — Project Documentation

Navigation index for the iGEM Claremont 2026 active-learning phage engineering project.

---

## `planning/` — Current planning & design docs

| File | Audience | Purpose |
|------|----------|---------|
| [`iGEM_2026_Project_Plan.md`](planning/iGEM_2026_Project_Plan.md) | PI · external reviewers | Full English plan (~5500 words). Background, system architecture, all dry-lab + wet-lab modules with literature backing, AL cycle, risks, timeline, deliverables. |
| [`iGEM_2026_Project_Plan_PI_Brief.md`](planning/iGEM_2026_Project_Plan_PI_Brief.md) | PI quick read | Abridged English brief covering only Sections 3 (dry lab), 4 (wet lab), 5 (AL cycle), 9 (resources), 11 (PI questions). With status indicators. |
| [`PI_briefing_2026-05-11.md`](planning/PI_briefing_2026-05-11.md) | PI · team | **Current status** (updated 2026-05-12): all modules 00–07 done, Boltz-2 results, BALD first run, validation strategy (4 tiers), PI decisions, dry-wet interface. Bilingual EN + ZH. |
| [`iGEM_2026_项目大纲_中文版.md`](planning/iGEM_2026_项目大纲_中文版.md) | Team (Chinese-preferring) | Chinese plan, deepest version. Same scope as PI plan + extended dry-lab module mechanics + dedicated **6-layer data-scarcity strategy** chapter. |
| [`project_pivot_summary_for_team.md`](planning/project_pivot_summary_for_team.md) | Team sync | Chinese narrative explaining the May 2026 pivot from the original 6-factor approach to active learning. Includes Path A/B/C strain-acquisition decision tree. |
| [`iGEM_2026_Project_Flow.excalidraw`](planning/iGEM_2026_Project_Flow.excalidraw) | All | Editable flow chart source (open at https://excalidraw.com). |
| [`iGEM_2026_Project_Flow.png`](planning/iGEM_2026_Project_Flow.png) | Wiki / slides | Rendered flow chart. |

## `reference/` — Lookups & technical references

| File | Purpose |
|------|---------|
| [`glossary.md`](reference/glossary.md) | All technical vocabulary used across docs (ELISA, BALD, ESM-2, MALDI-TOF, etc.) — alphabetical index + 12 categorized sections. |
| [`prodigal_manual.md`](reference/prodigal_manual.md) | Bilingual EN/ZH guide to Prodigal (bacterial gene caller used in Module 02). |
| [`prodigal_manual.pdf`](reference/prodigal_manual.pdf) | PDF version. |

## `protocols/` — Wet lab Benchling protocols

Sourced from team Benchling. Renamed from the original `· Benchling.pdf` filenames for shell-friendliness.

| File | Used in module |
|------|---------------|
| [`Xanthomonas_Cultivation.pdf`](protocols/Xanthomonas_Cultivation.pdf) | 4.1 strain isolation, 4.5 ELISA cell prep, 4.6 receptor KO host growth |
| [`Xanthomonas_Transformation.pdf`](protocols/Xanthomonas_Transformation.pdf) | 4.6 receptor KO (electroporation of pK18mobsacB) |
| [`Phage_Plaque_Assay_SmallScale.pdf`](protocols/Phage_Plaque_Assay_SmallScale.pdf) | 4.2 phage isolation, 4.7 plaque assay validation |
| [`Phage_Liquid_Amplification.pdf`](protocols/Phage_Liquid_Amplification.pdf) | 4.2 phage stock preparation |
| [`Phage_Infection_Curves.pdf`](protocols/Phage_Infection_Curves.pdf) | ⚠️ Current version measures bacterial density only — needs rewrite to add phage challenge step (see Project Plan §8.2 Flag 1). |

## `archive/` — Pre-pivot artifacts

Kept for reference. Describe the pre-May-2026 6-factor pipeline; superseded by the active-learning plan in `planning/`.

| File | What it was |
|------|-------------|
| [`integrated_pipeline_guide.md`](archive/integrated_pipeline_guide.md) | Bilingual walkthrough of the original 7-step pipeline (Step 3 = 6-factor weighting). |
| [`integrated_pipeline_guide.pdf`](archive/integrated_pipeline_guide.pdf) | PDF version. |
| [`workflow_chart.jpeg`](archive/workflow_chart.jpeg) | Original 7-step workflow diagram referenced in the pre-pivot CLAUDE.md. |

---

## Adding new docs

- **Planning / design decisions** → `planning/`
- **Lookup tables / tool manuals** → `reference/`
- **Wet lab SOPs from Benchling** → `protocols/`
- **Anything no longer current** → `archive/` (don't delete; keep for history)
- Update this README's tables when you add anything.
