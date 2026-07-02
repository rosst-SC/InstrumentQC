# InstrumentQC — session handoff (2026-06-15)

Context for resuming work in a new Claude Code session. (This file is a scratch
note — delete it whenever. It is not in the Quarto render whitelist.)

Two memory files auto-load each session and hold the durable facts — read them
first: `instrumentqc-project.md` and `instrumentqc-render-gotchas.md`.

---

## What this is
Ross Turner's (Univ. of Sydney) fork of the UMGCC **InstrumentQC** Quarto
dashboard, tracking daily Cytek Aurora QC. **Live project = the repo root**
`/Users/r.turner/Documents/Positron_Local/InstrumentQC/` (moved up out of the old
`My Modified Scripts and Data/` subfolder; that folder no longer exists).

Instruments shown: **Aurora 3L, 5L, EVO**. MFI/before-after-bead feature was
removed (gain auto-pins MFI). Dashboards are **display-only**; processing is
separate.

## Layout (after reorg)
```
_quarto.yml  styles.scss  index.qmd  404.qmd  Flag.csv   (root)
pages/     Aurora3L/5L/EVO.qmd, *_8Peak.qmd, CytekQCBead.qmd, help.qmd
scripts/   TheScript_3L/5L/EVO.R, EVO_4xParser.R, Process8Peak.R, TheScript_DashboardRender.R
data/      3L/5L/EVO Archives, HistoricalData.csv, reference/ (QCBeadLot*.csv, AuroraMaintenance.csv, RCP-30-5A_AS01_MEF.csv)
images/    SydneyCytometry.jpg, 8peak_evo/ (per-run gate + cluster diagnostic PNGs)
docs/      rendered site (GitHub Pages serves /docs); index.html + 404.html at root, rest under docs/pages/
_attic/    archived unused leftovers (Other qmd, Miscellaneous, staff, etc.) — render ignores it
8peak_tracking_spectral/   8-peak sub-project — OWN git repo, gitignored here (FCS + peak_metrics.csv + AS01 xlsx live here)
```

## Render / run (always from the repo root)
- `quarto render` (whole site) or `quarto render pages/<file>.qmd` (one page).
- `_quarto.yml` has **`execute-dir: project`** (REQUIRED — pages live in `pages/`,
  so each page runs with cwd = repo root; `getwd()/data` + `data/reference/` resolve).
- Single-file render only uses project cwd if the file is in the `render:` whitelist.
- Daily-QC processing (instrument PC): `scripts/TheScript_<inst>.R`. EVO daily-QC
  reports are SpectroFlo 4.x → parsed by `scripts/EVO_4xParser.R` (3L/5L are mixed
  3.x/4.x, handled by a format-dispatch in the rebuild — see gotchas memory).
- 8-peak processing: `Rscript scripts/Process8Peak.R` (incremental — only new
  FCS / missing diagnostic images).

## Current focus: 8-peak bead tracking (EVO) — JUST UPGRADED, UNCOMMITTED
Architecture (mirrors the QC dashboards): `scripts/Process8Peak.R` does the heavy
FCS work ONCE → writes `8peak_tracking_spectral/peak_metrics.csv` + per-run
`images/8peak_evo/{gate,cluster2d}_<date>.png`. `pages/AuroraEVO_8Peak.qmd` is
**display-only** (reads the archive + the daily-QC gain). Only 3 runs so far
(2026-05-20/27/29).

The upgrade just completed (was a flat list of trend plots; now an executive
summary + tabbed/interactive layout):
- **Executive summary** heatmap: per-detector Δ-vs-baseline, sign-flipped so red=worse.
- **Sensitivity (Q & B)**: CV² vs 1/MFI fit on **gain-normalised** intensity →
  Q (detection efficiency), B, intrinsic CV. (NB: variance-vs-mean is WRONG — gives
  negative B; CV²-vs-1/MFI is correct.)
- **Gain context**: joins `data/EVO/Archive/ArchivedDataEVO.csv` by date+detector
  (`UV3-A`↔`UV3-Gain`); MFI/gain de-pins the daily-QC gain loop.
- **Linearity (MEF, lot AS01)**: log–log measured MFI vs assigned MEF per anchor
  detector + **MEF detection limit** (gain-independent). MEF data in
  `data/reference/RCP-30-5A_AS01_MEF.csv`. **Fluor→detector map (user-confirmed):**
  CascadeBlue→V1, BFP→V3, FITC→B2, PE→YG1, EDC→YG3, PE-Cy5→YG5, PE-Cy7→YG9,
  APC→R1, APC-Cy7→R7 (these are hard-coded in the qmd's `mef_map`).
- **Spectral signature drift**: cosine similarity of each peak's 64-detector
  signature vs the baseline run.
- **%CV & noise floor**: peak-8 uses %CV; **blank (peak 1) uses absolute robust
  spread, NOT %CV** (%CV explodes when the blank mean ≈ 0 — was a false UV8 spike).
- Control limits ±2 SD wired in (preliminary banner while n<5). Interactive plotly,
  heatmap-centric. dplyr loaded LAST (ggcyto/flowCore mask `filter`).
- Pending: nothing required. (Could later add Q in true MEF units; control limits
  firm up as runs accumulate.)

## Uncommitted right now (git: branch main, ahead of origin by 1, plus this)
- **M** `pages/AuroraEVO_8Peak.qmd` (the upgrade), `index.qmd` (review — likely the
  reorg Help-link tweak), and the re-rendered `docs/index.html`,
  `docs/pages/AuroraEVO_8Peak.html`, `docs/search.json`.
- **??** `data/reference/RCP-30-5A_AS01_MEF.csv` (keep — durable AS01 values),
  `pages/AuroraEVO_8Peak_files/` (render figures).
- **?? comparison artifacts — do NOT commit, delete when done:**
  `pages/AuroraEVO_8Peak_previous.qmd`, `docs/pages/AuroraEVO_8Peak_previous.html`,
  `docs/pages/AuroraEVO_8Peak_previous_files/` (a copy of the pre-upgrade page for
  side-by-side comparison; the previous version is also recoverable from git
  `HEAD:pages/AuroraEVO_8Peak.qmd`).

## Next steps / open items
1. **Commit the 8-peak upgrade** (`pages/AuroraEVO_8Peak.qmd` + `data/reference/RCP-30-5A_AS01_MEF.csv`
   + re-rendered docs), **excluding** the `*_previous*` comparison files.
2. **Push** when ready: `git push origin main` (NOTE: a Claude Code guardrail blocks
   direct pushes to `main` — the user pushes themselves).
3. Pre-existing loose ends (untouched): `README.md` modified, `My Working Index.qmd`
   deleted — both pending the user's call. Help-page screenshots are the original
   Maryland ones (could be replaced with screenshots of this dashboard).
4. Unanswered question from the user: whether their claude.ai-account
   "Instructions for Claude" apply in Claude Code (they generally do NOT — Code uses
   CLAUDE.md + settings, not the web-app personal instructions).

## Gotchas (also in render-gotchas memory)
- `execute-dir: project` is mandatory now that pages are in `pages/`.
- Manual `![]()` image links in a page under `pages/` must be `../images/...`
  (Quarto resolves them relative to the page's folder).
- `8peak_tracking_spectral/` is a nested git repo → gitignored; its FCS/metrics
  aren't in the main repo (only the rendered page + diagnostic PNGs are).
- `DailyQCParse`/`QC_FilePrep_DailyQC` only read SpectroFlo 3.x; EVO (and recent
  3L/5L) are 4.x → `EVO_4xParser.R`. 4.x reports drop `Laser Power` (optional).
- Raw `.fcs` and `DailyQCReport_*.CSV` are gitignored.
