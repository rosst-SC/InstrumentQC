# InstrumentQC

Ross Turner's (University of Sydney) fork of the UMGCC FCSS **InstrumentQC** Quarto
dashboard (orig. David Rach), tracking daily QC for Cytek Aurora flow cytometers:
**3L, 5L, EVO**. Dashboards are **display-only** — processing raw instrument data
is a separate, deliberately non-render-time step (see Architecture below).

## Layout

```
_quarto.yml  styles.scss  index.qmd  404.qmd  Flag.csv     (root)
pages/       Aurora3L/5L/EVO.qmd, *_8Peak.qmd, AuroraEVO_6PeakUltra.qmd,
             CytekQCBead.qmd, help.qmd, _lj_helpers.R (shared LJ/Westgard helpers)
scripts/     TheScript_3L/5L/EVO.R (instrument-PC cron wrappers), EVO_4xParser.R,
             Process8Peak.R, Process6PeakUltra.R, TheScript_DashboardRender.R
data/        3L/5L/EVO Archive dirs, HistoricalData.csv, reference/
             (QCBeadLot*.csv, AuroraMaintenance.csv, RCP-30-5A_AS01_MEF.csv,
             URCP-38-2K_MEF.csv)
images/      SydneyCytometry.jpg, 8peak_evo/, 6peak_ultra_evo/ (per-run diagnostic PNGs)
docs/        rendered site — GitHub Pages serves /docs directly (no CI build step;
             rendered HTML is committed alongside source)
_attic/      archived unused leftovers — render ignores it
8peak_tracking_spectral/       8-peak sub-project — OWN git repo, gitignored here
6peak Ultra_tracking_spectral/ 6-peak sub-project — same pattern
```

Dropped from the original Maryland dashboard: BD CantoII/LSRII/AriaII/Fortessa
pages (stale duplicate data), the shinylive "History" widget (hardwired to
Maryland's GitHub raw URL), and the whole MFI / Before-After-QC-bead feature
(this lab doesn't run Before/After beads — gain auto-pins MFI).

## Render / run (always from the repo root)

- `quarto render` (whole site) or `quarto render pages/<file>.qmd` (one page).
- `_quarto.yml` has **`execute-dir: project`** — REQUIRED. Pages live in `pages/`;
  this makes every page run with cwd = repo root so `getwd()/data` and
  `data/reference/` resolve. A single-file render only gets project cwd if the
  file is in the `render:` whitelist in `_quarto.yml`.
- Manual `![]()` image links inside a `pages/*.qmd` file must be `../images/...`
  (Quarto resolves them relative to the page's own folder, not the project root).
- Tooling: quarto 1.8.27, R 4.5.2, Luciernaga 0.99.7 (system lib).

## Data processing architecture

Processing is split from display so renders stay reproducible:

- **Daily QC** (Gain + %rCV): `scripts/TheScript_<inst>.R` runs on the instrument
  PC via cron — it `git pull`s, parses new reports, appends to
  `data/<inst>/Archive/ArchivedData<inst>.csv` and `data/HistoricalData.csv`, then
  **auto-commits + pushes** with a PAT. Don't run this manually to backfill data.
- **SpectroFlo format split**: `Luciernaga::DailyQCParse` / `QC_FilePrep_DailyQC`
  only parse SpectroFlo **3.x** reports. EVO (and recent 3L/5L) reports are **4.x**
  and need `scripts/EVO_4xParser.R` instead — it reverse-engineers
  `QC_FilePrep_DailyQC`'s output to handle the 4.x diffs (BOM on line 1, no
  `Detector Settings` sub-header, renamed columns, `Flow Rate`+`Temperature` on
  one line) and produces the exact 305-col 5L schema. It's a re-runnable driver:
  rebuilds `ArchivedDataEVO.csv` from **all** `data/EVO/DailyQCReport_*.CSV`
  present and replaces the EVO rows in `HistoricalData.csv` (back up the archive
  first — see Gotchas). Coerce `Date` with `ymd()` before `bind_rows`-ing into
  `HistoricalData.csv` or it throws a vctrs type error (Date vs character).
- **`DailyQCParse` has no dedup** — appends via `bind_rows`; re-running it on
  already-archived days duplicates rows. `QCBeadParse` dedups by `DATE`+`TIME`
  and is safe to re-run, but also `file.remove`s the source `.fcs`.
- **8-peak bead tracking** (`RCP-30-5A` lot): `scripts/Process8Peak.R` does the
  heavy FCS work once (flowClust clustering of 8 populations) → writes
  `8peak_tracking_spectral/peak_metrics.csv` + per-run `images/8peak_evo/`
  gate/cluster PNGs. Incremental — only processes new dates / missing images.
  `pages/AuroraEVO_8Peak.qmd` is display-only (reads the archive + PNGs, no
  flowClust). Load `dplyr` **last** in these pages — `ggcyto`/`flowCore` mask
  `filter`.
- **6-peak Ultra Rainbow tracking** (`URCP-38-2K` lot): mirrors the 8-peak
  layout exactly. `scripts/Process6PeakUltra.R` writes
  `6peak Ultra_tracking_spectral/peak_metrics.csv` + `images/6peak_ultra_evo/`
  PNGs, keyed on acquisition **datetime** (not date — multiple runs/day happen).
  `pages/AuroraEVO_6PeakUltra.qmd` is the display-only dashboard.

## Gotchas

- **`VisualQCSummary` detector filter:** its default `detectorType="-A"` filters
  by the `-A` area suffix, but the Gain archive (`ArchivedData<inst>.csv`) names
  detectors *without* it (`V3`, `B3`, `FSC` — only the Bead/MFI csv uses `-A`).
  Left at default → 0-row summary → `SmallTable` crashes in `gt::data_color`.
  Pass detectors explicitly, e.g.
  `detectorType = "^FSC$|^SSC|V3$|B3$|R3$|YG3$"` (the Cytek QC-fail trigger set,
  per the Help page).
- **`index.qmd` instrument-count assumptions:** code like `Data[2:5] <- ...`
  assumes 4 instruments (a leftover from the original 4-instrument Maryland
  dashboard); with fewer instruments the frame has fewer columns and it throws
  "undefined columns selected". Use `mutate(across(-Date, ~ na_if(., "Unknown")))`
  instead of positional indexing.
- **Empty-archive corruption:** a bad run can overwrite `data/<inst>/Archive/*.csv`
  with header + all-empty rows. Symptom: `CurrentData(type="Gain")` errors with
  `missing value where TRUE/FALSE needed` (an un-`na.rm`'d `any(str_detect(...))`).
  Recover from the last good git commit.
- **Fragile per-tab plot enumerations:** hand-indexed blocks like
  `ggplotly(List[[1]])` … `[[16]]` break silently if a detector count changes
  and invite copy-paste bugs. Prefer
  `htmltools::tagList(lapply(List, ggplotly))` — robust to any list length.
- **No-data guard for a not-yet-populated instrument:** a dashboard `.qmd` can
  render cleanly before its Archive exists via (1) a placeholder
  `#| content: valuebox` with `#| eval: !expr '!HasData'`, then (2)
  `if (!HasData) knitr::knit_exit()` in an `include:false` chunk before any
  `## column {.tabset}` markdown. The home page guards similarly and builds its
  instrument list via `intersect(c("3L","5L","EVO"), unique(ShinyData$Instrument))`
  so a new instrument auto-joins once its Archive has rows.
- **Nested git repos:** `8peak_tracking_spectral/` and
  `6peak Ultra_tracking_spectral/` are their own git repos and gitignored here —
  their FCS/metrics CSVs aren't in this repo, only the rendered page + diagnostic
  PNGs under `images/`.
- Raw `.fcs` and `DailyQCReport_*.CSV` are gitignored.

## Q/B sensitivity model

Before touching `qb_fit()`, the `stn17` block, or `mef_map` in
`pages/AuroraEVO_8Peak.qmd` / `pages/AuroraEVO_6PeakUltra.qmd`, read
`Sensitivity Articles/QB_sensitivity_reference.md` — it ties every modeling
choice to a source (Hoffman & Wood 2007, Chase 1998, STN-17, BD CS&T Qr/Br).
Key point: the model is `CV² = CV_intr² + (1/Q)(1/S) + (B/Q)(1/S²)` fit as
**CV² vs 1/MFI** on gain-normalised intensity — fitting variance-vs-mean instead
gives a wrong (negative) B.

## Git / deploy

- GitHub Pages serves `/docs` directly — rendered HTML is committed alongside
  source, so content-change commits include large generated diffs. This is
  expected, not a bug.
- Direct pushes to `main` are blocked by a guardrail for Claude Code — the user
  pushes manually.
