# InstrumentQC

Ross Turner's (University of Sydney) fork of the UMGCC FCSS **InstrumentQC** Quarto
dashboard (orig. David Rach), tracking daily QC for Cytek Aurora flow cytometers:
**3L, 5L, EVO**. Dashboards are **display-only** — processing raw instrument data
is a separate, deliberately non-render-time step (see Architecture below).

## Layout

```
_quarto.yml  styles.scss  index.qmd  404.qmd  Flag.csv     (root)
pages/       Aurora3L/5L/EVO.qmd, AuroraEVO_6PeakUltra*.qmd,
             CytekQCBead.qmd, help.qmd, _lj_helpers.R (shared LJ/Westgard helpers)
scripts/     TheScript_3L/5L/EVO.R (instrument-PC cron wrappers), EVO_4xParser.R,
             Process6PeakUltra.R, TheScript_DashboardRender.R
data/        3L/5L/EVO Archive dirs, HistoricalData.csv, reference/
             (QCBeadLot*.csv, AuroraMaintenance.csv, URCP-38-2K_MEF.csv)
images/      SydneyCytometry.jpg, 6peak_ultra_evo/ (per-run diagnostic PNGs)
docs/        rendered site — GitHub Pages serves /docs directly (no CI build step;
             rendered HTML is committed alongside source)
_attic/      archived unused leftovers — render ignores it. Includes
             8peak_2026-11-13/ (retired 8-peak pages, script, PNGs, MEF csv and
             the nested 8peak_tracking_spectral/ repo; see its README)
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
- **8-peak bead tracking** (`RCP-30-5A`) was retired on 2026-11-13 and archived to
  `_attic/8peak_2026-11-13/`.
- **6-peak Ultra Rainbow tracking** (`URCP-38-2K` lot) is the live bead tracker.
  `scripts/Process6PeakUltra.R` does the heavy FCS work once (flowClust) and
  writes `6peak Ultra_tracking_spectral/peak_metrics.csv` +
  `images/6peak_ultra_evo/` PNGs, keyed on acquisition **datetime** (not date —
  multiple runs/day happen). It is incremental. `pages/AuroraEVO_6PeakUltra*.qmd`
  are display-only. Load `dplyr` **last** in these pages — `ggcyto`/`flowCore`
  mask `filter`.

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
- **Nested git repos:** `6peak Ultra_tracking_spectral/` (and the archived
  `_attic/8peak_2026-11-13/8peak_tracking_spectral/`) are their own git repos
  and gitignored here —
  their FCS/metrics CSVs aren't in this repo, only the rendered page + diagnostic
  PNGs under `images/`.
- Raw `.fcs` and `DailyQCReport_*.CSV` are gitignored.

## Q/B sensitivity model

Two methods exist; only the 6-peak one is live (8-peak retired 2026-11-13):

- **EVO 6-peak Ultra (current method)**: `pages/_qb_parks.R`, sourced by the shared
  loader `pages/_6peak_load.R`. That loader is used by both `AuroraEVO_6PeakUltra.qmd`
  and `_Dashboard.qmd`; edit the loader, not the page setup chunk. The method is
  Parks 2017 weighted quadratic `V = c0 + c1·M + c2·M²`, fit on all non-saturated
  peaks **including the blank**, then the Parks 2018 **2+2SD limit of resolution**.
  - Headline: `lor_mfi = 4·(c1 + blank SD)` for all 64 detectors. Daily QC pins the
    MFI scale, so it tracks the MEF limit; it's background-limited and stable at
    about 5–10%. Q itself has about 30% fit SE per run.
  - MEF-unit values exist for the 18 `mef_map` detectors. Only some pass the
    log–log slope check (0.9–1.1); the rest are flagged via `mef_ok`.
  - Needs the archive's raw `sd` (blank = right-half SD) and `saturated` columns
    from `Process6PeakUltra.R`. Never use `rcv × median`.
  - Read `Sensitivity Articles/QB_2plus2SD_method_summary.md` first.
  - The previous version is archived in `_attic/6peak_ultra_v1_2026-11-11/`.
- **EVO 8-peak (old method, archived in `_attic/8peak_2026-11-13/`)**: `qb_fit()`
  (CV² vs 1/MFI, brightest-CV subtraction) plus the STN-17 closed form, documented in
  `Sensitivity Articles/QB_sensitivity_reference.md`. The method summary lists why
  these were superseded. To port it, reuse `_qb_parks.R`; the 8-peak archive needs
  the same `sd`/`saturated` columns first.

**Render locale:** run `quarto render` with `LANG=en_AU.UTF-8 LC_ALL=en_AU.UTF-8`.
Without a UTF-8 locale, R writes Unicode in plot titles and captions as `<U+2014>` etc.

## Git / deploy

- GitHub Pages serves `/docs` directly — rendered HTML is committed alongside
  source, so content-change commits include large generated diffs. This is
  expected, not a bug.
- Direct pushes to `main` are blocked by a guardrail for Claude Code — the user
  pushes manually.
