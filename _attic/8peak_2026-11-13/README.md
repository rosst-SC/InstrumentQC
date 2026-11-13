# 8-peak Rainbow (RCP-30-5A) tracking — archived 2026-11-13

Retired from the live site. URCP-38-2K 6-peak Ultra tracking
(`pages/AuroraEVO_6PeakUltra*.qmd`) replaces it. Paths below mirror where each
file used to live.

| Here | Was |
|---|---|
| `pages/AuroraEVO_8Peak.qmd` | EVO 8-peak dashboard (display-only; `qb_fit()` + STN-17 Q/B) |
| `pages/Aurora3L_8Peak.qmd`, `pages/Aurora5L_8Peak.qmd` | placeholder pages (never populated) |
| `scripts/Process8Peak.R` | flowClust processing → `peak_metrics.csv` + diagnostic PNGs |
| `images/8peak_evo/` | per-run gate / cluster PNGs (7 runs, 2026-05-20 → 2026-07-01) |
| `data/reference/RCP-30-5A_AS01_MEF.csv` | lot AS01 MEF assignments |
| `docs/AuroraEVO_8Peak.html` | last rendered page (a snapshot — its image links no longer resolve) |
| `8peak_tracking_spectral/` | nested sub-project git repo (FCS, `peak_metrics.csv`) — still gitignored |

The method is documented in `Sensitivity Articles/QB_sensitivity_reference.md`.
`Sensitivity Articles/QB_2plus2SD_method_summary.md` explains why it was superseded.

To restore it, move the files back to the paths above. Then re-add the three
pages to the `render:` list and the "8 Peak Beads" navbar menu in `_quarto.yml`.
