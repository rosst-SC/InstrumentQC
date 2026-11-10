# 6-peak Ultra Rainbow page — v1 snapshot (2026-11-11)

This is the EVO URCP-38-2K page as it was before the Q/B/2+2SD method from the
Sensitivity project (Parks 2017 / 2018) replaced it.
It covers 27 runs, from 2026-06-26 to 2026-11-10.

## v1 method

These files used the old sensitivity method:
- `qb_fit()`: fits CV² vs 1/x on the resolved fluorescent peaks (peak > 1) only.
  - It subtracts the brightest peak's CV as "intrinsic" and uses 1/CV⁴ weights.
  - Q = 1/slope, on gain-normalised MFI or on the URCP MEF values.
- STN-17 closed form: Q_stn and B_stn, from the blank, the dim peak and the brightest peak.
- MEF detection limit: blank + 2 × (rcv × median), inverted through the free-slope log-log calibration.
- Blank SD was taken as rcv × median. That gives the wrong sign when the blank median is ≤ 0.

`QB_2plus2SD_method_summary.md`, §5, explains why these were replaced.

## Files

| File | Restore to |
|---|---|
| `AuroraEVO_6PeakUltra.qmd` | `pages/` |
| `AuroraEVO_6PeakUltra_Dashboard.qmd` | `pages/` |
| `_6peak_load.R` | `pages/` |
| `Process6PeakUltra.R` | `scripts/` |
| `peak_metrics.csv` | `6peak Ultra_tracking_spectral/` (v1 schema, with no `sd` / `saturated` columns) |
| `*_v1.html` | Rendered pages, for reference only. Their CSS and JS point to `../site_libs`, so they show unstyled when opened from here. |

## Restoring v1

1. Copy the files back to the locations in the table.
2. From the repo root, run:
   ```
   quarto render pages/AuroraEVO_6PeakUltra.qmd
   quarto render pages/AuroraEVO_6PeakUltra_Dashboard.qmd
   ```

The v1 `Process6PeakUltra.R` only processes runs that are not yet in the archive.
If you restore the script but keep the v2 `peak_metrics.csv`, it will still work, because v1 ignores the extra columns.
