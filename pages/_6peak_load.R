# Shared data-load + metric computation for the 6-peak Ultra Rainbow (URCP-38-2K)
# pages. Sourced by pages/AuroraEVO_6PeakUltra.qmd and
# pages/AuroraEVO_6PeakUltra_Dashboard.qmd so the metric math (Q/B fit, MEF
# calibration, resolution, SI, spectral signature) lives in one place. The Q/B/
# 2+2SD estimator is pages/_qb_parks.R (Parks 2017/2018); see Sensitivity
# Articles/QB_2plus2SD_method_summary.md before touching it or mef_map.
# v1 (CV^2-vs-1/x qb_fit + STN-17) is archived in _attic/6peak_ultra_v1_2026-11-11/.
library(ggcyto)   # scale_*_flowjo_biexp; pulls in flowCore (masks dplyr verbs)
library(ggplot2)
library(tidyr)
library(readr)
library(plotly)
library(dplyr)    # LAST so dplyr verbs win the namespace masking

results_file <- if (file.exists("6peak Ultra_tracking_spectral/peak_metrics.csv")) {
  "6peak Ultra_tracking_spectral/peak_metrics.csv"
} else "peak_metrics.csv"
pm <- read_csv(results_file, show_col_types = FALSE)
pm$datetime <- as.POSIXct(pm$datetime, tz = "UTC")
pm$date     <- as.Date(pm$datetime)

# URCP-38-2K resolves into 6 populations: blank particles (peak 1) + 5 fluorescent
# rainbow levels (peaks 2-6). Peak 6 is the brightest.
n_peaks  <- 6
biex_max <- 4194304; biex_pos <- 5.62; biex_neg <- 0; biex_width_basis <- -250

# Detector acquisition order (UV -> V -> B -> YG -> R) from the archive.
laser_levels <- c("UV", "V", "B", "YG", "R")
laser_of <- function(d) sub("[0-9]+-A$", "", as.character(d))
det_order <- pm |> distinct(detector) |>
  mutate(laser = laser_of(detector),
         num   = as.integer(sub("-A$", "", sub("^(UV|V|B|YG|R)", "", detector)))) |>
  arrange(match(laser, laser_levels), num)
detector_channels <- det_order$detector
pm$detector <- factor(pm$detector, levels = detector_channels)

# Run axis. Unlike the 8-peak set (one run/day) URCP can be acquired more than
# once a day, so the tracking unit is the acquisition TIMESTAMP. run_label is an
# ORDERED factor (chronological) used as the categorical x on every panel; it
# behaves like `date` did on the 8-peak page (arrange/first/last/min/max all work).
run_ord  <- pm |> distinct(run, datetime) |> arrange(datetime)
run_levs <- format(run_ord$datetime, "%Y-%m-%d %H:%M")
pm$run_label <- factor(format(pm$datetime, "%Y-%m-%d %H:%M"), levels = run_levs, ordered = TRUE)
runs   <- levels(pm$run_label); n_runs <- length(runs)
base_run <- runs[1]

# ---- GAIN cross-reference: the daily QC gain applied on each run ------------
# On the Aurora, daily QC pins bead MFI to a target by adjusting APD gain, so raw
# MFI is largely held constant by design. Dividing MFI by gain de-pins it, giving
# a ~gain-independent signal; CV-based sensitivity (Q/B) is gain-invariant. Gain
# is recorded once per day, so it is joined on the calendar date.
gfile <- "data/EVO/Archive/ArchivedDataEVO.csv"
has_gain <- file.exists(gfile)
if (has_gain) {
  g <- read.csv(gfile, check.names = FALSE)
  g$Date <- as.Date(lubridate::ymd_hms(g$DateTime))
  glong <- g |>
    select(Date, ends_with("-Gain")) |>
    pivot_longer(-Date, names_to = "detg", values_to = "gain") |>
    mutate(detector = paste0(sub("-Gain$", "", detg), "-A")) |>
    group_by(Date, detector) |> summarise(gain = dplyr::first(gain), .groups = "drop")
  pm <- pm |> left_join(glong, by = c("date" = "Date", "detector"))
} else pm$gain <- NA_real_
# Gain de-pins MFI. When SOME runs have recorded gain but a newer run does NOT yet
# (its daily-QC gain has not been exported to the archive), do NOT substitute raw
# MFI for that run: mixing raw and gain-normalised values across runs makes a
# spurious order-of-magnitude jump on the gain panels. Instead drop the gainless
# run from gain-normalised views (mfi_norm = NA). Only when NO run has any gain do
# we fall back to raw MFI globally (a consistent, if relative, scale).
any_gain <- has_gain && any(is.finite(pm$gain))
pm <- pm |> mutate(mfi_norm = dplyr::case_when(
  !is.na(gain) & gain > 0 ~ median / gain,
  any_gain                ~ NA_real_,
  TRUE                    ~ median))
has_gain <- any_gain
runs_no_gain <- pm |> filter(peak == n_peaks) |>
  group_by(run_label) |> summarise(g = any(is.finite(gain)), .groups = "drop") |>
  filter(!g) |> pull(run_label) |> as.character()

# ---- per-detector x per-run metrics ----
# Sensitivity (Parks 2017 weighted quadratic, V = c0 + c1*M + c2*M^2 on all
# non-saturated peaks incl. the blank) -> 2+2SD limit of resolution (Parks 2018).
# Headline `lor_mfi` (MFI, all 64 detectors); Q/B in Spe as relative views. Fits
# with chi2/df > 50 are suppressed (qflag says why). See pages/_qb_parks.R.
source("pages/_qb_parks.R")
qb <- pm |> group_by(run_label, detector) |> group_modify(~ qb_spe(.x)) |> ungroup()

# Resolution, brightest-/blank-peak %CV, brightest gain-normalised intensity
res   <- pm |> group_by(run_label, detector) |> summarise(n_resolved = sum(resolved), .groups = "drop")
cvtop <- pm |> filter(peak == n_peaks, in_range) |> select(run_label, detector, rcv) |> rename(cvt = rcv)
cvp1  <- pm |> filter(peak == 1, in_range) |> mutate(sd1 = sd) |>
  select(run_label, detector, sd1)   # blank NOISE FLOOR = right-half robust SD (raw MFI; never rcv*median)
ptop  <- pm |> filter(peak == n_peaks, in_range) |> select(run_label, detector, median, gain, mfi_norm)

# Staining index (lowest resolved peak vs blank peak 1)
si_one <- function(g) {
  cand <- g$peak[g$peak > 1 & g$resolved & g$in_range]
  if (!length(cand)) return(tibble(si = NA_real_, pos_peak = NA_integer_))
  p <- min(cand); m1 <- g$median[g$peak == 1]; rsd1 <- g$sd[g$peak == 1]   # right-half blank SD
  tibble(si = (g$median[g$peak == p] - m1) / (2 * rsd1), pos_peak = as.integer(p))
}
si <- pm |> group_by(run_label, detector) |> group_modify(~ si_one(.x)) |> ungroup()

# Linearity (rank-based, all 64): measured MFI vs peak rank, per detector-run. Blank
# (peak 1) excluded -- it is the blank, not a fluorescence-ladder level (matches calib).
lin <- pm |> filter(resolved, in_range, peak > 1) |> group_by(run_label, detector) |>
  filter(n() >= 3) |>
  summarise(r2 = summary(lm(log10(median) ~ peak))$r.squared, .groups = "drop")

# ---- MEF calibration (SPHERO URCP-38-2K certificate) for anchor detectors ---
# fluor -> detector map (standard Cytek Aurora 5L peak channels; every pairing was
# validated against these beads by a log-log median-vs-MEF fit, R^2 >= 0.96). The
# two blue-laser cross-excited tandems (PE-TexasRed/PE-Cy5 at 488ex) are excluded.
mef_map <- tibble::tribble(
  ~fluor,              ~detector,
  "BUV395",            "UV2-A",
  "BUV737",            "UV14-A",
  "BV421",             "V1-A",
  "BV510",             "V5-A",
  "BV605",             "V10-A",
  "BV650",             "V11-A",
  "BV711",             "V13-A",
  "BV786",             "V15-A",
  "FITC",              "B2-A",
  "PerCP-Cy5.5",       "B9-A",
  "PE",                "YG1-A",
  "PE-TexasRed(561)",  "YG3-A",
  "PE-Cy5(561)",       "YG5-A",
  "PE-Cy5.5",          "YG6-A",
  "PE-Cy7",            "YG9-A",
  "APC",               "R1-A",
  "Alexa700",          "R4-A",
  "APC-Cy7",           "R7-A")
mef_file <- "data/reference/URCP-38-2K_MEF.csv"
has_mef  <- file.exists(mef_file)
if (has_mef) {
  mefv <- readr::read_csv(mef_file, show_col_types = FALSE) |> inner_join(mef_map, by = "fluor")
  # PEAK ALIGNMENT: data peak 1 = blank particles (no assigned MEF; the
  # certificate's peak 1 is also blank/NA and absent from the CSV). The 5
  # fluorescent levels are data peaks 2-6 and line up DIRECTLY with certificate
  # peaks 2-6, so the join is 1:1 on peak (R^2 ~0.999; the log-log slope varies by
  # detector, 0.79-1.36 -- see mef_ok below).
  lin_mef <- pm |> mutate(detector = as.character(detector)) |>
    inner_join(mefv, by = c("detector", "peak")) |>
    filter(resolved, in_range, median > 0, MEF > 0)
  # Sensitivity in MEF units (Parks 2017/2018): k from a slope-1 log-log fit on
  # peaks 2-6; Q = k/c1, B = c0/(c1 k), BSD = blank SD / k, lor = 4(1/Q + BSD).
  # mef_ok = free log-log slope 0.9-1.1 & R^2 > 0.98; when FALSE the certificate's
  # MEF scale does not match this detector (most UV/V channels, slope 1.15-1.36),
  # so the MEF values are shown but flagged lower-confidence (mef_flag).
  mef_pm <- pm |> mutate(detector = as.character(detector)) |>
    inner_join(mef_map, by = "detector") |>
    left_join(mefv |> select(detector, peak, MEF), by = c("detector", "peak"))
  qb_mef <- mef_pm |> group_by(run_label, detector, fluor) |>
    group_modify(~ qb_mef_one(.x, .y$detector)) |> ungroup()
  calib <- qb_mef |> transmute(run_label, detector, fluor, b = ll_slope, r2 = ll_r2, k)
}

# Spectral signature drift: cosine similarity of the 64-detector signature of each
# peak vs the baseline run (gain-normalised). 1 = unchanged; drops = optical drift.
cos_sim <- function(a, b) sum(a * b) / (sqrt(sum(a * a)) * sqrt(sum(b * b)))
sig <- pm |> filter(in_range) |> select(run_label, peak, detector, mfi_norm) |>
  arrange(peak, detector)
spec <- sig |> group_by(peak) |> group_modify(function(g, k) {
  bv <- g |> filter(run_label == base_run)
  g |> group_by(run_label) |> summarise(
    cosine = cos_sim(mfi_norm, bv$mfi_norm[match(detector, bv$detector)]), .groups = "drop")
}) |> ungroup()

# Combined per-detector-run table + baseline deltas (for the overview heatmap)
metrics <- res |>
  left_join(qb |> select(run_label, detector, lor_mfi, Q, B, CV0, chi2_df, qflag),
            by = c("run_label","detector")) |>
  left_join(cvtop, by = c("run_label","detector")) |>
  left_join(si,    by = c("run_label","detector")) |>
  left_join(ptop |> select(run_label, detector, mfi_norm), by = c("run_label","detector"))
# Baseline = each detector's first run with a finite value (a suppressed first run,
# e.g. Q unresolved, would otherwise blank that detector's whole Δ series).
add_base <- function(df, v) df |> group_by(detector) |> arrange(run_label) |>
  mutate(.b = .data[[v]][which(is.finite(.data[[v]]))[1]],
         "{v}_d" := (.data[[v]] - .b) / abs(.b) * 100) |> select(-.b) |> ungroup()
for (v in c("lor_mfi","Q","B","cvt","si","mfi_norm","n_resolved")) metrics <- add_base(metrics, v)

# Control limits (mean +/- 2 SD) across runs, per detector, for the trend bands.
clim <- function(df, v) df |> group_by(detector) |>
  summarise(mu = mean(.data[[v]], na.rm = TRUE), sd = sd(.data[[v]], na.rm = TRUE), .groups = "drop") |>
  mutate(lcl = mu - 2*sd, ucl = mu + 2*sd)
cv_y_scale <- scale_y_continuous(expand = expansion(mult = c(0.06, 0.18)))

# Helper: per-detector x run heatmap (rows = detector in acquisition order).
# .heat1 / .heat_fc1 draw one panel for the detectors in `dets`; the public heat() /
# heat_fc() (below) split by laser so each laser gets its own colour scale.
.heat1 <- function(df, fill, title, sub = NULL, diverging = FALSE, limitsym = FALSE,
                   dets = detector_channels) {
  p <- ggplot(df, aes(run_label, detector, fill = .data[[fill]],
                      text = sprintf("%s | %s<br>%s: %.3g", detector, run_label, fill, .data[[fill]]))) +
    geom_tile(colour = "white", linewidth = 0.1) +
    scale_y_discrete(limits = rev(dets)) +
    labs(x = NULL, y = NULL, title = title, subtitle = sub, fill = "") +
    theme_bw(base_size = 8) + theme(axis.text.y = element_text(size = 7),
                                    axis.text.x = element_text(angle = 30, hjust = 1),
                                    plot.subtitle = element_text(size = 8, colour = "grey40"))
  if (diverging) {
    lim <- if (limitsym) max(abs(df[[fill]]), na.rm = TRUE) else NA
    p <- p + scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0,
                                  limits = if (is.na(lim)) NULL else c(-lim, lim))
  } else p <- p + scale_fill_viridis_c(option = "viridis")
  p
}
# One heatmap per laser, each with its own colour scale, so small changes within a
# laser are not washed out by the range of the other lasers. Returns a list of
# ggplots (class "laser_heat"); gp() renders it as stacked plotly panels.
heat_rows_px <- 20     # px per detector row in a per-laser panel
heat_pad_px  <- 120    # title, subtitle and x-axis labels
by_laser <- function(df, build, title) {
  df  <- df |> filter(!is.na(detector))
  las <- laser_of(as.character(df$detector))
  out <- list()
  for (L in laser_levels) {
    dets <- detector_channels[laser_of(detector_channels) == L]
    d    <- df[las == L, , drop = FALSE]
    if (!nrow(d)) next
    p <- build(d, sprintf("%s — %s laser", title, L), dets)
    attr(p, "n_det") <- length(dets)
    out[[L]] <- p
  }
  structure(out, class = "laser_heat")
}
heat <- function(df, fill, title, sub = NULL, diverging = FALSE, limitsym = FALSE)
  by_laser(df, function(d, t, dets) .heat1(d, fill, t, sub, diverging, limitsym, dets), title)

gp <- function(p, h = 620) {
  if (inherits(p, "laser_heat")) {
    # Each panel sits in a fixed-height, non-flexing div: on the dashboard, cards
    # otherwise squeeze the 5 stacked widgets to fit, so the card scrolls instead.
    return(do.call(htmltools::tagList, lapply(p, function(pp) {
      ph <- heat_pad_px + heat_rows_px * attr(pp, "n_det")
      htmltools::div(style = sprintf("height:%dpx; min-height:%dpx; flex:none;", ph, ph),
        ggplotly(pp, tooltip = "text", height = ph) |> plotly::config(displayModeBar = FALSE))
    })))
  }
  ggplotly(p, tooltip = "text", height = h) |> plotly::config(displayModeBar = FALSE)
}

# ---- Per-detector fold-change reference for the FC maps --------------------
# Absolute heatmaps span orders of magnitude across detectors, so the companion
# "FC" maps show each detector relative to its OWN reference. Two modes:
#   * NULL  -> DYNAMIC / preliminary: reference = each detector's own best run so far.
#   * c("YYYY-MM-DD HH:MM","...") -> LOCKED: per-detector robust centre (median) + SD
#     over a qualification window, frozen to fc_ref_file (two-sided deviation view).
fc_ref_window <- NULL                                   # <- set a run-label range to LOCK
fc_ref_file   <- "data/reference/EVO_6peakUltra_reference.csv"
fc_locked     <- !is.null(fc_ref_window)
fc_ref        <- NULL
if (fc_locked) {
  if (file.exists(fc_ref_file)) {
    fc_ref <- readr::read_csv(fc_ref_file, show_col_types = FALSE)
  } else {
    fc_metrics <- bind_rows(
      qb    |> transmute(run_label, detector, metric = "lor_mfi",  value = lor_mfi),
      qb    |> transmute(run_label, detector, metric = "Q",        value = Q),
      si    |> transmute(run_label, detector, metric = "si",       value = si),
      cvp1  |> transmute(run_label, detector, metric = "sd1",      value = sd1),
      ptop  |> transmute(run_label, detector, metric = "mfi_norm", value = mfi_norm))
    fc_ref <- fc_metrics |>
      filter(run_label >= fc_ref_window[1], run_label <= fc_ref_window[2], is.finite(value)) |>
      group_by(detector, metric) |>
      summarise(ref = median(value), ref_sd = sd(value), n_ref = dplyr::n(), .groups = "drop")
    readr::write_csv(fc_ref, fc_ref_file)
  }
}

# DYNAMIC reference: each detector vs its own best (one-sided). deg >= 0 (0 = best).
best_fc <- function(df, v, higher_better = TRUE) {
  df |> group_by(detector) |> mutate(
    .val  = ifelse(is.finite(.data[[v]]) & .data[[v]] > 0, .data[[v]], NA_real_),  # log2 needs > 0
    .best = { x <- .val[is.finite(.val)]
              if (!length(x)) NA_real_ else if (higher_better) max(x) else min(x) },
    deg   = if (higher_better) log2(.best / .val) else log2(.val / .best),
    fold  = 2^deg) |> ungroup()
}
.heat_fc1 <- function(df, v, title, sub = NULL, higher_better = TRUE, dets = detector_channels) {
  if (fc_locked) {
    r <- fc_ref |> filter(metric == v) |> transmute(detector = as.character(detector), ref, ref_sd)
    d <- df |> mutate(detector = as.character(detector)) |> left_join(r, by = "detector") |>
      mutate(value = .data[[v]],
             lfc   = log2(value / ref)        * (if (higher_better) -1 else 1),
             sddev = (value - ref) / ref_sd   * (if (higher_better) -1 else 1))
    hi <- stats::quantile(abs(d$lfc), 0.95, na.rm = TRUE); if (!is.finite(hi) || hi <= 0) hi <- 1
    ggplot(d, aes(run_label, detector, fill = pmax(pmin(lfc, hi), -hi),
          text = sprintf("%s | %s<br>%s = %.3g<br>%+.2f log₂ vs normal | %+.1f SD",
                         detector, run_label, v, value, lfc, sddev))) +
      geom_tile(colour = "white", linewidth = 0.1) +
      scale_y_discrete(limits = rev(dets)) +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                           limits = c(-hi, hi), name = "log₂ vs\nqualified\nnormal") +
      labs(x = NULL, y = NULL, title = title,
           subtitle = "deviation from locked reference (red = worse)") +
      theme_bw(base_size = 8) + theme(axis.text.y = element_text(size = 7),
                                      axis.text.x = element_text(angle = 30, hjust = 1),
                                      plot.subtitle = element_text(size = 8, colour = "grey40"))
  } else {
    d  <- best_fc(df, v, higher_better)
    hi <- stats::quantile(d$deg, 0.95, na.rm = TRUE); if (!is.finite(hi) || hi <= 0) hi <- 1
    ggplot(d, aes(run_label, detector, fill = pmin(deg, hi),
          text = sprintf("%s | %s<br>%s = %.3g<br>%.2f× its own best", detector, run_label, v, .data[[v]], fold))) +
      geom_tile(colour = "white", linewidth = 0.1) +
      scale_y_discrete(limits = rev(dets)) +
      scale_fill_gradient(low = "white", high = "#b2182b", limits = c(0, hi),
                          name = "log₂ fold\nvs best", na.value = "grey92") +
      labs(x = NULL, y = NULL, title = title,
           subtitle = if (is.null(sub)) "preliminary: vs each detector's best so far" else sub) +
      theme_bw(base_size = 8) + theme(axis.text.y = element_text(size = 7),
                                      axis.text.x = element_text(angle = 30, hjust = 1),
                                      plot.subtitle = element_text(size = 8, colour = "grey40"))
  }
}
heat_fc <- function(df, v, title, sub = NULL, higher_better = TRUE)
  by_laser(df, function(d, t, dets) .heat_fc1(d, v, t, sub, higher_better, dets), title)
