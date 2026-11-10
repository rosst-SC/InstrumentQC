# Shared data-load + metric computation for the 6-peak Ultra Rainbow (URCP-38-2K)
# pages. Sourced by pages/AuroraEVO_6PeakUltra.qmd and
# pages/AuroraEVO_6PeakUltra_Dashboard.qmd so the metric math (Q/B fit, MEF
# calibration, resolution, SI, spectral signature) lives in one place. See
# Sensitivity Articles/QB_sensitivity_reference.md before touching qb_fit()/
# stn17/mef_map.
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
# Sensitivity (Q & B): CV^2 = intrinsic_CV^2 + (1/Q)*(1/x) + (B/Q)*(1/x^2),
# where x is the per-event intensity. Two flavours of x:
#   - xvar = "mfi_norm" (default): gain-normalised MFI -> Q/B in arbitrary but
#     gain-stabilised units; available for all 64 detectors (relative drift).
#   - xvar = "MEF": URCP-38-2K assigned MEF -> Q/B in true MEF units (STN-17),
#     gain-independent; anchor detectors only (absolute sensitivity).
# GUARD: Q = 1/slope explodes when the slope -> 0, so a noisy run can emit a phantom
# Q spike. We suppress Q/B (-> NA, with a `qflag` reason) when the fit is ill-
# conditioned: near-flat/negative slope, poor R^2, a non-monotonic CV profile
# (CV^2 should rise with 1/x), or an inflated brightest-peak CV.
qb_fit <- function(d, xvar = "mfi_norm", true_b = FALSE,
                   min_r2 = 0.8, min_mono = 0.3, max_bright_cv = 15) {
  d0 <- d                                                              # full peak set (incl. blank) for SD-ratio gate
  # The blank (peak 1) is used ONLY as the blank (background/SD-ratio below, and the
  # STN-17 closed form) -- NOT as a regression bead (STN-17/Hoffman); so the fit uses
  # the resolved FLUORESCENT peaks only (peak > 1).
  d  <- d[which(d$resolved & d$in_range & d$peak > 1 & !is.na(d[[xvar]]) & d[[xvar]] > 0), ]
  if (nrow(d) < 3)
    return(tibble(intr_cv = NA_real_, Q = NA_real_, B = NA_real_, npts = nrow(d), qflag = "n<3", qadv = NA_character_))
  o    <- order(d[[xvar]])                                             # dim -> bright
  cv2  <- (d$rcv / 100)^2; inv <- 1 / d[[xvar]]
  # Intrinsic CV is MEASURED, not extrapolated: the brightest in-range peak's rCV is
  # the high-signal floor (its photon term ~0; Hoffman & Wood 2007). Subtract it and
  # fit the intrinsic-corrected CV^2 through the origin, so Q (linear coef) and B/Q
  # (quadratic coef) come from a 2-param fit rather than a 3-param fit with a free,
  # extrapolated intercept. INVERSE-VARIANCE WEIGHTED (Var(CV^2) ~ CV^4 at const n),
  # since Hoffman warns unweighted fits are dominated by one extreme point and the
  # error in B is up to 10x that in Q.
  intr <- d$rcv[o[length(o)]] / 100                                    # measured intrinsic CV (brightest)
  cv2c <- pmax(cv2 - intr^2, 0)                                        # intrinsic-corrected CV^2
  w    <- 1 / pmax(cv2, 1e-6)^2                                        # inverse-variance weights
  s1   <- coef(lm(cv2c ~ 0 + inv, weights = w))[1]                     # slope = 1/Q (weighted, through origin)
  r2   <- summary(lm(cv2 ~ inv))$r.squared                            # linearity diagnostic (unweighted)
  mono <- suppressWarnings(stats::cor(inv, cv2, method = "spearman"))  # CV^2 vs 1/x: want > 0
  # Validity-gate reference peaks (Chase 1998; Perfetto 2014): "Dim-2" = the dimmest
  # fit peak (the blank is already excluded from d, so o[1] is the dimmest fluorescent)
  # vs the blank (peak 1, from d0) for the SD ratio. Applied ADVISORY below, not as a
  # hard gate -- spectral rainbow beads often miss these thresholds.
  dim_i   <- o[1]                                                     # dimmest fluorescent = Chase "Dim-2"
  rcv_dim <- d$rcv[dim_i]; rcv_brt <- d$rcv[o[length(o)]]
  sd_dim  <- d$rcv[dim_i] / 100 * d$median[dim_i]
  bl      <- d0[which(d0$peak == 1), ]                               # blank absent in the MEF fit (lin_mef
  sd_blk  <- if (nrow(bl)) bl$rcv[1] / 100 * bl$median[1] else NA_real_  # has no blank) -> SDratio then skipped
  # HARD guards: the fit is ill-conditioned -> suppress Q/B (blank tile, NOT a gain).
  flags <- c(if (!is.finite(s1) || s1 <= 0)                    "flat slope"              else NULL,
             if (is.finite(r2)   && r2   < min_r2)              sprintf("R2=%.2f", r2)    else NULL,
             if (is.finite(mono) && mono < min_mono)            "non-monotonic"          else NULL,
             if (is.finite(intr) && intr * 100 > max_bright_cv) sprintf("brightCV=%.0f%%", intr * 100) else NULL)
  # ADVISORY photon-domination criteria (Chase 1998; Hoffman allows 2x once the CV is
  # intrinsic-corrected, which it is here): dim-peak rCV >= 2x brightest, dim-peak rSD
  # >= 3x blank. These do NOT suppress -- on spectral rainbow beads they are frequently
  # unmet (in-range peaks sit in the bright/intrinsic-dominated regime and the "blank"
  # is itself a dim bead), which is the caveat to flag, not a fit failure.
  adv <- c(if (is.finite(rcv_dim) && is.finite(rcv_brt) && rcv_dim < 2 * rcv_brt) sprintf("CVratio=%.1f", rcv_dim / rcv_brt) else NULL,
           if (is.finite(sd_blk)  && sd_dim < 3 * sd_blk)       sprintf("SDratio=%.1f", sd_dim / sd_blk) else NULL)
  qadv <- if (length(adv)) paste(adv, collapse = "; ") else NA_character_
  if (length(flags))
    return(tibble(intr_cv = intr, Q = NA_real_, B = NA_real_, npts = nrow(d),
                  qflag = paste(flags, collapse = "; "), qadv = qadv))
  s2 <- if (nrow(d) >= 4) coef(lm(cv2c ~ 0 + inv + I(inv^2), weights = w))[2] else NA_real_
  B  <- if (true_b) s2 / s1 else s2
  tibble(intr_cv = intr, Q = 1 / s1, B = B, npts = nrow(d), qflag = NA_character_, qadv = qadv)
}
qb <- pm |> group_by(run_label, detector) |> group_modify(~ qb_fit(.x)) |> ungroup()

# Resolution, brightest-/blank-peak %CV, brightest gain-normalised intensity
res   <- pm |> group_by(run_label, detector) |> summarise(n_resolved = sum(resolved), .groups = "drop")
cvtop <- pm |> filter(peak == n_peaks, in_range) |> select(run_label, detector, rcv) |> rename(cvt = rcv)
cvp1  <- pm |> filter(peak == 1, in_range) |> mutate(sd1 = rcv / 100 * median) |>
  select(run_label, detector, sd1)   # blank NOISE FLOOR = absolute robust SD; %CV unstable near 0 mean
ptop  <- pm |> filter(peak == n_peaks, in_range) |> select(run_label, detector, median, gain, mfi_norm)

# Staining index (lowest resolved peak vs blank peak 1)
si_one <- function(g) {
  cand <- g$peak[g$peak > 1 & g$resolved & g$in_range]
  if (!length(cand)) return(tibble(si = NA_real_, pos_peak = NA_integer_))
  p <- min(cand); m1 <- g$median[g$peak == 1]; rsd1 <- abs(g$rcv[g$peak == 1] * m1 / 100)
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
  # peaks 2-6, so the join is 1:1 on peak (median R^2 ~0.9999, slope ~1.17).
  lin_mef <- pm |> mutate(detector = as.character(detector)) |>
    inner_join(mefv, by = c("detector", "peak")) |>
    filter(resolved, in_range, median > 0, MEF > 0)
  # log-log calibration per anchor x run (slope ~1 if linear; intercept holds gain)
  calib <- lin_mef |> group_by(run_label, detector, fluor) |> filter(n() >= 3) |>
    summarise(b = coef(lm(log10(median) ~ log10(MEF)))[2],
              a = coef(lm(log10(median) ~ log10(MEF)))[1],
              r2 = summary(lm(log10(median) ~ log10(MEF)))$r.squared, .groups = "drop")
  # Detection limit in MEF: blank (peak 1) + 2 robust SD, inverted through the
  # calibration. The applied gain cancels, so MEF_LOD is gain-independent.
  blank <- pm |> filter(peak == 1) |> mutate(detector = as.character(detector)) |>
    semi_join(mef_map, by = "detector") |>
    transmute(run_label, detector, thr = median * (1 + 2 * rcv / 100))
  sens <- calib |> left_join(blank, by = c("run_label", "detector")) |>
    mutate(MEF_LOD = 10^((log10(thr) - a) / b))
  # MEF-unit Q & B (STN-17): same CV^2-vs-1/x regression as `qb`, x = assigned MEF
  # -> Q/B in true molecules-of-equivalent-fluorophore units, gain-independent.
  qb_mef <- lin_mef |> group_by(run_label, detector, fluor) |>
    group_modify(~ qb_fit(.x, xvar = "MEF", true_b = TRUE)) |> ungroup()
  # ---- Literal STN-17 closed form (independent cross-check) ----
  # All three populations selected PER DETECTOR-RUN from lin_mef (already resolved,
  # in-range, MEF-assigned): blank (Dim Bead 1 = peak 1, extrapolated via calibration);
  # dim reference = the LOWEST fluorescent peak resolved from the blank; brightest = the
  # HIGHEST resolved in-range peak (not a fixed peak 6, which may be saturated) so its
  # CV^2 matches the high-signal floor the regression uses for intrinsic.
  p_cv <- function(p) pm |> filter(peak == p, median > 0) |>
    mutate(detector = as.character(detector)) |> semi_join(mef_map, by = "detector") |>
    transmute(run_label, detector, m = median, cv = rcv / 100)
  b1   <- p_cv(1) |> rename(m1 = m, cv1 = cv)
  bdim <- lin_mef |> group_by(run_label, detector) |> filter(peak == min(peak)) |> ungroup() |>
    transmute(run_label, detector, cvdim = rcv / 100, MEFdim = MEF, dimpeak = as.integer(peak))
  btop <- lin_mef |> group_by(run_label, detector) |> filter(peak == max(peak)) |> ungroup() |>
    transmute(run_label, detector, cvtop = rcv / 100, toppeak = as.integer(peak))
  stn17 <- calib |> select(run_label, detector, fluor, a, b) |>
    left_join(b1,   by = c("run_label", "detector")) |>
    left_join(bdim, by = c("run_label", "detector")) |>
    left_join(btop, by = c("run_label", "detector")) |>
    mutate(
      MEF_blank = 10^((log10(m1) - a) / b),                  # blank MEF via calib
      # GUARD: like qb_fit, the closed form is invalid on a noisy run. It needs the
      # blank/dim peaks BROADER than the brightest (cv1,cvdim > cvtop) or the SD^2
      # terms go negative, and a finite positive extrapolated MEF_blank; else -> NA.
      ok_blank  = is.finite(MEF_blank) & MEF_blank > 0 & is.finite(cv1) & is.finite(cvtop) & cv1 > cvtop,
      ok_dim    = is.finite(cvdim) & is.finite(cvtop) & cvdim > cvtop,
      Q_stn     = ifelse(ok_blank, 1 / (cv1^2 * MEF_blank), NA_real_),
      SD2_blank = ifelse(ok_blank, MEF_blank^2 * (cv1^2   - cvtop^2), NA_real_),
      SD2_dim   = ifelse(ok_dim,   MEFdim^2    * (cvdim^2 - cvtop^2), NA_real_),
      B_stn     = ifelse(ok_blank & ok_dim & is.finite(SD2_dim) & SD2_dim > 0,
                         (SD2_blank / SD2_dim) * MEF_blank, NA_real_))
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
  left_join(qb,    by = c("run_label","detector")) |>
  left_join(cvtop, by = c("run_label","detector")) |>
  left_join(si,    by = c("run_label","detector")) |>
  left_join(ptop |> select(run_label, detector, mfi_norm), by = c("run_label","detector"))
add_base <- function(df, v) df |> group_by(detector) |> arrange(run_label) |>
  mutate("{v}_d" := (.data[[v]] - first(.data[[v]])) / abs(first(.data[[v]])) * 100) |> ungroup()
for (v in c("Q","B","cvt","si","mfi_norm","n_resolved")) metrics <- add_base(metrics, v)

# Control limits (mean +/- 2 SD) across runs, per detector, for the trend bands.
clim <- function(df, v) df |> group_by(detector) |>
  summarise(mu = mean(.data[[v]], na.rm = TRUE), sd = sd(.data[[v]], na.rm = TRUE), .groups = "drop") |>
  mutate(lcl = mu - 2*sd, ucl = mu + 2*sd)
cv_y_scale <- scale_y_continuous(expand = expansion(mult = c(0.06, 0.18)))

# Helper: per-detector x run heatmap (rows = detector in acquisition order)
heat <- function(df, fill, title, sub = NULL, diverging = FALSE, limitsym = FALSE) {
  p <- ggplot(df, aes(run_label, detector, fill = .data[[fill]],
                      text = sprintf("%s | %s<br>%s: %.3g", detector, run_label, fill, .data[[fill]]))) +
    geom_tile(colour = "white", linewidth = 0.1) +
    scale_y_discrete(limits = rev(detector_channels)) +
    labs(x = NULL, y = NULL, title = title, subtitle = sub, fill = "") +
    theme_bw(base_size = 8) + theme(axis.text.y = element_text(size = 5),
                                    axis.text.x = element_text(angle = 30, hjust = 1),
                                    plot.subtitle = element_text(size = 8, colour = "grey40"))
  if (diverging) {
    lim <- if (limitsym) max(abs(df[[fill]]), na.rm = TRUE) else NA
    p <- p + scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0,
                                  limits = if (is.na(lim)) NULL else c(-lim, lim))
  } else p <- p + scale_fill_viridis_c(option = "viridis")
  p
}
gp <- function(p, h = 620) ggplotly(p, tooltip = "text", height = h) |>
  plotly::config(displayModeBar = FALSE)

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
heat_fc <- function(df, v, title, sub = NULL, higher_better = TRUE) {
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
      scale_y_discrete(limits = rev(detector_channels)) +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                           limits = c(-hi, hi), name = "log₂ vs\nqualified\nnormal") +
      labs(x = NULL, y = NULL, title = title,
           subtitle = "deviation from locked reference (red = worse)") +
      theme_bw(base_size = 8) + theme(axis.text.y = element_text(size = 5),
                                      axis.text.x = element_text(angle = 30, hjust = 1),
                                      plot.subtitle = element_text(size = 8, colour = "grey40"))
  } else {
    d  <- best_fc(df, v, higher_better)
    hi <- stats::quantile(d$deg, 0.95, na.rm = TRUE); if (!is.finite(hi) || hi <= 0) hi <- 1
    ggplot(d, aes(run_label, detector, fill = pmin(deg, hi),
          text = sprintf("%s | %s<br>%s = %.3g<br>%.2f× its own best", detector, run_label, v, .data[[v]], fold))) +
      geom_tile(colour = "white", linewidth = 0.1) +
      scale_y_discrete(limits = rev(detector_channels)) +
      scale_fill_gradient(low = "white", high = "#b2182b", limits = c(0, hi),
                          name = "log₂ fold\nvs best", na.value = "grey92") +
      labs(x = NULL, y = NULL, title = title,
           subtitle = if (is.null(sub)) "preliminary: vs each detector's best so far" else sub) +
      theme_bw(base_size = 8) + theme(axis.text.y = element_text(size = 5),
                                      axis.text.x = element_text(angle = 30, hjust = 1),
                                      plot.subtitle = element_text(size = 8, colour = "grey40"))
  }
}
