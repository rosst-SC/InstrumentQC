# =============================================================================
# Process6PeakUltra.R - 6-peak Ultra Rainbow bead processing (Cytek Aurora EVO)
#
# Reads SPHERO URCP-38-2K Ultra Rainbow FCS, gates beads, clusters the 6
# populations (blank + 5 fluorescent levels), and
#   1. APPENDS per-detector x per-peak metrics to an archive CSV (peak_metrics.csv)
#   2. SAVES two per-run diagnostic images (bead gate + 2D cluster assignment)
#      to images/6peak_ultra_evo/ so the website page can show them without FCS.
#
# Mirrors scripts/Process8Peak.R, with the URCP-specific differences baked into
# the cfg block:
#   * n_peaks = 6 (blank particles = peak 1; 5 fluorescent levels = peaks 2-6).
#     The two dimmest events are ONE population (the blank); at K=7 flowClust
#     over-split it in two, so K is fixed at 6 to keep the blank intact.
#   * Anchors are a non-saturating B/V/R set. They must NOT include UV (UV reads
#     fine sub-structure in the blank, which made flowClust split it in two) and
#     must NOT include detectors that SATURATE on the brightest bead (YG/red peg
#     at the 4.19M ceiling and merge the top populations). This set keeps every
#     peak on-scale and gives a stable 6-cluster fit across runs.
#   * Runs are keyed by acquisition TIMESTAMP (YYYYMMDD_HHMM), not date, because
#     more than one URCP run can be acquired in a day.
#
# INCREMENTAL: a run is (re)processed only if its metrics are not yet in the
# archive OR its diagnostic PNGs are missing. Re-running is cheap; the website
# page (AuroraEVO_6PeakUltra.qmd) only reads the archive + saved PNGs.
#
# Run from the InstrumentQC project root whenever new FCS arrive:
#   Rscript scripts/Process6PeakUltra.R
# =============================================================================

suppressMessages({
  library(flowCore); library(MASS); library(flowClust)
  library(ggcyto); library(ggplot2); library(scales)
  library(tidyr); library(readr)
  library(dplyr)   # load LAST so dplyr's verbs win over flowCore's masking
})

## ---- Instrument configuration (EVO, URCP-38-2K) -----------------------------
cfg_EVO6 <- list(
  data_dir = if (dir.exists("6peak Ultra_tracking_spectral/data/EVO_6peak Ultra"))
               "6peak Ultra_tracking_spectral/data/EVO_6peak Ultra" else "data/EVO_6peak Ultra",
  results_file = if (dir.exists("6peak Ultra_tracking_spectral"))
               "6peak Ultra_tracking_spectral/peak_metrics.csv" else "peak_metrics.csv",
  diag_dir = "images/6peak_ultra_evo",           # per-run diagnostic PNGs for the website
  detector_pattern = "^(UV|V|B|YG|R)[0-9]+-A$",
  scatter_x = "FSC-A", scatter_y = "SSC-A", gate_distance = 4,
  singlet_x = "FSC-A", singlet_y = "FSC-H", singlet_distance = 4,
  biex_max = 4194304, biex_pos = 5.62, biex_neg = 0, biex_width_basis = -250,
  n_peaks = 6, cluster_cofactor = 150,
  anchor_detectors = c("B4-A", "B7-A", "B8-A", "B9-A", "V3-A", "V8-A", "R2-A", "R5-A"),
  diag_detectors_2d = c("V8-A", "B9-A"),         # non-saturating; show 6 separated blobs
  resolution_k = 2
)

## ---- pipeline (extracted from the analysis qmd) -----------------------------
n_events <- function(ff) nrow(exprs(ff))

# Unique run id = the "YYYYMMDD_HHMM" stamp (falls back to date at 00:00). This is
# the de-dup / tracking key so two runs on the same day stay distinct.
run_id <- function(fname) {
  m <- regmatches(fname, regexpr("[0-9]{8}_[0-9]{4}", fname))
  if (length(m) == 1) return(m)
  d <- regmatches(fname, regexpr("[0-9]{8}", fname))
  if (length(d) == 1) return(paste0(d, "_0000"))
  NA_character_
}
parse_run_datetime <- function(fname) {
  r <- run_id(fname)
  if (is.na(r)) return(as.POSIXct(NA))
  as.POSIXct(r, format = "%Y%m%d_%H%M", tz = "UTC")
}

make_bead_gate <- function(ff, cfg) {
  m  <- exprs(ff)[, c(cfg$scatter_x, cfg$scatter_y), drop = FALSE]
  rc <- MASS::cov.rob(m)
  ellipsoidGate(filterId = "beads", .gate = rc$cov, mean = rc$center, distance = cfg$gate_distance)
}
make_singlet_gate <- function(ff, cfg) {
  m  <- exprs(ff)[, c(cfg$singlet_x, cfg$singlet_y), drop = FALSE]
  rc <- MASS::cov.rob(m)
  ellipsoidGate(filterId = "singlets", .gate = rc$cov, mean = rc$center, distance = cfg$singlet_distance)
}

robust_cv <- function(x) 100 * (IQR(x) / 1.349) / abs(median(x))

assign_peaks <- function(ff, cfg) {
  te  <- asinh(exprs(ff)[, cfg$anchor_detectors, drop = FALSE] / cfg$cluster_cofactor)
  set.seed(1)   # flowClust uses random EM starts; seed for reproducible peak labels/rCVs
  fit <- flowClust(flowFrame(te), varNames = cfg$anchor_detectors,
                   K = cfg$n_peaks, B = 500, trans = 0)
  lab <- flowClust::Map(fit, rm.outliers = FALSE)
  brightness <- tapply(rowMeans(te), lab, median)
  rank_map   <- setNames(rank(brightness), names(brightness))
  as.integer(rank_map[as.character(lab)])
}

peak_metrics_one <- function(ff, lab, name, detector_channels, cfg) {
  run <- run_id(name); dt <- parse_run_datetime(name)
  e <- exprs(ff); np <- cfg$n_peaks
  bind_rows(lapply(detector_channels, function(d) {
    a     <- asinh(e[, d] / cfg$cluster_cofactor)
    med_a <- tapply(a, lab, median)[as.character(1:np)]
    sd_a  <- tapply(a, lab, function(z) IQR(z) / 1.349)[as.character(1:np)]
    gap   <- pmin(c(Inf, diff(med_a)), c(diff(med_a), Inf))
    data.frame(
      file = name, run = run, datetime = dt, date = as.Date(dt),
      detector = d, peak = 1:np,
      median   = tapply(e[, d], lab, median)[as.character(1:np)],
      rcv      = tapply(e[, d], lab, robust_cv)[as.character(1:np)],
      n        = as.integer(table(factor(lab, levels = 1:np))),
      in_range = tapply(e[, d], lab, median)[as.character(1:np)] < cfg$biex_max,
      resolved = gap > cfg$resolution_k * sd_a
    )
  }))
}

## ---- per-run diagnostic plots (saved as PNGs for the website) ---------------
gate_plot_one <- function(ff, bead_res, run, cfg) {
  raw <- exprs(ff)[, c(cfg$scatter_x, cfg$scatter_y), drop = FALSE]
  df  <- data.frame(x = raw[, cfg$scatter_x], y = raw[, cfg$scatter_y], bead = bead_res@subSet)
  ggplot(df, aes(x = x, y = y, colour = bead)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick"),
                        labels = c("excluded", "bead"), name = NULL) +
    scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    labs(x = cfg$scatter_x, y = cfg$scatter_y, title = paste("Bead population gate —", run)) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme_bw()
}

cluster2d_plot_one <- function(gated, lab, run, cfg) {
  e  <- exprs(gated)
  df <- data.frame(x = e[, cfg$diag_detectors_2d[1]], y = e[, cfg$diag_detectors_2d[2]],
                   peak = factor(lab, levels = 1:cfg$n_peaks))
  ggplot(df, aes(x = x, y = y, colour = peak)) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_x_flowjo_biexp(maxValue = cfg$biex_max, pos = cfg$biex_pos,
                         neg = cfg$biex_neg, widthBasis = cfg$biex_width_basis) +
    scale_y_flowjo_biexp(maxValue = cfg$biex_max, pos = cfg$biex_pos,
                         neg = cfg$biex_neg, widthBasis = cfg$biex_width_basis) +
    labs(x = sprintf("%s (biexponential scale)", cfg$diag_detectors_2d[1]),
         y = sprintf("%s (biexponential scale)", cfg$diag_detectors_2d[2]),
         colour = "Peak", title = paste("Cluster assignment (2D) —", run)) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme_bw()
}

## ---- incremental driver -----------------------------------------------------
process_new <- function(cfg) {
  files <- list.files(cfg$data_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) stop("No FCS files found in: ", cfg$data_dir)
  file_runs <- vapply(basename(files), run_id, character(1))
  dir.create(cfg$diag_dir, showWarnings = FALSE, recursive = TRUE)

  done <- character(0)
  if (file.exists(cfg$results_file)) {
    done <- as.character(unique(read_csv(cfg$results_file, show_col_types = FALSE)$run))
  }
  gate_png <- function(r) file.path(cfg$diag_dir, paste0("gate_",      r, ".png"))
  clus_png <- function(r) file.path(cfg$diag_dir, paste0("cluster2d_", r, ".png"))

  # (Re)process a run if its metrics are missing OR a diagnostic image is missing.
  idx <- which(vapply(seq_along(files), function(i) {
    r <- file_runs[i]
    !(r %in% done) || !file.exists(gate_png(r)) || !file.exists(clus_png(r))
  }, logical(1)))

  if (length(idx) == 0) {
    message("Nothing to do - archive and diagnostics already current (", length(done), " run(s)).")
    return(invisible(FALSE))
  }
  message("Processing ", length(idx), " run(s) (new or missing diagnostics): ",
          paste(basename(files[idx]), collapse = ", "))

  set.seed(1)  # cov.rob (gating) and flowClust draw random subsamples
  new_metrics <- bind_rows(lapply(idx, function(i) {
    r        <- file_runs[i]
    ff       <- read.FCS(files[i], transformation = FALSE)
    dch      <- grep(cfg$detector_pattern, colnames(ff), value = TRUE)
    bead_res <- flowCore::filter(ff, make_bead_gate(ff, cfg))
    beads    <- Subset(ff, bead_res)
    gated    <- Subset(beads, flowCore::filter(beads, make_singlet_gate(beads, cfg)))
    lab      <- assign_peaks(gated, cfg)
    ggsave(gate_png(r), gate_plot_one(ff, bead_res, r, cfg), width = 6, height = 5, dpi = 110)
    ggsave(clus_png(r), cluster2d_plot_one(gated, lab, r, cfg), width = 6, height = 5, dpi = 110)
    peak_metrics_one(gated, lab, basename(files[i]), dch, cfg)
  }))

  key <- c("run", "detector", "peak")
  combined <- new_metrics
  if (file.exists(cfg$results_file)) {
    prior <- read_csv(cfg$results_file, show_col_types = FALSE)
    combined <- bind_rows(anti_join(prior, new_metrics, by = key), new_metrics)
  }
  combined <- arrange(combined, datetime, detector, peak)
  write_csv(combined, cfg$results_file)
  message("Wrote ", cfg$results_file, ": ", nrow(combined), " rows across ",
          length(unique(combined$run)), " run(s); diagnostics in ", cfg$diag_dir, "/")
  invisible(TRUE)
}

## ---- run --------------------------------------------------------------------
process_new(cfg_EVO6)
