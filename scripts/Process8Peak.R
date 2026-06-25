# =============================================================================
# Process8Peak.R - 8 peak bead processing (Cytek Aurora EVO)
#
# Reads RCP-30-5A 8-peak FCS, gates beads, clusters the 8 populations, and
#   1. APPENDS per-detector x per-peak metrics to an archive CSV (peak_metrics.csv)
#   2. SAVES two per-run diagnostic images (bead gate + 2D cluster assignment)
#      to images/8peak_evo/ so the website page can show them without the FCS.
#
# INCREMENTAL: a run is (re)processed only if its metrics are not yet in the
# archive OR its diagnostic PNGs are missing. So re-running is cheap, adding one
# new run only costs that run, and existing runs are backfilled with images the
# first time this updated script is run. The website page (AuroraEVO_8Peak.qmd)
# only reads the archive + the saved PNGs, so the site render stays light.
#
# Run from the InstrumentQC project root whenever new FCS arrive:
#   Rscript scripts/Process8Peak.R
#
# To add 3L/5L later: copy the EVO `cfg` block (anchor detectors differ by laser
# count), point data_dir/results_file/diag_dir at that instrument, call process_new().
# =============================================================================

suppressMessages({
  library(flowCore); library(MASS); library(flowClust)
  library(ggcyto); library(ggplot2); library(scales)
  library(tidyr); library(readr)
  library(dplyr)   # load LAST so dplyr's verbs win over flowCore's masking
})

## ---- Instrument configuration (EVO) -----------------------------------------
cfg_EVO <- list(
  data_dir = if (dir.exists("8peak_tracking_spectral/data/EVO_8peak"))
               "8peak_tracking_spectral/data/EVO_8peak" else "data/EVO_8peak",
  results_file = if (dir.exists("8peak_tracking_spectral"))
               "8peak_tracking_spectral/peak_metrics.csv" else "peak_metrics.csv",
  diag_dir = "images/8peak_evo",                 # per-run diagnostic PNGs for the website
  detector_pattern = "^(UV|V|B|YG|R)[0-9]+-A$",
  scatter_x = "FSC-A", scatter_y = "SSC-A", gate_distance = 4,
  singlet_x = "FSC-A", singlet_y = "FSC-H", singlet_distance = 4,
  biex_max = 4194304, biex_pos = 5.62, biex_neg = 0, biex_width_basis = -250,
  n_peaks = 8, cluster_cofactor = 150,
  anchor_detectors = c("UV9-A", "B4-A", "B7-A", "B8-A", "B9-A", "V8-A", "R2-A", "YG3-A"),
  diag_detectors_2d = c("V1-A", "YG1-A"),
  resolution_k = 2
)

## ---- pipeline (extracted from the analysis qmd) -----------------------------
n_events <- function(ff) nrow(exprs(ff))

parse_run_date <- function(fname) {
  m <- regmatches(fname, regexpr("[0-9]{8}", fname))
  if (length(m) == 1) as.Date(m, format = "%Y%m%d") else as.Date(NA)
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
  date <- parse_run_date(name); e <- exprs(ff); np <- cfg$n_peaks
  bind_rows(lapply(detector_channels, function(d) {
    a     <- asinh(e[, d] / cfg$cluster_cofactor)
    med_a <- tapply(a, lab, median)[as.character(1:np)]
    sd_a  <- tapply(a, lab, function(z) IQR(z) / 1.349)[as.character(1:np)]
    gap   <- pmin(c(Inf, diff(med_a)), c(diff(med_a), Inf))
    data.frame(
      file = name, date = date, detector = d, peak = 1:np,
      median   = tapply(e[, d], lab, median)[as.character(1:np)],
      rcv      = tapply(e[, d], lab, robust_cv)[as.character(1:np)],
      n        = as.integer(table(factor(lab, levels = 1:np))),
      in_range = tapply(e[, d], lab, median)[as.character(1:np)] < cfg$biex_max,
      resolved = gap > cfg$resolution_k * sd_a
    )
  }))
}

## ---- per-run diagnostic plots (saved as PNGs for the website) ---------------
gate_plot_one <- function(ff, bead_res, date, cfg) {
  raw <- exprs(ff)[, c(cfg$scatter_x, cfg$scatter_y), drop = FALSE]
  df  <- data.frame(x = raw[, cfg$scatter_x], y = raw[, cfg$scatter_y], bead = bead_res@subSet)
  ggplot(df, aes(x = x, y = y, colour = bead)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick"),
                        labels = c("excluded", "bead"), name = NULL) +
    scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    labs(x = cfg$scatter_x, y = cfg$scatter_y, title = paste("Bead population gate —", date)) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme_bw()
}

cluster2d_plot_one <- function(gated, lab, date, cfg) {
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
         colour = "Peak", title = paste("Cluster assignment (2D) —", date)) +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme_bw()
}

## ---- incremental driver -----------------------------------------------------
process_new <- function(cfg) {
  files <- list.files(cfg$data_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) stop("No FCS files found in: ", cfg$data_dir)
  file_dates <- vapply(basename(files), function(f) as.character(parse_run_date(f)), character(1))
  dir.create(cfg$diag_dir, showWarnings = FALSE, recursive = TRUE)

  done <- character(0)
  if (file.exists(cfg$results_file)) {
    done <- as.character(unique(read_csv(cfg$results_file, show_col_types = FALSE)$date))
  }
  gate_png <- function(d) file.path(cfg$diag_dir, paste0("gate_",      d, ".png"))
  clus_png <- function(d) file.path(cfg$diag_dir, paste0("cluster2d_", d, ".png"))

  # (Re)process a run if its metrics are missing OR a diagnostic image is missing.
  idx <- which(vapply(seq_along(files), function(i) {
    d <- file_dates[i]
    !(d %in% done) || !file.exists(gate_png(d)) || !file.exists(clus_png(d))
  }, logical(1)))

  if (length(idx) == 0) {
    message("Nothing to do - archive and diagnostics already current (", length(done), " date(s)).")
    return(invisible(FALSE))
  }
  message("Processing ", length(idx), " run(s) (new or missing diagnostics): ",
          paste(basename(files[idx]), collapse = ", "))

  set.seed(1)  # cov.rob (gating) and flowClust draw random subsamples
  new_metrics <- bind_rows(lapply(idx, function(i) {
    d        <- file_dates[i]
    ff       <- read.FCS(files[i], transformation = FALSE)
    dch      <- grep(cfg$detector_pattern, colnames(ff), value = TRUE)
    bead_res <- flowCore::filter(ff, make_bead_gate(ff, cfg))
    beads    <- Subset(ff, bead_res)
    gated    <- Subset(beads, flowCore::filter(beads, make_singlet_gate(beads, cfg)))
    lab      <- assign_peaks(gated, cfg)
    ggsave(gate_png(d), gate_plot_one(ff, bead_res, d, cfg), width = 6, height = 5, dpi = 110)
    ggsave(clus_png(d), cluster2d_plot_one(gated, lab, d, cfg), width = 6, height = 5, dpi = 110)
    peak_metrics_one(gated, lab, basename(files[i]), dch, cfg)
  }))

  key <- c("date", "detector", "peak")
  combined <- new_metrics
  if (file.exists(cfg$results_file)) {
    prior <- read_csv(cfg$results_file, show_col_types = FALSE)
    combined <- bind_rows(anti_join(prior, new_metrics, by = key), new_metrics)
  }
  combined <- arrange(combined, date, detector, peak)
  write_csv(combined, cfg$results_file)
  message("Wrote ", cfg$results_file, ": ", nrow(combined), " rows across ",
          length(unique(combined$date)), " date(s); diagnostics in ", cfg$diag_dir, "/")
  invisible(TRUE)
}

## ---- run --------------------------------------------------------------------
process_new(cfg_EVO)
