# Shared Levey-Jennings + Westgard helpers for the Aurora daily-QC dashboards
# (EVO / 5L / 3L). Sourced from each page's setup chunk via
# source("pages/_lj_helpers.R") -- _quarto.yml sets execute-dir: project, so the
# path resolves from the project root regardless of which page is rendering.
# The "_" filename prefix keeps Quarto from treating this as a renderable page.
#
# Callers must already have loaded ggplot2, plotly, htmltools and knitr.

# Levey-Jennings control bands. QC_Plots() returns a list of trend ggplots
# (one per detector); this layers robust control limits *behind* the existing
# trace, QC-fail squares and engineer-visit verticals. Center/spread are robust
# (median +/- MAD-scaled SD); QC-fail points are excluded from the baseline but
# still plot. Limits are recomputed from the visible history at render time.
add_lj_bands <- function(p, k = c(2, 3), min_n = 8,
                         exclude_flagged = TRUE, clamp0 = FALSE) {
  ycol <- p$labels$title
  d <- p$data
  if (is.null(ycol) || !ycol %in% names(d) || !"DateTime" %in% names(d)) return(p)
  vals <- suppressWarnings(as.numeric(d[[ycol]]))
  keep <- is.finite(vals)
  fcol <- paste0("Flag-", ycol)
  if (exclude_flagged && fcol %in% names(d)) {
    fl <- suppressWarnings(as.logical(d[[fcol]]))
    keep <- keep & !(fl %in% TRUE)
  }
  base <- vals[keep]
  if (length(base) < min_n) return(p)                       # preliminary: no bands
  center <- stats::median(base, na.rm = TRUE)
  s <- stats::mad(base, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(p)                    # degenerate / flat
  k <- sort(unique(k)); inner <- k[1]; outer <- k[length(k)]
  xr <- range(d$DateTime, na.rm = TRUE)
  lo <- function(m) { v <- center - m * s; if (clamp0) max(v, 0) else v }
  hi <- function(m) center + m * s
  rect <- function(ymin, ymax, fillcol)
    geom_rect(data = data.frame(xmin = xr[1], xmax = xr[2], ymin = ymin, ymax = ymax),
              inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = fillcol, alpha = 0.12)
  band_layers <- list(
    rect(lo(inner), hi(inner), "#2ca02c"),   # in-control (green)
    rect(lo(outer), lo(inner), "#ff7f0e"),   # warning below (amber)
    rect(hi(inner), hi(outer), "#ff7f0e"),   # warning above (amber)
    geom_hline(yintercept = center, color = "grey40", linewidth = 0.4),
    geom_hline(yintercept = c(lo(outer), hi(outer)),
               color = "#d62728", linetype = "dashed", linewidth = 0.3)
  )
  p$layers <- c(band_layers, p$layers)
  p
}

# Thin wrapper: add bands to every plot in a list, then make them interactive.
lj_render <- function(plotlist, clamp0 = FALSE)
  htmltools::tagList(lapply(lapply(plotlist, add_lj_bands, clamp0 = clamp0), ggplotly))

# ---- 1-month trailing bands (Gain) -------------------------------------------
# Gain is a CONTROLLED output (daily QC sets the voltage to park the bead peak at
# target), so it drifts up as detectors age -- a full-history band is the wrong
# reference. Gain plots instead get a band that MOVES: each run is judged against
# the mean +/- SD of the same detector's runs in the trailing month BEFORE it
# (current run excluded; needs >= min_n prior runs and a non-zero SD, else no band
# at that run). This is exactly the home-page scoring in pages/_home_status.R
# (HOME_TRAIL_DAYS / HOME_MIN_N), so a point in the amber zone here is a home-page
# Caution and a point past the red dashed line is a Fail.
LJ_TRAIL_DAYS <- 30
LJ_MIN_N      <- 8

# Per-run trailing center / SD for one series (vectors aligned with dt).
trailing_stats <- function(dt, x, trail_days = LJ_TRAIL_DAYS, min_n = LJ_MIN_N) {
  t <- as.numeric(as.POSIXct(dt)); win <- trail_days * 86400
  n <- length(x); m <- rep(NA_real_, n); s <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    if (!is.finite(t[i])) next
    sel <- is.finite(t) & t >= t[i] - win & t < t[i] & is.finite(x)
    if (sum(sel) < min_n) next
    sd_i <- stats::sd(x[sel])
    if (!is.finite(sd_i) || sd_i <= 0) next
    m[i] <- mean(x[sel]); s[i] <- sd_i
  }
  list(center = m, sd = s)
}

# Layer the trailing band behind one Gain plot's existing layers. The band is
# drawn through each run's own limits, broken wherever a run has no band (start
# of history, or a gap with < min_n runs in the prior month).
add_lj_bands_trailing <- function(p, k = c(2, 3), trail_days = LJ_TRAIL_DAYS,
                                  min_n = LJ_MIN_N) {
  ycol <- p$labels$title
  d <- p$data
  if (is.null(ycol) || !ycol %in% names(d) || !"DateTime" %in% names(d)) return(p)
  d <- d[order(d$DateTime), , drop = FALSE]
  x <- suppressWarnings(as.numeric(d[[ycol]]))
  st <- trailing_stats(d$DateTime, x, trail_days, min_n)
  ok <- is.finite(st$center)
  if (!any(ok)) return(p)                                   # preliminary: no bands
  k <- sort(unique(k)); inner <- k[1]; outer <- k[length(k)]
  b <- data.frame(DateTime = d$DateTime, c = st$center, s = st$sd,
                  seg = cumsum(!ok))[ok, , drop = FALSE]
  b$lo2 <- b$c - inner * b$s; b$hi2 <- b$c + inner * b$s
  b$lo3 <- b$c - outer * b$s; b$hi3 <- b$c + outer * b$s
  rib <- function(ymin, ymax, fillcol)
    geom_ribbon(data = b, inherit.aes = FALSE,
                aes(x = DateTime, ymin = .data[[ymin]], ymax = .data[[ymax]], group = seg),
                fill = fillcol, alpha = 0.12)
  ln <- function(y, ...)
    geom_line(data = b, inherit.aes = FALSE,
              aes(x = DateTime, y = .data[[y]], group = seg), ...)
  band_layers <- list(
    rib("lo2", "hi2", "#2ca02c"),                           # in-control (green)
    rib("lo3", "lo2", "#ff7f0e"),                           # warning below (amber)
    rib("hi2", "hi3", "#ff7f0e"),                           # warning above (amber)
    ln("c",   color = "grey40", linewidth = 0.4),
    ln("lo3", color = "#d62728", linetype = "dashed", linewidth = 0.3),
    ln("hi3", color = "#d62728", linetype = "dashed", linewidth = 0.3)
  )
  p$layers <- c(band_layers, p$layers)
  p
}

# ---- Westgard-style rule flagging (Gain only) --------------------------------
# Reduced, drift-aware ruleset: 1_3s (reject) and 1_2s (warn) only. The
# drift-counting rules (10x, 4_1s) are dropped because gain legitimately drifts as
# detectors age, R_4s needs two control levels (we have one/day), and 2_2s was
# removed. Each point is scored against the same 1-month trailing band drawn on
# the plot (z = (x - trailing mean) / trailing SD), so markers line up with the
# amber / red zones. Returns the plot data with z / loc_center / wg_rules /
# wg_sev / wg_label, or NULL.
westgard_eval <- function(p, trail_days = LJ_TRAIL_DAYS, min_n = LJ_MIN_N) {
  ycol <- p$labels$title
  d <- p$data
  if (is.null(ycol) || !ycol %in% names(d) || !"DateTime" %in% names(d)) return(NULL)
  d <- d[order(d$DateTime), , drop = FALSE]
  x  <- suppressWarnings(as.numeric(d[[ycol]]))
  dt <- as.POSIXct(d$DateTime)
  st <- trailing_stats(dt, x, trail_days, min_n)
  z  <- (x - st$center) / st$sd
  rule <- ifelse(abs(z) > 3, "1_3s", ifelse(abs(z) > 2, "1_2s", NA_character_))
  sev  <- ifelse(is.na(rule), NA_character_, ifelse(rule == "1_3s", "reject", "warn"))
  d$z <- z; d$loc_center <- st$center; d$wg_rules <- rule; d$wg_sev <- sev
  d$wg_label <- ifelse(is.na(rule), NA_character_,
    sprintf("%s | %s | %.1f | z=%.2f | %s", ycol,
            format(as.Date(dt), "%Y-%m-%d"), x, z, rule))
  d
}

# Overlay Westgard violation markers ON TOP of the trailing bands for one plot.
add_westgard <- function(p) {
  p <- add_lj_bands_trailing(p)
  e <- westgard_eval(p)
  if (is.null(e)) return(p)
  v <- e[!is.na(e$wg_rules), , drop = FALSE]
  if (nrow(v) == 0) return(p)
  ycol <- p$labels$title
  v$.sev <- factor(v$wg_sev, levels = c("warn", "reject"))
  p + geom_point(data = v, inherit.aes = FALSE,
        aes(x = DateTime, y = .data[[ycol]], color = .sev, text = wg_label),
        shape = 1, stroke = 1.1, size = 4) +
      scale_color_manual(values = c(warn = "#ff7f0e", reject = "#d62728"),
                         drop = FALSE, guide = "none")
}

# Render Gain plot list with trailing bands + Westgard markers (interactive).
lj_render_wg <- function(plotlist)
  htmltools::tagList(lapply(lapply(plotlist, add_westgard), ggplotly,
                            tooltip = c("x", "y", "text")))

# Restrict a single named detector's plot to data on/after `cutoff` (a Date or
# "YYYY-MM-DD" string), e.g. after a deliberate gain step, so its bands / Westgard
# baseline and x-axis reflect only the new operating level. The trace layers
# inherit p$data (waiver) and the band/Westgard helpers read p$data directly, so
# replacing the matched plot's $data trims everything at once. Other plots in the
# list are returned untouched.
trim_detector_before <- function(plotlist, title, cutoff) {
  cutoff <- as.Date(cutoff)
  for (i in seq_along(plotlist)) {
    p <- plotlist[[i]]
    if (!is.null(p$labels$title) && p$labels$title == title && "DateTime" %in% names(p$data)) {
      d <- p$data
      plotlist[[i]]$data <- d[as.Date(d$DateTime) >= cutoff, , drop = FALSE]
    }
  }
  plotlist
}

# ---- Laser scaling factors: fixed tolerance ---------------------------------
# Area scaling factors should stay CONSISTENT (they sit near 1.0-1.3 and barely
# move day to day), so an SD band is the wrong check: a rock-stable factor gets a
# hair-trigger band. Instead each factor is held to a fixed absolute tolerance
# around its own long-run median (QC-fail runs excluded from the median). Drift
# of more than SCALING_TOL suggests a suspect laser or a damaged telescope.
SCALING_TOL <- 0.5

scaling_eval <- function(p, tol = SCALING_TOL) {
  ycol <- p$labels$title
  d <- p$data
  if (is.null(ycol) || !ycol %in% names(d) || !"DateTime" %in% names(d)) return(NULL)
  d <- d[order(d$DateTime), , drop = FALSE]
  x <- suppressWarnings(as.numeric(d[[ycol]]))
  fcol <- paste0("Flag-", ycol)
  flagged <- if (fcol %in% names(d)) (suppressWarnings(as.logical(d[[fcol]])) %in% TRUE)
             else rep(FALSE, nrow(d))
  base <- x[is.finite(x) & !flagged]
  if (length(base) == 0) return(NULL)
  center <- stats::median(base)
  dev <- x - center
  rule <- ifelse(is.finite(dev) & abs(dev) > tol, sprintf("drift>%g", tol), NA_character_)
  d$z <- NA_real_; d$dev <- dev; d$loc_center <- center
  d$wg_rules <- rule; d$wg_sev <- ifelse(is.na(rule), NA_character_, "reject")
  d$wg_label <- ifelse(is.na(rule), NA_character_,
    sprintf("%s | %s | %.3f | %+.3f from median %.3f | %s", ycol,
            format(as.Date(d$DateTime), "%Y-%m-%d"), x, dev, center, rule))
  d
}

# Fixed tolerance band (median +/- tol) behind the trace, plus red markers on
# runs outside it.
add_scaling_band <- function(p, tol = SCALING_TOL) {
  e <- scaling_eval(p, tol)
  if (is.null(e)) return(p)
  center <- e$loc_center[1]
  xr <- range(e$DateTime, na.rm = TRUE)
  band_layers <- list(
    geom_rect(data = data.frame(xmin = xr[1], xmax = xr[2],
                                ymin = center - tol, ymax = center + tol),
              inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#2ca02c", alpha = 0.12),
    geom_hline(yintercept = center, color = "grey40", linewidth = 0.4),
    geom_hline(yintercept = c(center - tol, center + tol),
               color = "#d62728", linetype = "dashed", linewidth = 0.3)
  )
  p$layers <- c(band_layers, p$layers)
  v <- e[!is.na(e$wg_rules), , drop = FALSE]
  if (nrow(v) == 0) return(p)
  ycol <- p$labels$title
  p + geom_point(data = v, inherit.aes = FALSE,
                 aes(x = DateTime, y = .data[[ycol]], text = wg_label),
                 color = "#d62728", shape = 1, stroke = 1.1, size = 4)
}

lj_render_scaling <- function(plotlist, tol = SCALING_TOL)
  htmltools::tagList(lapply(lapply(plotlist, add_scaling_band, tol = tol), ggplotly,
                            tooltip = c("x", "y", "text")))

# Collect violations across a set of detector plots into one recent-first table.
# `plotlist` is scored with the Westgard rules (Gain); `scaling` (optional) with
# the fixed scaling tolerance -- those rows have no z, and Value is the factor.
westgard_summary <- function(plotlist, scaling = list()) {
  one <- function(p, fun, digits) {
    e <- fun(p); if (is.null(e)) return(NULL)
    v <- e[!is.na(e$wg_rules), , drop = FALSE]
    if (nrow(v) == 0) return(NULL)
    data.frame(Detector = p$labels$title, Date = as.Date(v$DateTime),
               Rule = v$wg_rules,
               Value = round(suppressWarnings(as.numeric(v[[p$labels$title]])), digits),
               z = round(v$z, 2), Severity = v$wg_sev, stringsAsFactors = FALSE)
  }
  rows <- c(lapply(plotlist, one, fun = westgard_eval, digits = 1),
            lapply(scaling,  one, fun = scaling_eval,  digits = 3))
  out <- do.call(rbind, rows)
  if (is.null(out)) return(data.frame(Detector = character(), Date = as.Date(character()),
                                      Rule = character(), Value = numeric(),
                                      z = numeric(), Severity = character()))
  out[order(out$Date, decreasing = TRUE), , drop = FALSE]
}
