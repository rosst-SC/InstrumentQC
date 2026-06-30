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

# ---- Westgard-style rule flagging (Gain only) --------------------------------
# Reduced, drift-aware ruleset. EVO Gain is a CONTROLLED output (the daily-QC
# routine sets the voltage to park the bead peak at target), so it drifts up as
# detectors age -- that is correct, not a fault. So we DROP the drift-counting
# rules (10x, 4_1s) and R_4s (needs two control levels; we have one/day), and keep
# 1_3s (reject), 1_2s (warn) and 2_2s (sustained step). Each point is scored
# against a ROLLING ~45-day baseline BEFORE it (not the full-window band stats):
# z = (x - local median) / s_eff, where s_eff is floored at the detector's global
# robust SD AND a practical-significance fraction of the level (rel_floor), so a
# trivial gain nudge on a rock-stable detector isn't called "2 SD". Returns the
# plot data with z / loc_center / wg_rules / wg_sev / wg_label, or NULL.
westgard_eval <- function(p, trail_days = 45, min_n = 8, rel_floor = 0.025) {
  ycol <- p$labels$title
  d <- p$data
  if (is.null(ycol) || !ycol %in% names(d) || !"DateTime" %in% names(d)) return(NULL)
  d <- d[order(d$DateTime), , drop = FALSE]
  x  <- suppressWarnings(as.numeric(d[[ycol]]))
  dt <- as.POSIXct(d$DateTime)
  fcol <- paste0("Flag-", ycol)
  flagged <- if (fcol %in% names(d)) (suppressWarnings(as.logical(d[[fcol]])) %in% TRUE)
             else rep(FALSE, nrow(d))
  n <- nrow(d)
  z <- rep(NA_real_, n); loc <- rep(NA_real_, n)
  gmad <- stats::mad(x[!flagged & is.finite(x)], constant = 1.4826)   # characteristic noise
  if (is.finite(gmad) && gmad > 0) {
    win <- trail_days * 86400
    for (i in seq_len(n)) {
      if (!is.finite(x[i])) next
      sel <- which(dt >= (dt[i] - win) & dt < dt[i] & !flagged & is.finite(x))
      if (length(sel) < min_n) next
      m <- stats::median(x[sel]); s <- stats::mad(x[sel], constant = 1.4826)
      s_eff <- max(s, gmad, rel_floor * abs(m))
      loc[i] <- m; z[i] <- (x[i] - m) / s_eff
    }
  }
  rule <- rep(NA_character_, n); sev <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    if (is.na(z[i])) next
    rs <- character(0)
    if (abs(z[i]) > 3) rs <- c(rs, "1_3s") else if (abs(z[i]) > 2) rs <- c(rs, "1_2s")
    prev <- tail(which(!is.na(z[seq_len(i - 1)])), 1)                  # last scored point
    if (length(prev) == 1 &&
        ((z[i] > 2 && z[prev] > 2) || (z[i] < -2 && z[prev] < -2)))
      rs <- c(rs, "2_2s")
    if (length(rs) > 0) {
      rule[i] <- paste(rs, collapse = "+")
      sev[i]  <- if (any(rs %in% c("1_3s", "2_2s"))) "reject" else "warn"
    }
  }
  d$z <- z; d$loc_center <- loc; d$wg_rules <- rule; d$wg_sev <- sev
  d$wg_label <- ifelse(is.na(rule), NA_character_,
    sprintf("%s | %s | %.1f | z=%.2f | %s", ycol,
            format(as.Date(dt), "%Y-%m-%d"), x, z, rule))
  d
}

# Overlay Westgard violation markers ON TOP of the LJ bands for one plot.
add_westgard <- function(p, clamp0 = FALSE) {
  p <- add_lj_bands(p, clamp0 = clamp0)
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

# Render Gain plot list with bands + Westgard markers (interactive).
lj_render_wg <- function(plotlist)
  htmltools::tagList(lapply(lapply(plotlist, add_westgard), ggplotly,
                            tooltip = c("x", "y", "text")))

# Collect violations across a set of detector plots into one recent-first table.
westgard_summary <- function(plotlist) {
  rows <- lapply(plotlist, function(p) {
    e <- westgard_eval(p); if (is.null(e)) return(NULL)
    v <- e[!is.na(e$wg_rules), , drop = FALSE]
    if (nrow(v) == 0) return(NULL)
    data.frame(Detector = p$labels$title, Date = as.Date(v$DateTime),
               Rule = v$wg_rules, Value = round(suppressWarnings(as.numeric(v[[p$labels$title]])), 1),
               z = round(v$z, 2), Severity = v$wg_sev, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) return(data.frame(Detector = character(), Date = as.Date(character()),
                                      Rule = character(), Value = numeric(),
                                      z = numeric(), Severity = character()))
  out[order(out$Date, decreasing = TRUE), , drop = FALSE]
}
