# Home-page QC status rules (index.qmd). Sourced from the project root
# (_quarto.yml sets execute-dir: project). Callers must have loaded dplyr.
#
# Each daily-QC run is scored over EVERY detector in the Gain archive
# (data/<inst>/Archive/ArchivedData<inst>.csv; FSC, SSC, SSC-B and all
# fluorescence detectors -- laser power/delay/scaling are not scored):
#
#   Gain : z = (gain - mean) / SD of the same detector's runs in the trailing
#          month BEFORE this run (current run excluded). Needs >= min_n prior
#          runs and a non-zero SD; otherwise that detector's gain is unscored.
#   rCV  : absolute %rCV.
#
#   Pass    : every |z| <= 2  and every rCV <= 6%
#   Caution : every |z| <= 3  (but some |z| > 2) and every rCV <= 6%
#   Fail    : any   |z| >  3  or  any   rCV >  6%
#
# rCV is pass/fail only: it never produces a Caution.
#
# Status codes use the Luciernaga colour vocabulary (Green / Yellow / Red, Gray =
# unscored) so SmallTable / SmallTableGlobal / InstrumentText / InstrumentColor
# render them unchanged.

HOME_TRAIL_DAYS <- 30
HOME_MIN_N      <- 8
HOME_Z_PASS     <- 2
HOME_Z_CAUTION  <- 3
HOME_RCV_LIMIT  <- 6

.worst <- function(s) {
  if (any(s == "Red", na.rm = TRUE)) "Red"
  else if (any(s == "Yellow", na.rm = TRUE)) "Yellow"
  else if (any(s == "Green", na.rm = TRUE)) "Green"
  else "Gray"
}

# Long per-run x per-detector scores for one instrument.
home_qc_scores <- function(inst, MainFolder, trail_days = HOME_TRAIL_DAYS,
                           min_n = HOME_MIN_N) {
  d <- Luciernaga:::CurrentData(x = inst, MainFolder = MainFolder, type = "Gain")
  d <- d[!is.na(d$DateTime), , drop = FALSE]
  d <- d[order(d$DateTime), , drop = FALSE]
  gcols <- grep("-Gain$", names(d), value = TRUE)
  gcols <- gcols[!startsWith(gcols, "Flag")]
  dets  <- sub("-Gain$", "", gcols)
  dets  <- dets[paste0(dets, "-% rCV") %in% names(d)]
  dt  <- as.numeric(d$DateTime)
  win <- trail_days * 86400

  out <- lapply(dets, function(det) {
    g <- suppressWarnings(as.numeric(d[[paste0(det, "-Gain")]]))
    r <- suppressWarnings(as.numeric(d[[paste0(det, "-% rCV")]]))
    z <- rep(NA_real_, length(g))
    for (i in seq_along(g)) {
      if (!is.finite(g[i])) next
      sel <- dt >= dt[i] - win & dt < dt[i] & is.finite(g)
      if (sum(sel) < min_n) next
      s <- stats::sd(g[sel])
      if (!is.finite(s) || s <= 0) next
      z[i] <- (g[i] - mean(g[sel])) / s
    }
    data.frame(
      DateTime = d$DateTime, Detector = det,
      GainValue = g, z = z, rCVValue = r,
      Gain = dplyr::case_when(is.na(z) ~ "Gray",
                              abs(z) <= HOME_Z_PASS ~ "Green",
                              abs(z) <= HOME_Z_CAUTION ~ "Yellow",
                              TRUE ~ "Red"),
      rCV = dplyr::case_when(is.na(r) ~ "Gray",
                             r <= HOME_RCV_LIMIT ~ "Green",
                             TRUE ~ "Red"),
      stringsAsFactors = FALSE)
  })
  dplyr::bind_rows(out)
}

# One status per run.
home_run_status <- function(scores) {
  scores |>
    dplyr::group_by(DateTime) |>
    dplyr::summarise(QCStatus = .worst(c(Gain, rCV)), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(DateTime))
}

# Per-detector table for the latest run: detectors that are not Pass first
# (worst first), then the Cytek indicator set, for SmallTable().
home_latest_table <- function(scores, indicators = "^FSC$|^SSC|V3$|B3$|R3$|YG3$") {
  latest <- scores[scores$DateTime == max(scores$DateTime), , drop = FALSE]
  rank <- c(Red = 3, Yellow = 2, Gray = 1, Green = 0)
  latest$sev <- pmax(rank[latest$Gain], rank[latest$rCV])
  keep <- latest$sev >= 2 | grepl(indicators, latest$Detector)
  latest[keep, , drop = FALSE] |>
    dplyr::arrange(dplyr::desc(sev), !grepl(indicators, Detector)) |>
    dplyr::transmute(Detector,
                     GainValue = round(GainValue),
                     `Gain z` = round(z, 1),
                     Gain,
                     rCVValue = round(rCVValue, 2),
                     rCV)
}

# Wide Date x Instrument history (latest run per day), newest first.
home_history <- function(run_status_list) {
  wide <- NULL
  for (inst in names(run_status_list)) {
    s <- run_status_list[[inst]] |>
      dplyr::mutate(Date = as.Date(DateTime)) |>
      dplyr::group_by(Date) |>
      dplyr::slice_max(DateTime, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::select(Date, !!inst := QCStatus)
    wide <- if (is.null(wide)) s else dplyr::full_join(wide, s, by = "Date")
  }
  wide |>
    dplyr::mutate(dplyr::across(-Date, ~ dplyr::na_if(., "Gray"))) |>
    dplyr::arrange(dplyr::desc(Date))
}

# SmallTable() variant with the gain z-score column under the Gain spanner.
home_small_table <- function(data) {
  pal <- function(x) dplyr::case_when(x == "Green" ~ "#0B6623", x == "Yellow" ~ "#BA8E23",
                                      x == "Red" ~ "#C80815", x == "Gray" ~ "#D3D3D3",
                                      TRUE ~ NA_character_)
  data |>
    gt::gt() |>
    gt::data_color(columns = c(Gain, rCV), fn = pal) |>
    gt::sub_values(values = "Green",  replacement = "Pass") |>
    gt::sub_values(values = "Yellow", replacement = "Caution") |>
    gt::sub_values(values = "Red",    replacement = "Fail") |>
    gt::sub_values(values = "Gray",   replacement = "n/a") |>
    gt::sub_missing(missing_text = "") |>
    gt::opt_table_font(font = "Montserrat") |>
    gt::cols_align(align = "center") |>
    gt::tab_spanner(label = "Gain ", columns = c(GainValue, `Gain z`, Gain)) |>
    gt::tab_spanner(label = "%RCV ", columns = c(rCVValue, rCV)) |>
    gt::cols_label(GainValue = "Value", `Gain z` = "z (1 mo)", Gain = "Status",
                   rCVValue = "Value", rCV = "Status")
}

# ---- What triggered each non-Pass run ---------------------------------------
# One row per Caution/Fail run: which check fired (Gain, rCV or both), how many
# gain detectors moved up vs down (a coordinated shift across many detectors is
# the signature to cross-check against the 6-peak MFI data -- bead degradation
# can fake a gain change), and the offending detectors, worst first.
home_triggers <- function(scores, inst, max_list = 6) {
  fmt_list <- function(x) {
    if (length(x) == 0) return("")
    if (length(x) > max_list)
      x <- c(x[seq_len(max_list)], sprintf("+%d more", length(x) - max_list))
    paste(x, collapse = ", ")
  }
  scores |>
    dplyr::group_by(DateTime) |>
    dplyr::filter(any(Gain %in% c("Yellow", "Red") | rCV == "Red")) |>
    dplyr::group_modify(function(s, key) {
      g  <- s[s$Gain %in% c("Yellow", "Red"), , drop = FALSE]
      g  <- g[order(-abs(g$z)), , drop = FALSE]
      r  <- s[s$rCV == "Red", , drop = FALSE]
      r  <- r[order(-r$rCVValue), , drop = FALSE]
      data.frame(
        Status  = .worst(c(s$Gain, s$rCV)),
        Trigger = paste(c(if (nrow(g)) "Gain", if (nrow(r)) "rCV"), collapse = " + "),
        Shift   = if (nrow(g)) sprintf("%d up / %d down", sum(g$z > 0), sum(g$z < 0)) else "",
        `Gain detectors` = fmt_list(sprintf("%s %+.2f SD", g$Detector, g$z)),
        `rCV detectors`  = fmt_list(sprintf("%s %.1f%%", r$Detector, r$rCVValue)),
        check.names = FALSE, stringsAsFactors = FALSE)
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(Date = as.Date(DateTime), Time = format(DateTime, "%H:%M"),
                  Instrument = inst, .before = 1) |>
    dplyr::select(-DateTime)
}

# gt table of triggers across instruments, newest first. Instrument links to its
# Levey-Jennings page so the flagged detectors can be inspected.
home_trigger_table <- function(trig) {
  pages <- c("3L" = "pages/Aurora3L.html", "5L" = "pages/Aurora5L.html",
             "EVO" = "pages/AuroraEVO.html")
  trig |>
    dplyr::arrange(dplyr::desc(Date), dplyr::desc(Time), Instrument) |>
    dplyr::mutate(Instrument = sprintf("[%s](%s)", Instrument, pages[Instrument])) |>
    gt::gt() |>
    gt::fmt_markdown(columns = Instrument) |>
    gt::data_color(columns = Status, fn = function(x) dplyr::case_when(
      x == "Yellow" ~ "#BA8E23", x == "Red" ~ "#C80815", TRUE ~ NA_character_)) |>
    gt::sub_values(columns = Status, values = "Yellow", replacement = "Caution") |>
    gt::sub_values(columns = Status, values = "Red",    replacement = "Fail") |>
    gt::opt_table_font(font = "Montserrat") |>
    gt::cols_align(align = "center", columns = c(Date, Time, Instrument, Status, Trigger, Shift)) |>
    gt::cols_label(Shift = "Gain shift") |>
    gt::cols_width(`Gain detectors` ~ gt::px(420), `rCV detectors` ~ gt::px(220),
                   Status ~ gt::px(90), Trigger ~ gt::px(100)) |>
    gt::opt_interactive(use_search = TRUE, use_filters = TRUE, use_pagination = TRUE,
                        page_size_default = 25, use_compact_mode = TRUE)
}
