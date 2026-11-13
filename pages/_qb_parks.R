# =============================================================================
# _qb_parks.R - Q, B, BSD, CV0 and the 2+2SD limit of resolution from multi-peak
# beads (Parks et al. 2017 weighted quadratic + Parks et al. 2018 2+2SD).
#
# Ported from the Sensitivity project (~/Documents/Positron_Local/Sensitvity/
# sensitivity_functions.R: parks_fit(), mef_scale()). Every choice is explained in
# Sensitivity Articles/QB_2plus2SD_method_summary.md, which supersedes the
# conflicting parts of QB_sensitivity_reference.md (the STN-17 closed form,
# brightest-CV subtraction, and 1/CV^4 weights are all retired).
#
# Model, per detector x run, in raw MFI:   V = c0 + c1*M + c2*M^2
#   fit on ALL non-saturated peaks INCLUDING the blank, weights n/(2*Vhat^2),
#   iterated. V = sd^2 where sd is IQR/1.349 (blank: right-half P84-P50).
#
# Headline for ALL 64 detectors -- 2+2SD in as-acquired MFI (no MEF needed):
#   lor_mfi = 4*(c1 + blank_sd)      (= LoR_MEF * k; c1 = MFI per Spe)
# On the Aurora, daily QC pins bead MFI to fixed targets, so the MFI scale is held
# constant (URCP k varies only ~2% across all runs, incl. the v2-gain switch) and
# lor_mfi tracks the MEF limit. It is background-limited (4/Q share ~2-5%), so it
# is stable run-to-run (~5-10% CV), whereas Q itself carries ~30% fit SE per run.
# A detector losing sensitivity shows as more blank noise at the pinned scale.
#
# Photoelectron (Spe) units -- gain-invariant in principle, but they divide by the
# noisy c1, so they are secondary/relative views:
#   Q_I = 1/c1 Spe per MFI (gain-dependent -> Q = Q_I * gain, relative);
#   B = B_Spe = c0/c1^2; BSD_Spe = sqrt(c0)/c1; lor_spe = 4*(1 + BSD_Spe); CV0 = sqrt(c2)
# MEF units (detectors with a URCP certificate assignment), k = MFI per MEF:
#   Q = k/c1,  B = c0/(c1*k),  BSD = blank sd / k,  LoR = 4*(1/Q + BSD)
# =============================================================================

## ---- Parks 2017 weighted quadratic (variance space) --------------------------
# d: one detector x run, columns peak, median, sd, n, saturated.
parks_fit <- function(d, max_iter = 20) {
  f <- d[which(!d$saturated & is.finite(d$sd) & d$sd > 0), ]
  M <- f$median; V <- f$sd^2; n <- f$n
  na_out <- function(flag) tibble(c0 = NA_real_, c1 = NA_real_, c2 = NA_real_,
                                  c1_se = NA_real_, vc = list(NULL), chi2_df = NA_real_,
                                  npts = length(M), flag = flag)
  if (length(M) < 4) return(na_out("too few points"))
  fit_once <- function(quadratic) {
    w <- n / (2 * V^2)
    for (i in seq_len(max_iter)) {
      fit <- if (quadratic) lm(V ~ M + I(M^2), weights = w) else lm(V ~ M, weights = w)
      Vh  <- pmax(fitted(fit), 1e-12 * max(V))
      w_new <- n / (2 * Vh^2)
      if (max(abs(w_new - w) / w) < 1e-6) { w <- w_new; break }
      w <- w_new
    }
    fit
  }
  fit <- fit_once(TRUE); flag <- character(0)
  cc  <- coef(fit)
  if (!is.finite(cc[3]) || cc[3] < 0) {
    fit <- fit_once(FALSE); flag <- c(flag, "c2<0: CV0 fixed at 0"); cc <- c(coef(fit), 0)
  }
  vc <- vcov(fit); if (ncol(vc) == 2) vc <- rbind(cbind(vc, 0), 0)
  if (!is.finite(cc[2]) || cc[2] <= 0) return(na_out("c1<=0: Q undefined"))
  if (cc[1] < 0) flag <- c(flag, "c0<0: B undefined")
  Vf <- fitted(fit)
  # Goodness of fit vs pure sampling Var(V) = 2V^2/n; IQR-based SDs are less
  # efficient, so ~2-7 is normal; > 50 means the model does not describe the data.
  w_res   <- (V - Vf) / (Vf * sqrt(2 / n))
  chi2_df <- sum(w_res^2) / max(length(M) - length(coef(fit)), 1)
  tibble(c0 = unname(cc[1]), c1 = unname(cc[2]), c2 = unname(cc[3]),
         c1_se = sqrt(vc[2, 2]), vc = list(vc), chi2_df = chi2_df, npts = length(M),
         flag = if (length(flag)) paste(flag, collapse = "; ") else NA_character_)
}

## ---- MEF scale k (MFI per MEF): slope-1 log-log fit, weighted by n -----------
# Free-slope fit is the linearity diagnostic: valid MEF units need slope 0.9-1.1
# and R^2 > 0.98, otherwise the lot's MEF scale does not fit this detector.
mef_scale <- function(d) {
  f <- d[which(!is.na(d$MEF) & !d$saturated & d$median > 0), ]
  if (nrow(f) < 3) return(tibble(k = NA_real_, ll_slope = NA_real_, ll_r2 = NA_real_, n_cal = nrow(f)))
  lk   <- weighted.mean(log10(f$median) - log10(f$MEF), f$n)
  free <- lm(log10(median) ~ log10(MEF), f)
  tibble(k = 10^lk, ll_slope = unname(coef(free)[2]),
         ll_r2 = summary(free)$r.squared, n_cal = nrow(f))
}

qb_chi2_max <- 50     # above this the fit is not valid -> Q/B/LoR suppressed
# Q resolvability: when the photon term c1 is swamped by the blank (c0) and bead
# CV0 (c2) terms on the dim peaks, c1 collapses toward 0 and Q = 1/c1 explodes
# (e.g. V4-A 2026-07-06 13:54 and 2026-11-05: Q x60-130, SE 3400-5900%). The fit
# itself is fine (chi2/df ~3-5) and 2+2SD is ~95% blank-driven, so only Q and the
# Q-derived B are suppressed; 2+2SD is kept.
qb_q_rse_max <- 1     # Q fit SE above 100% -> Q (and B) not resolved
mef_slope_ok <- c(0.9, 1.1); mef_r2_ok <- 0.98

## ---- Sensitivity for any detector (MFI headline + Spe) ----------------------
# Returns one row; Q (relative, gain-normalised) and B (= B_Spe) keep the column
# names the page's heatmap/FC/baseline helpers already use.
qb_spe <- function(d) {
  pf  <- parks_fit(d)
  nb  <- parks_fit(d[d$peak > 1, ])                     # robustness: blank excluded
  g   <- suppressWarnings(max(d$gain, na.rm = TRUE)); if (!is.finite(g)) g <- NA_real_
  ok  <- is.finite(pf$c1) && is.finite(pf$chi2_df) && pf$chi2_df <= qb_chi2_max
  c0p <- is.finite(pf$c0) && pf$c0 >= 0
  Q_rse <- if (ok) pf$c1_se / pf$c1 else NA_real_
  q_ok  <- ok && is.finite(Q_rse) && Q_rse <= qb_q_rse_max            # Q resolved?
  qflag <- c(if (!is.na(pf$flag)) pf$flag else NULL,
             if (is.finite(pf$chi2_df) && pf$chi2_df > qb_chi2_max) sprintf("chi2/df=%.0f", pf$chi2_df) else NULL,
             if (ok && !q_ok) sprintf("Q unresolved (SE %.0f%%): Q/B suppressed", 100 * Q_rse) else NULL)
  Q_I     <- if (q_ok) 1 / pf$c1 else NA_real_
  BSD_Spe <- if (q_ok && c0p) sqrt(pf$c0) / pf$c1 else NA_real_
  bsd_mfi <- { b <- d$sd[d$peak == 1]; if (length(b)) b[1] else NA_real_ }
  tibble(lor_mfi = if (ok) 4 * (pf$c1 + bsd_mfi) else NA_real_,                 # headline 2+2SD, MFI (kept when Q unresolved)
         bsd_mfi = bsd_mfi,
         photon_share = if (ok) pf$c1 / (pf$c1 + bsd_mfi) else NA_real_,         # 4/Q share
         Q_I = Q_I, Q_rse = Q_rse,
         Q   = Q_I * g,                                    # relative Q: Spe per gain-normalised MFI
         Q_noblank_ratio = if (q_ok && is.finite(nb$c1)) pf$c1 / nb$c1 else NA_real_,  # ~1 = robust
         B   = if (q_ok && c0p) pf$c0 / pf$c1^2 else NA_real_,                         # B_Spe
         BSD_Spe = BSD_Spe,
         lor_spe = 4 * (1 + BSD_Spe),                                               # 2+2SD, Spe
         CV0 = if (ok) sqrt(max(pf$c2, 0)) else NA_real_,
         chi2_df = pf$chi2_df, npts = pf$npts,
         qflag = if (length(qflag)) paste(qflag, collapse = "; ") else NA_character_)
}

## ---- MEF-unit sensitivity, MEF-mapped detectors -------------------------------
# d additionally carries MEF (NA for the blank); det = detector name (group_modify
# drops the grouping columns, so pass .y$detector).
qb_mef_one <- function(d, det = NA_character_) {
  sc <- mef_scale(d); pf <- parks_fit(d); k <- sc$k
  bl <- d[d$peak == 1, ]
  ok <- is.finite(k) && is.finite(pf$c1) && is.finite(pf$chi2_df) && pf$chi2_df <= qb_chi2_max
  mef_ok <- is.finite(sc$ll_slope) && sc$ll_slope >= mef_slope_ok[1] && sc$ll_slope <= mef_slope_ok[2] &&
            is.finite(sc$ll_r2) && sc$ll_r2 > mef_r2_ok
  q_ok <- ok && pf$c1_se / pf$c1 <= qb_q_rse_max                 # Q resolved?
  Qfit <- if (ok) k / pf$c1 else NA_real_                         # used for 2+2SD even if unresolved
  Q    <- if (q_ok) Qfit else NA_real_
  Q_se <- if (q_ok) k / pf$c1^2 * pf$c1_se else NA_real_
  c0p  <- q_ok && is.finite(pf$c0) && pf$c0 >= 0
  B    <- if (c0p) pf$c0 / (pf$c1 * k) else NA_real_
  BSD  <- if (ok && nrow(bl)) bl$sd[1] / k else NA_real_          # measured blank SD (MEF)
  BSDf <- if (ok && is.finite(pf$c0) && pf$c0 >= 0) sqrt(pf$c0) / k else NA_real_  # fitted cross-check
  lor  <- 4 * (1 / Qfit + BSD)
  fl <- c(if (!is.na(pf$flag)) pf$flag else NULL,
          if (is.finite(pf$chi2_df) && pf$chi2_df > qb_chi2_max) sprintf("chi2/df=%.0f: fit invalid", pf$chi2_df) else NULL,
          if (ok && !q_ok) sprintf("Q unresolved (SE %.0f%%): Q/B suppressed", 100 * pf$c1_se / pf$c1) else NULL,
          if (!mef_ok && is.finite(sc$ll_slope)) sprintf("log-log slope %.2f: MEF scale mismatched", sc$ll_slope) else NULL,
          if (isTRUE(grepl("^R[0-9]", det))) "red detector: dim-bead intrinsic CV (Hoffman & Wood 2007)" else NULL)
  tibble(k = k, ll_slope = sc$ll_slope, ll_r2 = sc$ll_r2, n_cal = sc$n_cal,
         Q = Q, Q_se = Q_se, B = B, BSD = BSD, BSD_fit = BSDf,
         lor = lor, lor_se = 4 * Q_se / Q^2, photon_share = (4 / Qfit) / lor,
         CV0 = if (ok) sqrt(max(pf$c2, 0)) else NA_real_, chi2_df = pf$chi2_df,
         mef_ok = mef_ok, mef_flag = if (length(fl)) paste(fl, collapse = "; ") else NA_character_)
}
