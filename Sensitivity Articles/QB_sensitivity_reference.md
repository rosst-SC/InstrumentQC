# Q/B Fluorescence-Sensitivity Reference for the InstrumentQC Aurora Dashboard

**Purpose.** A project reference that ties every methodological choice in the Q/B
sensitivity sections of `pages/AuroraEVO_8Peak.qmd` (SPHERO RCP-30-5A, 8-peak) and
`pages/AuroraEVO_6PeakUltra.qmd` (SPHERO URCP-38-2K, 6-peak) to a traceable source,
resolves the spectral/APD/auto-gain questions the classical PMT papers do not
address, and yields implementable validity rules and a prioritized list of
dashboard changes. Companion to `QB_rainbow_bead_sensitivity_bibliography.md`.

**Tags.** `[SRC]` = directly supported by a cited source (verified this session).
`[SYNTH]` = my synthesis across sources. `[PROJECT]` = a project-specific decision
or open question not settled by the literature. Verification marks (✔/⚠/✖) carry
over from the bibliography.

**Verified this session (read/mined):** Chase & Hoffman 1998; Hoffman & Wood 2007
(Unit 1.20); Wood & Hoffman 1998 overview; STN-17; BD Qr/Br app note (all as PDFs);
Parks 2017 (PMC5483398 full text); Cytek QBSURE TDS + patent abstracts; and DOIs/PMCs
for every other bibliography entry (see §7).

---

## 1. Unified model (one notation)

### 1.1 Symbols

| Symbol | Meaning | Units in this doc |
|---|---|---|
| `S` (or `f`, `MI`) | per-event signal intensity, **linear scale** | linear fluorescence, MEF/MESF, or MFI/gain |
| `CV` | robust coefficient of variation of a bead population | fraction (or %) |
| `SD` | robust standard deviation = `CV · S` | signal units |
| `Q` | detection efficiency = photoelectrons per unit signal | photoelectrons / MEF (or /S-unit) |
| `B` | optical background as an equivalent signal | MEF (or S-units) |
| `CV_intr` (`CV0`) | intrinsic CV: bead + illumination floor at high signal | fraction |
| `MEF` / `MESF` / `ERF` | molecules of equivalent (soluble) fluorophore / equiv. reference fluorophore | molecules |

### 1.2 Master equation `[SRC]`

The photoelectron-statistics model (Poisson conversion of photons to
photoelectrons) gives measured spread as a sum of an intrinsic floor, a
signal-limited (shot-noise) term, and a background term. Hoffman & Wood 2007
(Unit 1.20, p. 1.20.12–13) state it as:

```
CV_measured = sqrt[ (f + B) / (Q · f²) + CV_intr² ]
```

which expands to the exact form the dashboard fits `[SRC: Hoffman & Wood 2007]`:

```
CV² = CV_intr² + (1/Q)·(1/S) + (B/Q)·(1/S²)                     … (M)
```

Multiplying (M) by `S²` gives the **variance form** (Wood 1998; Chase 1998 Eq. 4;
Hoffman & Wood 2007):

```
SD² = CV_intr²·S² + (1/Q)·S + B/Q                                … (V)
```

If the intrinsic term is removed first (subtract the bright-bead CV in quadrature),
(V) becomes a straight line in `S`:

```
SD²_corrected = (1/Q)·S + (B/Q)      →  slope = 1/Q,  intercept = B/Q   … (L)
```

Chase 1998 Eq. 8/10 give the single-/two-bead closed forms: `Q = 1/(CV²·F)` and
`B = (SD_background / SD_particle)² · MESF_particle`.

### 1.3 The three equivalent forms actually in use

All three below are algebraically the same model; they differ in axis, weighting,
and how the intrinsic term is handled.

1. **Dashboard "CV²-vs-1/S"** (form M, intrinsic pre-subtracted, through origin):
   `qb_fit()` measures `CV_intr` as the brightest resolved peak's rCV, subtracts
   `CV_intr²`, then fits `cv2c ~ 0 + inv + I(inv^2)` with `inv = 1/S`. Linear coef
   = `1/Q`; quadratic coef = `B/Q`. `[SRC: Hoffman & Wood 2007 for the model + bright-bead intrinsic correction]`
2. **BD / Hoffman-Wood "SD²-vs-S"** (form L): linear regression of intrinsic-corrected
   `SD²` on `S`; `Q = 1/slope`, `B = intercept/slope = Q·intercept`.
   `[SRC: BD app note; Hoffman & Wood 2007 p. 1.20.13]`
3. **STN-17 three-point closed form**: blank + one dim bead + brightest bead;
   `Q = 1/(CV²_blank · MEF_blank)`, `SD² = MEF²·(CV² − CV²_bright)`,
   `B = (SD²_blank / SD²_dim)·MEF_blank`. `[SRC: STN-17 Fig. 2]`

**Verification that forms 1 and 2 are identical:** `SD² = S²·CV²`, so dividing (L)
by `S²` returns (M). The handoff note "variance-vs-mean is WRONG — gives negative B"
is therefore **not** a property of form (V)/(L) itself — Hoffman & Wood, BD, Parks
and STN-17 all fit the variance form successfully. It is a symptom of **unweighted,
few-point** fitting, which Hoffman & Wood 2007 (p. 1.20.16) explicitly warn is
"heavily influenced by the data highest on the scale" and inflates the error in B.
See §3.4 and §6. `[SYNTH]`

**Dashboard implication.** The dashboard's chosen form (M) is a legitimate,
source-backed rearrangement of the canonical model. The one substantive divergence
from the canonical estimator is that the dashboard *measures and pre-subtracts*
`CV_intr` and fits **two** free parameters, whereas Parks 2017 fits **three**
(`Q, B, CV0`) simultaneously; that trade-off is analysed in §2 and §3.4.

---

## 2. Method-comparison matrix

| Feature | Dashboard **relative** (`xvar=mfi_norm`) | Dashboard **MEF** (`xvar=MEF`) | STN-17 closed form | BD Qr/Br (CS&T) | Parks 2017 (weighted quadratic) | Patrone 2025 (multi-gain) | Cytek QBSURE (SpectroFlo) |
|---|---|---|---|---|---|---|---|
| Beads | RCP-30-5A / URCP-38-2K, all peaks | same, anchor detectors only | RCP-30-5A: blank + dim + bright | CS&T beads (bright/mid/dim) | multi-level beads and/or LED pulses | multi-level, swept across gains | QBSURE beads |
| x-axis unit | MFI/gain (**arbitrary, gain-stabilised**) | assigned **MEF** | assigned/extrapolated **MEF** | **ABD** units (≈ anti-CD4 Ab) | arbitrary intensity **or** MEF | intensity across gains | **MEFL** |
| # populations fitted | resolved fluorescent peaks (2–8 / 2–6) | same | 3 | 3 | ≥3 (or many LED levels) | many (levels × gains) | bright + mid + dim + laser-off blank |
| Intrinsic CV | **measured** (brightest peak), pre-subtracted | same | subtract `CV²_bright` in SD² | corrected from bead-lot file | **free 3rd parameter (CV0)** | modelled | corrected from lot |
| Fit / estimator | weighted `CV²`-vs-`1/S`, through origin | same | algebraic 3-point | `SD²`-vs-median linear | **simultaneous** weighted quadratic, iterated | global fit across gains | `Q=1/(MEFL·CV²)`; `b` from SD ratio |
| Weighting | `w = 1/CV²²` (∝ `1/CV⁴`) `[PROJECT]` | same | none (algebraic) | statistical (rSD/rCV, lot-corrected) | **weights ∝ 1/point-variance, iterated** `[SRC]` | proper, gain-aware | vendor |
| B source | quadratic coef `B/Q` (curvature) | `B/Q`, `true_b` → `B` in MEF | blank/dim SD ratio → **B in MEF** | intercept → `B = Q·intercept`, ABD | intercept term `c0 = B/Q²` | **explicit gain-independent B** | laser-off blank SD |
| Gain treatment | divide MFI by daily-QC gain `[PROJECT]` | MEF is gain-independent | MEF is gain-independent | fixed PMT V per config | **exploits multiple gains** | **core idea: multi-gain** | fixed per QC target |
| Validity check | hard guards + advisory gates (§3) | same + STN-17 cross-check | needs `CV_blank,CV_dim > CV_bright` | R²; lot QC | **std errors + per-peak weighted residuals** `[SRC]` | stable B across gains | pass/fail vs spec |
| Portability | single-detector trend only | cross-day (per detector) | cross-instrument (if MEF valid) | within-config; not cross-instrument | cross-site (demonstrated, Parks 2018) | designed for inter-instrument | within Cytek config |

**Dashboard implication.** The dashboard already implements three columns of this
matrix (relative, MEF, STN-17). The two columns it does **not** implement — Parks
2017 (simultaneous weighted quadratic with residual diagnostics) and Patrone 2025
(multi-gain) — are the two most defensible upgrades and are exactly where the Aurora's
auto-gain behaviour (§5) could be turned from a nuisance into an asset.

---

## 3. Validity & decision rules (implementable)

This section is written as a checklist the dashboard can enforce. Thresholds and
their sources are explicit; where the dashboard already implements a rule, the code
symbol is named.

### 3.1 Preconditions (source-mandated) `[SRC: Hoffman & Wood 2007; Chase 1998]`

- **Linear scale only.** The model assumes photoelectron statistics on a linearly
  amplified signal. Log-amp data "are usually not accurate enough to estimate B
  accurately" (Hoffman & Wood 2007, p. 1.20.6). Aurora MFI + robust CV are linear →
  **satisfied.** Keep any future gain/intensity handling on the linear scale.
- **Robust statistics.** Hoffman & Wood 2007 (p. 1.20.17) explicitly recommend
  percentile-based robust CV/SD. Dashboard uses robust CV (`rcv`) and IQR/1.349 for
  the blank → **satisfied.**
- **Singlets only, adequate events.** 1000–5000 events/population, singlet-gated,
  compensation/unmixing **off** for the raw-channel CV. The dashboard reads
  per-detector raw peak metrics from `Process*.R`; confirm the upstream gate excludes
  multiplets (a multiplet inflates CV → distorts B). `[SRC: Hoffman & Wood 2007 p.1.20.8]`

### 3.2 Hard suppression rules (fit is ill-conditioned → Q/B = NA) `[PROJECT, motivated by SRC]`

`Q = 1/slope` blows up as slope → 0, so a noisy run can emit a phantom Q spike (the
documented 2026-05-29 red-laser fault, `R4-A` read Q≈35). `qb_fit()` suppresses when
**any** of:

| Guard | Threshold | Basis |
|---|---|---|
| flat/negative slope | `s1 ≤ 0` or non-finite | `Q=1/slope` undefined `[SYNTH]` |
| poor linearity | `R² < 0.8` (`min_r2`) | model demands linear `CV²` vs `1/S` `[SYNTH]` |
| non-monotonic CV | Spearman(`1/S`,`CV²`) `< 0.3` (`min_mono`) | `CV²` must rise with `1/S` `[SYNTH]` |
| inflated bright CV | brightest-peak CV `> 15%` (`max_bright_cv`) | bright bead is the intrinsic floor; if it is that broad, the run is faulty `[SYNTH]` |

These are engineering guards, not literature thresholds; they are defensible but the
specific cutoffs (0.8, 0.3, 15%) are `[PROJECT]` choices. **Recommendation:** keep,
but expose them as documented parameters and note they are project-tuned.

### 3.3 Advisory (photon-domination) gates — flag, do not suppress `[SRC: Chase 1998; Hoffman & Wood 2007]`

For the shot-noise term to dominate (so `Q = 1/(CV²·F)` is valid), the dim bead's
CV must be dominated by photoelectron statistics, not background or intrinsic spread.
The literature criteria are:

- **SD ratio:** dim-bead `SD ≥ 3× blank SD` (background broadening negligible).
  Chase 1998, p. 272; Hoffman & Wood 2007. Dashboard: `sd_dim < 3*sd_blk` → advisory.
- **CV ratio:** dim-bead `CV ≥ 3× brightest-bead CV` **if uncorrected**, or **≥ 2×**
  **if the CV is intrinsic-corrected**. Chase 1998 uses 3×; Hoffman & Wood 2007
  (p. 1.20.13) state the relaxed factor: intrinsic CV must be "a factor of three
  lower… if no correction is used, or a factor of two lower if the measured CV is
  corrected with the intrinsic CV." The dashboard intrinsic-corrects, so its `2×`
  (`rcv_dim < 2*rcv_brt`) is **correct and source-backed.** `[SRC]`

**Why advisory not hard:** on spectral rainbow beads the in-range peaks often sit in
the bright/intrinsic-dominated regime and the "blank" is itself a dim bead, so these
are routinely unmet. That is a caveat to surface, not a fit failure. `[SYNTH]` This
matches Hoffman & Wood 2007's own warning that Rainbow dim beads have too high an
intrinsic CV to give valid Q/B on **red** detectors (§3.5).

### 3.4 The weighting / negative-B rule `[SRC: Hoffman & Wood 2007; Parks 2017]`

Hoffman & Wood 2007 (p. 1.20.16, "Sources of error"): *"The effect of an error in
estimated CV on the error in the calculation of B can be as much as ten times as
great as the effect on Q… particularly true when simple, unweighted linear fits to
the data are used."* This is the literature basis for the dashboard's inverse-variance
weighting (`w = 1/CV²²`, since `Var(CV²) ∝ CV⁴`). Parks 2017 formalises it: weights
∝ 1/point-variance, **iterated** to convergence, fitting `Q, B, CV0` simultaneously.

**Decision rule:** unweighted fits over few (3–5) points are the documented cause of
unstable/negative B. The dashboard's weighting mitigates this; adopting Parks' full
simultaneous 3-parameter weighted fit with reported standard errors would resolve it
properly (§6, Rec 3).

### 3.5 Per-channel validity (which detectors can yield valid absolute Q/B) `[SRC: Hoffman & Wood 2007; bibliography §"Known limitation"]`

- **Red / far-red detectors:** RCP-30-5A dim beads carry too high an intrinsic CV for
  valid Q/B on red-laser channels (Hoffman & Wood 2007 explicitly: works for FITC/PE,
  fails for red). This gates the **APC (R1-A)** and **APC-Cy7 (R7-A)** MEF anchors on
  the 8-peak page and the red anchors (`R1-A`, `R4-A`, `R7-A`) on the 6-peak page.
- **URCP-38-2K** spans a wider intensity range than RCP-30-5A, so it is *less* affected,
  but the same spectral-mismatch caution applies to hard-dyed beads vs real fluorophores.
- **Cross-excited tandems** (e.g. PE-TexasRed/PE-Cy5 at 488 ex) are already excluded from
  the 6-peak `mef_map` — correct.

**Rule:** mark red/far-red absolute Q/B (MEF & STN-17) as **lower-confidence** in the
UI; keep them out of any pass/fail gate. Relative single-detector trending is still fine
on those channels (it does not claim absolute accuracy).

### 3.6 Levey-Jennings control-limit & bead-lot rules `[SRC: BD app note; Wang & Hoffman 2017]`

- Relative Q/B is valid for **single-detector L-J trending only**; a `±2 SD`
  control band per detector is appropriate once ~8–10 qualified runs exist
  (dashboard's `clim()` + the locked `fc_ref_window`). `[SRC: BD; STN-17 "Levey-Jennings … monitor performance over time"]`
- **Bead-lot change:** expect a small step in relative Q on lot change, because
  non-normalised units are not lot-corrected — "a small change in the Levey-Jennings
  Qr tracking chart may be seen when a new bead lot is used" (BD app note). **Rule:**
  annotate lot changes on L-J charts and re-baseline the locked reference after a lot
  change. `[SRC: BD]`
- **A 50% change in non-normalised Q = a 50% change in normalised Q** (BD): so
  percentage-change aggregation across detectors is valid even though absolute
  relative-Q values are not comparable across detectors — this is exactly why the
  dashboard aggregates `Q_d` (% vs baseline) per laser rather than median absolute Q.

**Dashboard implication.** §3.1–3.3 are already implemented and now fully sourced.
The actionable gaps are: (a) surface red/far-red as lower-confidence (§3.5), (b) move
from 2-param pre-subtracted to Parks' 3-param weighted fit (§3.4), (c) formalise
lot-change annotation (§3.6).

---

## 4. Vendor ↔ academic reconciliation (Cytek-centred)

### 4.1 Cytek SpectroFlo / QBSURE — what the Aurora actually reports `[SRC: Cytek QBSURE TDS B9-50002; patents US 11,131,618 / 11,879,826 / 12,366,516]` ⚠ *(vendor/patent, not peer-reviewed)*

Per **raw detector** (per APD channel), for each dye parameter marked on the QC bead
for the given configuration:

```
Q = 1 / (MEFL · CV²)                     … detection efficiency, photoelectrons/MEFL
b = (SD_blank / SD_ModBright)² · MEFL     … background, MEFL
R                                         … resolution limit (derived from Q, b)
```

- `MEFL` = molecules of equivalent fluorescein (the NIST fluorescein/ERF scale).
- `SD_blank` = SD with **laser off** (electronic/optical dark background);
  `SD_ModBright` = SD with **modulated laser** on the bright bead.
- **These map 1:1 onto the dashboard's STN-17 closed form:** Cytek's
  `Q = 1/(MEFL·CV²)` is Chase Eq. 8 / STN-17 `Q = 1/(CV²_blank·MEF_blank)`; Cytek's
  `b = (SD_blank/SD_bright)²·MEFL` is Chase Eq. 10 / STN-17
  `B = (SD²_blank/SD²_dim)·MEF_blank`. `[SYNTH]` The only structural difference is
  Cytek's *laser-off* blank vs the dashboard's *blank-bead* population as the
  background reference — see §4.4.

### 4.2 BD Qr/Br (secondary reference) `[SRC: BD app note 2012]`

- `Q = 1/(Mean[MESF] · CV²)`; `Q = 1/slope` of `SD²`-vs-median, `B = Q · intercept`.
- Reported in **ABD units** (Assigned BD): 1 ABD ≈ fluorescence of one anti-CD4 Ab
  (≈40,000 Ab/CD4⁺ lymphocyte anchor). `Qr, Br` = Q, B in ABD units.
- **Non-qualified configs** (SORP/custom filters) are **not** normalised to ABD and
  **not** lot-corrected, but "these non-normalised Qr values are still valid and serve
  as a powerful metric for tracking changes in detector performance over time" —
  the direct precedent for the dashboard's relative Q.
- BD caveat: `Q`/`Qr` for equivalent detectors "can vary as much as two- to three-fold
  between two optimized flow cytometers" (also Chase 1998: PMT photocathode sensitivity
  varies ~3×). → absolute Q is **not** cross-instrument portable without spectrally
  matched standards.

### 4.3 Quantitative mapping to MEF-unit academic Q/B `[SYNTH]`

- **Cytek Q (MEFL) ↔ academic Q (MESF):** same definition (photoelectrons per
  equivalent-fluorophore molecule); the unit differs only in the reference fluorophore
  (fluorescein/ERF vs the target dye's MESF). To compare a dashboard MEF-anchor Q
  (e.g. FITC on B2-A) with SpectroFlo's Q, both must use the *same* reference scale;
  the dashboard's per-anchor MEF (lot AS01 / URCP certificate) is dye-specific MEF,
  Cytek's is MEFL — a fixed per-dye ratio, not 1:1. **Flag:** do not compare the
  dashboard's PE-MEF Q against SpectroFlo's MEFL Q numerically without converting.
- **BD ABD ↔ MESF:** `Qr` (ABD) and `Q` (MESF) differ by the per-dye ABD↔MESF factor
  in the CS&T lot file; only *relative changes* (%) are directly comparable across the
  two (BD's "50% = 50%" rule).
- **ERF (NIST) ↔ MESF:** ERF is the cross-instrument-portable unit (Wang & Gaigalas
  2011; Hoffman/NIST-ISAC 2012); MESF is dye-specific. URCP-38-2K lots are increasingly
  ERF-certified. `[SRC]`

### 4.4 Definition mismatches that make numbers differ `[SYNTH]`

1. **Background reference:** Cytek uses **laser-off** SD (dark background);
   dashboard/STN-17 use the **blank-bead** population (which includes scatter/Raman/
   autofluorescence of the bead). These are *different B's* — dashboard B ≥ Cytek b in
   general. Do not expect the dashboard's blank-derived B to equal SpectroFlo's b.
2. **Unit of x:** MEFL (Cytek) vs dye-MEF (dashboard) vs ABD (BD) vs arbitrary MFI/gain
   (dashboard relative). Only same-unit comparisons are numeric; cross-unit comparisons
   are relative-only.
3. **Which term:** the dashboard's relative-fit "B" heatmap is the raw `I(inv²)`
   coefficient = **B/Q**, *not* B. A change there can be background *or* efficiency.
   The blank-derived **B in MEF** (STN-17 panel) is the canonical B. The dashboard
   already labels this correctly — keep that distinction prominent.
4. **Simultaneous vs pre-subtracted intrinsic:** vendor tools and Parks fit CV0/intrinsic
   differently than the dashboard's measured-and-subtracted approach (§3.4).

**Dashboard implication.** The dashboard's absolute (MEF/STN-17) Q/B is
*definitionally* the Cytek/Hoffman quantity, so co-trending against SpectroFlo's
per-detector Q/B/R is a valid cross-check — **but only as trends/percentages**, not
absolute values, because of the unit and background-reference mismatches above. Add a
note to the MEF panels stating this explicitly.

---

## 5. Spectral / APD / auto-gain open questions

### 5.1 APD vs PMT photoelectron statistics — is "Q" interpretable identically? `[SRC: Lawrence 2008; Hoffman & Wood 2007; de Rond 2021]`

- The Q/B model only assumes **Poisson photon→photoelectron conversion on a linear
  scale** (Chase 1998 Eq. 1–4; Wood 1998). It is detector-agnostic in principle:
  de Rond 2021 transported the identical Q/B/R framework to a **scatter** detector,
  proving the model travels beyond PMT fluorescence channels. `[SRC]`
- **But APDs differ from PMTs** in ways that change the *value* and *interpretation*
  of Q: APDs have much higher quantum efficiency in the **red/far-red/NIR** (Lawrence
  2008 — the reason Aurora/ID7000 use APD arrays), and APDs carry **excess noise
  (multiplication/excess-noise factor F)** plus higher dark current than PMTs. The
  simple Poisson `SDe = √n` (Chase Eq. 1) becomes `√(F·n)` for an APD. `[SYNTH]`
- **Consequence for the dashboard:** a per-APD "Q" is still a valid *relative*
  detection-efficiency index (higher = more signal per event), and trends are
  meaningful. But the *absolute* photoelectron interpretation (`Q` = photoelectrons/MEF)
  is only approximate for APDs unless the excess-noise factor is accounted for — which
  neither the dashboard nor STN-17 does. **Open question `[PROJECT]`:** whether to
  present absolute Q on APDs at all, or only relative Q + MEF detection limit.

### 5.2 The "no Q for an unmixed parameter" structural problem `[SRC: bibliography §4; SYNTH]`

On a spectral instrument, Q and B are naturally defined **per raw detector (per APD
channel)**, because unmixing is a linear recombination that mixes the per-channel
photoelectron scales. A single "Q of an unmixed parameter" is not well defined the way
it is for a filter-based PMT channel. This is why:

- The dashboard correctly computes Q/B **per raw detector** (all 64), never per unmixed
  fluorophore — this is the right design and should be stated as a deliberate choice.
- Unmixed-space resolution is characterised by **stain index, SSM/SSE** (Nguyen 2013)
  and gain-independent spread (SQI, Bhowmick 2021, validated on the Aurora), **not** by
  a per-parameter Q/B. The dashboard's SI panels cover part of this; SSM/SQI are a
  possible future addition. `[SRC]`

### 5.3 Can the Aurora's daily gain variation be exploited as a Patrone-style multi-gain series? `[SRC: Patrone 2025; Parks 2017; Hoffman & Wood 2007 Fig 1.20.2]`

**The idea.** Patrone 2025 proposes estimating Q/B by combining measurements across
**multiple gains**, explicitly modelling the (often dominant) gain-independent
background, which stabilises B. The Aurora *already* varies each APD's gain daily (to
pin bead MFI to target). Over many QC days we therefore have, per detector, a series of
(gain, MFI, CV) points at **different gains** — potentially a free multi-gain dataset.

**Feasibility assessment `[SYNTH]`:**

- *Supporting physics:* Hoffman & Wood 2007 (Fig 1.20.2) show `Q ∝ laser power` and
  that varying illumination while re-pinning gain changes CV predictably — the same
  lever Patrone formalises. Q and B are gain-invariant in CV-space, so a multi-gain
  series is analysable.
- *What we'd gain:* a gain-independent B estimate (the dashboard's weakest quantity),
  without needing MESF beads for the *relative* series — Parks 2017 confirms the
  quadratic model "operates in arbitrary intensity units without requiring MESF
  calibration." `[SRC]`
- *Obstacles specific to this dashboard:*
  1. **Confounded, not designed, gain steps.** Daily gain changes are small and driven
     by QC targeting, not a deliberate sweep — limited gain range → weak leverage.
  2. **One gain value per detector per day** (from `ArchivedDataEVO.csv`), and gain is
     recorded once/day; the multi-gain series is *across days*, so detector aging/optical
     drift confounds the gain axis. Patrone assumes same-session multi-gain.
  3. **Sparse runs** (3–4 so far) — far short of what a stable global fit needs.
- *Verdict:* **not feasible now as a primary method**, but worth prototyping later as a
  *secondary* B-stabilisation once ≥10–15 runs exist, treating gain as a covariate with
  an explicit aging term. A *deliberate* single-session multi-gain acquisition (run the
  bead at 3–4 gains in one sitting) would be the clean way to realise Patrone's method
  on the Aurora. `[PROJECT]`

**Dashboard implication.** Keep per-raw-detector Q/B (§5.2). Present APD absolute Q
cautiously (§5.1). Log the multi-gain approach as a documented future experiment
(§6, Rec 5), not a current feature.

---

## 6. Recommendations for the dashboard (prioritized)

Documentation only — no code changed this pass. Effort = rough dev size; Impact =
effect on QC decision quality.

| # | Recommendation | Basis | Effort | Impact |
|---|---|---|---|---|
| 1 | **Make relative single-detector Q the primary displayed metric; demote absolute MEF/STN-17 Q/B to a gated "lower-confidence" panel** (already largely done — finish the framing). Relative Q is the only cross-day-valid quantity given daily re-gaining. | BD app note (non-normalised Qr valid for tracking); §3.6 | Low | High |
| 2 | **Flag red/far-red anchors (APC R1-A, APC-Cy7 R7-A, Alexa700 R4-A) as lower-confidence** on every absolute Q/B/MEF panel; exclude from any pass/fail gate. | Hoffman & Wood 2007 (Rainbow dim beads invalid on red); §3.5 | Low | High |
| 3 | **Replace the 2-param pre-subtracted fit with Parks 2017's simultaneous weighted quadratic** (`Q, B, CV0` fitted together, iterated weights, report standard errors + per-peak weighted residuals). Resolves the negative/unstable-B problem properly and gives an objective validity diagnostic. | Parks 2017 (PMC5483398); Hoffman & Wood 2007 p.1.20.16 | Med | High |
| 4 | **State the per-raw-detector design explicitly** and add a note that no per-unmixed-parameter Q is computed (by design); optionally add SSM/SQI for unmixed-space spread. | Bibliography §4; Nguyen 2013; Bhowmick 2021 | Low (note) / Med (SSM) | Med |
| 5 | **Log a future multi-gain B-stabilisation experiment** (Patrone-style): deliberate single-session 3–4-gain bead acquisition; do **not** rely on across-day QC gain steps. | Patrone 2025; §5.3 | Med (experiment) | Med |
| 6 | **Annotate bead-lot changes on all L-J/Q trend charts** and re-baseline the locked FC reference after a lot change. | BD app note (lot step); §3.6 | Low | Med |
| 7 | **Add a units/definition banner to the MEF & STN-17 panels** clarifying that dashboard blank-derived B ≠ SpectroFlo laser-off b, and cross-vendor comparisons are trend-only. | §4.4 | Low | Med |
| 8 | **Adopt ERF-certified bead lots where available** (URCP) to gain cross-instrument-portable absolute units; keep MESF for dye-specific work. | Wang & Gaigalas 2011; Hoffman/NIST-ISAC 2012 | Low (procurement) | Med |
| 9 | **Expose the hard-guard thresholds** (`min_r2=0.8`, `min_mono=0.3`, `max_bright_cv=15`) as documented, tunable parameters and label them project-tuned, not literature values. | §3.2 | Low | Low |

**What beads would actually enable valid absolute Q/B here `[SYNTH]`:** for red/far-red
channels, spectrally appropriate dim references (surface-labelled APC standards or
ERF-certified Ultra Rainbow lots) or an LED-pulse source (Parks 2017) — **not** the
8-peak RCP-30-5A dim beads. This is the single biggest limiter on absolute-sensitivity
validity on the Aurora's red detectors.

---

## 7. Bibliography (each: the specific claim it authoritatively supports)

**Foundational model**
- **Wood & Hoffman 1998, Cytometry 33:256–259** ✔ (PDF) — sensitivity is jointly governed
  by Q and B; peak-channel/CV specs alone are inadequate.
- **Wood 1998, Cytometry 33:260–266** ✔ (doi:10.1002/(SICI)1097-0320(19981001)33:2<260) —
  the photoelectron derivation; `CV=√(f+B)/√(Q·f)`; a performance level = a family of (Q,B) pairs.
- **Chase & Hoffman 1998, Cytometry 33:267–279** ✔ (PDF) — dim-bead protocol; Eq. 8
  `Q=1/(CV²F)`, Eq. 10 `B=(SD_bg/SD_particle)²·MESF`; **3× SD-ratio & 3× CV-ratio**
  photon-domination criteria; quadrature intrinsic subtraction; linear-scale preference.
- **Steen 1992, Cytometry 13:822–830** ⚠ (doi:10.1002/cyto.990130804; **PMID 1458999** —
  bibliography's 1451336 is wrong) — LED-pulse photoelectron characterisation; precursor.
- **Hoffman & Wood 2007, Curr Protoc Cytom Unit 1.20** ✔ (PDF) — canonical protocol;
  `SD²=f/Q+B/Q` (slope 1/Q, intercept B/Q); `SD²=MESF²(CV²_dim−CV²_bright)`; **2× relaxed
  criterion when intrinsic-corrected**; **B error ≤10× Q error when unweighted**; robust
  CV recommended; **Rainbow dim beads invalid on red detectors**; `Q ∝ laser power`.
- **Wood 2021, Cytometry A 99:664–667** ✔ (doi:10.1002/cyto.a.24281) — Poisson (<100
  photons) vs Gaussian regimes; modern restatement of detection efficiency.

**Bead calibration / units**
- **Hoffman/NIST-ISAC 2012, Cytometry A 81A:785–796** ✔ (doi:10.1002/cyto.a.22086) —
  MESF/ERF assignment varies by manufacturer; hard-dyed beads not spectrally matched →
  verify per instrument model.
- **Wang & Gaigalas 2011, J Res NIST 116:671–683** ✔ (nvlpubs V116.N03.A03) — defines ERF
  as the cross-instrument-portable unit tied to NIST fluorescein SRM 1932.
- **Schwartz et al. 2002, J Res NIST 107:83–91** ✔ (PMID 27446720) — formal MESF definition
  and assignment methodology.

**Validity / improved estimators**
- **Parks et al. 2017, Cytometry A 91A:232–249** ✔ (doi:10.1002/cyto.a.23052; PMC5483398) —
  automated **simultaneous weighted-quadratic** fit of `Q,B,CV0`; `V=B+M+CV0²M²`; iterated
  1/variance weights; std errors + per-peak weighted residuals; works in **arbitrary units**.
- **Parks et al. 2018, Cytometry A 93:1087–1091** ✔ (doi:10.1002/cyto.a.23605) — 23-instrument
  multisite validation; the bead/LED Q-scale method is portable across sites.
- **Patrone et al. 2025, Cytometry A** ✔ (doi:10.1002/cyto.a.24955; PDF) — "beyond Q and B":
  **multi-gain global fit** modelling gain-independent background → stable B.
- **Nguyen et al. 2013, Cytometry A 83A:306–315** ✔ (doi:10.1002/cyto.a.22251) — spillover
  spreading matrix (SSM); post-compensation spread ∝ √fluorescence; the experimental-background
  term feeding B in real panels.
- **Bhowmick et al. 2021, Sci Rep 11:20553** ✔ (doi:10.1038/s41598-021-99831-7) — gain/dynamic-range
  **independent** spread index (SQI); validated on the **Cytek Aurora**; raw SSM is sensitivity-dependent.

**Spectral / APD**
- **Nolan et al. 2013, Cytometry A 83A:253–264** ✔ (doi:10.1002/cyto.a.22241; PMC3594514) —
  bead-based quantitative performance assessment transfers to a spectral platform;
  virtual-bandpass vs classic-least-squares unmixing.
- **Lawrence et al. 2008, Cytometry A 73A:767–776** ✔ (doi:10.1002/cyto.a.20595; PDF) —
  APD vs PMT: APDs' higher QE in red/far-red raises Q where PMT Q collapses (the Aurora APD rationale).
- **de Rond et al. 2021, Cytometry A 99:671–679** ✔ (doi:10.1002/cyto.a.24243; PMC8359315) —
  derives Q/B/R for a **scatter** detector via the same photoelectron framework → model travels
  to non-canonical channels.
- **Ortyn et al. 2006, Cytometry A 69A:852–862** ✔ (doi:10.1002/cyto.a.20306) — sensitivity &
  compensation in a spectral **imaging** cytometer.

**Standardization / QC**
- **Perfetto et al. 2014, Cytometry A 85A:1037–1048** ✔ (PDF) — total effective background =
  instrument (optical+electronic) + experiment (Raman, spillover spreading); why Q/B matters
  for panel design & inter-instrument transfer.
- **Wang & Hoffman 2017, Curr Protoc Cytom Unit 1.3** ✔ (doi:10.1002/cpcy.14) — standardization/
  calibration/control; surface-labelled vs hard-dyed beads; MESF/ERF; L-J tracking.

**Vendor / standards documents**
- **STN-17 (Spherotech)** ✔ (PDF) — RCP-30-5A blank + dim + bright → `Q=1/(CV²·MEF_blank)`,
  `B=(SD²_blank/SD²_mid)·MEF_blank`; the 8-peak method the dashboard implements.
- **STN-8 / STN-9 / STN-14 (Spherotech)** ✔ (spherotech.com) — RCP/URCP identities, 256-relative-
  channel calibration, MEF/MEFL/MEPE assignment and cross-calibration; L-J tracking of %CV/gain.
- **Ultra Rainbow Quantitative Particle Kit (Spherotech)** ✔ — URCP MEF/ERF provenance for the 6-peak page.
- **BD "Qr and Br in BD FACSDiva" app note 2012** ✔ (PDF) — `Q=1/(Mean[MESF]·CV²)`; ABD units;
  non-normalised Qr valid for tracking ("50% = 50%"); lot-change step; `Q` varies 2–3× across identical instruments.
- **Cytek QBSURE / SpectroFlo (TDS B9-50002; patents US 11,131,618 / 11,879,826 / 12,366,516)** ⚠ —
  per-raw-detector `Q=1/(MEFL·CV²)`, `b=(SD_blank/SD_ModBright)²·MEFL`, resolution limit R.

---

*Prepared for InstrumentQC, July 2026. Cross-checked against `qb_fit()`, the `stn17`
block, and `mef_map` in `pages/AuroraEVO_8Peak.qmd` and `pages/AuroraEVO_6PeakUltra.qmd`.*
