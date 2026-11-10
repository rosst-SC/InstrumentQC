> **Copied 2026-11-11 from the Sensitivity project** (`~/Documents/Positron_Local/Sensitvity/`).
> This is the method now used by `pages/_qb_parks.R` (EVO 6-peak Ultra page + dashboard).
> Where it disagrees with `QB_sensitivity_reference.md`, this file wins.
> The "Implementation in this folder" section describes the Sensitivity project, not this repo.

# Calculating Q, B and the 2+2SD limit of resolution from multi-peak beads

A practical summary for future projects. It distils what was learned while comparing four cytometers
(A8, Aurora Evo, NovoCyte Opteon, Attune Xenith) with SPHERO RCP-30-5A 8-peak beads (lot AS01), September–October 2026.
The working implementation is in this folder (see [Implementation](#7-implementation-in-this-folder)).

Where this summary disagrees with `QB_sensitivity_reference.md`, this summary is the corrected version. The disagreements are:
- the claimed STN-17 = Chase = Cytek equivalence;
- pre-subtracting the brightest-peak CV;
- the wording about red detectors.

---

## 1. What the three numbers mean

| Quantity | Meaning | Units | Better |
|---|---|---|---|
| **Q** | Detection efficiency: statistical photoelectrons (Spe) collected per unit of dye | Spe per MEF | higher |
| **B** | Background: optical plus electronic, expressed as an equivalent amount of dye | MEF | lower |
| **BSD** | Background SD, the noise floor a dim stain must rise above. In the model, BSD = √(B/Q) | MEF | lower |
| **CV₀** | Intensity-proportional spread: bead intrinsic CV plus illumination and fluidics | % | lower |
| **2+2SD** | **Limit of resolution:** the smallest dye amount that separates cleanly from background. It combines Q and BSD | MEF | lower |

**Use 2+2SD as the headline metric.** Parks et al. (2018) chose it as the most informative single number for comparing instruments.

- Q and B on their own can mislead. In this study the Xenith's blank SD was second-best on 6 channels, yet its very low Q put it last on 2+2SD.
- Every value is in the MEF units of one fluorochrome. Compare the same fluorochrome across instruments, never across fluorochromes.

---

## 2. The model

Photoelectron statistics (Wood 1998) give each bead population a variance, in instrument units (MFI), of

```
V = c0 + c1·M + c2·M²
```

where `M` is the population's median.

| Term | Source | Grows with |
|---|---|---|
| `c0` | Background variance (optical + electronic) | constant |
| `c1·M` | Poisson photon (shot) noise | signal |
| `c2·M²` | Bead intrinsic CV + illumination (CV₀) | signal² |

### Converting coefficients to physical quantities

`k` = MFI per MEF, from the bead lot's MEF values (step 5 below).

```
Q   (Spe per MEF)  = k / c1
B   (MEF)          = c0 / (c1 · k)
BSD (MEF)          = √c0 / k                = √(B / Q)
CV0                = √c2
```

The same quantities in photoelectron units, as flowQB reports them, are:

```
Q_I   (Spe per MFI)  = 1 / c1
B_Spe (Spe)          = c0 / c1²
```

Both conventions are correct; they differ only in the units of B.
- With B in MEF (signal units): `c0 = B/Q`, using Q per MFI.
- With B in Spe: `c0 = B_Spe/Q_I²`.

### The 2+2SD limit of resolution (Parks et al. 2018, from Maecker et al. 2004 and Chase & Hoffman 1998)

Let the background SD be σ, in Spe. Adding S Spe of dye gives a population with SD √(σ² + S). It is "resolved" when its −2 SD point reaches the background's +2 SD point:

```
S − 2·√(σ² + S) = 2σ   ⇒   S = 4·(1 + σ)          (in Spe)

2+2SD (MEF) = 4 · (1/Q + BSD_MEF)                  (divide by Q)
```

- When Q is high, the limit approaches 4 × BSD, so it is **background-limited**.
- When Q is low, the 4/Q term dominates, so it is **photon-limited**. Report the 4/Q share to show which applies.
- CV₀ is ignored, as in Parks et al. At these low intensities its contribution is negligible.

---

## 3. Acquisition requirements

- **Linear data only.** Log-amplified data are not accurate enough (Hoffman & Wood 2007).
- **Singlet beads**, gated on scatter. Multiplets inflate SDs and so distort B.
- **Enough events per peak.** About 600–1,400 per peak worked here; Hoffman & Wood suggest 1,000–5,000.
- **No saturated peaks.** Keep the brightest peak below about 90% of full scale.
- **Record the settings**: gains or voltages, flow rate, and whether they are the instrument's QC defaults. In this study every instrument ran at its QC defaults, which shows routine performance; Parks et al. 2018 did the same.
- **Watch scale placement.** Detectors whose brightest bead sits below about 10% of full scale picked up extra, intensity-dependent noise. The Xenith, at 4–13% of scale and 409–673 V, showed this: its model fits failed (χ²/df ≈ 100). A voltage titration shows how much sensitivity the default settings leave unused.

---

## 4. Procedure, step by step

1. **Read the FCS files** with `read.FCS(transformation = FALSE, truncate_max_range = FALSE)`. Keep one raw `-A` detector per dye.
2. **Assign peaks on several detectors at once.**
   - Use `flowClust` with K = 8 on asinh-transformed data, giving each detector its own cofactor (about the SD of its dim events).
   - Set `min.count = max.count = −1`. Otherwise flowClust silently drops events piled at the scale limits, such as the Xenith's values clipped at 0.
   - Check the result: about 12.5% of events per peak, and medians that rise from peak 1 to 8 in every detector. Fall back to k-means on the summed asinh if needed.
3. **Calculate per-peak statistics**: n, median, and SD = IQR/1.349.
   - **For the blank, use the right half: SD = P84.13 − P50** (Chase & Hoffman 1998). This is robust when the scale is clipped at 0.
   - Exclude saturated peaks (more than 0.5% of events at ≥ 98% of the range).
   - **Do not apply a "resolved" filter.** Because the peaks were assigned on several detectors, dim peaks that overlap in one detector still have valid statistics, and they carry most of the information about B.
4. **Use SD directly, never rCV × median.** Blank medians can be negative, and the product then has the wrong sign. (This bug was found in the InstrumentQC dashboard.)
5. **Calculate the MEF scale k.**
   - Fit `log10(M) = log10(k) + log10(MEF)` on peaks 2–8 with the slope fixed at 1, weighted by n.
   - **Check the free-slope fit:** the slope should be within about 0.9–1.1 and R² above 0.98. If not, the lot's MEF scale does not fit that detector, and Q and B are not valid in MEF units.
   - In this study, BFP (slope 0.83) and PE-Cy7 (0.86) failed on every instrument.
6. **Fit the Parks model** (Parks et al. 2017) on all non-saturated peaks, **including the blank**:
   - weighted least squares on `V ~ M + I(M^2)` with weights `n / (2·V̂²)`;
   - start from the observed V, recompute the weights from the fitted V̂, and repeat until stable (≤ 20 iterations);
   - standard errors by the delta method from `vcov()`, e.g. `SE(Q) = k·SE(c1)/c1²`.
7. **Check validity:**
   - `c1 ≤ 0` → Q is undefined.
   - `c0 < 0` → B is undefined.
   - `c2 < 0` → refit with `V ~ M` and set CV₀ = 0.
   - **χ²/df of the weighted residuals:** about 2–7 is normal with IQR-based SDs. Above about 50, the model does not describe the data, so Q, B and 2+2SD are not valid. Report the model-free blank SD and separation index instead.
   - Plot the residuals per peak. They should show no trend.
8. **Calculate 2+2SD** = 4·(1/Q + BSD).
   - Take BSD from the measured blank SD divided by k, and check it against the fitted √c0/k. They agreed within 4% here, which shows the blank bead's own fluorescence adds little.
   - SE(2+2SD) ≈ 4·SE(Q)/Q².
9. **Cross-check:**
   - Refit without the blank: Q should barely change.
   - Run flowQB (`flowqb_crosscheck.R`). Its Q agreed within ±15% on 28 of 32 detectors.
   - The Chase closed form, *only* when its criteria are met (see §5). With RCP-30-5A they rarely are.
10. **Report** 2+2SD as the headline, with Q, B, BSD, CV₀, χ²/df, the log–log slope, scale placement, and caveats.

---

## 5. Pitfalls found and corrected (do not repeat)

| Pitfall | Why it is wrong | Effect seen |
|---|---|---|
| Q from a linear-only fit and B from a separate quadratic fit (old InstrumentQC `qb_fit()`) | The linear fit absorbs the background term into the slope; Q and B then come from different models | Q 25–60% too low |
| **STN-17 closed form**: `Q = 1/(CV²_blank·MEF_blank)`, `B = (SD²_blank/SD²_dim)·MEF_blank` | Chase & Hoffman 1998 Eq. 8/10 apply to a *moderately bright* bead (SD ≥ 3× blank, CV ≥ 3× bright) and scale B by **that bead's** MEF. The blank's variance is background-dominated | Q 3–6× too low |
| Subtracting the brightest peak's CV as "intrinsic" | The brightest peak still contains photon noise (up to about 30% of its CV² here) | Biased Q; fit CV₀ as a third parameter instead |
| Weights `1/CV⁴` that ignore n, on peaks 2–8 | Puts nearly all the weight on the brightest peaks, which carry little photon information | Unstable Q and B |
| Excluding "unresolved" dim peaks | Throws away the peaks that determine B | Unstable B |
| Ranking on blank SD alone | Ignores Q. A low-Q instrument can have a small background SD | Xenith appeared 2nd; it is last on 2+2SD |
| Using the lot's MEF where the detector doesn't match its assignment | Hard-dyed beads are not spectrally matched (Hoffman et al. 2012) | Log–log slope ≠ 1; Q and B in wrong units |

---

## 6. Interpretation caveats

- **Q is an effective (statistical-photoelectron) Q.** It includes each detector's excess-noise factor (PMT dynode statistics or APD multiplication noise). That is the right quantity for comparing how well dim populations resolve, but it is not a quantum efficiency. On APDs it can also change with gain.
- **Check the MEF definitions on the lot certificate.** Spherotech assigns each MEF unit at a specific excitation wavelength and filter. A secondary source (Cytoflow documentation) gives MECSB as 405 → 450/50, MEBFP as 405 → 530/40, and the PE family at 488 nm. Detectors far from these definitions give log–log slopes ≠ 1.
- **Red detectors** (APC, APC-Cy7): Hoffman & Wood (2007) found the Rainbow dim beads' *intrinsic CV* too large for accurate Q/B there. The beads do have signal in the red; the limitation is their spread. Treat these rows as lower-confidence.
- **Spectral instruments: this measures one raw detector per dye.** That is best-single-detector performance, not sensitivity after unmixing.
  - Bin widths differ between instruments, so per-detector Q partly reflects how much of the emission one bin catches.
  - A fairer comparison sums the bins inside each MEF filter band into a "pseudo-filter" (Parks 2020). This is valid for Rainbow beads, because all 8 peaks share one dye mix.
  - Sensitivity in unmixed space needs single-stain controls (Parks 2020; Feher et al. 2016, 2020).
- **One acquisition per instrument.** The SEs cover fit uncertainty only, not day-to-day variation.

### Typical values from this study (QC-default settings, 2+2SD in MEF)

| Fluor | A8 | Aurora Evo | Opteon | Xenith |
|---|---|---|---|---|
| Cascade Blue | 48 | 38 | **24** | 99 |
| FITC | 1,100 | 750 | **260** | — |
| PE | 220 | 190 | **59** | — |
| ECD | 96 | 76 | **25** | — |
| PE-Cy5 | 620 | 520 | **160** | — |
| PE-Cy7 | 1,800 | 1,600 | **450** | 2,900 |
| APC | 67 | 140 | **53** | — |
| APC-Cy7 | 88 | 97 | **52** | — |

"—" means the model fit was poor (χ²/df > 50), so there is no valid limit. On the A8, Aurora and Opteon the limit was background-limited: 4/Q was a median of about 5% of the total.

---

## Worked example: Opteon, PE detector (Y582-A)

Peak statistics (MFI):

| Peak | n | Median | SD | MEPE |
|---|---|---|---|---|
| 1 (blank) | 1214 | 104 | 591 (right half) | — |
| 2 | 1172 | 18,766 | 1,181 | 594 |
| 3 | 1234 | 58,424 | 1,868 | 1,448 |
| 4 | 1225 | 188,573 | 4,653 | 4,136 |
| 5 | 1370 | 582,990 | 13,071 | 12,439 |
| 6 | 1263 | 1,730,996 | 39,547 | 38,053 |
| 7 | 1182 | 5,483,212 | 131,836 | 112,904 |
| 8 | 1243 | 12,471,450 | 280,188 | 300,942 |

Results:

```
k   = 42.6 MFI per MEPE                       (slope-1 log–log fit, peaks 2–8)
fit:  c0 = 3.63e5   c1 = 33.7   c2 = 4.92e-4  (all 8 peaks, weighted, iterated; χ²/df = 12)

Q    = k/c1        = 42.6/33.7            = 1.26 Spe per MEPE   (SE 0.26)
B    = c0/(c1·k)   = 3.63e5/(33.7 × 42.6) = 253 MEPE
BSD  = √c0/k       = 602/42.6             = 14.1 MEPE  (measured blank: 591/42.6 = 13.9)
CV0  = √c2                                = 2.2 %
2+2SD = 4·(1/Q + BSD) = 4·(0.79 + 13.9)   = 59 MEPE
4/Q share                                 = 5 %   → background-limited
```

---

## 7. Implementation in this folder

| File | Role |
|---|---|
| `sensitivity_functions.R` | All of the calculations: `read_param_map()`, `read_mef()`, `assign_peaks()`, `peak_stats()`, `mef_scale()`, `parks_fit()`, `chase_qb()`, `analyse_channel()`, `caveat_flags()` |
| `Sensitivity_Comparison.qmd` | Report. Computes 2+2SD and its SE and the 4/Q share in the setup chunk, writes `results/`, and renders the HTML (`quarto render`) |
| `flowqb_crosscheck.R` | Independent check with flowQB 2.11.1, installed into `.flowqb-lib/` and patched in memory for current R. Writes `results/flowqb_crosscheck.csv` |
| `results/sensitivity_summary.xlsx` | Per detector: 2+2SD, Q, B, CV₀ with SEs, plus fit coefficients, per-peak statistics and residuals |
| `Parameter Assignment Cytometers.xlsx` | The detector used for each fluorochrome on each instrument. Entries in FlowJo's `name :: desc` form are matched on the text before ` :: ` |
| `8 peak linearity QC record AS01.xlsx` | Lot MEF values: one sheet per fluor, peaks 2–8 in column 3 |

**To reuse with new data:**
1. Put the singlet-gated `export_*.fcs` files in one folder per instrument.
2. Update the `instruments` table at the top of `sensitivity_functions.R`, the parameter-assignment sheet, and the MEF workbook for the new bead lot.
3. Render the report.
4. Re-run `flowqb_crosscheck.R` if you want the flowQB comparison updated.

---

## References

- Chase ES, Hoffman RA. Resolution of dimly fluorescent particles: a practical measure of fluorescence sensitivity. *Cytometry* 1998;33:267–279. doi:10.1002/(SICI)1097-0320(19981001)33:2<267::AID-CYTO24>3.0.CO;2-R
- Feher K et al. Multispectral flow cytometry: the consequences of increased light collection. *Cytometry A* 2016;89A:681–689 (corrigendum 2019;95A:803–804); and Response to Parks. *Cytometry A* 2020;97A:1278. doi:10.1002/cyto.a.24208
- Hoffman RA, Wang L, Bigos M, Nolan JP. NIST/ISAC standardization study: variability in assignment of intensity values to fluorescence standard beads… *Cytometry A* 2012;81A:785–796. doi:10.1002/cyto.a.22086
- Hoffman RA, Wood JCS. Characterization of flow cytometer instrument sensitivity. *Curr Protoc Cytom* 2007;Unit 1.20. doi:10.1002/0471142956.cy0120s40
- Maecker HT, Frey T, Nomura L, Trotter J. Selecting fluorochrome conjugates for maximum sensitivity. *Cytometry A* 2004;62A:169–173.
- Parks DR et al. Evaluating flow cytometer performance with weighted quadratic least squares analysis of LED and multi-level bead data. *Cytometry A* 2017;91A:232–249. doi:10.1002/cyto.a.23052
- Parks DR et al. Methodology for evaluating and comparing flow cytometers: a multisite study of 23 instruments. *Cytometry A* 2018;93A:1087–1091. doi:10.1002/cyto.a.23605
- Parks DR. Multispectral flow cytometry: unaddressed issues and recommendations for improvement. *Cytometry A* 2020;97A:1276–1277. doi:10.1002/cyto.a.24166
- Spidlen J et al. flowQB: automated quadratic characterization of flow cytometer instrument sensitivity. Bioconductor 2.11.1, 2019 (removed in Bioconductor 3.10).
- Wood JCS. Fundamental flow cytometer properties governing sensitivity and resolution. *Cytometry* 1998;33:260–266. doi:10.1002/(SICI)1097-0320(19981001)33:2<260::AID-CYTO23>3.0.CO;2-R
