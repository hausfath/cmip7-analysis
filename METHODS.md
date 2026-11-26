# CMIP7 analysis: warming and climate sensitivity (methods)

Status: snapshot of ESGF holdings taken **2026-11-25** (all 8,727 manifest files acquired). CMIP7 publication
is ongoing daily; every number here is specific to that snapshot and will change.

## 1. Data

### 1.1 CMIP7 model output (ESGF)
- Discovery: ESGF CMIP7 STAC indexes (`discovery.east.esgf.io`, `discovery.west.esgf.io`) via
  `intake-esgf` 2026.9.28, per the [CMIP7 user guidance](https://wcrp-cmip.github.io/cmip7-guidance/docs/CMIP7/Guidance_for_users/#1-accessing-cmip7-data).
  Script: `scripts/build_manifest.py` → `inventory/manifest_2026-11-25.csv` (file-level: URLs,
  sizes, sha256 checksums, versions; latest version per dataset only). Both indexes responded when this
  manifest was built, and no dataset in the previous manifest is missing from this one.
  The manifest was rebuilt the same day to add the radiation files for §2.10 (4,829 files); the 3,898
  tas and ECS radiation files in the first build are unchanged (same versions and checksums).
- Branded variables: `tas_tavg-h2m-hxy-u` (monthly, region `glb`); `rsdt`, `rsut`, `rlut`
  `_tavg-u-hxy-u` (monthly) for piControl and abrupt-4xCO2 (ECS) and for historical, esm-hist,
  esm-piControl and the scen7/esm-scen7 files that start by 2024 (energy imbalance, §2.10). Sub-monthly output is not used:
  CAS-ESM2-1-CSM `scen7-vl` publishes `tas` only as 3-hourly fields, so it is absent here.
- Experiments: `historical`, `esm-hist`, `piControl`, `esm-piControl`, `abrupt-4xCO2`, `1pctCO2`,
  and all ScenarioMIP `scen7-*` / `esm-scen7-*` Tier-1 scenarios (only h, m, vl published so far).
  Experiment definitions (WCRP-universe CV, saved in `reference/cmip7_cv/`): `historical` and
  `esm-hist` run 1850–2021; `scen7-*` run 2022–2100. `scen7-h` = "High emission scenario … expected
  to result in forcings below SSP5-8.5"; `scen7-m` = "Medium … extending current policies";
  `scen7-vl` = "Very Low … limiting warming to 1.5C by the end of the century".
- Acquisition: `scripts/acquire_global_means.py`. Each file is downloaded over HTTPS (dropped
  connections resume with an HTTP Range request where the node supports it), sha256-verified
  against the index checksum (8,714 of 8,727 files pass; see stale-index note below), reduced to an area-weighted global mean per month, and deleted. No raw
  field is retained (peak raw footprint ≈ workers × one file). Outputs: `derived/monthly/<dataset>/`,
  global attributes in `derived/attrs/`, per-file log in `derived/acquire_log.csv`.
- TLS: the CAS node `esg-earthlab.iap.ac.cn` served an expired certificate (notAfter 2026-11-22), so
  its 18 FGOALS-g3-5 files were fetched with `--insecure-host esg-earthlab.iap.ac.cn`, which turns off
  certificate checking for that host only. Such files must match the index sha256 exactly (the
  stale-index exception does not apply to them); all 18 do, and their attributes record
  `checksum_status: sha256-ok, tls-unverified`.
- Area weights: exact spherical band areas from `lat_bnds`, |sin φ₂ − sin φ₁|. All 13 models acquired
  are rectilinear lon-lat grids with bounds. Check vs cos(lat) on the first time step of every file:
  |difference| ≤ 0.011 (native units; K or W m⁻²) for every model except FGOALS-g3-5 (1.23 K). Its
  GAMIL grid has unevenly spaced latitudes (2.0° in the tropics to 5.2° near the poles), so
  centre-point cos(lat) weights are wrong for it. Its bounds are contiguous, span ±90° and give band
  weights summing to 2; cos(lat) × Δlat weights agree with the bounds weights to 0.04 K (historical
  r1, first month). `areacella` not used.
- Holdings snapshot (members with `tas_tavg-h2m-hxy-u`, i.e. what was acquired; the raw inventory
  `inventory/cmip7_mon_member_counts.csv` also counts `tmaxavg`/`tminavg`/`tmax` brandings, which are
  monthly means of daily extremes and are **not** used):

| Model | Institution | hist | esm-hist | scen7 (h/m/vl) | esm-scen7 (h/m/vl) | piCtl+4xCO2 rad | 1pctCO2 |
|---|---|---|---|---|---|---|---|
| CanESM6-0-MR | CCCma | 25 | 25 | 25/25/25 | 25/25/25 | ✓ | ✓ |
| CanESM5-1 | CCCma | 16 | 25 | – | 25/25/25 | ✓ | ✓ |
| UKESM1-3-LL | UKNCSP | 5 (f2) | 8 | 3/–/3 | 8/7/8 | ✓ | ✓ |
| UKCM2-0-LL | UKNCSP | 10 | – | –/–/1 | – | ✓ | ✓ |
| UKCM2a-0-HH | UKNCSP | 3 | – | 2/–/2 | – | ✓ | ✓ |
| ACCESS-ESM1-6 | ACCESS-Consortium | 10 | 30 | 10/–/10 | 30/–/30 | ✓ | ✓ |
| CAS-ESM2-1-CSM | CAS | 1 | – | – (vl `tas` 3-hourly only) | – | piCtl only | – |
| EC-Earth3-ESM-1-1 | EC-Earth-Consortium | – | 2 | – | 2/–/2 | – | – |
| CNRM-ESM2-1e | CNRM-CERFACS | 1 (r2) | – (r6 has tmaxavg/tminavg only) | – | – | – | – |
| KACE-2-0-G | NIMS-KMA | 2 | – | 2/–/2 | – | – | – |
| FGOALS-g3-5 | CAS | 6 | – | 6/–/6 | – | – | – |
| FGOALS-f3-L | CAS | – | – | 3/–/3 | – | – | – |
| GFDL-ESM4p5 | NOAA-GFDL | – | – | – | – | piCtl `tas` only | – |

- Data-quality notes:
  - **CanESM5-1** datasets carry version `v20190429` (impossible for scen7 experiments). Global
    attributes confirm genuine CMIP7 output (`mip_era: CMIP7`, `creation_date: 2026-06-30`, CCCma runid
    `v51p2-c7a4xco2-r1`, parent `piControl` in CMIP7). Treated as a metadata error; included. Physics
    variant p2.
  - **UKESM1-3-LL** `historical` is forcing index f2 while its piControl/abrupt/1pct are f1.
  - **FGOALS-f3-L** has only `scen7-h`/`scen7-vl` (r1–r3). Each declares parent `historical r1i1p1f1`,
    which is not published, and there is no control run, so the model is listed but not analysed.
  - **GFDL-ESM4p5** has only a 1,000-year `piControl` (`tas`), so it is listed but not analysed.
  - **FGOALS-g3-5 scenario metadata.** All 12 `scen7-h`/`scen7-vl` runs (r1–r6) are 79 years long, but
    their `time_units` are `days since` 1972 (r2), 2022 (r3), 2122 (r1), 2172 (r4), 2222 (r5) and
    2272 (r6). Every run declares parent `historical r1i1p1f1` with `branch_time_in_parent` 1850 days
    after 1850-01-01 (1855). The values clearly run 2022–2100: each starts near the model's 2017–2021
    level, `scen7-h` ends about 1.6 K warmer and `scen7-vl` cooler. The parent cannot be identified from
    the data. Monthly Dec 2021 → Jan 2022 jumps to the six historical members (0.00–0.40 K) are no
    smaller than ordinary month-to-month changes (median 0.07 K, 90th percentile 0.20 K). No
    `scen7-h` run tracks any `scen7-vl` run over its first six months (smallest RMS difference 0.05 K),
    as two runs branched from the same state would. Following a user decision, `gmst.py` relabels each
    run to 2022–2100 and joins it to the historical member with the same number (rule
    `relabel-same-label` in `splice_check.csv`). That pairing is an assumption. Every metric is linear in the
    member values, so any other one-to-one pairing gives the same ensemble means. It does change the
    member ranges and the member-spread σ_iv of the spliced metrics (2015–2024 warming, 1980–2024 trend),
    and so their z. The 2012–2021 test uses only years up to 2021 and is unaffected. The joins look ordinary: 2022–2026 minus parent 2017–2021 is +0.15 °C (SD 0.05)
    against a median of +0.20 °C (SD 0.11) over all 354 joins. Worth reporting to CAS / ESGF errata.
  - **KACE-2-0-G scenario branch time.** `scen7-h` and `scen7-vl` (r1 and r7) give `branch_time_in_parent` = 0
    (1850-01-01), apparently a placeholder. Their time axes start 2022-01 and each names the historical
    run with its own member label (r1 or r7) as parent; 2021→2022 jumps are 0.05 and 0.02 K for r1 and
    0.09 and 0.02 K for r7. Each is joined to its named parent (rule `placeholder-branch`: branch time
    exactly 0, run starts in 2022, declared parent published).
  - **EC-Earth3-ESM-1-1 `esm-scen7-h` and `esm-scen7-vl` r4** have consistent metadata (parent
    `esm-hist r4i1p1f1`, branch 2022-01-01), but their 2022–2026 means are 0.65 °C (`esm-scen7-h`) and
    0.61 °C (`esm-scen7-vl`) below the parent's 2017–2021 mean (median over all joins +0.20 °C, SD 0.11),
    and their first years are 0.93 and 0.97 °C below the parent's 2021. Both are kept under the metadata
    rule, as is any other r4 scenario join. The model is not tested (2 members, no control), but these
    joins lower its 2015–2024 value (which continues r4 with `esm-scen7-h`) and the low end of both
    scenario ranges (2081–2100 vs 1850–1900: r4 3.24 vs r1 3.44 °C under `esm-scen7-h`, r4 1.51 vs r1
    1.73 °C under `esm-scen7-vl`).
    `splice_check.csv` has `jump_2021_2022` and `offset_5yr` for every join; the report flags any
    `offset_5yr` more than 0.5 °C from the median.
  - **Stale index checksums.** CanESM5-1 piControl: 8 files (all 4 `tas` files; 2 each of `rsut`, `rlut`) on the
    CCCma node were rewritten on 2026-08-19 under the same dataset version. The served files are
    complete (size equals the server `Content-Length`; e.g. 26,589,888 vs 26,612,556 bytes in the index),
    readable, and contiguous (1,200 months per file), but neither size nor sha256 matches the ESGF index.
    Accepted with `--accept-stale` and logged as `ok-stale` in `derived/acquire_log.csv`. CanESM5-1 ECS
    therefore rests on the server copy of its piControl, which cannot be checksum-verified. The
    KACE-2-0-G `historical` r1 `tas` file and its four `scen7-vl` r1 files (`tas`, `rsdt`, `rsut`, `rlut`) are in
    the same state (size equals the server `Content-Length`, sha256 differs from the index); their global means
    are complete and contiguous (2,064 months, 1850–2021; 948 months, 2022–2100). Worth
    reporting to CCCma, NIMS-KMA / [ESGF errata](https://errata.esgf.io/).
  - **Scenario parents are matched by metadata, not variant label.** UKESM1-3-LL `esm-scen7-*` r7/r8
    and `scen7-h`/`scen7-vl` r2 are `…f1` but declare an f2 parent (`esm-hist r7i1p1f2`/`r8i1p1f2`,
    `historical r2i1p1f2`), which is published, so they are spliced. UKESM1-3-LL `scen7-h`/`scen7-vl`
    r1 and r3 declare parent `historical r1i1p1f1`/`r3i1p1f1`, and `esm-scen7-m` r8 declares
    `esm-hist r8i1p1f1`; none of these is published (the historical members are f2, esm-hist r8 is
    f2), so they are not spliced (`outputs/splice_check.csv`). Every splice used under the metadata rule branches at 2022-01-01 in
    the parent's time axis; the two exceptions are the FGOALS-g3-5 and KACE-2-0-G cases above.

### 1.2 Observations
- GSAT ensemble (10,000 members, 1850–2025, anomalies vs 1850–1900 by construction), the build of
  2026-09-24 (HadCRUT5 5.2.0.0). The full ensemble is not public. `data/obs/` vendors the only obs
  statistics this analysis uses: the annual median, 5 %, 95 % and mean (`annual_summary.csv`), and each
  metric's median, 5 %, 95 % and SD over the members (`metric_summary.csv`). `scripts/gmst.py --obs <file>`
  regenerates both from the full ensemble. Run with the vendored files, it reproduces every output
  exactly (checked on 2026-11-02). The 2026-11-02 ensemble update (Berkeley Earth → BEST-HR beta) changes
  the annual median by at most 0.008 °C (1852), and by ≤ 0.001 °C for 2011–2020, 2012–2021 and each year
  2022–2025. It is not used here.
- Sanity check: obs 2011–2020 vs 1850–1900 = 1.088 °C [5–95 %: 0.979–1.230] vs AR6 SPM A.1.2
  1.09 [0.95–1.20] °C.
- The ensemble is anomalised so that the *cross-member* mean over 1850–1900 is zero; individual
  members are not, so per-member warming metrics include baseline uncertainty (intended).
- Naming: following the chapter convention from the SOD onward, the observational ensemble is called the
  GSAT ensemble and treated as equivalent to GSAT (Ch. 2 finds the GSAT-product median matches the
  GMST-product median over 1950–2021 and 1980–2021). It is compared directly with model `tas`, with
  **no** blended-tas/masking correction. Direction of the omitted effect: to the extent the observations
  reflect SST rather than air temperature over the ocean, and model blended warming is smaller than model
  `tas` warming, comparing `tas` with the observations makes models look *warmer* relative to obs, so "too warm" findings
  are, if anything, overstated and "too cool" findings understated. Rather than assume a size for the
  effect, `historical_metrics.csv` reports `<metric>_reduction_to_consistent`: the fractional reduction
  of model warming that would bring a too-warm model to the one-sided 5 % consistency threshold.

- Pangeo CMIP6 validation data (§2.5) have no latitude bounds in the zarr stores; band edges are
  taken as midpoints between grid latitudes (poles at ±90°), unlike the `lat_bnds` path used for CMIP7.
  The effect is bounded by the validation itself: ECS, ERF₂ₓ and λ reproduce Zelinka et al. (2020) to
  within the tolerances in §2.5, so the weighting difference is below that level.
- Pangeo time units can differ between variables of the same run (CanESM5 r1i1p1f1 abrupt-4xCO2 `tas`
  uses `hours since 1850-01-16 12:00`; its radiation stores use `days since 1850-01-01`). Units are
  therefore stored per variable (`time_units_<var>`).
- **Energy imbalance observations** (`data/obs/eei/`, built by `scripts/obs_eei.py` from the downloaded
  sources; raw files not kept):
  - Wu, Gregory, Zanna & Khatiwala (2025, PNAS, [10.1073/pnas.2408839122](https://doi.org/10.1073/pnas.2408839122)),
    Zenodo [15395211](https://doi.org/10.5281/zenodo.15395211) (CC BY 4.0; the latest version of the record
    11107298): annual 0–2000 m global ocean heat content, 72 members, 1870–2020 (`gohc`, ZJ), and the
    authors' 20-yr energy imbalance `Nm`/`Nd` (validation only). The earlier record's `dohu_dt` is a
    20-yr running trend of the 0–2000 m layer, not an annual full-depth value, so it is not used.
  - Minière, von Schuckmann & Speich (2026, ESSD Discussions,
    [10.5194/essd-2026-106](https://doi.org/10.5194/essd-2026-106), preprint), Zenodo
    [18485246](https://doi.org/10.5281/zenodo.18485246) (CC BY 4.0): monthly 0–2000 m global ocean heat
    content from 13 gridded products (EN4 once; its four XBT-correction variants are not used), annual
    means of complete years only. Primary: each product's native ocean domain (`no_mask`), the basis of
    the paper's 2005–2020 headline; sensitivity: their reference mask (60°S–60°N, deeper than 300 m),
    which excludes heat gain the models include.
  - Hakuba, Frederikse & Landerer (2021, GRL, [10.1029/2021GL093624](https://doi.org/10.1029/2021GL093624)),
    Zenodo [5104970](https://doi.org/10.5281/zenodo.5104970): published energy imbalance and 5–95 % for
    2005–2019 from hydrography (0.76) and from sea level minus ocean mass (0.94 W m⁻², the abstract value),
    `trends_global.xlsx`. Their non-ocean and deep terms are already included.
  - CERES EBAF-TOA Ed4.2.1 global monthly net flux (ceres-tool.larc.nasa.gov), annual means 2001–2025.
    Its net flux was adjusted once so that July 2005–June 2015 matches an in-situ ocean heating value of
    0.71 W m⁻² ([Loeb et al. 2018](https://doi.org/10.1175/JCLI-D-17-0208.1)), kept in Ed4.2 (Data
    Quality Summary V6). Its level is therefore not independent of the ocean heat estimates; its
    changes are.

### 1.3 CMIP6 reference values
- Primary: Hausfather, Marvel, Schmidt, Nielsen-Gammon & Zelinka (2022), *Nature* 605, 26–29,
  [doi:10.1038/d41586-022-01192-2](https://doi.org/10.1038/d41586-022-01192-2); supporting data
  [Zenodo 6476375](https://zenodo.org/records/6476375), sheet ECSTCR ("Summary of Mark Zelinka and IPCC
  calculations (Table 7.SM5) for ECS and TCR"). Columns ECS150, TCR. Extracted to
  `reference/hausfather2022_cmip6_ecs_tcr.csv`. Two rows that copy another model's values
  (CanESM5-CanOE, E3SM-1-1) are excluded → 53 models with ECS150, 52 with TCR.
  Mean ± SD: ECS 3.75 ± 1.07 °C, TCR 1.99 ± 0.41 °C.
- Cross-checks: Zelinka et al. (2020) v2.2 JSON ([GitHub](https://github.com/mzelinka/cmip56_forcing_feedback_ecs));
  IPCC AR6 WGI Ch. 7 SM Table 7.SM.5 (transcribed from the PDF page images to
  `reference/ar6_table7sm5_cmip6_ecs_tcr.csv`; transcription reproduces the printed CMIP6 mean/SD
  to within 0.02, residual from the table's own rounding/model set).
- AR6 assessed ranges (AR6 WGI Ch. 7, Table 7.13): ECS best 3.0, likely 2.5–4.0, very likely
  2–5 °C; TCR best 1.8, likely 1.4–2.2, very likely 1.2–2.4 °C.

## 2. Methods

### 2.1 Annual means
Calendar-aware: monthly global means weighted by days in month for each model calendar
(360_day, 365_day, proleptic_gregorian/standard). Only complete years (12 months) are used.

### 2.2 Branch alignment
Each abrupt-4xCO2 / 1pctCO2 run is aligned to its parent piControl using `branch_time_in_parent`
(converted to a parent model year with `parent_time_units` and the parent calendar) and rounded to the
nearest 1 January. Child year *k* (k = 1 … N) is paired with parent year (branch year + k − 1). Checks
(`cmiplib.anomalies`, all fail loudly): the rounding offset is recorded (0 days for every run used);
the child's first year must equal `branch_time_in_child` decoded in the child's own units; child years
must be contiguous; the parent must cover every contemporaneous year. Monthly timestamps must fall
mid-month (day 10–20) at both ends of each series, so end-of-interval stamps cannot shift years.

### 2.3 Effective climate sensitivity (Gregory)
- Anomalies: Δx(k) = x_4x(k) − x̃_pi(k), where x̃_pi is the **21-yr centred running mean** of the
  contemporaneous piControl (Zelinka et al. 2020, following Caldwell et al. 2016: "we compute
  anomalies relative to a 21-year running average from piControl"). At segment ends where the
  centred window would leave the published piControl, the window is truncated (sensitivity below).
- N = rsdt − rsut − rlut (net downward TOA).
- OLS ΔN = F + λ ΔT over years 1–150. ERF₂ₓ = F/2, ECS = −F/(2λ) (Gregory et al. 2004; AR6 Ch. 7 SM 7.SM.4).
- Sensitivity tests: (a) anomalies vs a linear fit to the contemporaneous 150-yr piControl segment;
  (b) vs the contemporaneous 150-yr piControl mean (no drift removal); (c) parent alignment shifted by
  ±1 yr (−1 is not computable where the branch is at the control's first year); (d) fit over years
  21–150 and 1–20 only (pattern-effect curvature).
- Uncertainty on the fit: fixed-design circular moving-block bootstrap of the regression *residuals*
  (ΔT held fixed; residual blocks of 10 yr wrap around the series end, so early high-leverage years
  are not undersampled; 2,000 resamples) → 5–95 % range of ECS. Describes regression sampling within one
  realisation only: it excludes realisation-to-realisation spread in ECS and structural uncertainty.

### 2.4 TCR
**Primary: linear convention (b) below**, chosen because it reproduces AR6 Table 7.SM.5 (§2.5).
1pctCO2 minus contemporaneous piControl, mean of years 61–80 (AR6 7.SM.4: "years 60 to 79",
zero-indexed). Two reference conventions are reported: (a) the 21-yr running mean (as for ECS), and
(b) an OLS line fitted to the contemporaneous piControl over the full 1pctCO2 length (≤150 yr), as in
ESMValTool `climate_metrics/tcr.py`, which underlies AR6 Table 7.SM.5 via Meehl et al. (2020). ±1 yr
branch-shift sensitivity as for ECS.

### 2.5 Validation (before applying to CMIP7)
Reproduce the published CMIP6 ECS/TCR for CanESM5 (r1i1p1f1), UKESM1-0-LL (r1i1p1f2), and
ACCESS-ESM1-5 (r1i1p1f1) from Pangeo `gs://cmip6` output (`scripts/acquire_cmip6_validation.py`,
identical downstream code). Pre-registered pass criteria: |ECS − Zelinka v2.2| ≤ 0.05 °C,
|ERF₂ₓ| ≤ 0.05 W m⁻², |λ| ≤ 0.02 W m⁻² K⁻¹, and |TCR − published| ≤ 0.05 °C, checked before any CMIP7 number
was computed.

**Result (`outputs/cmip6_sensitivity.csv`): passes.**
- ECS (rm21) matches Zelinka v2.2 to ≤ 0.014 °C on the three validation models (ERF₂ₓ ≤ 0.003 W m⁻²,
  λ ≤ 0.001 W m⁻² K⁻¹), and to ≤ 0.020 °C across all seven recomputed CMIP6 runs (including the
  predecessors CanESM5 p2, HadGEM3-GC31-LL/MM and CNRM-ESM2-1).
- TCR with the linear control fit matches AR6 Table 7.SM.5 to ≤ 0.004 °C on the validation models and to
  ≤ 0.007 °C on all six p1 runs with an AR6 value; it matches H22 to ≤ 0.027 °C on the
  validation models. H22 differs from AR6 for HadGEM3-GC31-LL (2.49 vs 2.55 °C); our value agrees with AR6.
- TCR with the 21-yr running-mean control is up to 0.067 °C higher than H22 and fails the criterion, so the
  linear convention is the headline TCR and rm21 is reported as a sensitivity.
- The ±1 yr branch shifts move CMIP7 ECS by ≤ 0.01 °C and TCR by ≤ 0.004 °C (`outputs/cmip7_sensitivity.csv`).

### 2.6 Historical warming (`scripts/gmst.py`)
- Per member: annual tas → anomaly vs that member's own 1850–1900 mean. Historical runs are truncated
  at 2021 (ACCESS-ESM1-6 and one UKCM2a-0-HH member extend to December 2022, CAS-ESM2-1-CSM to
  November 2022).
- Per model: the value is the **mean of per-member metrics** over members covering the window
  (`<metric>_n` columns), never the metric of a skipna ensemble-mean series (which would change
  composition mid-window). One value per model × track in all cross-model statements.
  `historical` (concentration-driven) and `esm-hist` (emissions-driven) are separate tracks, never
  pooled.
- Extension past 2021: each scenario run is joined to the historical run named in its own
  `parent_experiment_id`/`parent_variant_label`, required to branch at 2022-01-01 (exceptions for
  FGOALS-g3-5 and KACE-2-0-G, §1.1). Preference scen7-m > scen7-h > scen7-vl (esm- analogues for
  esm-hist).
- Splice sensitivity (`outputs/splice_sensitivity.csv`): ensemble-mean metric with each scenario over
  the members continued in all of them. Max range across scenarios: 0.029 °C (2015–2024 warming) and
  0.008 °C/decade (1980–2024 trend), both EC-Earth3-ESM-1-1 (n = 2); UKCM2a-0-HH (n = 2) 0.021 °C and
  0.006 °C/decade; ≤ 0.012 °C and ≤ 0.003 °C/decade for the others. Per-member differences are larger (chaotic divergence after the branch) and are not a splice
  error.
- Metrics vs GSAT ensemble:
  1. Warming 2012–2021 vs 1850–1900 (every model; no splice needed). Also 2011–2020 for AR6 comparison.
  2. Warming 2015–2024 vs 1850–1900 (spliced models only).
  3. OLS trends 1980–2021, 1950–2021 (Ch. 2 GMST-vs-GSAT comparison periods) and 1980–2024 (spliced).
- Obs uncertainty: distribution of each metric across the 10,000 GSAT-ensemble members (median,
  5–95 %, SD).
- **Consistency test** (`<metric>_z`): z = (model value − obs median) / σ, with
  σ² = σ_obs² + σ_iv² (1 + 1/n). The model value is an n-member mean; obs is one realisation, hence
  the (1 + 1/n). σ_iv is the SD of the metric across members if ≥ 5 members, else from the detrended
  parent control (σ of overlapping 10-yr and 51-yr means combined in quadrature for warming metrics,
  treating the two windows as independent; σ of overlapping 42/45/72-yr OLS trends for trends).
  "Warmer than observed" = z > 1.645 (one-sided 5 %). Models with fewer than 5 members and no control
  run (EC-Earth3-ESM-1-1, CNRM-ESM2-1e, KACE-2-0-G) are not tested. FGOALS-g3-5 (n = 6, no control
  run) is tested with its member spread (σ_iv = 0.067 °C, z = 3.58 for 2012–2021). The σ_iv source
  matters for UKESM1-3-LL `historical` (n = 5, so the member spread is used: σ_iv = 0.164 °C, z = 0.90 for 2012–2021);
  the piControl σ_iv (0.127 °C) would give z = 1.09, also below the threshold. Headline counts are conservative toward the
  "too warm" claim: a model counts only if every published track exceeds the threshold.
- **Control drift** (`outputs/drift_check.csv`, sensitivity only). Two corrections, each converted
  to its spurious contribution to 2012–2021 minus 1850–1900 (slope × 141.5 yr) and tested with the same σ:
  - `<metric>_driftcorr`: removes the parent control's **full-record** OLS slope, i.e. the drift itself.
  - `<metric>_driftcorr_seg`: removes each member's **contemporaneous 1850–2021 segment** slope
    (branch-aligned as §2.2; full-record where the segment is incomplete, as for 8 of the 16 CanESM5-1
    `historical` members, which branch at control years 6731–6841 while the published piControl `tas`
    covers 6500–6900). Over 172 yr this mixes drift
    with centennial variability, which the member spread already samples, so it is the weaker test.
  Headline count of models warmer than observed over 2012–2021: 5 of 8 raw, 4 of 7 with the full-record
  correction, 3 of 7 with the segment correction (UKCM2-0-LL drops out). FGOALS-g3-5 has no control
  run, so it is in the raw count only. Notable cases (raw / full-record /
  segment, °C): UKCM2-0-LL 1.43 / 1.41 / 1.25, whose control has strong centennial variability
  (century means 12.95–13.47 °C; segment slopes −0.008 to 0.033 °C/decade for branches at control years
  2100–2280, full-record 0.002); UKCM2a-0-HH 1.48 / 1.63 / 1.63, a genuine cooling drift
  in a 250-yr control (−0.011 °C/decade), so correction *raises* its warming; CAS-ESM2-1-CSM
  1.13 / 0.99 / 0.94. Control runs are not published for EC-Earth3-ESM-1-1, CNRM-ESM2-1e,
  KACE-2-0-G and FGOALS-g3-5.

### 2.7 Projections
Per model ensemble mean (and member range) of ΔGSAT vs 1850–1900 for each scenario; 2041–2060 and
2081–2100 means. Too few models for a meaningful multi-model distribution; reported model by model.
ACCESS-ESM1-6 `scen7-*` and `esm-scen7-*` files run to 2150; only years to 2100 enter the figures and
metrics.

### 2.8 Comparison with CMIP6
Place each CMIP7 ECS/TCR on the CMIP6 distribution (H22), with AR6 assessed ranges. Lineage
predecessors present in the CMIP6 table are marked for context (CanESM5 → CanESM5-1, CanESM6-0-MR;
UKESM1-0-LL → UKESM1-3-LL; HadGEM3-GC31-LL/MM → UKCM2-0-LL / UKCM2a-0-HH; ACCESS-ESM1-5 → ACCESS-ESM1-6;
CNRM-ESM2-1 → CNRM-ESM2-1e). With ≤ 6 CMIP7 models, no inferential test of a CMIP6→CMIP7 shift is
attempted.
- Pairings are our judgement by model family, not an official lineage. The least certain is
  UKCM2a-0-HH ← HadGEM3-GC31-MM: HH is the high-resolution configuration, but CMIP6 HadGEM3-GC31-HH has
  no abrupt-4xCO2, so MM is the nearest relative with an ECS. CanESM6-0-MR is a substantially new model
  (CanAM6/GEM atmosphere, NEMO4.2) paired with CanESM5 p1 for context only.
- CanESM5-1 r1i1p2f1 is CanESM5.1-p2 ([Sigmond et al. 2023](https://doi.org/10.5194/gmd-16-6553-2023)),
  whose atmosphere was retuned to lower its effective climate sensitivity to 4.09 K (their value). The model
  is in the CMIP6 CVs as source_id `CanESM5-1` (release year 2022), and ESGF has CMIP6 output for it, including
  p2 abrupt-4xCO2. It is not in H22. So its lower ECS relative to CanESM5 dates from that CMIP6-era retuning,
  not from development for CMIP7. CanESM5 stays its predecessor here because CanESM5 is the CCCma model
  in the CMIP6 ensembles (H22, AR6). Our CMIP7 value agrees with Sigmond et al. to within 0.05 °C (see
  `outputs/lineage_pairs.csv`). The pairing uses CanESM5 r1i1p2f1, which is CanESM5.0-p2 and unrelated to
  CanESM5.1-p2 despite the shared label. Its ECS is within 0.03 °C of r1i1p1f1 (`outputs/cmip6_sensitivity.csv`).
- Predecessor values are recomputed with the same code (`outputs/lineage_pairs.csv`).
- "Lower ECS" (headline) requires the two models' 5–95 % regression ranges not to overlap; the
  looser test (predecessor's central value outside the CMIP7 range) is the sensitivity. A missing
  predecessor is "untested". Both ranges are within-realisation regression ranges and exclude
  realisation-to-realisation ECS spread, so they understate the uncertainty of each difference.

### 2.9 CMIP6 historical warming (`scripts/cmip6_historical.py`, fig02b)
Descriptive counterpart to fig02 for CMIP6, from the H22 supporting data (Zenodo 6476375, sheets
`Historical` 1850–2014 and `SSP126/245/370/585` 2015–2100; 37 models, annual means).
- **What the series are.** Neither the Comment nor the Zenodo record states the member or the
  variable ("global mean surface temperature anomalies"). Each column is a single realisation, not an
  ensemble mean: interannual noise (SD of 1850–1950 first differences / √2) is 0.054–0.128 °C across
  models, including CanESM5 (0.080) and MIROC6 (0.115), which have ~50 members each; a 50-member mean
  would be ~0.012 °C (`outputs/cmip6_h22_checks.txt`).
  Which member is unknown. Treated as tas, consistent with the no-blending convention (§1.2).
- **Baseline.** H22 anomalies are relative to 1850–1899 (column means 1.4e-10). Re-baselined to
  1850–1900 to match the obs and CMIP7 rows; this shifts models by −0.0033 to +0.0074 °C.
- **Splice.** Historical ends in 2014, so both metrics need a scenario. Primary: SSP2-4.5 (the only
  sheet with all 37 models). Sensitivity: SSP1-2.6, SSP3-7.0, SSP5-8.5 where present (10 model × SSP
  gaps; 2–4 joins per model). The range across joins is the bar in fig02b. Spread in 2012–2021 warming:
  median 0.067 °C, max 0.150 °C (IPSL-CM6A-LR); in the 1980–2021 trend: median 0.020, max
  0.041 °C/decade. Forced differences between SSPs are small by 2021; most of this spread is chaotic
  divergence of the scenario branches after 2015, which decorrelates annual noise within a year or two
  even from a common initial state (the same effect seen in the CMIP7 per-member splices, §2.6). It
  therefore measures splice-choice sensitivity and says nothing about member identity. The 2014→2015 joins to SSP2-4.5 are within each
  model's own 1985–2014 year-to-year changes (|z| > 2 for 3 of 37, ~1.7 expected by chance), so a
  hist/SSP member mismatch is not detectable but cannot be excluded.
- **Read check.** The workbook's `CMIP6 MMM` "50th Percentile" columns equal the across-model *mean*
  of the SSP sheets to ≤ 7e-7 °C for 2015–2099 (SSP1-2.6, 2-4.5, 3-7.0; SSP5-8.5 differs by up to
  0.058, presumably a different model set). Its 5th/95th columns match no simple statistic of the 37
  columns (quantiles, min/max, mean ± 1.645 SD), so they are not used. This confirms the columns are
  read correctly; it is not an external validation of the series.
- **Descriptive only.** One member and no control run per model, so no σ_iv and no z-test or count
  of "too warm" models (user choice). Results (1850–1900 baseline): CMIP6 median 2012–2021 warming
  1.10 °C (5–95 % across models 0.85–1.71) vs obs 1.11 °C; median 1980–2021 trend 0.241 °C/decade
  (0.165–0.411) vs obs 0.191; 29 of 37 members exceed the observed median trend.
- **Predecessors.** CanESM5, UKESM1-0-LL, HadGEM3-GC31-LL and ACCESS-ESM1-5 are highlighted (§2.8
  pairings). HadGEM3-GC31-MM (UKCM2a-0-HH) is not in H22. Comparisons with CMIP7 mix one CMIP6
  member against CMIP7 member means, and CMIP6 vs CMIP7 historical forcing.
- Outputs: `outputs/cmip6_historical_metrics.csv` (obs row copied from `historical_metrics.csv`),
  `outputs/cmip6_splice_sensitivity.csv`, `outputs/cmip6_h22_checks.txt`.

### 2.10 Earth's energy imbalance (`scripts/eei.py`, fig07–fig08)
- Model EEI: net TOA flux N = rsdt − rsut − rlut (global, annual, W m⁻²) for every member of historical
  and esm-hist with all three variables, continued to 2024 with the scenario member used for GSAT
  (`splice_check.csv`, preference m > h > vl). Members without scenario radiation end in 2021.
- Control correction: the OLS line fitted to the full control N record (historical → piControl,
  esm-hist → esm-piControl; the member named as parent, else the first) is subtracted at each member's
  contemporaneous control years (§2.2 branch year). The correction is almost entirely a constant energy
  leak; per-model control means and slopes are in `outputs/eei_control.csv`. Models without published
  control radiation are not analysed. Sensitivities (`eei_sensitivity.csv`): control mean only, and a
  quadratic fit (relevant where the control N curves, e.g. CAS-ESM2-1-CSM). A leak that changes as the
  climate warms cannot be detected by any control-based correction.
- Observations as W m⁻² of Earth's surface (radius 6,371 km; Minière's own tables use a TOA area with
  +20 km, 0.6 % smaller values). Ocean heat records: rate = OLS slope of the annual 0–2000 m heat content
  over the window, plus 0.06 W m⁻² for the ocean below 2000 m from 1992 (zero before; von Schuckmann et
  al. 2023, ±0.03 taken as 1σ, the wider reading), divided by f_ocean = 0.89, the ocean share of Earth's
  heat gain for 1971–2020 and 2006–2020 in von Schuckmann et al. (2023). Sensitivity: 0.91 (AR6 SPM
  A.4.2). Wu: mean and SD over the 72 members; Minière: mean and SD over the products that cover the
  whole window (13 for 2005–2019, 7 for 2015–2024; the shading in fig07 is the product range).
- Estimators: models are compared with the ocean heat estimates through the OLS slope of their
  cumulative heat (annual-mean integral of N) over the same years, and with CERES through the window
  mean of N. The two differ by a few hundredths of a W m⁻² for 10–15-yr windows.
- Fixed windows, chosen before the comparison: 2005–2019 (the Hakuba window; Wu, Minière and Hakuba
  tested, CERES shown), 2015–2024 (Minière tested, CERES shown), and the change from 2005–2014 to
  2015–2024 (CERES and Minière tested). The change cancels both a constant model leak and the CERES
  anchor. fig07 shows centred 10-yr running rates for description only.
- Test: z = (model ensemble mean − estimate) / σ, σ² = σ_obs² + σ_iv²(1 + 1/n), with σ_iv the SD of the
  same statistic over overlapping windows of the control run's residual about its linear fit (all
  analysed models have a control by construction). σ_obs: member/product SD (Wu, Minière), the published
  5–95 % / 3.29 (Hakuba), zero for the CERES change (instrument stability uncertainty not included).
  A model is "lower" or "higher" only if |z| > 1.645 in the same direction against every tested
  estimate for that window; otherwise it is consistent with at least one.
  The report also gives the count against the CERES change alone, since CERES and Minière disagree on
  the 2005–2014 to 2015–2024 increase and the "every estimate" rule hides that.
- Validation: Minière OLS rates reproduce their published trend table (FIG4) to 0.0001 W m⁻² for 2005–2020
  and 1960–2020 under both masks (e.g. 0.5741 native domain, 0.5417 reference mask, TOA area). Wu: our
  20-yr rate converted with their conventions (+0.07, ÷ 0.90) is 0.772–0.780 W m⁻² depending on the window
  ends, against their 0.775 for 2000–2020. CERES 2001–2024 trend 0.447 W m⁻² per decade (Mauritsen et al.
  2025: 0.45). Model side: the control-corrected ensemble-mean N averages −0.11 to +0.02 W m⁻² over
  1850–1900 across the 11 model ensembles. UKESM1-3-LL, the one model with both controls, has similar
  piControl and esm-piControl means (0.098 and 0.087 W m⁻²).
- Members reaching 2024: only members with a scenario join have 2022–2024, so some 2015–2024 and change
  statistics rest on fewer members than the 2005–2019 ones (UKCM2-0-LL 1 of 10, UKESM1-3-LL historical
  1 of 5, UKCM2a-0-HH 2 of 3); the report table gives n per window where it differs.

## 3. Caveats
- Early and non-random sample: CMIP7 holdings are dominated by UK, Canadian and Australian centres;
  CMIP7 ensemble statistics are not representative of the eventual CMIP7 ensemble.
- 150-yr Gregory ECS is an *effective* sensitivity and underestimates the long-run equilibrium for
  models with a strong pattern effect.
- No blending/masking correction (see 1.2); the omitted effect would, if anything, reduce
  apparent model warm bias.
- The consistency test assumes model internal variability is realistic; using member spread from
  esm-hist includes carbon-cycle variability. With fewer than 5 members (UKCM2a-0-HH, CAS-ESM2-1-CSM),
  σ_iv comes from the control, and the test is sensitive to that choice; UKESM1-3-LL historical has
  exactly 5 (§2.6).
- CMIP7 historical forcing (CMIP7 forcing datasets, ending 2021) differs from CMIP6; differences
  between a CMIP7 model and its CMIP6 predecessor combine model and forcing changes.
- esm-hist runs simulate their own CO₂; part of any esm-hist vs historical difference is CO₂, not
  climate response.

## 4. Citation / data availability (CMIP7 requirements)
No CMIP7 DOIs had been minted as of the guidance page (2026-09-23). Datasets are cited by DRS id and
version (`inventory/manifest_2026-11-25.csv`); a source/institution/version table will be generated
into `outputs/cmip7_dataset_table.csv`. Licence: per-dataset global attribute (CC-BY-4.0 for all files
inspected). Acknowledgement text per the CMIP7 guidance to be included in any publication.
