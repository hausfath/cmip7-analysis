# CMIP7 analysis: warming and climate sensitivity

**Report: <https://hausfath.github.io/cmip7-analysis/>**

A first look at the CMIP7 output published on ESGF as of **2026-10-24**: historical and projected
global surface air temperature compared with the GMST ensemble, and effective climate sensitivity
(ECS) and transient climate response (TCR) compared with CMIP6. CMIP7 publication is ongoing, so every
number here belongs to that snapshot. Full conventions, checks and caveats are in [METHODS.md](METHODS.md).

## Pipeline

```
python3.13 -m venv .venv && .venv/bin/pip install -r requirements.txt
```
CMIP6 validation data come from the Pangeo cloud archive and need a separate venv with zarr v2
(`python3.13 -m venv .venv-acq && .venv-acq/bin/pip install -r requirements-acq.txt`), used only in step 3.

| Step | Script | Output |
|---|---|---|
| 1. File manifest from the ESGF CMIP7 STAC indexes | `scripts/build_manifest.py` | `inventory/manifest_<date>.csv` |
| 2. Stream each file, verify sha256, reduce to an area-weighted global monthly mean, delete raw | `scripts/acquire_global_means.py` (`--accept-stale` for the CanESM5-1 index issue) | `derived/monthly/`, `derived/attrs/`, `derived/acquire_log.csv` |
| 3. CMIP6 validation + predecessor runs (Pangeo) | `scripts/acquire_cmip6_validation.py` (then `--backfill-attrs`) | `derived/cmip6_validation/` |
| 4. ECS / TCR; CMIP6 first, as a check against published values | `scripts/sensitivity.py cmip6`, then `cmip7` | `outputs/cmip6_sensitivity.csv`, `outputs/cmip7_sensitivity.csv`, `derived/gregory/` |
| 5. Historical and projected GSAT vs the GMST ensemble | `scripts/gmst.py` (obs from `data/obs/`; `--obs <file>` rebuilds it from the full ensemble) | `derived/gmst/`, `outputs/historical_metrics.csv`, `projection_metrics.csv`, `splice_*.csv`, `drift_check.csv` |
| 5b. CMIP6 historical metrics from the H22 supporting data (descriptive) | `scripts/cmip6_historical.py` | `outputs/cmip6_historical_metrics.csv`, `cmip6_splice_sensitivity.csv`, `cmip6_h22_checks.txt` |
| 6. Figures | `scripts/figures_gmst.py`, `scripts/figures_sensitivity.py` | `figures/fig01–fig06`, `fig02b` (PNG + PDF), `outputs/lineage_pairs.csv` |
| 7. Dataset citation table | `scripts/dataset_table.py` | `outputs/cmip7_dataset_table.csv` |
| 8. Report (GitHub Pages) | `scripts/build_report.py` | `docs/index.html`, `docs/figures/` |
| One-off | `scripts/oneoff_esmhist_scen7m.py`: every emissions-driven run to 2026 vs observations (not in the report) | `figures/oneoff_esmhist_scen7m` |

Run every script from the repository root with `.venv/bin/python -I scripts/<name>.py`. Steps 4–8 need
only the derived files committed here (about 50 MB), so every number, figure and the report regenerate
without downloading anything. Steps 1–2 move about 130 GB over ESGF and are resumable. ESGF holdings
change daily, so a fresh step 1 will not reproduce the 2026-10-24 manifest; `inventory/manifest_2026-10-24.csv`
records exactly what was used.

fig02b and `cmip6_historical.py` are a descriptive side analysis (CMIP6 models from the
Hausfather et al. 2022 supporting data, one member each). They are not part of the report.

## Inputs and licences

- CMIP7 model output: ESGF, all 444 datasets CC-BY-4.0. Cite per the
  [CMIP7 terms of use](https://wcrp-cmip.github.io/cmip7-guidance/docs/CMIP7/Guidance_for_users/#2-terms-of-use-and-citations-requirements);
  dataset IDs, versions and tracking IDs are in `outputs/cmip7_dataset_table.csv`.
- CMIP6 output: Pangeo `gs://cmip6` (CMIP6 terms of use).
- Observations: summary statistics of a 10,000-member GMST ensemble (`data/obs/`; see METHODS §1.2).
  The full ensemble is not public; its component datasets keep their own licences.
- CMIP6 reference ECS/TCR and CMIP6 historical/SSP GMST series: Hausfather et al. (2022) supporting data,
  [Zenodo 6476375](https://zenodo.org/records/6476375); Zelinka et al. (2020) v2.2,
  [GitHub](https://github.com/mzelinka/cmip56_forcing_feedback_ecs); IPCC AR6 WGI Table 7.SM.5
  (values transcribed into `reference/*.csv`, source noted in each header). The H22 workbook is CC-BY-4.0.

Code: [MIT licence](LICENSE). Derived model data keep the licences of the CMIP6/CMIP7 output they come from.
