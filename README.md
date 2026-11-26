# CMIP7 analysis: warming and climate sensitivity

**Report: <https://hausfath.github.io/cmip7-analysis/>**

A first look at the CMIP7 output published on ESGF as of **2026-11-25**: historical and projected
global surface air temperature compared with the GSAT ensemble, and effective climate sensitivity
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
| 2. Stream each file, verify sha256, reduce to an area-weighted global monthly mean, delete raw | `scripts/acquire_global_means.py` (`--accept-stale` for files rewritten after indexing: CanESM5-1 piControl, KACE-2-0-G historical and scen7-vl r1; `--insecure-host esg-earthlab.iap.ac.cn` while that node's TLS certificate is expired, which still requires an exact sha256 match) | `derived/monthly/`, `derived/attrs/`, `derived/acquire_log.csv` |
| 3. CMIP6 validation + predecessor runs (Pangeo) | `scripts/acquire_cmip6_validation.py` (then `--backfill-attrs`) | `derived/cmip6_validation/` |
| 4. ECS / TCR; CMIP6 first, as a check against published values | `scripts/sensitivity.py cmip6`, then `cmip7` | `outputs/cmip6_sensitivity.csv`, `outputs/cmip7_sensitivity.csv`, `derived/gregory/` |
| 5. Historical and projected GSAT vs the GSAT ensemble | `scripts/gmst.py` (obs from `data/obs/`; `--obs <file>` rebuilds it from the full ensemble) | `derived/gmst/`, `outputs/historical_metrics.csv`, `projection_metrics.csv`, `splice_*.csv`, `drift_check.csv` |
| 5b. CMIP6 historical metrics from the H22 supporting data (descriptive) | `scripts/cmip6_historical.py` | `outputs/cmip6_historical_metrics.csv`, `cmip6_splice_sensitivity.csv`, `cmip6_h22_checks.txt` |
| 5c. Earth's energy imbalance vs observations (obs inputs vendored once by `scripts/obs_eei.py`) | `scripts/eei.py` | `derived/eei/`, `outputs/eei_metrics.csv`, `eei_obs.csv`, `eei_tests.csv`, `eei_control.csv`, `eei_sensitivity.csv` |
| 6. Figures | `scripts/figures_gmst.py`, `scripts/figures_sensitivity.py`, `scripts/figures_eei.py` | `figures/fig01–fig08`, `fig02b` (PNG + PDF), `outputs/lineage_pairs.csv` |
| 7. Dataset citation table | `scripts/dataset_table.py` | `outputs/cmip7_dataset_table.csv` |
| 8. Report (GitHub Pages) | `scripts/build_report.py` | `docs/index.html`, `docs/figures/` |
| One-off | `scripts/oneoff_esmhist_scen7m.py`: every emissions-driven run to 2026 vs observations (not in the report) | `figures/oneoff_esmhist_scen7m` |

Run every script from the repository root with `.venv/bin/python -I scripts/<name>.py`. Steps 4–8 need
only the derived files committed here (about 130 MB), so every number, figure and the report regenerate
without downloading anything. Steps 1–2 move about 361 GB over ESGF and are resumable. ESGF holdings
change daily, so a fresh step 1 will not reproduce the 2026-11-25 manifest; `inventory/manifest_2026-11-25.csv`
records exactly what was used. That file was rebuilt on 2026-11-25, after the GSAT, ECS and TCR results were first published,
to add 4,829 radiation files (rsdt, rsut, rlut) for the energy imbalance section; the 3,898 files used for
those first results are unchanged in it (same versions and checksums).

### Weekly refresh

`scripts/weekly/refresh.sh` runs every Sunday night (launchd) from a separate clone:

1. It rebuilds the manifest and compares it with the last committed one (`scripts/weekly/scan.py`).
   The scan also lists models that have abrupt-4xCO2 or 1pctCO2 on ESGF but not yet the variables
   needed for ECS/TCR. If nothing is new, the run stops there.
2. Otherwise it acquires the new files and reruns steps 4–8.
3. `scripts/weekly/check.py` then flags anything that needs a person to look at it:
   - a new model;
   - withdrawn or re-versioned datasets;
   - incomplete acquisition or a stale checksum;
   - a failed area-weighting check;
   - a scenario join that is not under the metadata rule, or has an unusual offset;
   - a large change in an existing model's metrics.
4. Claude, headless and limited to editing README.md and METHODS.md, updates the numbers in those
   two files and gives a verdict (prompt: `scripts/weekly/claude_prompt.md`).
5. A clean run is pushed to `main`, which updates the report. Anything else goes to a
   `weekly-refresh` pull request for review.

fig02b and `cmip6_historical.py` are a descriptive side analysis (CMIP6 models from the
Hausfather et al. 2022 supporting data, one member each). They are not part of the report.

## Inputs and licences

- CMIP7 model output: ESGF, all 2,192 datasets CC-BY-4.0. Cite per the
  [CMIP7 terms of use](https://wcrp-cmip.github.io/cmip7-guidance/docs/CMIP7/Guidance_for_users/#2-terms-of-use-and-citations-requirements);
  dataset IDs, versions and tracking IDs are in `outputs/cmip7_dataset_table.csv`.
- CMIP6 output: Pangeo `gs://cmip6`. The six models used (ACCESS-ESM1-5, CanESM5, CNRM-ESM2-1,
  HadGEM3-GC31-LL/MM, UKESM1-0-LL) are all CC BY 4.0 per the
  [CMIP6 source_id registry](https://github.com/WCRP-CMIP/CMIP6_CVs/blob/main/CMIP6_source_id.json).
- CMIP7 experiment definitions in `reference/cmip7_cv/`: WCRP CVs (CC BY 4.0).
- Observations: summary statistics of a 10,000-member GSAT ensemble (`data/obs/`; see METHODS §1.2).
  The full ensemble is not public; its component datasets keep their own licences.
- Energy imbalance observations (`data/obs/eei/`, annual series derived from the sources): Wu et al. 2025
  (Zenodo 15395211, CC BY 4.0), Minière et al. 2026 (Zenodo 18485246, CC BY 4.0), Hakuba et al. 2021
  (Zenodo 5104970, "other-open"; three published values), CERES EBAF-TOA Ed4.2.1 (NASA Langley ASDC,
  freely available with citation).
- CMIP6 reference ECS/TCR and CMIP6 historical/SSP global-mean temperature series: Hausfather et al. (2022) supporting data,
  [Zenodo 6476375](https://zenodo.org/records/6476375); Zelinka et al. (2020) v2.2,
  [GitHub](https://github.com/mzelinka/cmip56_forcing_feedback_ecs) (no licence file, so `sensitivity.py`
  downloads it by tag instead of shipping it); IPCC AR6 WGI Table 7.SM.5
  (values transcribed into `reference/*.csv`, source noted in each header). The H22 workbook is CC-BY-4.0.

Code: [MIT licence](LICENSE). Derived model data keep the licences of the CMIP6/CMIP7 output they come from.
