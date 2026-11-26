"""Vendor the observational Earth energy imbalance inputs as small annual CSVs in data/obs/eei/.

Sources (downloaded once by hand; raw files are not kept in the repository):
  Wu et al. (2025, PNAS)        Zenodo 15395211: GF_OHC_estimate.nc (gohc, 0-2000 m, 72 members),
                                global_energy_budget.nc (Nm, Nd)
  Miniere et al. (2026, ESSD-D) Zenodo 18485246: GOHC_ensemble_1960-2024_<mask>.nc (monthly, J)
  Hakuba et al. (2021, GRL)     Zenodo 5104970: Data/trends_global.xlsx
  CERES EBAF-TOA Ed4.2.1        ceres-tool.larc.nasa.gov subset: gtoa_net_all_mon (monthly, W m-2)

Usage: python -I scripts/obs_eei.py --wu DIR --miniere DIR --hakuba XLSX --ceres NC
       (--wu DIR may hold wu_gohc_0-2000m_ZJ.csv instead of the 800 MB GF_OHC_estimate.nc)
"""
import argparse, pathlib
import numpy as np
import pandas as pd
import xarray as xr

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "obs" / "eei"
TC = xr.coders.CFDatetimeCoder(use_cftime=True)
# products in the Miniere ensemble mean; EN4 counted once (the EN4_* XBT variants are sensitivity runs)
MINIERE_PRODUCTS = ["EN4", "ISHII", "IAP", "NOAA", "CORA5.2", "RFROM", "ARMOR3D", "JAMSTEC", "SCRIPPS",
                    "CSIO", "GDCSM", "IPRC", "ISAS"]
MASKS = {"no_mask": "ocean_domain_sensitivity_test/GOHC_ensemble_1960-2024_no_mask.nc",
         "GCOS_surface": "GOHC_ensemble_1960-2024_Miniere_etal_2026.nc"}


def wu(d):
    d = pathlib.Path(d)
    if (d / "GF_OHC_estimate.nc").exists():
        x = xr.open_dataset(d / "GF_OHC_estimate.nc", decode_times=TC, drop_variables=["theta_e"])
        g = pd.DataFrame(x.gohc.values.T, index=[t.year for t in x.time.values],
                         columns=[f"m{int(e):02d}" for e in x.ENSEMBLE.values])
    else:
        g = pd.read_csv(d / "wu_gohc_0-2000m_ZJ.csv", index_col=0)
    g.index.name = "year"
    g.round(4).to_csv(OUT / "wu2025_gohc_0-2000m_ZJ.csv")
    e = xr.open_dataset(d / "global_energy_budget.nc", decode_times=TC)
    pd.DataFrame({"year": [t.year for t in e.time.values], "Nm": e.Nm.values, "Nd": e.Nd.values}).round(5) \
        .to_csv(OUT / "wu2025_energy_imbalance_20yr.csv", index=False)
    print("Wu:", g.shape, "years", g.index.min(), g.index.max())


def miniere(d):
    for mask, rel in MASKS.items():
        x = xr.open_dataset(pathlib.Path(d) / "outputs" / "netCDF" / rel)
        assert x.attrs.get("mask") == mask, (rel, x.attrs.get("mask"))
        df = pd.DataFrame({p: x[f"GOHC_0-2000m_{p}"].to_series() for p in MINIERE_PRODUCTS + ["ensemble_mean"]})
        g = df.groupby(df.index.year)
        ann = g.mean().where(g.count() == 12)          # complete years only
        ann.index.name = "year"
        ann.to_csv(OUT / f"miniere2026_gohc_0-2000m_J_{mask}.csv", float_format="%.5e")
        print("Miniere", mask, ann.notna().sum().to_dict())


def hakuba(xlsx):
    rows = []
    for sheet in ["2005-2019", "2003-2019", "2005-2015", "2010-2018"]:
        t = pd.read_excel(xlsx, sheet_name=sheet).set_index("Unnamed: 0")
        for label, method in [("Hydrographic observations 0-2000m", "hydrography"),
                              ("GMSL - barystatic", "geodetic")]:
            r = t.loc[label]
            rows.append({"window": sheet, "method": method, "eei_p05": r["EEI (W m-2) 5th"],
                         "eei_mean": r["EEI (W m-2) mean"], "eei_p95": r["EEI (W m-2) 95th"]})
    pd.DataFrame(rows).round(4).to_csv(OUT / "hakuba2021_eei_windows.csv", index=False)
    print("Hakuba:", len(rows), "rows")


def ceres(nc):
    x = xr.open_dataset(nc)
    s = x["gtoa_net_all_mon"].to_series()
    g = s.groupby(s.index.year)
    ann = g.mean().where(g.count() == 12).dropna()
    ann.index.name = "year"
    ann.rename("net_toa_Wm2").round(4).to_csv(OUT / "ceres_ebaf_ed4.2.1_annual.csv")
    print("CERES:", ann.index.min(), ann.index.max())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--wu"), ap.add_argument("--miniere"), ap.add_argument("--hakuba"), ap.add_argument("--ceres")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    for f, arg in [(wu, a.wu), (miniere, a.miniere), (hakuba, a.hakuba), (ceres, a.ceres)]:
        if arg:
            f(arg)
