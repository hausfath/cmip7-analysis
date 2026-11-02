"""Gregory ECS (150 yr) and TCR for CMIP6 (validation + lineage predecessors) and CMIP7.

  python -I scripts/sensitivity.py cmip6   # CMIP6 runs from Pangeo, checked vs published values
  python -I scripts/sensitivity.py cmip7   # CMIP7 models

Conventions (METHODS.md §2.3–2.5):
  ECS primary  = rm21 (Zelinka et al. 2020); sensitivities: linear, mean, ±1 yr branch shift,
                 years 21–150 fit.
  TCR primary  = linear (OLS fit over the full contemporaneous 1pctCO2-length piControl
                 segment, as ESMValTool climate_metrics/tcr.py, used for AR6 Table 7.SM.5 via
                 Meehl et al. 2020; reproduces AR6 to <=0.004 degC on the validation set);
                 rm21 and mean as sensitivities.
Outputs: outputs/cmip6_sensitivity.csv, outputs/cmip7_sensitivity.csv,
         derived/gregory/<model>_<member>.csv (annual ΔT, ΔN pairs, rm21)
"""
import sys, json
import numpy as np
import pandas as pd
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))  # -I drops the script dir
import cmiplib as L  # noqa: E402
ZELINKA_URL = "https://raw.githubusercontent.com/mzelinka/cmip56_forcing_feedback_ecs/v2.2/cmip56_forcing_feedback_ecs.json"

OUT = L.ROOT / "outputs"
GDIR = L.ROOT / "derived" / "gregory"
METHODS = ["rm21", "linear", "mean"]

CMIP6 = [("CanESM5", "r1i1p1f1", "r1i1p1f1"), ("UKESM1-0-LL", "r1i1p1f2", "r1i1p1f2"),
         ("ACCESS-ESM1-5", "r1i1p1f1", "r1i1p1f1"), ("CanESM5", "r1i1p2f1", "r1i1p2f1"),
         ("HadGEM3-GC31-LL", "r1i1p1f3", "r1i1p1f1"), ("HadGEM3-GC31-MM", "r1i1p1f3", "r1i1p1f1"),
         ("CNRM-ESM2-1", "r1i1p1f2", "r1i1p1f2")]
VALIDATION = {("CanESM5", "r1i1p1f1"), ("UKESM1-0-LL", "r1i1p1f2"), ("ACCESS-ESM1-5", "r1i1p1f1")}
CMIP7 = [("CanESM6-0-MR", "r1i1p1f1", "r1i1p1f1"), ("CanESM5-1", "r1i1p2f1", "r1i1p2f1"),
         ("UKESM1-3-LL", "r1i1p1f1", "r1i1p1f1"), ("UKCM2-0-LL", "r1i1p1f1", "r1i1p1f1"),
         ("UKCM2a-0-HH", "r1i1p1f1", "r1i1p1f1"), ("ACCESS-ESM1-6", "r1i1p1f1", "r1i1p1f1")]


def net(load, *args):
    sd, a = load(*args, "rsdt")
    su, _ = load(*args, "rsut")
    lu, _ = load(*args, "rlut")
    return (sd - su - lu).dropna(), sd


def run(model, mem, pmem, load):
    r = {"model": model, "member": mem, "pi_member": pmem}
    piT, pa = load(model, "piControl", pmem, "tas")
    piN, piSD = net(load, model, "piControl", pmem)
    r.update(pi_first=int(piT.index.min()), pi_last=int(piT.index.max()),
             pi_tas_K=float(piT.mean()), pi_rsdt=float(piSD.mean()), pi_N=float(piN.mean()))
    try:
        T4, a4 = load(model, "abrupt-4xCO2", mem, "tas")
        N4, _ = net(load, model, "abrupt-4xCO2", mem)
        by, off, dec = L.branch_year(a4, pa)
        r.update(branch_4x=by, branch_4x_offset_days=off, branch_4x_decoded=dec,
                 abrupt_first=int(T4.index.min()), abrupt_years=len(T4))
        for meth in METHODS:
            dT, nt = L.anomalies(T4, piT, by, 150, meth, a4)
            dN, _ = L.anomalies(N4, piN, by, 150, meth)
            for k, v in L.gregory(dT, dN).items():
                r[f"{k}_{meth}"] = v
            if meth == "rm21":
                r["rm21_truncated_windows_4x"] = nt
                pd.DataFrame({"year": np.arange(1, 151), "dT": dT, "dN": dN}).to_csv(
                    GDIR / f"{model}_{mem}.csv", index=False)
                g = L.gregory(dT[20:], dN[20:], nboot=1)
                r["ECS_21_150_rm21"] = g["ECS"]
                g = L.gregory(dT[:20], dN[:20], nboot=1)
                r["ECS_1_20_rm21"] = g["ECS"]
        for sh in (-1, 1):
            try:
                dT, _ = L.anomalies(T4, piT, by, 150, "rm21", shift=sh)
                dN, _ = L.anomalies(N4, piN, by, 150, "rm21", shift=sh)
                r[f"ECS_rm21_shift{sh:+d}"] = L.gregory(dT, dN, nboot=1)["ECS"]
            except ValueError as e:
                r[f"ECS_rm21_shift{sh:+d}"] = np.nan
    except (FileNotFoundError, ValueError) as e:
        r["abrupt_note"] = str(e)
    try:
        T1, a1 = load(model, "1pctCO2", mem, "tas")
        by1, off1, dec1 = L.branch_year(a1, pa)
        n1 = min(len(T1), 150)
        r.update(branch_1pct=by1, branch_1pct_offset_days=off1, branch_1pct_decoded=dec1,
                 onepct_years=len(T1))
        for meth in METHODS:
            dT, nt = L.anomalies(T1, piT, by1, 80, meth, a1, fit_len=n1)
            r[f"TCR_{meth}"] = L.tcr(dT)
            if meth == "rm21":
                r["rm21_truncated_windows_1pct"] = nt
        for sh in (-1, 1):
            try:
                r[f"TCR_linear_shift{sh:+d}"] = L.tcr(L.anomalies(T1, piT, by1, 80, "linear", shift=sh,
                                                                  fit_len=n1)[0])
            except ValueError:   # branch at the first/last control year
                r[f"TCR_linear_shift{sh:+d}"] = np.nan
    except (FileNotFoundError, ValueError) as e:
        r["onepct_note"] = str(e)
    return r


def zelinka_json():
    """Zelinka et al. (2020) v2.2 table. The source repo has no licence file, so it is fetched by tag, not vendored."""
    f = L.ROOT / "reference" / "zelinka" / "cmip56_forcing_feedback_ecs.json"
    if not f.exists():
        import urllib.request
        f.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(ZELINKA_URL, f)
    return f


def main(mode):
    OUT.mkdir(exist_ok=True)
    GDIR.mkdir(parents=True, exist_ok=True)
    rows = []
    if mode == "cmip6":
        ref = pd.read_csv(L.ROOT / "reference" / "hausfather2022_cmip6_ecs_tcr.csv", comment="#").set_index("model")
        ar6 = pd.read_csv(L.ROOT / "reference" / "ar6_table7sm5_cmip6_ecs_tcr.csv", comment="#").set_index("model")
        z = json.load(open(zelinka_json()))["CMIP6"]
        for m, mem, pmem in CMIP6:
            try:
                r = run(m, mem, pmem, L.load_cmip6_val)
            except FileNotFoundError as e:
                print("skip", m, mem, e)
                continue
            r["validation_set"] = (m, mem) in VALIDATION
            r["ECS_zelinka"] = z[m][mem]["ECS"]
            r["F2x_zelinka"] = z[m][mem]["ERF2x"]
            r["lambda_zelinka"] = z[m][mem]["NET"]
            # Published TCR refers to the p1 configuration; no reference value exists for other
            # physics variants (CanESM5 p2). Non-validation rows are an extra cross-check.
            p1 = "p1f" in mem
            r["TCR_H22"] = ref.loc[m, "TCR"] if p1 else np.nan
            r["TCR_AR6"] = ar6.loc[m, "TCR"] if p1 and m in ar6.index else np.nan
            rows.append(r)
        df = pd.DataFrame(rows)
        for k in ("ECS", "F2x", "lambda"):
            df[f"d{k}_vs_zelinka"] = df[f"{k}_rm21"] - df[f"{k}_zelinka"]
        for meth in ("rm21", "linear"):
            df[f"dTCR_{meth}_vs_H22"] = df[f"TCR_{meth}"] - df.TCR_H22
            df[f"dTCR_{meth}_vs_AR6"] = df[f"TCR_{meth}"] - df.TCR_AR6
        df.to_csv(OUT / "cmip6_sensitivity.csv", index=False)
    else:
        for m, mem, pmem in CMIP7:
            try:
                rows.append(run(m, mem, pmem, L.load_cmip7))
            except (FileNotFoundError, ValueError) as e:
                rows.append({"model": m, "member": mem, "abrupt_note": f"piControl: {e}"})
        df = pd.DataFrame(rows)
        df.to_csv(OUT / "cmip7_sensitivity.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 50, "display.max_rows", 200):
        print(df.round(3).T.to_string())


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "cmip7")
