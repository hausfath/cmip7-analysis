"""Shared analysis code: annual means, branch alignment, Gregory ECS, TCR.

Used identically for the CMIP6 validation set and for CMIP7 (see METHODS.md §2).
Every alignment step fails loudly rather than silently shifting years.
"""
import json, pathlib
import numpy as np
import pandas as pd
import cftime

ROOT = pathlib.Path(__file__).resolve().parents[1]
MONTHLY = ROOT / "derived" / "monthly"
ATTRS = ROOT / "derived" / "attrs"
VAL = ROOT / "derived" / "cmip6_validation"

_NOLEAP = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def days_in_month(year, month, calendar):
    cal = calendar.lower()
    if cal == "360_day":
        return 30
    if cal in ("noleap", "365_day"):
        return _NOLEAP[month - 1]
    if cal in ("all_leap", "366_day"):
        return 29 if month == 2 else _NOLEAP[month - 1]
    leap = (year % 4 == 0 and year % 100 != 0) or year % 400 == 0
    if cal == "julian":
        leap = year % 4 == 0
    return 29 if (month == 2 and leap) else _NOLEAP[month - 1]


def check_midmonth(mon, units, calendar, label):
    """Month labels must come from mid-interval timestamps. End-of-interval
    stamps would shift every month by one and still pass the 12-month filter."""
    d = cftime.num2date(mon.time.to_numpy()[[0, -1]], units, calendar)
    for x in d:
        if not 10 <= x.day <= 20:
            raise ValueError(f"{label}: timestamp {x} is not mid-month")


def annual(mon, var, calendar, label=""):
    """Month-length-weighted annual means; complete years only, no NaN months."""
    dup = mon[mon.duplicated(["year", "month"], keep=False)]
    if len(dup):
        spread = dup.groupby(["year", "month"])[var].agg(lambda s: s.max() - s.min()).max()
        if spread > 1e-6:
            raise ValueError(f"{label}: overlapping records disagree (max diff {spread})")
    m = mon.drop_duplicates(["year", "month"]).copy()
    if m[var].isna().any():
        raise ValueError(f"{label}: {int(m[var].isna().sum())} NaN monthly values")
    m["w"] = [days_in_month(y, mo, calendar) for y, mo in zip(m.year, m.month)]
    g = m.groupby("year")
    n = g.month.nunique()
    a = (m[var] * m.w).groupby(m.year).sum() / g.w.sum()
    return a[n == 12]


# ---------------------------------------------------------------- CMIP7 I/O
def dskey(source, exp, member, var):
    hits = sorted(MONTHLY.glob(f"MIP-DRS7.CMIP7.*.*.{source}.{exp}.{member}.glb.mon.{var}.*"))
    hits = [h for h in hits if h.name.split(".")[10] in ("tavg-h2m-hxy-u", "tavg-u-hxy-u")]
    if len(hits) != 1:
        raise FileNotFoundError(f"{source} {exp} {member} {var}: {len(hits)} matches")
    return hits[0].name


def load_cmip7(source, exp, member, var):
    k = dskey(source, exp, member, var)
    files = sorted((MONTHLY / k).glob("*.csv"))
    mon = pd.concat([pd.read_csv(f) for f in files]).sort_values("time")
    at = json.loads((ATTRS / f"{k}.json").read_text())
    check_midmonth(mon, at["time_units"], at["calendar"], k)
    if int(at.get("masked_points", 0)):
        raise ValueError(f"{k}: masked points in a global field")
    return annual(mon, var, at["calendar"], k), at


def members_cmip7(source, exp, var="tas"):
    out = []
    for h in MONTHLY.glob(f"MIP-DRS7.CMIP7.*.*.{source}.{exp}.*.glb.mon.{var}.*"):
        p = h.name.split(".")
        if p[10] in ("tavg-h2m-hxy-u", "tavg-u-hxy-u"):
            out.append(p[6])
    return sorted(set(out), key=lambda s: (int(s[1:s.index("i")]), s))


def load_cmip6_val(model, exp, member, var):
    mon = pd.read_csv(VAL / f"{model}.{exp}.{member}.{var}.csv")
    at = json.loads((VAL / f"{model}.{exp}.{member}.attrs.json").read_text())
    check_midmonth(mon, at.get(f"time_units_{var}", at["time_units"]), at["calendar"], f"{model} {exp} {var}")
    return annual(mon, var, at["calendar"], f"{model} {exp} {var}"), at


# ------------------------------------------------------- branch alignment
def _nearest_jan1(bt, units, calendar, what):
    d = cftime.num2date(bt, units, calendar)
    jan_this = cftime.datetime(d.year, 1, 1, calendar=calendar)
    jan_next = cftime.datetime(d.year + 1, 1, 1, calendar=calendar)
    dt_this, dt_next = (d - jan_this).days, (jan_next - d).days
    y, off = (d.year, dt_this) if dt_this <= dt_next else (d.year + 1, -dt_next)
    return y, off, str(d)


def branch_year(child_attrs, parent_attrs):
    """Parent model year in which the child branched, rounded to the nearest
    1 January. Returns (year, offset_days, decoded_date). |offset| > 15 d is flagged
    by the caller."""
    units = child_attrs.get("parent_time_units") or parent_attrs["time_units"]
    return _nearest_jan1(float(child_attrs["branch_time_in_parent"]), units,
                         parent_attrs["calendar"], "parent")


def child_start_year(child_attrs):
    bt = child_attrs.get("branch_time_in_child")
    if bt in (None, "None"):
        return None
    return _nearest_jan1(float(bt), child_attrs["time_units"], child_attrs["calendar"], "child")[0]


def piC_reference(pi, years, method, fit_years=None):
    """Reference piControl value for each parent year in `years`.

    rm21  : 21-yr centred running mean (calendar years; windows truncated only at
            the ends of the published record, minimum 11 yr)
    linear: OLS line fitted over `fit_years` (default: the contemporaneous years)
    mean  : mean over `fit_years`
    Returns (values, n_truncated_windows).
    """
    years = np.asarray(years)
    full = pi.reindex(range(int(pi.index.min()), int(pi.index.max()) + 1))
    if method == "rm21":
        need = range(int(years.min()) - 10, int(years.max()) + 11)
        inside = [y for y in need if full.index.min() <= y <= full.index.max()]
        if full.loc[inside].isna().any():
            raise ValueError("gap in piControl inside a running-mean window")
        rm = full.rolling(21, center=True, min_periods=11).mean()
        cnt = full.rolling(21, center=True, min_periods=1).count().reindex(years)
        vals = rm.reindex(years).to_numpy()
        if np.isnan(vals).any():
            raise ValueError("running mean undefined for some contemporaneous years")
        return vals, int((cnt < 21).sum())
    fy = np.asarray(fit_years if fit_years is not None else years)
    seg = full.loc[fy.min():fy.max()]
    if seg.isna().any() or len(seg) != fy.max() - fy.min() + 1:
        raise ValueError("piControl incomplete over the fit segment")
    if method == "linear":
        b = np.polyfit(seg.index.to_numpy(float), seg.to_numpy(), 1)
        return np.polyval(b, years.astype(float)), 0
    if method == "mean":
        return np.full(len(years), seg.mean()), 0
    raise ValueError(method)


def anomalies(child, pi, by, n, method, child_attrs=None, fit_len=None, shift=0):
    """Anomalies for child years 1..n, paired with parent years by+shift .. by+shift+n-1.
    fit_len: for linear/mean, length of the contemporaneous parent segment to fit
    (default n). shift: deliberate misalignment for the ±1 yr sensitivity."""
    idx = child.index.to_numpy()
    if len(idx) < n:
        raise ValueError(f"child has {len(idx)} < {n} years")
    if not np.all(np.diff(idx[:n]) == 1):
        raise ValueError("child years not contiguous")
    if child_attrs is not None:
        cs = child_start_year(child_attrs)
        if cs is not None and cs != idx[0]:
            raise ValueError(f"child starts {idx[0]} but branch_time_in_child decodes to {cs}")
    c = child.iloc[:n].to_numpy()
    by = by + shift
    py = np.arange(by, by + n)
    missing = sorted(set(py) - set(pi.index))
    if missing:
        raise ValueError(f"piControl lacks {len(missing)} of the {n} contemporaneous years "
                         f"({missing[0]}–{missing[-1]})")
    fy = np.arange(by, by + (fit_len or n))
    ref, ntr = piC_reference(pi, py, method, fit_years=fy)
    return c - ref, ntr


# ---------------------------------------------------------- metrics
def _fit(dT, dN):
    b, a = np.polyfit(dT, dN, 1)
    return a, b


def gregory(dT, dN, nboot=2000, block=10, seed=0):
    """OLS Gregory fit + fixed-design circular moving-block residual bootstrap
    (ΔT kept fixed; residual blocks wrap around, so early high-leverage years are
    not undersampled). The range reflects regression sampling only."""
    a, b = _fit(dT, dN)
    out = {"F2x": a / 2, "lambda": b, "ECS": -a / (2 * b)}
    resid = dN - (a + b * dT)
    rng = np.random.default_rng(seed)
    n = len(dT)
    nb = int(np.ceil(n / block))
    ecs = np.empty(nboot)
    for i in range(nboot):
        starts = rng.integers(0, n, nb)
        idx = (np.concatenate([np.arange(s, s + block) for s in starts]) % n)[:n]
        aa, bb = _fit(dT, a + b * dT + resid[idx])
        ecs[i] = -aa / (2 * bb)
    out["ECS_p05"], out["ECS_p95"] = np.percentile(ecs, [5, 95])
    return out


def tcr(dT):
    """Mean of years 61–80 (AR6 7.SM.4 'years 60 to 79', zero-indexed)."""
    return float(np.mean(dT[60:80]))


def pi_variability(pi, lengths=(10,), trend_lengths=(42, 72)):
    """Internal-variability σ from a linearly detrended piControl: σ of overlapping
    k-yr means, and σ of overlapping k-yr OLS trends (°C/decade)."""
    x = pi.to_numpy()
    t = np.arange(len(x), dtype=float)
    x = x - np.polyval(np.polyfit(t, x, 1), t)
    s = pd.Series(x)
    out = {f"sd_mean{k}": float(s.rolling(k).mean().std()) for k in lengths}
    for k in trend_lengths:
        tt = np.arange(k, dtype=float) - (k - 1) / 2
        w = tt / (tt ** 2).sum()
        tr = np.convolve(x, w[::-1], mode="valid") * 10
        out[f"sd_trend{k}"] = float(np.std(tr, ddof=1))
    return out
