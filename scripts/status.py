"""Completeness of derived global means vs the manifest, per model × experiment × variable."""
import pathlib, pandas as pd
ROOT = pathlib.Path(__file__).resolve().parents[1]
m = pd.read_csv(ROOT / "inventory" / "manifest_latest.csv")
m["done"] = [(ROOT / "derived" / "monthly" / k / (f[:-3] + ".csv")).exists() for k, f in zip(m.dskey, m.filename)]
g = m.groupby(["source_id", "experiment_id", "variable_id"]).agg(files=("done", "size"), done=("done", "sum"))
g["pct"] = (100 * g.done / g.files).round(0).astype(int)
pd.set_option("display.max_rows", 300)
print(g[g.pct < 100].to_string() if len(g[g.pct < 100]) else "all complete")
print(f"\nTOTAL {m.done.sum()}/{len(m)} files; {m.loc[m.done,'size'].sum()/1e9:.1f}/{m['size'].sum()/1e9:.1f} GB")
