"""Query ESGF CMIP7 STAC indexes (east+west) for monthly global variables needed
for GSAT and climate-sensitivity assessment. Writes inventory/*.csv."""
import sys, pathlib
import pandas as pd
import intake_esgf
from intake_esgf import ESGFCatalog

OUT = pathlib.Path(__file__).resolve().parents[1] / "inventory"
VARS = ["tas", "rsdt", "rsut", "rlut"]

cat = ESGFCatalog()
cat.search(project="CMIP7", variable_id=VARS, frequency="mon")
df = cat.df.copy()
print(df.columns.tolist())
print(len(df))
df.drop(columns=[c for c in df.columns if c in ("id",)], errors="ignore").to_csv(OUT / "cmip7_mon_tas_rad_raw.csv", index=False)
g = df.groupby(["source_id", "experiment_id", "variable_id"]).variant_label.nunique().unstack("variable_id").fillna(0).astype(int)
g.to_csv(OUT / "cmip7_mon_member_counts.csv")
pd.set_option("display.width", 200); pd.set_option("display.max_rows", 500)
print(g)
