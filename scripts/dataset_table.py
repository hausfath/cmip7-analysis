"""Citation / provenance table for every CMIP7 dataset used (METHODS.md §4).

  python -I scripts/dataset_table.py  ->  outputs/cmip7_dataset_table.csv
"""
import json, pathlib
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
man = pd.read_csv(ROOT / "inventory" / "manifest_latest.csv")
log = pd.read_csv(ROOT / "derived" / "acquire_log.csv")
log = log.sort_values("utc").drop_duplicates(["dskey", "filename"], keep="last")
rows = []
for dk, g in man.groupby("dskey"):
    a = json.loads((ROOT / "derived" / "attrs" / f"{dk}.json").read_text())
    st = log[log.dskey == dk].status.value_counts().to_dict()
    p = dk.split(".")
    rows.append({"dataset_id": f"{dk}.v{g.version.iloc[0]}", "institution_id": p[3], "source_id": p[4],
                 "experiment_id": p[5], "variant_label": p[6], "variable": f"{p[9]}_{p[10]}", "grid": p[11],
                 "version": g.version.iloc[0], "n_files": len(g), "size_GB": round(g["size"].sum() / 1e9, 3),
                 "files_sha256_ok": st.get("ok", 0), "files_ok_stale_index": st.get("ok-stale", 0),
                 "license_id": a.get("license_id", ""), "tracking_id_first_file": a.get("tracking_id", ""),
                 "creation_date_first_file": a.get("creation_date", "")})
df = pd.DataFrame(rows)
df.to_csv(ROOT / "outputs" / "cmip7_dataset_table.csv", index=False)
print(len(df), "datasets;", df.license_id.value_counts().to_dict(), "; stale-index files:",
      int(df.files_ok_stale_index.sum()))
