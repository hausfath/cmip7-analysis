"""Build a file-level manifest (URLs, sizes, checksums) of the CMIP7 monthly
datasets used in this assessment, from the ESGF CMIP7 STAC indexes.

Output: inventory/manifest_<YYYY-MM-DD>.csv and inventory/manifest_latest.csv
"""
import datetime, json, pathlib
import pandas as pd
from intake_esgf import ESGFCatalog

ROOT = pathlib.Path(__file__).resolve().parents[1]
INV = ROOT / "inventory"

# Experiments needing only tas
TAS_EXPS = ["historical", "esm-hist", "piControl", "esm-piControl", "abrupt-4xCO2", "1pctCO2",
            "scen7-h", "scen7-hl", "scen7-m", "scen7-ml", "scen7-l", "scen7-vl", "scen7-ln",
            "esm-scen7-h", "esm-scen7-hl", "esm-scen7-m", "esm-scen7-ml", "esm-scen7-l",
            "esm-scen7-vl", "esm-scen7-ln"]
# TOA radiation needed for the Gregory regression
RAD_EXPS = ["piControl", "abrupt-4xCO2"]

KEYF = ["drs_specs", "project", "activity_id", "institution_id", "source_id", "experiment_id",
        "variant_label", "region", "frequency", "variable_id", "branding", "grid_label"]


def query(variable, branding, exps):
    cat = ESGFCatalog()
    cat.search(project="CMIP7", variable_id=variable, frequency="mon", region="glb",
               experiment_id=exps, variable_branding_suffix=branding)
    rows = []
    for i in cat._get_file_info(quiet=True):
        urls = sorted(set(i.get("HTTPServer", [])))
        row = dict(zip(KEYF, i["key"].split(".")))
        row.update(version=i["dataset_id"].split("|")[0].split(".")[-1], filename=pathlib.Path(i["path"]).name,
                   size=i.get("size"), checksum=i.get("checksum"), checksum_type=i.get("checksum_type"),
                   file_start=str(i.get("file_start")), file_end=str(i.get("file_end")),
                   urls=json.dumps(urls))
        rows.append(row)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    df = pd.concat([query("tas", "tavg-h2m-hxy-u", TAS_EXPS),
                    query(["rsdt", "rsut", "rlut"], "tavg-u-hxy-u", RAD_EXPS)], ignore_index=True)
    # Keep only the latest version of each dataset
    df["dskey"] = df[KEYF].agg(".".join, axis=1)
    latest = df.groupby("dskey").version.transform("max")
    df = df[df.version == latest].drop_duplicates(["dskey", "filename"]).sort_values(["dskey", "filename"])
    today = datetime.date.today().isoformat()
    df.to_csv(INV / f"manifest_{today}.csv", index=False)
    df.to_csv(INV / "manifest_latest.csv", index=False)
    print(len(df), "files;", df.dskey.nunique(), "datasets;", round(df["size"].sum() / 1e9, 1), "GB")
