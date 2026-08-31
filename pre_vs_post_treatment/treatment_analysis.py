import re
from pathlib import Path
import numpy as np
import pandas as pd


# load files

metadata_files = {
    "IMMU": "../data/clinical/IMMU-BC1_curated_export_20260423.xlsx",
    "SYNG": "../data/clinical/SYNG-BC1_curated_export_20260423.xlsx"}

sbs_files = {
    "IMMU": "../SBS_assignment/IMMUcan_after_filtering/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
    "SYNG": "../SBS_assignment/SYNERGY_after_filtering/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"}

outdir = Path("treatment_signature_analysis")
outdir.mkdir(exist_ok=True)


# helper functions

def core_id(x):
    x = str(x).strip()
    parts = x.split("_vs_")
    x = next((p for p in parts if "-FIXT-" in p), parts[0])
    return re.sub(r"-DNA-\d+$", "", x)


def load_metadata(path, cohort):
    metadata = pd.read_excel(path, sheet_name=None)
    samples, meds = metadata["Samples"].copy(), metadata["Medications"].copy()

    samples["cohort"] = cohort
    samples["sampling_day"] = pd.to_numeric(samples["sampling_day"], errors="coerce")
    samples = samples[samples["sample_id"].str.contains("FIXT", na=False)].copy()
    samples["core_sample_id"] = samples["sample_id"].map(core_id)

    start_cols = [c for c in ["drug_start_day", "medication_start_day"] if c in meds.columns]
    meds[start_cols] = meds[start_cols].apply(pd.to_numeric, errors="coerce")
    meds["treatment_start_day"] = meds[start_cols].min(axis=1)
    meds["cohort"] = cohort

    first_treatment = meds.groupby("patient_id")["treatment_start_day"].min().rename("first_treatment_day")
    samples = samples.merge(first_treatment, on="patient_id", how="left")
    samples["timepoint"] = np.select(
        [samples["sampling_day"] <= samples["first_treatment_day"],
         samples["sampling_day"] > samples["first_treatment_day"]],
        ["pre", "post"], default="unknown"
    )
    return samples, meds


def load_sbs_samples(path, cohort):
    df = pd.read_csv(path, sep="\t", usecols=[0])
    df = df.rename(columns={df.columns[0]: "activity_sample_id"})
    df["core_sample_id"] = df["activity_sample_id"].map(core_id)
    df["cohort"] = cohort
    return df[["cohort", "core_sample_id"]].drop_duplicates()


def select_pairs(samples):
    rows = []

    for (cohort, patient), group in samples.groupby(["cohort", "patient_id"]):
        pre, post = group[group["timepoint"] == "pre"], group[group["timepoint"] == "post"]
        if pre.empty or post.empty:
            continue

        pre, post = pre.sort_values("sampling_day").iloc[-1], post.sort_values("sampling_day").iloc[-1]
        rows.append({
            "cohort": cohort, "patient_id": patient, "first_treatment_day": pre["first_treatment_day"],
            "pre_sample": pre["core_sample_id"], "post_sample": post["core_sample_id"],
            "pre_sampling_day": pre["sampling_day"], "post_sampling_day": post["sampling_day"],
            "days_pre": pre["first_treatment_day"] - pre["sampling_day"],
            "days_post": post["sampling_day"] - post["first_treatment_day"]
        })

    return pd.DataFrame(rows)


# metadata and samples availability

sample_tables, med_tables = [], []

for cohort, path in metadata_files.items():
    samples, meds = load_metadata(path, cohort)
    sample_tables.append(samples)
    med_tables.append(meds)

samples = pd.concat(sample_tables, ignore_index=True)
meds = pd.concat(med_tables, ignore_index=True)
sbs_samples = pd.concat([load_sbs_samples(path, cohort) for cohort, path in sbs_files.items()], ignore_index=True)

eligible = samples.merge(sbs_samples, on=["cohort", "core_sample_id"], how="inner")
pairs = select_pairs(eligible)

if pairs.empty:
    raise ValueError("No pre/post pairs could be matched to SBS activities.")


# find treatments between pre and post samples

treatment_rows = []

for _, pair in pairs.iterrows():
    pm = meds[(meds["cohort"] == pair["cohort"]) & (meds["patient_id"] == pair["patient_id"])].copy()
    pm = pm[(pm["treatment_start_day"] > pair["pre_sampling_day"]) & (pm["treatment_start_day"] <= pair["post_sampling_day"])]

    treatment_rows.append({
        "cohort": pair["cohort"], "patient_id": pair["patient_id"],
        "drugs_between": "; ".join(sorted(pm["drug"].dropna().astype(str).unique())) if "drug" in pm else "",
        "drug_classes_between": "; ".join(sorted(pm["drug_class"].dropna().astype(str).unique())) if "drug_class" in pm else ""
    })

pairs = pairs.merge(pd.DataFrame(treatment_rows), on=["cohort", "patient_id"], how="left")
pairs = pairs.sort_values(["cohort", "patient_id"]).reset_index(drop=True)
pairs["pair_label"] = [f"Pair {i:02d}" for i in range(1, len(pairs) + 1)]

pairs.to_csv(outdir / "pair_mapping.csv", index=False)

print(f"\nSaved {len(pairs)} pre/post pairs to {outdir / 'pair_mapping.csv'}")
print(pairs.groupby("cohort")["patient_id"].nunique())
