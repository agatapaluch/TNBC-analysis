from pathlib import Path
import pandas as pd


# Input files

SYNG_SBS_FILE = Path("../SBS_assignment/SYNERGY_after_filtering/Assignment_Solution/Activities/Assignment_Solution_Activities.txt")
SYNG_ID_FILE = Path("../DBS_ID_assignment/ID/ID_SYNERGY_after_filtering/Assignment_Solution/Activities/Assignment_Solution_Activities.txt")
SYNG_METADATA_FILE = Path("../data/clinical/SYNG-BC1_curated_export_20260423.xlsx")

SBS3_OUTPUT_DIR = Path("SBS3")
HRD_OUTPUT_DIR = Path("SBS3_ID6")
MMR_OUTPUT_DIR = Path("SBS6_15_44")

SBS3_OUTPUT_DIR.mkdir(exist_ok=True)
HRD_OUTPUT_DIR.mkdir(exist_ok=True)
MMR_OUTPUT_DIR.mkdir(exist_ok=True)

SBS3_ALL_SAMPLES = SBS3_OUTPUT_DIR / "SBS3_all_samples.csv"
SBS3_PATIENT_LEVEL = SBS3_OUTPUT_DIR / "SBS3_patient_level.csv"

HRD_ALL_SAMPLES = HRD_OUTPUT_DIR / "SBS3_ID6_all_samples.csv"
HRD_PATIENT_LEVEL = HRD_OUTPUT_DIR / "SBS3_ID6_patient_level.csv"

MMR_ALL_SAMPLES = MMR_OUTPUT_DIR / "SBS6_15_44_all_samples.csv"
MMR_PATIENT_LEVEL = MMR_OUTPUT_DIR / "SBS6_15_44_patient_level.csv"

SBS3_THRESHOLD = 0


def clean_sample_id(sample):
    sample = str(sample).strip().split("_vs_", maxsplit=1)[0]
    return pd.Series([sample]).str.replace(r"-DNA-\d+$", "", regex=True).iloc[0]


def read_sbs_file(path):
    df = pd.read_csv(path, sep="\t")
    df["sample_id"] = df["Samples"].apply(clean_sample_id)

    signatures = ["SBS3", "SBS6", "SBS15", "SBS44"]

    for signature in signatures:
        if signature not in df.columns:
            print(f"{signature} is absent from {path.name}; setting activity to 0.")
            df[signature] = 0

        df[signature] = pd.to_numeric(df[signature], errors="coerce").fillna(0)

    for signature in signatures:
        duplicate_check = df.groupby("sample_id")[signature].nunique(dropna=False)
        conflicts = duplicate_check[duplicate_check > 1]

        if not conflicts.empty:
            print(f"\nWARNING: duplicated samples with different {signature} activities:")
            print(conflicts)

    return df[["sample_id", "SBS3", "SBS6", "SBS15", "SBS44"]].drop_duplicates("sample_id", keep="first")


def read_id_file(path):
    df = pd.read_csv(path, sep="\t")
    df["sample_id"] = df["Samples"].apply(clean_sample_id)

    if "ID6" not in df.columns:
        print(f"ID6 is absent from {path.name}; setting activity to 0.")
        df["ID6"] = 0

    df["ID6"] = pd.to_numeric(df["ID6"], errors="coerce").fillna(0)

    duplicate_check = df.groupby("sample_id")["ID6"].nunique(dropna=False)
    conflicts = duplicate_check[duplicate_check > 1]

    if not conflicts.empty:
        print("\nWARNING: duplicated samples with different ID6 activities:")
        print(conflicts)

    return df[["sample_id", "ID6"]].drop_duplicates("sample_id", keep="first")


def get_marker(row, signatures):
    markers = [signature for signature in signatures if row[signature] > 0]
    return " + ".join(markers) if markers else "none"


def collapse_to_patient_level(all_samples, status_col, output_label):
    known = all_samples[~all_samples[status_col].str.contains("unknown", case=False, na=False)].copy()

    sample_counts = known.groupby("patient_id").size().rename("n_samples").reset_index()
    multiple_samples = (sample_counts["n_samples"] > 1).sum()

    print(f"\n{output_label}:")
    print(f"Patients with known status: {known['patient_id'].nunique()}")
    print(f"Patients with multiple evaluable samples: {multiple_samples}")

    patient_level = (
        known.sort_values(["patient_id", "sampling_day", "sample_id"])
        .drop_duplicates("patient_id", keep="last")
        .reset_index(drop=True))

    print("Latest-sample status:")
    print(patient_level[status_col].value_counts())

    return patient_level


# Read data

sbs = read_sbs_file(SYNG_SBS_FILE)
ids = read_id_file(SYNG_ID_FILE)

metadata = pd.read_excel(SYNG_METADATA_FILE, sheet_name="Samples")
metadata = metadata[["patient_id", "sample_id", "sampling_day"]].copy()

metadata["patient_id"] = metadata["patient_id"].astype("string").str.strip()
metadata["sample_id"] = metadata["sample_id"].astype("string").str.strip()
metadata["sample_id"] = metadata["sample_id"].apply(clean_sample_id)
metadata["sampling_day"] = pd.to_numeric(metadata["sampling_day"], errors="coerce")
metadata = metadata.drop_duplicates("sample_id", keep="first")

all_samples = sbs.merge(metadata, on="sample_id", how="left", validate="one_to_one")

missing_metadata = all_samples[
    all_samples["patient_id"].isna() | all_samples["sampling_day"].isna()].copy()

if not missing_metadata.empty:
    print("\nSamples excluded because of missing metadata:")
    print(missing_metadata[["sample_id", "SBS3", "SBS6", "SBS15", "SBS44"]].to_string(index=False))

all_samples = all_samples[
    all_samples["patient_id"].notna() & all_samples["sampling_day"].notna()].copy()

all_samples["cohort"] = "Synergy"

print("\nSynergy SBS samples with metadata:", len(all_samples))
print("Synergy patients with SBS data:", all_samples["patient_id"].nunique())

# SBS3

sbs3_samples = all_samples[["cohort", "patient_id", "sample_id", "sampling_day", "SBS3"]].copy()

sbs3_samples["SBS3_status"] = "SBS3 negative"
sbs3_samples.loc[sbs3_samples["SBS3"] > SBS3_THRESHOLD, "SBS3_status"] = "SBS3 positive"

sbs3_patient_level = collapse_to_patient_level(
    sbs3_samples, status_col="SBS3_status", output_label="SBS3")

sbs3_samples.to_csv(SBS3_ALL_SAMPLES, index=False)
sbs3_patient_level.to_csv(SBS3_PATIENT_LEVEL, index=False)

print("\nSBS3 status - all samples:")
print(sbs3_samples["SBS3_status"].value_counts())

print("\nSBS3 status - patient level, latest sample:")
print(sbs3_patient_level["SBS3_status"].value_counts())


# SBS3 + ID6

hrd_samples = all_samples[
    ["cohort", "patient_id", "sample_id", "sampling_day", "SBS3"]].merge(ids, on="sample_id", how="left", validate="one_to_one")

hrd_samples["ID6_available"] = hrd_samples["ID6"].notna()

sbs3_positive = hrd_samples["SBS3"] > SBS3_THRESHOLD
id6_positive = hrd_samples["ID6"].fillna(0) > 0

hrd_samples["HRD_status"] = "HRD unknown"
hrd_samples.loc[sbs3_positive | id6_positive, "HRD_status"] = "HRD positive"
hrd_samples.loc[
    (~sbs3_positive) & hrd_samples["ID6_available"] & (~id6_positive), "HRD_status"] = "HRD negative"

hrd_samples["HRD_marker"] = hrd_samples.apply(
    lambda row: get_marker(row.fillna(0), ["SBS3", "ID6"]), axis=1)

hrd_patient_level = collapse_to_patient_level(
    hrd_samples, status_col="HRD_status", output_label="SBS3/ID6")

hrd_samples.to_csv(HRD_ALL_SAMPLES, index=False)
hrd_patient_level.to_csv(HRD_PATIENT_LEVEL, index=False)

print("\nSBS3/ID6 status - all samples:")
print(hrd_samples["HRD_status"].value_counts())

print("\nSBS3/ID6 status - patient level, latest evaluable sample:")
print(hrd_patient_level["HRD_status"].value_counts())

print("\nMarkers among SBS3/ID6-positive patients:")
print(hrd_patient_level.loc[hrd_patient_level["HRD_status"] == "HRD positive", "HRD_marker"].value_counts())


# SBS6 + SBS15 + SBS44

mmr_samples = all_samples[
    ["cohort", "patient_id", "sample_id", "sampling_day", "SBS6", "SBS15", "SBS44"]].copy()

mmr_positive = (
    (mmr_samples["SBS6"] > 0) |
    (mmr_samples["SBS15"] > 0) |
    (mmr_samples["SBS44"] > 0))

mmr_samples["MMR_status"] = "SBS6/15/44 negative"
mmr_samples.loc[mmr_positive, "MMR_status"] = "SBS6/15/44 positive"

mmr_samples["MMR_marker"] = mmr_samples.apply(
    lambda row: get_marker(row, ["SBS6", "SBS15", "SBS44"]), axis=1)

mmr_patient_level = collapse_to_patient_level(
    mmr_samples, status_col="MMR_status", output_label="SBS6/15/44")

mmr_samples.to_csv(MMR_ALL_SAMPLES, index=False)
mmr_patient_level.to_csv(MMR_PATIENT_LEVEL, index=False)

print("\nSBS6/15/44 status - all samples:")
print(mmr_samples["MMR_status"].value_counts())

print("\nSBS6/15/44 status - patient level, latest sample:")
print(mmr_patient_level["MMR_status"].value_counts())

print("\nMarkers among SBS6/15/44-positive patients:")
print(mmr_patient_level.loc[mmr_patient_level["MMR_status"] == "SBS6/15/44 positive", "MMR_marker"].value_counts())

# Output

print("\nSaved:")
print(SBS3_ALL_SAMPLES)
print(SBS3_PATIENT_LEVEL)
print(HRD_ALL_SAMPLES)
print(HRD_PATIENT_LEVEL)
print(MMR_ALL_SAMPLES)
print(MMR_PATIENT_LEVEL)