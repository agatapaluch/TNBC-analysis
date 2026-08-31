from pathlib import Path
import numpy as np
import pandas as pd


METADATA_FILE = Path("../../data/clinical/SYNG-BC1_curated_export_20260423.xlsx")

SBS_ASSIGNMENT_FILE = Path(
    "../../SBS_assignment/SYNERGY_after_filtering/"
    "Assignment_Solution/Activities/Assignment_Solution_Activities.txt")

DBS_ASSIGNMENT_FILE = Path(
    "../../DBS_ID_assignment/DBS/DBS_SYNERGY_after_filtering/"
    "Assignment_Solution/Activities/Assignment_Solution_Activities.txt")

ID_ASSIGNMENT_FILE = Path(
    "../../DBS_ID_assignment/ID/ID_SYNERGY_after_filtering/"
    "Assignment_Solution/Activities/Assignment_Solution_Activities.txt")

OUTPUT_DIR = Path("SYNERGY_outcome_groups")
OUTPUT_DIR.mkdir(exist_ok=True)

UNFAVORABLE_MAX_MONTHS = 24
FAVORABLE_MIN_MONTHS = 60

def normalize_sample_id(sample_id):
    sample_id = str(sample_id).strip()
    sample_id = sample_id.split("_vs_", maxsplit=1)[0]
    sample_id = pd.Series([sample_id]).str.replace(r"-DNA-\d+$", "", regex=True).iloc[0]
    return sample_id


# Read SBS samples

sbs = pd.read_csv(SBS_ASSIGNMENT_FILE, sep="\t")
sbs["sample_id_clean"] = sbs["Samples"].apply(normalize_sample_id)

dna_samples = set(sbs["sample_id_clean"])

print(f"\nSBS samples available: {len(dna_samples)}")


# Read metadata

samples = pd.read_excel(METADATA_FILE, sheet_name="Samples")
samples = samples[["patient_id", "sample_id", "sampling_day"]].copy()

samples["patient_id"] = samples["patient_id"].astype("string").str.strip()
samples["sample_id"] = samples["sample_id"].astype("string").str.strip()
samples["sample_id_clean"] = samples["sample_id"].apply(normalize_sample_id)
samples["sampling_day"] = pd.to_numeric(samples["sampling_day"], errors="coerce")

samples = samples.dropna(subset=["patient_id", "sampling_day"])


# Keep only samples with corresponding SBS data

dna_metadata = samples[
    samples["sample_id_clean"].isin(dna_samples)
].copy()

print(f"Metadata samples with SBS data: {len(dna_metadata)}")
print(f"Patients with DNA/SBS data: {dna_metadata['patient_id'].nunique()}")


# Select the latest sample

patients = (
    dna_metadata
    .sort_values(["patient_id", "sampling_day", "sample_id_clean"])
    .drop_duplicates("patient_id", keep="last")
    .reset_index(drop=True)
)

print(f"Patients after selecting latest DNA sample: {len(patients)}")


# Read survival data

survival = pd.read_excel(METADATA_FILE, sheet_name="Patients")
survival = survival[["patient_id", "last_FU_OS_days", "last_FU_OS_event"]].copy()

survival["patient_id"] = survival["patient_id"].astype("string").str.strip()
survival["last_FU_OS_days"] = pd.to_numeric(survival["last_FU_OS_days"], errors="coerce")

event = (survival["last_FU_OS_event"].astype("string").str.strip().str.lower())

survival["OS_EVENT"] = np.select(
    [event.isin(["1", "1.0", "dead"]), event.isin(["0", "0.0", "alive", "lost to follow-up", "lost to follow up"])],
    [1, 0],
    default=np.nan)

survival["OS_MONTHS"] = survival["last_FU_OS_days"] / 30.4375

survival = survival[["patient_id", "last_FU_OS_days", "OS_MONTHS", "OS_EVENT"]].copy()

if survival["patient_id"].duplicated().any():
    duplicated = survival.loc[survival["patient_id"].duplicated(keep=False), "patient_id"].unique()

    raise ValueError(
        f"Duplicated patient IDs in Patients sheet: {len(duplicated)}")

print(f"Patients in Patients sheet: {len(survival)}")


# Merge data

df = patients.merge(
    survival,
    on="patient_id",
    how="inner",
    validate="one_to_one")

valid = (
    df["OS_MONTHS"].notna() &
    (df["OS_MONTHS"] > 0) &
    df["OS_EVENT"].notna())

print(f"\nPatients before OS filtering: {len(df)}")
print(f"Patients with valid OS: {valid.sum()}")
print(f"Patients excluded: {(~valid).sum()}")

df = df[valid].copy()

# Define outcome groups

unfavorable = df[
    (df["OS_EVENT"] == 1) &
    (df["OS_MONTHS"] <= UNFAVORABLE_MAX_MONTHS)].copy()

favorable = df[
    df["OS_MONTHS"] > FAVORABLE_MIN_MONTHS].copy()

unfavorable["OUTCOME_GROUP"] = "Unfavorable"
favorable["OUTCOME_GROUP"] = "Favorable"


# Patients not meeting either definition

assigned_patients = (
    set(unfavorable["patient_id"]) |
    set(favorable["patient_id"])
)

not_assigned = df[
    ~df["patient_id"].isin(assigned_patients)
].copy()

not_assigned["OUTCOME_GROUP"] = "Not assigned"


# Save groups

columns = [
    "patient_id",
    "sample_id_clean",
    "sampling_day",
    "last_FU_OS_days",
    "OS_MONTHS",
    "OS_EVENT",
    "OUTCOME_GROUP"]

unfavorable = unfavorable[columns]
favorable = favorable[columns]
not_assigned = not_assigned[columns]

unfavorable.to_csv(OUTPUT_DIR / "SYNERGY_unfavorable.csv", index=False)
favorable.to_csv(OUTPUT_DIR / "SYNERGY_favorable.csv", index=False)
not_assigned.to_csv(OUTPUT_DIR / "SYNERGY_not_assigned.csv", index=False)


print("\nFINAL OUTCOME GROUPS")

print(
    f"Unfavorable: death <= {UNFAVORABLE_MAX_MONTHS} months"
    f"n = {len(unfavorable)}"
)

print(
    f"Favorable: OS > {FAVORABLE_MIN_MONTHS} months"
    f"n = {len(favorable)}"
)

print(
    f"Not assigned: "
    f"n = {len(not_assigned)}"
)

print("\nEvent status in favorable group:")
print(favorable["OS_EVENT"].value_counts().sort_index())

print("\nOS summary:")
print(
    pd.DataFrame({
        "Unfavorable": unfavorable["OS_MONTHS"].describe(),
        "Favorable": favorable["OS_MONTHS"].describe()}))


# Split SBS, DBS and ID files

favorable_samples = set(favorable["sample_id_clean"])
unfavorable_samples = set(unfavorable["sample_id_clean"])


def split_assignment_file(assignment_file, modality):
    activities = pd.read_csv(assignment_file, sep="\t")
    activities["sample_id_clean"] = activities["Samples"].apply(normalize_sample_id)

    favorable_activities = activities[
        activities["sample_id_clean"].isin(favorable_samples)].copy()

    unfavorable_activities = activities[
        activities["sample_id_clean"].isin(unfavorable_samples)].copy()

    favorable_activities = favorable_activities.drop(
        columns="sample_id_clean")

    unfavorable_activities = unfavorable_activities.drop(
        columns="sample_id_clean")

    favorable_file = OUTPUT_DIR / f"{modality}_favorable.txt"
    unfavorable_file = OUTPUT_DIR / f"{modality}_unfavorable.txt"

    favorable_activities.to_csv(
        favorable_file,
        sep="\t",
        index=False)

    unfavorable_activities.to_csv(
        unfavorable_file,
        sep="\t",
        index=False)

    print(f"\n{modality}")
    print(f"Favorable: {len(favorable_activities)} samples")
    print(f"Unfavorable: {len(unfavorable_activities)} samples")


split_assignment_file(SBS_ASSIGNMENT_FILE, "SBS")
split_assignment_file(DBS_ASSIGNMENT_FILE, "DBS")
split_assignment_file(ID_ASSIGNMENT_FILE, "ID")