import re
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator


# Input files
metadata_files = {
    "IMMU": "../data/clinical/IMMU-BC1_curated_export_20260423.xlsx",
    "SYNG": "../data/clinical/SYNG-BC1_curated_export_20260423.xlsx"}

# to characterize only samples for which we have DNA data
sbs_files = {
    "IMMU": "../SBS_assignment/IMMUcan_after_filtering/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
    "SYNG": "../SBS_assignment/SYNERGY_after_filtering/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"}

outdir = Path("cohort_characteristics")
outdir.mkdir(exist_ok=True)

def core_id(sample_id):
    """Extract the tumour sample identifier used for matching with metadata."""
    sample_id = str(sample_id).strip()
    parts = sample_id.split("_vs_")
    t_part = next((part for part in parts if "-FIXT-" in part.upper()), parts[0])

    return re.sub(r"-DNA-\d+$", "", t_part, flags=re.IGNORECASE)


def patient_key(sample_id):
    """Extract the patient identifier directly from the SBS profile identifier."""
    sample_id = core_id(sample_id)
    match = re.match(r"(.+?)-FIXT-", sample_id, flags=re.IGNORECASE)

    return match.group(1).upper() if match else pd.NA


def clean_category(value):
    if pd.isna(value) or str(value).strip() == "":
        return "Unknown"

    return re.sub(r"\s+", " ", str(value).strip())


def normalise_location(value, cohort):
    """Standardise sampling-site names without regard to letter case."""
    value = clean_category(value)
    key = re.sub(r"\s+", " ", value).casefold()

    if cohort == "IMMU" and "breast" in key:
        return "Breast"

    location_mapping = {}

    if cohort == "IMMU":
        location_mapping.update({
            "left cerebellar metastasis": "Brain",
            "right occipital brain metastasis": "Brain",
            "right pleura": "Pleura"})

    elif cohort == "SYNG":
        location_mapping.update({
            "distant skin": "Skin",
            "distant lymph nodes": "Lymph nodes",
            "pulmonary": "Lung",
            "hepatic": "Liver"})

    return location_mapping.get(key, key.capitalize())


def load_metadata(path, cohort):
    metadata = pd.read_excel(path, sheet_name="Samples")
    required_columns = ["sample_id", "sample_origin", "sampling_site"]
    metadata = metadata[required_columns].copy()
    metadata["cohort"] = cohort
    metadata["core_sample_id"] = metadata["sample_id"].map(core_id)
    metadata["sample_origin"] = metadata["sample_origin"].map(clean_category)
    metadata["sampling_site"] = metadata["sampling_site"].map(lambda value: normalise_location(value, cohort))

    return metadata.drop_duplicates(subset=["cohort", "core_sample_id"], keep="first")


def load_sbs_profiles(path, cohort):
    profiles = pd.read_csv(path, sep="\t", usecols=[0])
    profiles = profiles.rename(columns={profiles.columns[0]: "activity_sample_id"})
    profiles["activity_sample_id"] = (profiles["activity_sample_id"].astype(str).str.strip())
    profiles["cohort"] = cohort
    profiles["core_sample_id"] = profiles["activity_sample_id"].map(core_id)
    profiles["patient_key"] = profiles["activity_sample_id"].map(patient_key)

    return profiles.drop_duplicates(subset=["cohort", "activity_sample_id"])


def prepare_cohort(cohort):
    metadata = load_metadata(metadata_files[cohort], cohort)
    profiles = load_sbs_profiles(sbs_files[cohort], cohort)

    data = profiles.merge(metadata, on=["cohort", "core_sample_id"], how="left", validate="many_to_one")
    missing_metadata = data["sample_id"].isna()
    data.loc[missing_metadata, "sample_origin"] = "Missing metadata"
    data.loc[missing_metadata, "sampling_site"] = "Missing metadata"

    origin = data["sample_origin"].str.casefold()
    data["origin_group"] = "Other"
    data.loc[origin.str.contains("primary", na=False), "origin_group"] = "Primary"
    data.loc[origin.str.contains("metasta", na=False), "origin_group"] = "Metastatic"
    data.loc[(cohort == "IMMU") & origin.str.contains("recur", na=False), "origin_group"] = "Recurrent"
    data.loc[missing_metadata, "origin_group"] = "Missing metadata"

    return data

cohort_data = {
    "IMMU": prepare_cohort("IMMU"),
    "SYNG": prepare_cohort("SYNG")
}

all_locations = sorted({location
    for data in cohort_data.values()
    for location in data.loc[data["origin_group"].isin(["Primary", "Metastatic", "Recurrent"]),"sampling_site"].unique()})

location_colors = {
    "Breast": "#64A8ED",
    "Brain": "#DA8BC3",
    "Pleura": "#F5F187",
    "Chest wall": "#C44E52",
    "Skin": "#F5F187",
    "Lymph nodes": "#DD8452",
    "Liver": "#65A165",
    "Lung": "#755744",
    "Other": "#F5112E",
    "Hepatogastric": "#65A165"
}

missing_color = "#BDBDBD"

def plot_cohort(data, cohort, cohort_name, output_name):
    total_samples = data["activity_sample_id"].nunique()
    total_patients = data["patient_key"].nunique()

    def get_location_counts(origin_group):
        return (data.loc[data["origin_group"] == origin_group].groupby("sampling_site")["activity_sample_id"].nunique().sort_values(ascending=False))

    primary_by_location = get_location_counts("Primary")
    metastatic_by_location = get_location_counts("Metastatic")
    recurrent_by_location = get_location_counts("Recurrent")

    missing_count = data.loc[data["origin_group"] == "Missing metadata", "activity_sample_id"].nunique()

    categories = ["Primary", "Metastatic"]
    location_groups = [primary_by_location, metastatic_by_location]

    if not recurrent_by_location.empty:
        categories.append("Recurrent")
        location_groups.append(recurrent_by_location)

    if missing_count > 0:
        categories.append("Missing\nmetadata")
        location_groups.append(None)

    fig, ax = plt.subplots(figsize=(8, 6.5))

    def draw_stacked_bar(position, location_counts):
        bottom = 0
        for location, count in location_counts.sort_values(ascending=False).items():
            count = int(count)

            ax.bar(
                position,
                count,
                bottom=bottom,
                width=0.65,
                color=location_colors[location],
                edgecolor="black",
                linewidth=0.9)

            if len(location_counts) > 1 and ((cohort == 'SYNG' and count >= 3) or (cohort == 'IMMU' and count > 1)):
                ax.text(
                    position,
                    bottom + count / 2.3,
                    str(count),
                    ha="center",
                    va="center",
                    fontsize=12)

            bottom += count

        return bottom

    totals = []

    for position, location_counts in enumerate(location_groups):
        if location_counts is None:
            ax.bar(
                position,
                missing_count,
                width=0.65,
                color=missing_color,
                edgecolor="black",
                linewidth=0.9)
            totals.append(missing_count)
        else:
            totals.append(
                draw_stacked_bar(position, location_counts))

    max_height = max(totals)
    label_offset = max_height * 0.02

    for position, total in enumerate(totals):
        ax.text(
            position,
            total + label_offset,
            str(total),
            ha="center",
            va="bottom",
            fontsize=13,
            fontweight="bold")

    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(categories, fontsize=15)
    ax.set_ylabel("Number of samples", fontsize=15, labelpad=10)

    ax.set_title(
        f"{cohort_name}\n"
        f"(n = {total_samples} samples, {total_patients} patients)",
        fontsize=16,
        fontweight="bold",
        pad=14)

    ax.set_ylim(0, max_height + max(3, max_height * 0.12))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.tick_params(axis="y", labelsize=13)
    ax.grid(axis="y", alpha=0.2)
    ax.set_axisbelow(True)


    metastatic_order = (metastatic_by_location.sort_values(ascending=False).index.tolist())

    extra_locations = []

    for counts in [primary_by_location, recurrent_by_location,]:
        for location in counts.sort_values(ascending=False).index:
            if (location not in metastatic_order and location not in extra_locations):
                extra_locations.append(location)

    legend_order = list(reversed(metastatic_order)) + extra_locations

    legend_handles = [
        Patch(
            facecolor=location_colors[location],
            edgecolor="black",
            label=location)
        for location in legend_order]

    ax.legend(
        handles=legend_handles,
        title="Sampling site",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=True,
        fontsize=13,
        title_fontsize=14
    )

    fig.tight_layout(rect=[0, 0, 0.96, 1])
    plt.savefig(
        outdir / output_name,
        dpi=300,
        bbox_inches="tight"
    )
    #plt.show()

    print(
        f"\n{cohort_name}: "
        f"{total_samples} samples, {total_patients} patients")

    for label, counts in [
        ("Primary", primary_by_location),
        ("Metastatic", metastatic_by_location),
        ("Recurrent", recurrent_by_location)]:
        if not counts.empty:
            print(f"\n{label} locations:")
            print(counts.to_string())

    if missing_count > 0:
        print(f"\nMissing metadata: {missing_count}")

plot_cohort(
    cohort_data["IMMU"],
    cohort="IMMU",
    cohort_name="IMMUcan",
    output_name="immucan_sample_origin_and_locations.png")

plot_cohort(
    cohort_data["SYNG"],
    cohort="SYNG",
    cohort_name="Synergy",
    output_name="synergy_sample_origin_and_locations.png")