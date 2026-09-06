import re
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests


# Input files

subtype_files = {
    "IMMU": "../immucan_result.csv",
    "SYNG": "../synergy_result.csv"
}

sbs_files = {
    "IMMU": "../../SBS_assignment/IMMUcan_after_filtering/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
    "SYNG": "../../SBS_assignment/SYNERGY_after_filtering/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
}

dbs_files = {
    "IMMU": "../../DBS_ID_assignment/DBS/DBS_IMMUcan_after_filtering/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
    "SYNG": "../../DBS_ID_assignment/DBS/DBS_SYNERGY_after_filtering/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
}

id_files = {
    "IMMU": "../../DBS_ID_assignment/ID/ID_IMMUcan_after_filtering/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
    "SYNG": "../../DBS_ID_assignment/ID/ID_SYNERGY_after_filtering/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
}

outdir = Path("subtype_activity_boxplots")
outdir.mkdir(exist_ok=True)

# Settings

SUBTYPE_ORDER = ["BL1", "BL2", "IM", "M", "MSL", "LAR"]

POSITIVE_ONLY = True
LOG_SCALE = True

MIN_TOTAL_POSITIVE = 4
MIN_POSITIVE_PER_SUBTYPE = 2
MIN_SUBTYPES = 2
MIN_POSITIVE_FOR_BOX = 2
MIN_N_FOR_TEST = 2

ALPHA = 0.05

Y_MIN = None
Y_MAX = None

def core_id(sample_id):
    sample_id = str(sample_id).strip().replace(".", "-")
    parts = sample_id.split("_vs_")
    tumour_part = next((part for part in parts if "-FIXT-" in part.upper()), parts[0])
    tumour_part = re.sub(r"-(RNA|DNA)-\d+$", "", tumour_part, flags=re.IGNORECASE)
    return tumour_part


def read_subtypes(path, cohort):
    df = pd.read_csv(path)

    if "Unnamed: 0" in df.columns:
        sample_col = "Unnamed: 0"
    elif "sample" in df.columns:
        sample_col = "sample"
    elif "sample_id" in df.columns:
        sample_col = "sample_id"
    else:
        raise ValueError(
            f"Could not identify sample ID column in {path}.\n"
            f"Available columns: {list(df.columns)}")

    if "subtype" not in df.columns:
        raise ValueError(
            f"'subtype' column not found in {path}.\n"
            f"Available columns: {list(df.columns)}")

    df = df[[sample_col, "subtype"]].copy()
    df.columns = ["sample_original", "subtype"]

    df["sample_id"] = df["sample_original"].apply(core_id)
    df["subtype"] = df["subtype"].astype(str).str.strip()
    df["cohort"] = cohort

    df = df[df["subtype"].isin(SUBTYPE_ORDER)].copy()

    conflicts = df.groupby("sample_id")["subtype"].nunique()

    if (conflicts > 1).any():
        bad_samples = conflicts[conflicts > 1].index.tolist()
        raise ValueError(
            f"Conflicting subtype assignments in {cohort}:\n"
            + "\n".join(bad_samples))

    df = df.drop_duplicates("sample_id", keep="first")
    return df[["sample_id", "subtype", "cohort"]]


def read_activities(path, cohort, modality):
    df = pd.read_csv(path, sep="\t")
    sample_col = df.columns[0]

    df = df.rename(columns={sample_col: "sample_original"})
    df["sample_id"] = df["sample_original"].apply(core_id)

    signature_cols = [
        col for col in df.columns
        if re.fullmatch(rf"{modality}\d+", str(col), flags=re.IGNORECASE)]

    if not signature_cols:
        raise ValueError(
            f"No {modality} signature columns detected in {path}.\n"
            f"Available columns: {list(df.columns)}")

    print(f"\n{cohort} {modality}")
    print(f"Samples: {len(df)}")
    print(f"Signatures: {len(signature_cols)}")
    print(f"Detected signatures: {signature_cols}")

    duplicates = df[df["sample_id"].duplicated(keep=False)]

    if not duplicates.empty:
        print("\nTechnical duplicates detected:")
        for sample_id, group in duplicates.groupby("sample_id"):
            print(f"  {sample_id}: {len(group)} rows -> keeping first")

        df = df.drop_duplicates("sample_id", keep="first").copy()

    long_df = df.melt(
        id_vars=["sample_id"],
        value_vars=signature_cols,
        var_name="signature",
        value_name="activity")

    long_df["activity"] = pd.to_numeric(long_df["activity"], errors="coerce")
    long_df["cohort"] = cohort
    long_df["modality"] = modality

    return long_df

# Load subtype data
subtypes = pd.concat(
    [read_subtypes(path, cohort) for cohort, path in subtype_files.items()],
    ignore_index=True)

subtype_counts = (
    subtypes.groupby(["cohort", "subtype"])
    .size()
    .unstack(fill_value=0)
    .reindex(columns=SUBTYPE_ORDER, fill_value=0)
)

print(subtype_counts)
print(f"\nTotal samples with defined subtype:\n{len(subtypes)}")

# Load signature activities
activity_tables = []

for modality, files in [("SBS", sbs_files), ("DBS", dbs_files), ("ID", id_files)]:
    for cohort, path in files.items():
        activity_tables.append(read_activities(path, cohort, modality))

activities = pd.concat(activity_tables, ignore_index=True)


# Match activities with subtype assignments
data = activities.merge(subtypes, on=["sample_id", "cohort"], how="inner")

matched_counts = (
    data[["sample_id", "cohort", "modality", "subtype"]]
    .drop_duplicates()
    .groupby(["modality", "cohort", "subtype"])
    .size()
    .unstack(fill_value=0)
    .reindex(columns=SUBTYPE_ORDER, fill_value=0)
)

data = data.dropna(subset=["activity"]).copy()

if POSITIVE_ONLY:
    data = data[data["activity"] > 0].copy()

# Filter sparse signatures
signature_counts = (
    data.groupby(["modality", "signature", "subtype"])
    .size()
    .reset_index(name="n_positive"))

signature_summary = (
    signature_counts.groupby(["modality", "signature"])
    .agg(
        total_positive=("n_positive", "sum"),
        n_subtypes_sufficient=("n_positive", lambda x: (x >= MIN_POSITIVE_PER_SUBTYPE).sum())
    )
    .reset_index())

signature_summary["keep"] = (
    (signature_summary["total_positive"] >= MIN_TOTAL_POSITIVE)
    & (signature_summary["n_subtypes_sufficient"] >= MIN_SUBTYPES))


keep_signatures = signature_summary.loc[signature_summary["keep"], "signature"].tolist()
dropped_signatures = signature_summary.loc[~signature_summary["keep"], "signature"].tolist()

print("\nSignatures retained:")
print(keep_signatures)

print("\nSparse signatures omitted:")
print(dropped_signatures)

data_plot = data[data["signature"].isin(keep_signatures)].copy()

# Remove subtype-signature groups with fewer than two positive samples
group_counts = (
    data_plot.groupby(["modality", "signature", "subtype"])
    .size()
    .reset_index(name="n_positive"))

valid_groups = group_counts[
    group_counts["n_positive"] >= MIN_POSITIVE_FOR_BOX][["modality", "signature", "subtype"]]

data_plot = data_plot.merge(
    valid_groups,
    on=["modality", "signature", "subtype"],
    how="inner")

print("\nSubtype-signature groups retained for boxplots:")
print(
    group_counts[group_counts["n_positive"] >= MIN_POSITIVE_FOR_BOX]
    .to_string(index=False))

print(f"\nSubtype-signature groups omitted because n < {MIN_POSITIVE_FOR_BOX}:")
print(
    group_counts[group_counts["n_positive"] < MIN_POSITIVE_FOR_BOX]
    .to_string(index=False))

signature_summary.to_csv(
    outdir / "signature_visualization_filter.tsv",
    sep="\t",
    index=False)

plot_counts = (
    data_plot.groupby(["modality", "signature", "subtype"])
    .size()
    .reset_index(name="n_positive"))

plot_counts.to_csv(
    outdir / "signature_activity_boxplot_counts.tsv",
    sep="\t",
    index=False)


# Kruskal-Wallis tests
global_results = []

for (modality, signature), sig_data in data_plot.groupby(["modality", "signature"]):
    groups = []
    tested_subtypes = []
    group_sizes = []

    for subtype in SUBTYPE_ORDER:
        values = sig_data.loc[sig_data["subtype"] == subtype, "activity"].dropna().to_numpy()

        if len(values) >= MIN_N_FOR_TEST:
            groups.append(values)
            tested_subtypes.append(subtype)
            group_sizes.append(f"{subtype}:{len(values)}")

    if len(groups) < 2:
        continue

    try:
        statistic, p_value = kruskal(*groups)
    except ValueError as error:
        if "All numbers are identical" in str(error):
            statistic, p_value = 0.0, 1.0
        else:
            raise

    global_results.append({
        "modality": modality,
        "signature": signature,
        "n_subtypes_tested": len(groups),
        "subtypes_tested": ", ".join(tested_subtypes),
        "group_sizes": ", ".join(group_sizes),
        "n_total": sum(len(group) for group in groups),
        "H_statistic": statistic,
        "p_value": p_value,
        "nominal_significant": p_value < ALPHA})

global_results = pd.DataFrame(global_results)

if not global_results.empty:
    global_results["BH_adjusted_p"] = multipletests(
        global_results["p_value"],
        alpha=ALPHA,
        method="fdr_bh")[1]

    global_results["significant_BH"] = global_results["BH_adjusted_p"] < ALPHA
    global_results = global_results.sort_values(["p_value", "BH_adjusted_p"]).reset_index(drop=True)

print("KRUSKAL-WALLIS TESTS")

if global_results.empty:
    print(f"No signature had at least two subtype groups with n >= {MIN_N_FOR_TEST}.")
else:
    print(global_results.to_string(index=False))

global_results.to_csv(
    outdir / "signature_activity_kruskal_wallis.tsv",
    sep="\t",
    index=False)


# Pairwise Mann-Whitney U tests
pairwise_results = []

if not global_results.empty:
    nominal_global = global_results[global_results["p_value"] < ALPHA]
    significant_keys = set(zip(nominal_global["modality"], nominal_global["signature"]))

    for (modality, signature), sig_data in data_plot.groupby(["modality", "signature"]):
        if (modality, signature) not in significant_keys:
            continue

        valid_groups = {}

        for subtype in SUBTYPE_ORDER:
            values = sig_data.loc[sig_data["subtype"] == subtype, "activity"].dropna().to_numpy()

            if len(values) >= MIN_N_FOR_TEST:
                valid_groups[subtype] = values

        for group1, group2 in combinations(valid_groups, 2):
            values1 = valid_groups[group1]
            values2 = valid_groups[group2]

            statistic, p_value = mannwhitneyu(
                values1,
                values2,
                alternative="two-sided",
                method="auto")

            pairwise_results.append({
                "modality": modality,
                "signature": signature,
                "group1": group1,
                "group2": group2,
                "n1": len(values1),
                "n2": len(values2),
                "median1": np.median(values1),
                "median2": np.median(values2),
                "U_statistic": statistic,
                "p_value": p_value,
                "nominal_significant": p_value < ALPHA})

pairwise_results = pd.DataFrame(pairwise_results)

if not pairwise_results.empty:
    pairwise_results["BH_adjusted_p"] = np.nan
    pairwise_results["significant_BH"] = False

    for _, indices in pairwise_results.groupby(["modality", "signature"]).groups.items():
        indices = list(indices)
        adjusted_p = multipletests(
            pairwise_results.loc[indices, "p_value"],
            alpha=ALPHA,
            method="fdr_bh")[1]

        pairwise_results.loc[indices, "BH_adjusted_p"] = adjusted_p
        pairwise_results.loc[indices, "significant_BH"] = adjusted_p < ALPHA

    pairwise_results = pairwise_results.sort_values(["modality", "signature", "p_value"]).reset_index(drop=True)

print("PAIRWISE MANN-WHITNEY U TESTS")

if pairwise_results.empty:
    print("No pairwise tests were performed because no signature had a global nominal p < 0.05.")
else:
    print(pairwise_results.to_string(index=False))

pairwise_results.to_csv(
    outdir / "signature_activity_pairwise_mann_whitney.tsv",
    sep="\t",
    index=False)


# Signature order

def signature_number(signature):
    match = re.search(r"\d+", signature)
    return int(match.group()) if match else 9999


sbs_order = sorted(
    data_plot.loc[data_plot["modality"] == "SBS", "signature"].unique(),
    key=signature_number)

dbs_order = sorted(
    data_plot.loc[data_plot["modality"] == "DBS", "signature"].unique(),
    key=signature_number)

id_order = sorted(
    data_plot.loc[data_plot["modality"] == "ID", "signature"].unique(),
    key=signature_number)

dbs_id_order = dbs_order + id_order

data_plot["subtype"] = pd.Categorical(
    data_plot["subtype"],
    categories=SUBTYPE_ORDER,
    ordered=True)


palette = {
    "BL1": "#64A8ED",
    "BL2": "#E17C05",
    "IM": "#65A165",
    "M": "#F5112E",
    "MSL": "#DA8BC3",
    "LAR": "#F5F187",
    "UNS": "#B0B0B0",
    "Excluded": "#F2F2F2"
}


# Common y-axis

if data_plot.empty:
    raise ValueError("No signatures remained after filtering. Reduce filtering thresholds.")

if LOG_SCALE:
    positive_values = data_plot.loc[data_plot["activity"] > 0, "activity"]
    global_min = positive_values.min()
    global_max = positive_values.max()
    ymin = Y_MIN if Y_MIN is not None else global_min * 0.85
    ymax = Y_MAX if Y_MAX is not None else global_max * 1.35
else:
    global_max = data_plot["activity"].max()
    ymin = Y_MIN if Y_MIN is not None else 0
    ymax = Y_MAX if Y_MAX is not None else global_max * 1.08

def format_p(value):
    return "<0.001" if value < 0.001 else f"={value:.3f}"

def draw_compact_boxplots(ax, panel_data, signature_order, palette,
                          global_results, pairwise_results, log_scale=True):

    current_x = 0
    box_width = 0.72
    within_gap = 0.05
    between_signature_gap = 0.55

    signature_centers = []
    signature_labels = []
    signature_ranges = []
    signature_stats = []
    all_positions = []
    max_annotation_y = 0

    for signature in signature_order:
        sig_data = panel_data[panel_data["signature"] == signature]

        if sig_data.empty:
            continue

        modality = sig_data["modality"].iloc[0]
        groups = []

        for subtype in SUBTYPE_ORDER:
            values = sig_data.loc[sig_data["subtype"] == subtype, "activity"].dropna().to_numpy()

            if len(values) > 0:
                groups.append((subtype, values))

        if not groups:
            continue

        positions = []
        subtype_positions = {}
        subtype_values = {}

        for subtype, values in groups:
            pos = current_x

            positions.append(pos)
            all_positions.append(pos)
            subtype_positions[subtype] = pos
            subtype_values[subtype] = values

            bp = ax.boxplot(
                values,
                positions=[pos],
                widths=box_width,
                patch_artist=True,
                showfliers=False,
                manage_ticks=False,
                whis=1.5)

            for box in bp["boxes"]:
                box.set_facecolor(palette[subtype])
                box.set_edgecolor("0.35")
                box.set_linewidth(1)

            for median in bp["medians"]:
                median.set_color("0.25")
                median.set_linewidth(1.1)

            for whisker in bp["whiskers"]:
                whisker.set_color("0.4")
                whisker.set_linewidth(0.9)

            for cap in bp["caps"]:
                cap.set_color("0.4")
                cap.set_linewidth(0.9)

            upper_cap_y = max(cap.get_ydata()[0] for cap in bp["caps"])

            ax.annotate(
                f"n={len(values)}",
                xy=(pos, upper_cap_y),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=7,
                color="0.15",
                clip_on=False)

            current_x += box_width + within_gap

        signature_center = (positions[0] + positions[-1]) / 2
        signature_centers.append(signature_center)

        global_row = global_results[
            (global_results["modality"] == modality)
            & (global_results["signature"] == signature)]

        label = signature

        signature_labels.append(signature)

        if not global_row.empty:
            global_p = global_row.iloc[0]["p_value"]
            global_adjusted_p = global_row.iloc[0]["BH_adjusted_p"]

            if global_p < ALPHA:
                stats_text = f"KW p{format_p(global_p)}, BH-adj. p{format_p(global_adjusted_p)}"
                signature_stats.append((signature_center, stats_text))

        signature_ranges.append((positions[0], positions[-1]))

        if not pairwise_results.empty:
            sig_pairs = pairwise_results[
                (pairwise_results["modality"] == modality)
                & (pairwise_results["signature"] == signature)
                & (pairwise_results["p_value"] < ALPHA)
            ].copy()

            usable_pairs = []

            for _, row in sig_pairs.iterrows():
                group1 = row["group1"]
                group2 = row["group2"]

                if group1 not in subtype_positions or group2 not in subtype_positions:
                    continue

                x1 = subtype_positions[group1]
                x2 = subtype_positions[group2]

                usable_pairs.append({
                    "x1": min(x1, x2),
                    "x2": max(x1, x2),
                    "span": abs(x2 - x1),
                    "p_value": row["p_value"],
                    "BH_adjusted_p": row["BH_adjusted_p"]
                })

            usable_pairs.sort(key=lambda pair: (pair["span"], pair["p_value"]))

            levels = []

            for pair in usable_pairs:
                for level, occupied in enumerate(levels):
                    overlaps = any(
                        not (pair["x2"] < old_x1 or pair["x1"] > old_x2)
                        for old_x1, old_x2 in occupied)

                    if not overlaps:
                        pair["level"] = level
                        occupied.append((pair["x1"], pair["x2"]))
                        break
                else:
                    pair["level"] = len(levels)
                    levels.append([(pair["x1"], pair["x2"])])

            if usable_pairs:
                signature_max = max(np.max(values) for values in subtype_values.values())

                for pair in usable_pairs:
                    level = pair["level"]
                    x1 = pair["x1"]
                    x2 = pair["x2"]

                    if log_scale:
                        bracket_y = signature_max * 1.85 * (3.5 ** level)
                        bracket_bottom = bracket_y / 1.12
                        text_y = bracket_y * 1.08
                    else:
                        data_range = panel_data["activity"].max() - panel_data["activity"].min()
                        bracket_y = signature_max + 0.12 * data_range + level * 0.12 * data_range
                        bracket_bottom = bracket_y - 0.025 * data_range
                        text_y = bracket_y + 0.012 * data_range

                    ax.plot(
                        [x1, x1, x2, x2],
                        [bracket_bottom, bracket_y, bracket_y, bracket_bottom],
                        color="0.25",
                        linewidth=0.7,
                        clip_on=False)

                    ax.text(
                        (x1 + x2) / 2,
                        text_y,
                        f"p{format_p(pair['p_value'])}\n"
                        f"BH-adj. p{format_p(pair['BH_adjusted_p'])}",
                        ha="center",
                        va="bottom",
                        fontsize=8,
                        linespacing=1.2,
                        color="0.20",
                        clip_on=False)

                    max_annotation_y = max(max_annotation_y, text_y)

        current_x += between_signature_gap

    for i in range(len(signature_ranges) - 1):
        last_box_current = signature_ranges[i][1]
        first_box_next = signature_ranges[i + 1][0]
        separator_x = (last_box_current + first_box_next) / 2

        ax.axvline(separator_x, color="0.75", linewidth=0.5, alpha=0.80, zorder=0)

    ax.set_xticks(signature_centers)
    ax.set_xticklabels(signature_labels)

    for x, text in signature_stats:
        ax.annotate(
            text,
            xy=(x, 0),
            xycoords=("data", "axes fraction"),
            xytext=(0, -25),
            textcoords="offset points",
            ha="center",
            va="top",
            fontsize=8,
            color="0.25",
            clip_on=False)

    if all_positions:
        left = min(all_positions) - box_width / 2 - 0.08
        right = max(all_positions) + box_width / 2 + 0.08
        ax.set_xlim(left, right)

    return max_annotation_y


# Plot

mid = (len(sbs_order) + 1) // 2

sbs_order_1 = sbs_order[:mid]
sbs_order_2 = sbs_order[mid:]

sbs_data = data_plot[data_plot["modality"] == "SBS"].copy()
dbs_id_data = data_plot[data_plot["modality"].isin(["DBS", "ID"])].copy()

sns.set_theme(style="white")

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "legend.title_fontsize": 9})

fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharey=True)

plot_info = [
    ("SBS", sbs_data, sbs_order_1, axes[0]),
    ("SBS", sbs_data, sbs_order_2, axes[1]),
    ("DBS + ID", dbs_id_data, dbs_id_order, axes[2])]

annotation_maxima = []

for panel_title, panel_data, signature_order, ax in plot_info:
    annotation_max = draw_compact_boxplots(
        ax,
        panel_data,
        signature_order,
        palette,
        global_results,
        pairwise_results,
        LOG_SCALE)

    annotation_maxima.append(annotation_max)

    if LOG_SCALE:
        ax.set_yscale("log")

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(panel_title, loc="left", fontsize=11, fontweight="bold", pad=5)
    ax.grid(axis="y", linewidth=0.55, alpha=0.20)
    ax.grid(axis="x", visible=False)
    ax.tick_params(axis="x", labelsize=11, pad=4)
    ax.tick_params(axis="y", labelsize=11)
    sns.despine(ax=ax)


annotation_max = max(annotation_maxima)

if annotation_max > 0:
    final_ymax = max(ymax, annotation_max * 1.20 if LOG_SCALE else annotation_max * 1.05)
else:
    final_ymax = ymax

for ax in axes:
    ax.set_ylim(ymin, final_ymax)


fig.supylabel("Signature activity (mutations)", x=0.001, fontsize=12, fontweight='bold')

fig.suptitle(
    "Mutational signature activity across TNBC molecular subtypes",
    fontsize=12,
    y=0.995,
    x=0.55)

legend_handles = [
    Patch(
        facecolor=palette[subtype],
        edgecolor="0.35",
        linewidth=0.8,
        label=subtype)
    for subtype in SUBTYPE_ORDER]

legend = fig.legend(
    handles=legend_handles,
    #title="TNBC subtype",
    loc="upper center",
    bbox_to_anchor=(0.55, 0.975),
    ncol=6,
    frameon=True,
    columnspacing=1.0,
    handletextpad=0.4,
    fontsize=11)

legend.get_frame().set_edgecolor("0.35")
legend.get_frame().set_linewidth(0.8)

plt.subplots_adjust(top=0.89, bottom=0.01, left=0.10, right=0.98, hspace=0.32)

fig.savefig(
    outdir / "signature_activity_by_subtype_statistics_A4.pdf",
    bbox_inches="tight")

fig.savefig(
    outdir / "signature_activity_by_subtype_statistics_A4.png",
    dpi=600,
    bbox_inches="tight")

#plt.show()
