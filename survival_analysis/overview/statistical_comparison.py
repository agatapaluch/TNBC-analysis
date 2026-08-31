from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests


FAVORABLE_FILES = {
    "SBS": Path("SYNERGY_outcome_groups/SBS_favorable.txt"),
    "DBS": Path("SYNERGY_outcome_groups/DBS_favorable.txt"),
    "ID": Path("SYNERGY_outcome_groups/ID_favorable.txt")}

UNFAVORABLE_FILES = {
    "SBS": Path("SYNERGY_outcome_groups/SBS_unfavorable.txt"),
    "DBS": Path("SYNERGY_outcome_groups/DBS_unfavorable.txt"),
    "ID": Path("SYNERGY_outcome_groups/ID_unfavorable.txt")}

OUTPUT_DIR = Path("SYNERGY_outcome_groups/comparison_plots")
OUTPUT_DIR.mkdir(exist_ok=True)

SIGNATURE_ETIOLOGY = {
    "SBS1": "Clock-like", "SBS5": "Clock-like",

    "SBS4": "Tobacco", "SBS29": "Tobacco", "SBS92": "Tobacco",
    "SBS100": "Tobacco", "DBS2": "Tobacco/acetaldehyde",

    "SBS7a": "UV", "SBS7b": "UV", "SBS7c": "UV", "SBS7d": "UV",
    "SBS38": "UV", "DBS1": "UV", "ID13": "UV",

    "SBS3": "HRD",
    "SBS6": "dMMR", "SBS14": "dMMR", "SBS15": "dMMR",
    "SBS20": "dMMR", "SBS21": "dMMR", "SBS26": "dMMR", "SBS44": "dMMR",
    "SBS30": "dBER", "SBS36": "dBER",

    "DBS13": "HRD", "ID6": "HRD",
    "DBS10": "dMMR", "ID7": "dMMR",
    "ID8": "NHEJ",

    "SBS9": "Polymerase mutation",
    "SBS10a": "Polymerase mutation",
    "SBS10b": "Polymerase mutation",
    "SBS10c": "Polymerase mutation",
    "SBS10d": "Polymerase mutation",
    "DBS3": "Polymerase mutation",

    "ID1": "Replication slippage",
    "ID2": "Replication slippage",

    "SBS2": "APOBEC activity",
    "SBS13": "APOBEC activity",
    "SBS84": "AID activity",
    "SBS85": "AID activity",

    "SBS11": "Cancer treatment",
    "SBS25": "Cancer treatment",
    "SBS31": "Cancer treatment",
    "SBS32": "Cancer treatment",
    "SBS35": "Cancer treatment",
    "SBS86": "Cancer treatment",
    "SBS87": "Cancer treatment",
    "SBS90": "Cancer treatment",
    "SBS99": "Cancer treatment",
    "DBS5": "Cancer treatment",

    "SBS22a": "Aristolochic acid",
    "SBS22b": "Aristolochic acid",
    "DBS20": "Aristolochic acid",
    "ID23": "Aristolochic acid",

    "SBS24": "Aflatoxin",

    "SBS88": "Colibactin",
    "ID18": "Colibactin",

    "SBS42": "Haloalkane exposure",
    "SBS18": "ROS",

    "ID4": "TOP mutation",
    "ID17": "TOP mutation",

    "SBS27": "Sequencing artefacts",
    "SBS43": "Sequencing artefacts",
    "SBS45": "Sequencing artefacts",
    "SBS46": "Sequencing artefacts",
    "SBS47": "Sequencing artefacts",
    "SBS48": "Sequencing artefacts",
    "SBS49": "Sequencing artefacts",
    "SBS50": "Sequencing artefacts",
    "SBS51": "Sequencing artefacts",
    "SBS52": "Sequencing artefacts",
    "SBS53": "Sequencing artefacts",
    "SBS54": "Sequencing artefacts",
    "SBS55": "Sequencing artefacts",
    "SBS56": "Sequencing artefacts",
    "SBS57": "Sequencing artefacts",
    "SBS58": "Sequencing artefacts",
    "SBS59": "Sequencing artefacts",
    "SBS60": "Sequencing artefacts",
    "SBS95": "Sequencing artefacts",

    "SBS8": "Unknown", "SBS12": "Unknown", "SBS16": "Unknown",
    "SBS17a": "Unknown", "SBS17b": "Unknown", "SBS19": "Unknown",
    "SBS23": "Unknown", "SBS28": "Unknown", "SBS33": "Unknown",
    "SBS34": "Unknown", "SBS37": "Unknown", "SBS39": "Unknown",
    "SBS40a": "Unknown", "SBS40b": "Unknown", "SBS40c": "Unknown",
    "SBS41": "Unknown", "SBS89": "Unknown", "SBS91": "Unknown",
    "SBS93": "Unknown", "SBS94": "Unknown", "SBS96": "Unknown",
    "SBS97": "Unknown", "SBS98": "Unknown",

    "DBS4": "Unknown", "DBS6": "Unknown", "DBS8": "Unknown",
    "DBS9": "Unknown", "DBS11": "Unknown", "DBS12": "Unknown",
    "DBS14": "Unknown", "DBS15": "Unknown", "DBS16": "Unknown",
    "DBS17": "Unknown", "DBS18": "Unknown", "DBS19": "Unknown",

    "ID5": "Unknown", "ID9": "Unknown", "ID10": "Unknown",
    "ID11": "Unknown", "ID12": "Unknown", "ID14": "Unknown",
    "ID15": "Unknown", "ID16": "Unknown", "ID19": "Unknown",
    "ID20": "Unknown", "ID21": "Unknown", "ID22": "Unknown"
}


FAVORABLE_COLOR = "#4C72B0"
UNFAVORABLE_COLOR = "#C44E52"

TITLE_SIZE = 17
AXIS_SIZE = 15
TICK_SIZE = 14
LEGEND_SIZE = 13
VALUE_SIZE = 12


def normalize_sample_id(sample_id):
    sample_id = str(sample_id).strip()
    sample_id = sample_id.split("_vs_", maxsplit=1)[0]
    return pd.Series([sample_id]).str.replace(r"-DNA-\d+$", "", regex=True).iloc[0]


def load_group(files):
    data_by_type = {}
    sample_sizes = {}

    for modality, path in files.items():
        if not path.exists():
            raise FileNotFoundError(f"{modality} file not found: {path}")

        data = pd.read_csv(path, sep="\t", index_col=0)
        data = data.apply(pd.to_numeric, errors="coerce").fillna(0)

        data.index = [normalize_sample_id(sample) for sample in data.index]
        data = data[~data.index.duplicated(keep="first")]

        columns = [col for col in data.columns if col.startswith(modality)]

        if not columns:
            raise ValueError(f"No {modality} signature columns found in {path}")

        data_by_type[modality] = data[columns]
        sample_sizes[modality] = len(data)

    return data_by_type, sample_sizes


def get_etiologies(data_by_type):
    etiologies = set()

    for data in data_by_type.values():
        for signature in data.columns:
            etiologies.add(SIGNATURE_ETIOLOGY.get(signature, "Unknown"))

    return etiologies


def get_etiology_columns(data_by_type, etiology):
    columns = {}

    for modality, data in data_by_type.items():
        columns[modality] = [
            signature for signature in data.columns
            if SIGNATURE_ETIOLOGY.get(signature, "Unknown") == etiology
        ]

    return {modality: sigs for modality, sigs in columns.items() if sigs}


def get_prevalence_status(data_by_type, etiology):
    columns = get_etiology_columns(data_by_type, etiology)

    if not columns:
        return pd.Series(dtype=float)

    samples = data_by_type["SBS"].index
    status = pd.Series(np.nan, index=samples, dtype=float)

    for sample in samples:
        positive = False
        all_available = True

        for modality, signatures in columns.items():
            data = data_by_type[modality]

            if sample not in data.index:
                all_available = False
                continue

            if (data.loc[sample, signatures] > 0).any():
                positive = True

        if positive:
            status.loc[sample] = 1
        elif all_available:
            status.loc[sample] = 0

    return status


def get_complete_activity(data_by_type, etiology):
    columns = get_etiology_columns(data_by_type, etiology)

    if not columns:
        return pd.Series(dtype=float)

    samples = data_by_type["SBS"].index
    activity = pd.Series(np.nan, index=samples, dtype=float)

    for sample in samples:
        total = 0
        complete = True

        for modality, signatures in columns.items():
            data = data_by_type[modality]

            if sample not in data.index:
                complete = False
                break

            total += data.loc[sample, signatures].sum()

        if complete:
            activity.loc[sample] = total

    return activity


def calculate_activity_share(data_by_type):
    totals = {}

    for data in data_by_type.values():
        for signature in data.columns:
            total = data[signature].sum()

            if total <= 0:
                continue

            etiology = SIGNATURE_ETIOLOGY.get(signature, "Unknown")
            totals[etiology] = totals.get(etiology, 0) + float(total)

    grand_total = sum(totals.values())

    return {
        etiology: value / grand_total * 100
        for etiology, value in totals.items()
    }


def calculate_prevalence_percent(data_by_type, etiology):
    status = get_prevalence_status(data_by_type, etiology).dropna()

    if len(status) == 0:
        return np.nan

    return (status == 1).mean() * 100


def adjust_fdr(p_values):
    p_values = np.asarray(p_values, dtype=float)
    adjusted = np.full(len(p_values), np.nan)

    valid = np.isfinite(p_values)

    if valid.any():
        adjusted[valid] = multipletests(
            p_values[valid],
            method="fdr_bh"
        )[1]

    return adjusted


def significance_stars(q):
    if pd.isna(q):
        return ""
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return ""


def prepare_results(favorable, unfavorable):
    etiologies = sorted(
        get_etiologies(favorable) |
        get_etiologies(unfavorable)
    )

    favorable_activity_share = calculate_activity_share(favorable)
    unfavorable_activity_share = calculate_activity_share(unfavorable)

    rows = []

    for etiology in etiologies:
        fav_status = get_prevalence_status(favorable, etiology).dropna()
        unfav_status = get_prevalence_status(unfavorable, etiology).dropna()

        fav_positive = int((fav_status == 1).sum())
        fav_negative = int((fav_status == 0).sum())
        unfav_positive = int((unfav_status == 1).sum())
        unfav_negative = int((unfav_status == 0).sum())

        if len(fav_status) > 0 and len(unfav_status) > 0:
            _, prevalence_p = fisher_exact([
                [fav_positive, fav_negative],
                [unfav_positive, unfav_negative]
            ])
        else:
            prevalence_p = np.nan

        fav_activity = get_complete_activity(favorable, etiology).dropna()
        unfav_activity = get_complete_activity(unfavorable, etiology).dropna()

        if len(fav_activity) > 0 and len(unfav_activity) > 0:
            if (
                fav_activity.nunique() == 1 and
                unfav_activity.nunique() == 1 and
                fav_activity.iloc[0] == unfav_activity.iloc[0]
            ):
                activity_p = 1.0
            else:
                activity_p = mannwhitneyu(
                    fav_activity,
                    unfav_activity,
                    alternative="two-sided"
                ).pvalue
        else:
            activity_p = np.nan

        fav_prevalence = (
            fav_positive / len(fav_status) * 100
            if len(fav_status) else np.nan
        )

        unfav_prevalence = (
            unfav_positive / len(unfav_status) * 100
            if len(unfav_status) else np.nan
        )

        fav_share = favorable_activity_share.get(etiology, 0)
        unfav_share = unfavorable_activity_share.get(etiology, 0)

        rows.append({
            "Etiology": etiology,

            "Favorable prevalence": fav_prevalence,
            "Unfavorable prevalence": unfav_prevalence,
            "Prevalence difference": unfav_prevalence - fav_prevalence,

            "Favorable prevalence n": len(fav_status),
            "Unfavorable prevalence n": len(unfav_status),
            "prevalence_p": prevalence_p,

            "Favorable activity": fav_share,
            "Unfavorable activity": unfav_share,
            "Activity difference": unfav_share - fav_share,

            "Favorable activity n": len(fav_activity),
            "Unfavorable activity n": len(unfav_activity),
            "activity_p": activity_p
        })

    results = pd.DataFrame(rows).set_index("Etiology")

    results["prevalence_FDR"] = adjust_fdr(results["prevalence_p"])
    results["activity_FDR"] = adjust_fdr(results["activity_p"])

    results["max_difference"] = results[
        ["Prevalence difference", "Activity difference"]
    ].abs().max(axis=1)

    results = results.sort_values("max_difference", ascending=False)
    results = results.drop(columns="max_difference")

    return results


def create_difference_plot(results):
    fig, axes = plt.subplots(
        1, 2,
        figsize=(15, 8),
        sharey=True
    )

    y = np.arange(len(results))

    panels = [
        ("Prevalence difference", "prevalence_FDR", "Etiology Prevalence"),
        ("Activity difference", "activity_FDR", "Etiology Activity")
    ]

    for ax, (column, fdr_column, title) in zip(axes, panels):
        values = results[column].to_numpy(dtype=float)
        fdr_values = results[fdr_column].to_numpy(dtype=float)

        colors = [
            UNFAVORABLE_COLOR if value > 0 else FAVORABLE_COLOR
            for value in values]

        bars = ax.barh(
            y,
            values,
            color=colors,
            edgecolor="black",
            linewidth=0.4,
            zorder=3)

        ax.axvline(0, color="black", linewidth=1)

        finite_values = values[np.isfinite(values)]
        max_abs = max(abs(finite_values.min()), abs(finite_values.max())) if len(finite_values) else 1
        max_abs = max(max_abs, 1)

        ax.set_xlim(-max_abs * 1.35, max_abs * 1.35)

        for bar, value, q in zip(bars, values, fdr_values):
            if not np.isfinite(value) or abs(value) < 0.05:
                continue

            stars = significance_stars(q)
            label = f"{value:+.1f}"

            if stars:
                label += f" {stars}"

            offset = max_abs * 0.025

            if value > 0:
                x = value + offset
                ha = "left"
            else:
                x = value - offset
                ha = "right"

            ax.text(
                x,
                bar.get_y() + bar.get_height() / 2,
                label,
                ha=ha,
                va="center",
                fontsize=VALUE_SIZE
            )

        ax.set_title(
            title,
            fontsize=TITLE_SIZE,
            fontweight="bold",
            pad=14)

        ax.tick_params(
            axis="both",
            labelsize=TICK_SIZE,
            width=1.2,
            length=5)

        ax.grid(axis="x", alpha=0.2, zorder=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0].set_yticks(y)
    axes[0].set_yticklabels(results.index, fontsize=TICK_SIZE)
    axes[0].invert_yaxis()

    handles = [
        Line2D(
            [0], [0],
            marker="s",
            linestyle="none",
            markerfacecolor=FAVORABLE_COLOR,
            markeredgecolor="none",
            markersize=10,
            label="Higher in favorable (OS > 5 years)"
        ),
        Line2D(
            [0], [0],
            marker="s",
            linestyle="none",
            markerfacecolor=UNFAVORABLE_COLOR,
            markeredgecolor="none",
            markersize=10,
            label="Higher in unfavorable (death ≤ 2 years)")]

    legend = fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.55, 1.0),
        ncol=2,
        frameon=True,
        fancybox=False,
        framealpha=1,
        edgecolor="black",
        fontsize=LEGEND_SIZE)

    legend.get_frame().set_linewidth(0.8)

    fig.supxlabel(
        "Difference: Unfavorable − Favorable (percentage points)",
        fontsize=AXIS_SIZE,
        fontweight="bold",
        x=0.55,
        y=0.0)

    all_non_significant = (
            (results["prevalence_FDR"].dropna() >= 0.05).all() and
            (results["activity_FDR"].dropna() >= 0.05).all())

    if all_non_significant:
        fig.text(
            0.25,
            0.002,
            "All comparisons: FDR ≥ 0.05",
            ha="right",
            va="bottom",
            fontsize=14)

    plt.subplots_adjust(
        left=0.20,
        right=0.98,
        bottom=0.12,
        top=0.84,
        wspace=0.15)

    output_file = OUTPUT_DIR / "etiology_difference_comparison.png"

    plt.savefig(output_file, dpi=300,bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {output_file}")


def create_hrd_activity_boxplot(favorable, unfavorable):
    fav = get_complete_activity(favorable, "HRD").dropna()
    unfav = get_complete_activity(unfavorable, "HRD").dropna()

    p = mannwhitneyu(fav, unfav, alternative="two-sided").pvalue

    print("\nHRD patient-level activity:")
    print(f"Favorable: n={len(fav)}, median={fav.median():.1f}, mean={fav.mean():.1f}, max={fav.max():.1f}")
    print(f"Unfavorable: n={len(unfav)}, median={unfav.median():.1f}, mean={unfav.mean():.1f}, max={unfav.max():.1f}")
    print(f"Mann-Whitney U p = {p:.6f}")

    fig, ax = plt.subplots(figsize=(7, 10))

    bp = ax.boxplot(
        [fav, unfav],
        positions=[1, 2],
        widths=0.5,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.8},
        boxprops={"edgecolor": "black", "linewidth": 1.0},
        whiskerprops={"color": "black", "linewidth": 1.0},
        capprops={"color": "black", "linewidth": 1.0}
    )

    bp["boxes"][0].set_facecolor(FAVORABLE_COLOR)
    bp["boxes"][1].set_facecolor(UNFAVORABLE_COLOR)
    bp["boxes"][0].set_alpha(0.45)
    bp["boxes"][1].set_alpha(0.45)

    rng = np.random.default_rng(42)

    ax.scatter(
        rng.normal(1, 0.06, len(fav)), fav,
        s=65, color=FAVORABLE_COLOR, edgecolor="black",
        linewidth=0.4, alpha=0.8, zorder=3
    )

    ax.scatter(
        rng.normal(2, 0.06, len(unfav)), unfav,
        s=65, color=UNFAVORABLE_COLOR, edgecolor="black",
        linewidth=0.4, alpha=0.8, zorder=3
    )

    ax.set_xticks([1, 2])
    ax.set_xticklabels([
        f"Favorable\n(n = {len(fav)})",
        f"Unfavorable\n(n = {len(unfav)})"
    ], fontsize=18)


    ax.set_ylabel("Patient-level HRD activity", fontsize=21, fontweight="bold", labelpad=12)
    ax.set_title("HRD Activity by clinical outcome", fontsize=21, fontweight="bold", pad=14)

    ax.text(
        1.5, 0.98,
        f"Mann–Whitney U p = {p:.3f}",
        transform=ax.get_xaxis_transform(),
        ha="center", va="top", fontsize=20
    )

    ax.tick_params(axis="both", labelsize=16, width=1.2, length=5)
    ax.grid(axis="y", alpha=0.2, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    output_file = OUTPUT_DIR / "HRD_patient_level_activity_boxplot.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    #plt.show()
    plt.close(fig)

    print(f"Saved: {output_file}")


favorable, favorable_sizes = load_group(FAVORABLE_FILES)
unfavorable, unfavorable_sizes = load_group(UNFAVORABLE_FILES)

print("\nFavorable:")
for modality, n in favorable_sizes.items():
    print(f"{modality}: n = {n}")

print("\nUnfavorable:")
for modality, n in unfavorable_sizes.items():
    print(f"{modality}: n = {n}")

fav_hrd = get_complete_activity(favorable, "HRD").dropna()
unfav_hrd = get_complete_activity(unfavorable, "HRD").dropna()

# Original test
p_original = mannwhitneyu(fav_hrd, unfav_hrd, alternative="two-sided").pvalue

# Remove the single highest HRD activity from the unfavorable group
outlier_sample = unfav_hrd.idxmax()
outlier_value = unfav_hrd.max()
unfav_hrd_no_outlier = unfav_hrd.drop(index=outlier_sample)

p_no_outlier = mannwhitneyu(
    fav_hrd,
    unfav_hrd_no_outlier,
    alternative="two-sided"
).pvalue

print("\nHRD sensitivity analysis")
print(f"Original:")
print(f"Favorable: n={len(fav_hrd)}, median={fav_hrd.median():.1f}")
print(f"Unfavorable: n={len(unfav_hrd)}, median={unfav_hrd.median():.1f}")
print(f"Mann-Whitney p={p_original:.6f}")

print(f"\nAfter removing highest unfavorable HRD activity ({outlier_value:.1f}):")
print(f"Favorable: n={len(fav_hrd)}, median={fav_hrd.median():.1f}")
print(f"Unfavorable: n={len(unfav_hrd_no_outlier)}, median={unfav_hrd_no_outlier.median():.1f}")
print(f"Mann-Whitney p={p_no_outlier:.6f}")

def remove_sample(data_by_type, sample):
    return {
        modality: data.drop(index=sample, errors="ignore").copy()
        for modality, data in data_by_type.items()
    }

unfavorable_no_outlier = remove_sample(unfavorable, outlier_sample)

fav_share = calculate_activity_share(favorable)
unfav_share = calculate_activity_share(unfavorable)
unfav_share_no_outlier = calculate_activity_share(unfavorable_no_outlier)

print("\nRelative HRD activity:")
print(f"Favorable: {fav_share.get('HRD', 0):.1f}%")
print(f"Unfavorable, original: {unfav_share.get('HRD', 0):.1f}%")
print(f"Unfavorable, without outlier: {unfav_share_no_outlier.get('HRD', 0):.1f}%")

print("\nHRD activity difference:")
print(f"Original: "
    f"{unfav_share.get('HRD', 0) - fav_share.get('HRD', 0):+.1f} pp")
print(
    f"Without outlier: "
    f"{unfav_share_no_outlier.get('HRD', 0) - fav_share.get('HRD', 0):+.1f} pp")

results = prepare_results(favorable, unfavorable)

results.to_csv(
    OUTPUT_DIR / "etiology_comparison_statistics.csv"
)

print("\nStatistical results:")
print(
    results[
        [
            "Favorable prevalence",
            "Unfavorable prevalence",
            "prevalence_p",
            "prevalence_FDR",
            "Favorable activity",
            "Unfavorable activity",
            "activity_p",
            "activity_FDR"
        ]
    ].round(4).to_string()
)

print("\nSignificant prevalence differences (FDR < 0.05):")
significant_prevalence = results[results["prevalence_FDR"] < 0.05]

if significant_prevalence.empty:
    print("None")
else:
    print(
        significant_prevalence[
            [
                "Favorable prevalence",
                "Unfavorable prevalence",
                "prevalence_FDR"
            ]
        ].round(4).to_string()
    )

print("\nSignificant activity differences (FDR < 0.05):")
significant_activity = results[results["activity_FDR"] < 0.05]

if significant_activity.empty:
    print("None")
else:
    print(
        significant_activity[
            [
                "Favorable activity",
                "Unfavorable activity",
                "activity_FDR"
            ]
        ].round(4).to_string()
    )

print("\nPrevalence:")
print(results[["prevalence_p", "prevalence_FDR"]].sort_values("prevalence_p"))

print("\nActivity:")
print(results[["activity_p", "activity_FDR"]].sort_values("activity_p"))

create_difference_plot(results)
create_hrd_activity_boxplot(favorable, unfavorable)