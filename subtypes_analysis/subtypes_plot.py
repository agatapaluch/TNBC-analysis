import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# SETTINGS

IMMUCAN_FILE = "immucan_result.csv"
SYNERGY_FILE = "synergy_result.csv"

SYNERGY_SAMPLES_TO_REMOVE = ["SYNG.BC1.1522.FIXT.01.RNA.01", "SYNG.BC1.1570.FIXT.02.RNA.01"]

IMMUCAN_EXCLUDED = 13
SYNERGY_EXCLUDED = 24

OUTPUT_DIR = "."
OUTPUT_FILE = "tnbc_subtypes_combined_donuts.png"

TITLE_FONTSIZE = 19
PIE_TITLE_FONTSIZE = 16
LEGEND_TITLE_FONTSIZE = 14
LEGEND_FONTSIZE = 12
PERCENT_FONTSIZE = 14

SUBTYPE_ORDER = ["BL1", "BL2", "IM", "M", "MSL", "LAR", "UNS", "Excluded"]

SUBTYPE_COLORS = {
    "BL1": "#64A8ED",
    "BL2": "#E17C05",
    "IM": "#65A165",
    "M": "#F5112E",
    "MSL": "#DA8BC3",
    "LAR": "#F5F187",
    "UNS": "#B0B0B0",
    "Excluded": "#F2F2F2"}

def make_autopct(labels, min_pct=4):
    counter = {"i": 0}

    def _autopct(pct):
        label = labels[counter["i"]]
        counter["i"] += 1

        if pct >= min_pct or label == "Excluded":
            return rf"$\mathbf{{{label}}}$" + f"\n{pct:.1f}%"

        return ""

    return _autopct


def style_pie_percentages(autotexts):
    for text in autotexts:
        text.set_fontsize(PERCENT_FONTSIZE)
        #text.set_fontweight("bold")
        text.set_color("black")


def prepare_subtype_counts(input_file, excluded_n=0, samples_to_remove=None):
    df = pd.read_csv(input_file, index_col=0)

    if samples_to_remove:
        found = [s for s in samples_to_remove if s in df.index]
        not_found = [s for s in samples_to_remove if s not in df.index]

        print(f"\nSamples removed from {input_file}: {len(found)}")
        print(found)

        if not_found:
            print("WARNING - samples not found:")
            print(not_found)

        df = df[~df.index.isin(samples_to_remove)].copy()

    counts = (
        df["subtype"]
        .value_counts()
        .reindex(SUBTYPE_ORDER, fill_value=0)
    )

    counts["Excluded"] = excluded_n
    counts = counts[counts > 0]

    total_n = counts.sum()

    return counts, total_n

def format_count_unit(count):
    return "sample" if count == 1 else "samples"

def create_combined_tnbc_subtype_donuts(immucan_file, synergy_file, output_dir, output_file):
    immu_counts, immu_n = prepare_subtype_counts(immucan_file, excluded_n=IMMUCAN_EXCLUDED)
    syng_counts, syng_n = prepare_subtype_counts(synergy_file, samples_to_remove=SYNERGY_SAMPLES_TO_REMOVE, excluded_n=SYNERGY_EXCLUDED)

    all_subtypes = [
        subtype for subtype in SUBTYPE_ORDER
        if subtype in immu_counts.index or subtype in syng_counts.index
    ]

    immu_values = [immu_counts.get(subtype, 0) for subtype in all_subtypes]
    syng_values = [syng_counts.get(subtype, 0) for subtype in all_subtypes]
    colors = [SUBTYPE_COLORS.get(subtype, "#CCCCCC") for subtype in all_subtypes]

    fig, axes = plt.subplots(1, 2, figsize=(14, 9))

    wedges_immu, _, autotexts_immu = axes[0].pie(
        immu_values,
        startangle=90,
        colors=colors,
        radius=1.35,
        wedgeprops=dict(width=0.65, linewidth=0.5, edgecolor="#777777"),
        autopct=make_autopct(all_subtypes),
        pctdistance=0.78
    )
    style_pie_percentages(autotexts_immu)

    for text, subtype in zip(autotexts_immu, all_subtypes):
        if subtype == "Excluded":
            x, y = text.get_position()
            text.set_position((x - 0.08, y))

    for wedge, subtype in zip(wedges_immu, all_subtypes):
        if subtype == "Excluded":
            #wedge.set_hatch("///")
            wedge.set_edgecolor("#777777")
            wedge.set_linewidth(0.5)

    axes[0].set_title(
        rf"$\mathbf{{IMMUcan}}$" + f"\n" + f"n = {immu_n}",
        fontsize=PIE_TITLE_FONTSIZE,
        pad=25
    )
    axes[0].set_position([0.05, 0.30, 0.40, 0.52])

    wedges_syng, _, autotexts_syng = axes[1].pie(
        syng_values,
        startangle=90,
        colors=colors,
        radius=1.35,
        wedgeprops=dict(width=0.65, linewidth=0.5, edgecolor="#777777"),
        autopct=make_autopct(all_subtypes),
        pctdistance=0.78
    )
    style_pie_percentages(autotexts_syng)

    for text, subtype in zip(autotexts_syng, all_subtypes):
        if subtype == "Excluded":
            x, y = text.get_position()
            text.set_position((x - 0.08, y))

    for wedge, subtype in zip(wedges_syng, all_subtypes):
        if subtype == "Excluded":
            #wedge.set_hatch("///")
            wedge.set_edgecolor("#777777")
            wedge.set_linewidth(0.5)

    axes[1].set_title(
        rf"$\mathbf{{Synergy}}$" + f"\n" + f"n = {syng_n}",
        fontsize=PIE_TITLE_FONTSIZE,
        pad=25
    )
    axes[1].set_position([0.55, 0.30, 0.40, 0.52])

    legend_handles = [
        Patch(facecolor=SUBTYPE_COLORS.get(subtype, "#CCCCCC"), edgecolor="none")
        for subtype in all_subtypes]

    legend_labels = []
    for subtype in all_subtypes:
        immu_count = immu_counts.get(subtype, 0)
        syng_count = syng_counts.get(subtype, 0)

        immu_pct = (immu_count / immu_n * 100) if immu_n else 0
        syng_pct = (syng_count / syng_n * 100) if syng_n else 0

        legend_labels.append(
            f"{subtype} | "
            f"IMMUcan: {immu_count} {format_count_unit(immu_count)} ({immu_pct:.1f}%) | "
            f"Synergy: {syng_count} {format_count_unit(syng_count)} ({syng_pct:.1f}%)"
        )

    legend = fig.legend(
        legend_handles,
        legend_labels,
        title="TNBC MOLECULAR SUBTYPES",
        fontsize=LEGEND_FONTSIZE,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.01),
        frameon=False,
        ncol=1,
        handlelength=1.2,
        handletextpad=0.7,
    )

    legend._legend_box.sep = 13
    legend.get_title().set_fontweight("bold")
    legend.get_title().set_fontsize(LEGEND_TITLE_FONTSIZE)

    fig.suptitle(
        "Distribution of TNBC molecular subtypes",
        fontsize=TITLE_FONTSIZE,
        fontweight="bold",
        y=0.99
    )

    plt.subplots_adjust(left=0.04, right=0.96, top=0.82, bottom=0.34, wspace=0.05)

    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_file)

    plt.savefig(output_path, dpi=300, bbox_inches="tight", pad_inches=0.2)
    print(f"Combined TNBC subtype donut charts saved: {output_path}")

    plt.show()
    plt.close()

create_combined_tnbc_subtype_donuts(
    immucan_file=IMMUCAN_FILE,
    synergy_file=SYNERGY_FILE,
    output_dir=OUTPUT_DIR,
    output_file=OUTPUT_FILE)