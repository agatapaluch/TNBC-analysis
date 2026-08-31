import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys
import os
from pathlib import Path
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import warnings

warnings.filterwarnings("ignore")

SIGNATURE_ETIOLOGY = {
    # Aging
    'SBS1': 'Clock-like', 'SBS5': 'Clock-like',
    # Smoking
    'SBS4': 'Tobacco', 'SBS29': 'Tobacco', 'SBS92': 'Tobacco', 'SBS100': 'Tobacco',
    'DBS2': 'Tobacco/acetaldehyde',
    # UV exposure
    'SBS7a': 'UV', 'SBS7b': 'UV', 'SBS7c': 'UV', 'SBS7d': 'UV', 'SBS38': 'UV',
    'DBS1': 'UV', 'ID13': 'UV',
    # DNA repair deficiency
    'SBS3': 'HRD',
    'SBS6': 'dMMR', 'SBS14': 'dMMR', 'SBS15': 'dMMR', 'SBS20': 'dMMR', 'SBS21': 'dMMR', 'SBS26': 'dMMR', 'SBS44': 'dMMR',
    'SBS30': 'dBER', 'SBS36': 'dBER',
    'DBS13': 'HRD', 'ID6': 'HRD',
    'DBS10': 'dMMR', 'ID7': 'dMMR',
    'ID8': 'NHEJ',
    # Polymerase
    'SBS9': 'Polymerase mutation', 'SBS10a': 'Polymerase mutation', 'SBS10b': 'Polymerase mutation',
    'SBS10c': 'Polymerase mutation', 'SBS10d': 'Polymerase mutation',
    'DBS3': 'Polymerase mutation',
    # Slippage
    'ID1': 'Replication slippage', 'ID2': 'Replication slippage',
    # AID/APOBEC
    'SBS2': 'APOBEC activity', 'SBS13': 'APOBEC activity', 'SBS84': 'AID activity',
    'SBS85': 'AID activity',
    # Cancer treatment
    'SBS11': 'Cancer treatment', 'SBS25': 'Cancer treatment', 'SBS31': 'Cancer treatment', 'SBS32': 'Cancer treatment',
    'SBS35': 'Cancer treatment', 'SBS86': 'Cancer treatment', 'SBS87': 'Cancer treatment', 'SBS90': 'Cancer treatment',
    'SBS99': 'Cancer treatment',
    'DBS5': 'Cancer treatment',
    # Aristolochic acid
    'SBS22a': 'Aristolochic acid', 'SBS22b': 'Aristolochic acid',
    'DBS20': 'Aristolochic acid',
    'ID23': 'Aristolochic acid',
    # Aflatoxin
    'SBS24': 'Aflatoxin',
    # Colibactin
    'SBS88': 'Colibactin',
    'ID18': 'Colibactin',
    # Haloalkanes
    'SBS42': 'Haloalkane exposure',
    # ROS
    'SBS18': 'ROS',
    # Topoisomerase activity
    'ID4': 'TOP mutation', 'ID17': 'TOP mutation',
    # Sequencing artefacts
    'SBS27': 'Sequencing artefacts', 'SBS43': 'Sequencing artefacts', 'SBS45': 'Sequencing artefacts',
    'SBS46': 'Sequencing artefacts', 'SBS47': 'Sequencing artefacts', 'SBS48': 'Sequencing artefacts',
    'SBS49': 'Sequencing artefacts', 'SBS50': 'Sequencing artefacts', 'SBS51': 'Sequencing artefacts',
    'SBS52': 'Sequencing artefacts', 'SBS53': 'Sequencing artefacts', 'SBS54': 'Sequencing artefacts',
    'SBS55': 'Sequencing artefacts', 'SBS56': 'Sequencing artefacts', 'SBS57': 'Sequencing artefacts',
    'SBS58': 'Sequencing artefacts', 'SBS59': 'Sequencing artefacts', 'SBS60': 'Sequencing artefacts',
    'SBS95': 'Sequencing artefacts',
    # Unknown
    'SBS8': 'Unknown', 'SBS12': 'Unknown', 'SBS16': 'Unknown', 'SBS17a': 'Unknown', 'SBS17b': 'Unknown',
    'SBS19': 'Unknown', 'SBS23': 'Unknown', 'SBS28': 'Unknown',
    'SBS33': 'Unknown', 'SBS34': 'Unknown', 'SBS37': 'Unknown', 'SBS39': 'Unknown',
    'SBS40a': 'Unknown', 'SBS40b': 'Unknown', 'SBS40c': 'Unknown', 'SBS41': 'Unknown',
    'SBS89': 'Unknown', 'SBS91': 'Unknown', 'SBS93': 'Unknown',
    'SBS94': 'Unknown', 'SBS96': 'Unknown', 'SBS97': 'Unknown',
    'SBS98': 'Unknown',
    'DBS4': 'Unknown', 'DBS6': 'Unknown', 'DBS8': 'Unknown',
    'DBS9': 'Unknown', 'DBS11': 'Unknown', 'DBS12': 'Unknown', 'DBS14': 'Unknown',
    'DBS15': 'Unknown', 'DBS16': 'Unknown', 'DBS17': 'Unknown', 'DBS18': 'Unknown',
    'DBS19': 'Unknown', 'ID5': 'Unknown', 'ID9': 'Unknown', 'ID10': 'Unknown',
    'ID11': 'Unknown', 'ID12': 'Unknown', 'ID14': 'Unknown', 'ID15': 'Unknown',
    'ID16': 'Unknown', 'ID19': 'Unknown', 'ID20': 'Unknown', 'ID21': 'Unknown',
    'ID22': 'Unknown',
}

ETIOLOGY_COLORS = {
    'Clock-like': '#64A8ED',
    'Tobacco': '#755744',
    'Tobacco/acetaldehyde': '#91725e',
    'UV': '#8172B2',
    'HRD': '#F5112E',
    'dMMR': '#65A165',
    'dBER': '#D1EDAF',
    'NHEJ': '#AB88F2',
    'Polymerase mutation': '#F5F187',
    'Replication slippage': '#64B5CD',
    'APOBEC activity': '#DE8CB8',
    'AID activity': '#FCAEF8',
    'Cancer treatment': '#E17C05',
    'Aristolochic acid': '#F7E4EC',
    'Aflatoxin': '#DAC2F2',
    'Colibactin': '#632600',
    'Haloalkane exposure': '#DCFADC',
    'ROS': '#DBF5FF',
    'TOP mutation': '#FF9DA7',
    'Sequencing artefacts': '#F5E5BF',
    'Unknown': '#D0D0D0'
}


SIGNATURE_TYPE_COLORS = {
    'SBS': '#4C72B0',
    'DBS': '#AF8CE6',
    'ID': '#F59FD2'
}

def find_input_file(input_path):
    """ Find the Assignment_Solution_Activities.txt file """
    if os.path.isfile(input_path):
        return input_path, extract_dataset_name(input_path)

    activities_file = os.path.join(input_path, "Assignment_Solution", "Activities", "Assignment_Solution_Activities.txt")

    if os.path.exists(activities_file):
        dataset_name = extract_dataset_name(input_path)
        return activities_file, dataset_name
    else:
        print(f"Error: Could not find Assignment_Solution_Activities.txt in {input_path}")
        print(f"Expected path: {activities_file}")
        return None, None


def extract_dataset_name(path):
    """ Extract dataset name from path """
    path_parts = Path(path).parts
    return path_parts[-1]


def load_data(input_file):
    """ Load mutational signatures data """
    try:
        data = pd.read_csv(input_file, sep='\t', index_col=0)
        data = data.apply(pd.to_numeric, errors='coerce').fillna(0)
        print(f"Successfully loaded data: {data.shape[0]} samples, {data.shape[1]} signatures")
        return data
    except Exception as e:
        print(f"Error loading data: {e}")
        return None

TITLE_FONTSIZE = 18
SUBTITLE_FONTSIZE = 16
AXIS_LABEL_FONTSIZE = 16
TICK_LABEL_FONTSIZE = 16
VALUE_LABEL_FONTSIZE = 14
PIE_TITLE_FONTSIZE = 16
PIE_PERCENT_FONTSIZE = 16
LEGEND_FONTSIZE = 15
LEGEND_COMBINED_FONTSIZE = 15
LEGEND_TITLE_FONTSIZE = 16

def calculate_figure_width(n_signatures, min_width=12, max_width=50):
    """ Calculate optimal figure width based on number of signatures """
    width = max(min_width, min(max_width, n_signatures * 0.5))
    return width

def format_count_unit(count, singular, plural):
    return singular if count == 1 else plural

def style_standard_axes(ax, x_rotation=90):
    ax.tick_params(axis='x', labelsize=TICK_LABEL_FONTSIZE, rotation=x_rotation, width=1.4, length=6)
    ax.tick_params(axis='y', labelsize=TICK_LABEL_FONTSIZE, width=1.4, length=6)

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)


def get_signature_type(signature):
    if signature.startswith('SBS'):
        return 'SBS'
    if signature.startswith('DBS'):
        return 'DBS'
    if signature.startswith('ID'):
        return 'ID'
    return 'Other'


def get_signature_type_sample_sizes(data):
    """Return the number of samples actually available for SBS, DBS and ID."""
    sizes = {}
    for signature_type in ['SBS', 'DBS', 'ID']:
        columns = [col for col in data.columns if get_signature_type(col) == signature_type]
        if columns:
            sizes[signature_type] = int(data[columns].notna().any(axis=1).sum())
    return sizes


def format_signature_type_sample_sizes(sample_sizes):
    return ', '.join(f'{signature_type}: n = {n}' for signature_type, n in sample_sizes.items())


def add_signature_type_annotations(ax, signatures):
    """Add SBS/DBS/ID labels and separators to combined barplots and boxplots."""
    groups = []
    for signature_type in ['SBS', 'DBS', 'ID']:
        positions = [i for i, signature in enumerate(signatures) if get_signature_type(signature) == signature_type]
        if positions:
            groups.append((signature_type, min(positions), max(positions)))

    for signature_type, start, end in groups:
        middle = (start + end) / 2
        ax.text(middle, 1.015, signature_type, transform=ax.get_xaxis_transform(), ha='center', va='bottom',
                fontsize=SUBTITLE_FONTSIZE, fontweight='bold', clip_on=False)

def add_signature_type_legend(ax, sample_sizes=None):
    handles = []

    for signature_type in ['SBS', 'DBS', 'ID']:
        if sample_sizes and signature_type not in sample_sizes:
            continue

        label = (
            f"{signature_type} (n = {sample_sizes[signature_type]})" if sample_sizes else signature_type)

        handles.append(
            Line2D(
                [0], [0],
                marker='o',
                linestyle='None',
                markerfacecolor=SIGNATURE_TYPE_COLORS[signature_type],
                markeredgecolor='none',
                markersize=10,
                label=label))

    legend = ax.legend(
        handles=handles,
        loc='upper right',
        bbox_to_anchor=(0.98, 0.98),
        frameon=True,
        fancybox=False,
        framealpha=1,
        edgecolor='black',
        fontsize=15,
        handlelength=1,
        handletextpad=0.6,
        labelspacing=0.6,
        borderpad=0.6)

    legend.get_frame().set_linewidth(0.8)


def load_combined_data(input_paths):
    """Load and combine SBS, DBS and ID Assignment activity matrices."""
    data_by_type = {}
    dataset_names = []

    for signature_type, input_path in zip(['SBS', 'DBS', 'ID'], input_paths):
        input_file, dataset_name = find_input_file(input_path)
        if input_file is None:
            return None, None, None

        print(f'\nLoading {signature_type} data from: {input_file}')
        signature_data = load_data(input_file)
        if signature_data is None:
            return None, None, None

        signature_data = signature_data[~signature_data.index.duplicated(keep='first')]

        columns = [col for col in signature_data.columns if col.startswith(signature_type)]
        if not columns:
            print(f'Error: no {signature_type} signature columns found in {input_file}')
            return None, None, None

        data_by_type[signature_type] = signature_data[columns]
        dataset_names.append(dataset_name)

    combined_data = pd.concat([data_by_type['SBS'], data_by_type['DBS'], data_by_type['ID']], axis=1, join='outer')

    if combined_data.columns.duplicated().any():
        duplicated = combined_data.columns[combined_data.columns.duplicated()].tolist()
        print(f'Error: duplicated signature columns found: {duplicated}')
        return None, None, None

    sample_sizes = {signature_type: data_by_type[signature_type].shape[0] for signature_type in ['SBS', 'DBS', 'ID']}
    common_prefix = os.path.commonprefix(dataset_names).rstrip('._- ')
    dataset_name = f'{common_prefix}_SBS_DBS_ID' if common_prefix else 'SBS_DBS_ID_combined'

    print('\nCombined dataset:')
    for signature_type in ['SBS', 'DBS', 'ID']:
        print(f'{signature_type}: {sample_sizes[signature_type]} samples, {data_by_type[signature_type].shape[1]} signatures')
    print(f'Unique samples overall: {combined_data.shape[0]}')
    print(f'Total signatures: {combined_data.shape[1]}')

    return combined_data, dataset_name, sample_sizes

def style_pie_percentages(autotexts):
    for autotext in autotexts:
        autotext.set_fontsize(PIE_PERCENT_FONTSIZE)
        autotext.set_color('black')
        autotext.set_fontweight('bold')

def autopct_visible(pct):
    return f'{pct:.1f}%' if pct >= 5 else ''

def calculate_etiology_prevalence(data, sample_id=None):
    """
    Calculate etiology prevalence.

    Cohort-level: number of samples with at least one active signature assigned to a given etiology.
    Sample-level: number of active signatures assigned to a given etiology in one sample.
    """
    etiology_counts = {}

    if sample_id:
        patient_data = data.loc[sample_id]
        active_signatures = patient_data[patient_data > 0].index

        for sig in active_signatures:
            etiology = SIGNATURE_ETIOLOGY.get(sig, 'Unknown')
            etiology_counts[etiology] = etiology_counts.get(etiology, 0) + 1

        return etiology_counts

    for etiology_name in set(SIGNATURE_ETIOLOGY.values()) | {'Unknown'}:
        etiology_signatures = [sig for sig, etiol in SIGNATURE_ETIOLOGY.items() if etiol == etiology_name]

        if etiology_name == 'Unknown':
            all_mapped_sigs = set(SIGNATURE_ETIOLOGY.keys())
            unknown_sigs = set(data.columns) - all_mapped_sigs
            etiology_signatures.extend(unknown_sigs)

        etiology_signatures_in_data = [sig for sig in etiology_signatures if sig in data.columns]

        if etiology_signatures_in_data:
            samples_with_etiology = (data[etiology_signatures_in_data] > 0).any(axis=1).sum()
            if samples_with_etiology > 0:
                etiology_counts[etiology_name] = int(samples_with_etiology)

    return etiology_counts

def calculate_etiology_activity(data, sample_id=None):
    """
    Calculate summed activity per etiology.

    Cohort-level: sum of all signature activities assigned to each etiology across all samples.
    Sample-level: sum of active signature activities assigned to each etiology in one sample.
    """
    etiology_activity = {}

    if sample_id:
        patient_data = data.loc[sample_id]
        active_signatures = patient_data[patient_data > 0]

        for sig, activity in active_signatures.items():
            etiology = SIGNATURE_ETIOLOGY.get(sig, 'Unknown')
            etiology_activity[etiology] = etiology_activity.get(etiology, 0) + float(activity)

        return etiology_activity

    for sig in data.columns:
        total_activity = data[sig].sum()
        if total_activity > 0:
            etiology = SIGNATURE_ETIOLOGY.get(sig, 'Unknown')
            etiology_activity[etiology] = etiology_activity.get(etiology, 0) + float(total_activity)

    return etiology_activity


def sorted_metric_dict(metric_dict):
    return sorted(metric_dict.items(), key=lambda x: x[1], reverse=True)


def create_boxplot(data, output_dir, dataset_name, show_outliers=True, show_only=False, y_max=None,
                   outlier_threshold=None, combined=False, sample_sizes=None):
    """ Create boxplot of signature activities """

    if outlier_threshold is not None:
        data = data.copy()
        data[data > outlier_threshold] = outlier_threshold
        outlier_suffix = f"_filtered_{outlier_threshold}"
    else:
        outlier_suffix = ""

    if combined:
        active_data = data.loc[:, (data.fillna(0) != 0).any(axis=0)]
    else:
        active_data = data.loc[:, (data != 0).any(axis=0)]
    n_signatures = active_data.shape[1]

    if n_signatures == 0:
        print("No active signatures found for boxplot.")
        return

    width = calculate_figure_width(n_signatures)
    fig, ax = plt.subplots(figsize=(width, 8))
    data_melted = active_data.melt(var_name='Signature', value_name='Activity')

    if combined:
        palette = {sig: SIGNATURE_TYPE_COLORS.get(get_signature_type(sig), '#999999') for sig in active_data.columns}
        sns.boxplot(data=data_melted, x='Signature', y='Activity', order=active_data.columns,
                    palette=palette, showfliers=show_outliers, ax=ax)
    else:
        sns.boxplot(data=data_melted, x='Signature', y='Activity', showfliers=show_outliers, ax=ax)

    if show_outliers:
        if y_max is not None:
            ax.set_ylim(0 - (y_max * 0.04), y_max)
        if outlier_threshold is None:
            title_suffix = "with outliers" if combined else f"with outliers, n = {data.shape[0]} samples"
        else:
            title_suffix = f"filtered at {outlier_threshold}"
        filename = f"{dataset_name}_boxplot{outlier_suffix}_with_outliers.png"
    else:
        title_suffix = "without outliers" if combined else f"without outliers, n = {data.shape[0]} samples"
        filename = f"{dataset_name}_boxplot_no_outliers.png"

    if combined:
        sample_sizes = sample_sizes or get_signature_type_sample_sizes(data)
        title = f"Mutational Signature Activities ({title_suffix})"
        title_pad = 10
    else:
        title = f"Mutational Signature Activities ({title_suffix})"
        title_pad = 10

    ax.set_title(title, fontweight='bold', fontsize=TITLE_FONTSIZE, pad=title_pad)
    ax.set_xlabel('Mutational Signatures', labelpad=20, fontsize=AXIS_LABEL_FONTSIZE, fontweight='bold')
    ax.set_ylabel('Activity Level', labelpad=20, fontsize=AXIS_LABEL_FONTSIZE, fontweight='bold')
    style_standard_axes(ax, x_rotation=90)
    ax.grid(axis='y', alpha=0.25, zorder=0)

    if combined:
        add_signature_type_legend(ax, sample_sizes=sample_sizes)

    if show_only:
        plt.subplots_adjust(bottom=0.2)
        plt.show()
    else:
        plt.tight_layout()
        output_path = os.path.join(output_dir, filename)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Boxplot saved: {output_path}")
    plt.close()

def create_sample_activity_plot(data, output_dir, dataset_name, sample_id, show_only=False):
    """ Create activity plot for a single sample showing activity levels as points """
    if sample_id not in data.index:
        print(f"Error: Sample {sample_id} not found in data")
        return

    sample_data = data.loc[sample_id]
    active_data = sample_data[sample_data > 0]

    if len(active_data) == 0:
        print(f"No active signatures found for sample {sample_id}.")
        return

    width = calculate_figure_width(len(active_data))
    fig, ax = plt.subplots(figsize=(width, 7))

    x_pos = range(len(active_data))
    ax.scatter(x_pos, active_data.values, s=130, color='steelblue', zorder=3)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(active_data.index, rotation=90)
    ax.set_xlabel('Mutational Signatures', labelpad=20, fontsize=AXIS_LABEL_FONTSIZE, fontweight='bold')
    ax.set_ylabel('Activity Level', labelpad=20, fontsize=AXIS_LABEL_FONTSIZE, fontweight='bold')
    ax.set_title(
        f'Mutational Signature Activities\nSample: {sample_id}',
        fontweight='bold',
        fontsize=TITLE_FONTSIZE,
        pad=18)
    style_standard_axes(ax, x_rotation=90)
    ax.grid(axis='y', alpha=0.3, zorder=0)

    if show_only:
        plt.subplots_adjust(bottom=0.2)
        plt.show()
    else:
        plt.tight_layout()
        output_path = os.path.join(output_dir, f"{dataset_name}_activity_{sample_id}.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Activity plot saved: {output_path}")
    plt.close()


def create_barplot_active_signatures(data, output_dir, dataset_name, show_only=False, combined=False, sample_sizes=None):
    """ Create barplot showing number of samples with active signatures """
    active_counts = (data > 0).sum()
    active_counts = active_counts[active_counts > 0]

    if len(active_counts) == 0:
        print("No active signatures found for barplot.")
        return

    width = calculate_figure_width(len(active_counts))
    fig, ax = plt.subplots(figsize=(width, 8))

    if combined:
        bar_colors = [SIGNATURE_TYPE_COLORS.get(get_signature_type(sig), '#999999') for sig in active_counts.index]
        bars = ax.bar(range(len(active_counts)), active_counts.values, color=bar_colors, zorder=3, edgecolor='black', linewidth=0.5)
        sample_sizes = sample_sizes or get_signature_type_sample_sizes(data)
        title = f"Number of Samples with Active Signatures"
        title_pad = 10
        y_limit = active_counts.max() * 1.20 + 0.5
    else:
        bars = ax.bar(range(len(active_counts)), active_counts.values, zorder=3, edgecolor='black', linewidth=0.2)
        title = f"Number of Samples with Active Signatures (cohort size: n = {data.shape[0]} samples)"
        title_pad = 10
        y_limit = max(data.shape[0], active_counts.max()) * 1.12

    ax.set_title(title, fontweight='bold', fontsize=TITLE_FONTSIZE, pad=title_pad)
    ax.set_xlabel('Mutational Signatures', labelpad=22, fontsize=AXIS_LABEL_FONTSIZE, fontweight='bold')
    ax.set_ylabel('Number of Samples (with activity > 0)', labelpad=22, fontsize=AXIS_LABEL_FONTSIZE, fontweight='bold')
    ax.set_xticks(range(len(active_counts)))
    ax.set_xticklabels(active_counts.index, rotation=90)
    ax.set_ylim(0, y_limit)
    style_standard_axes(ax, x_rotation=90)
    ax.grid(axis='y', alpha=0.25, zorder=0)

    for i, v in enumerate(active_counts.values):
        ax.text(i, v + 0.2, str(v), ha='center', va='bottom', fontsize=VALUE_LABEL_FONTSIZE)

    if combined:
        add_signature_type_legend(ax, sample_sizes=sample_sizes)

    if show_only:
        plt.subplots_adjust(bottom=0.2)
        plt.show()
    else:
        plt.tight_layout()
        output_path = os.path.join(output_dir, f"{dataset_name}_active_signatures_barplot.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Active signatures barplot saved: {output_path}")
    plt.close()

def create_etiology_piechart(data, output_dir, dataset_name, sample_id=None, show_only=False, sample_sizes=None):
    """ Create pie chart of signature etiologies based on prevalence or signature count """
    if sample_id and sample_id not in data.index:
        print(f"Error: Sample {sample_id} not found in data.")
        return

    etiology_counts = calculate_etiology_prevalence(data, sample_id=sample_id)

    if sample_id:
        # For single sample: count active signatures per etiology
        title_text = f'Active Signatures Etiologies\nSample {sample_id}'
        filename = f"{dataset_name}_etiology_piechart_{sample_id}.png"
        singular_unit = "signature"
        plural_unit = "signatures"
    else:
        title_text = (f"Etiology Prevalence" if sample_sizes
                      else f'Etiology Prevalence (n = {data.shape[0]} samples)')
        filename = f"{dataset_name}_etiology_piechart_overall.png"
        singular_unit = "sample"
        plural_unit = "samples"

    if not etiology_counts:
        print("No active signatures found for pie chart.")
        return

    # Sorted legend
    sorted_etiologies = sorted(etiology_counts.items(), key=lambda x: x[1], reverse=True)
    sorted_values = [item[1] for item in sorted_etiologies]
    total = sum(sorted_values)

    fig, ax = plt.subplots(figsize=(14, 7))
    ax.set_position([0.1, 0.1, 0.5, 0.8])

    colors = [ETIOLOGY_COLORS.get(e, '#CCCCCC') for e, _ in sorted_etiologies]

    wedges, texts, autotexts = ax.pie(sorted_values, startangle=90, colors=colors, wedgeprops=dict(width=0.4), autopct=autopct_visible, pctdistance=0.8)
    style_pie_percentages(autotexts)

    legend_labels = []
    for etiology, count in sorted_etiologies:
        unit = format_count_unit(count, singular_unit, plural_unit)
        legend_labels.append(f"{etiology} ({count} {unit}, {count / total * 100:.1f}%)")

    legend = ax.legend(wedges, legend_labels, title="SIGNATURES ETIOLOGY\n", fontsize=LEGEND_FONTSIZE, loc="center left",
                       bbox_to_anchor=(1, 0, 0.5, 1.0), frameon=False)

    legend.get_title().set_fontweight('bold')
    legend.get_title().set_fontsize(LEGEND_TITLE_FONTSIZE)

    ax.set_title(title_text, fontsize=PIE_TITLE_FONTSIZE, fontweight='bold')

    if show_only:
        plt.show()
    else:
        plt.tight_layout()
        output_path = os.path.join(output_dir, filename)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Etiology pie chart saved: {output_path}")
    plt.close()


def create_etiology_activity_piechart(data, output_dir, dataset_name, sample_id=None, show_only=False, sample_sizes=None):
    """ Create pie chart of signature etiologies based on total activity """

    if sample_id and sample_id not in data.index:
        print(f"Error: Sample {sample_id} not found in data.")
        return

    etiology_activity = calculate_etiology_activity(data, sample_id=sample_id)

    if sample_id:
        # Single sample
        title_text = f'Etiology Activity\nSample: {sample_id}'
        filename = f"{dataset_name}_etiology_activity_piechart_{sample_id}.png"
    else:
        title_text = (f"Etiology Activity" if sample_sizes
                      else f'Etiology Activity (n = {data.shape[0]} samples)')
        filename = f'{dataset_name}_etiology_activity_piechart_overall.png'

    if not etiology_activity:
        print("No active signatures found for activity-based pie chart.")
        return

    sorted_items = sorted_metric_dict(etiology_activity)
    values = [value for _, value in sorted_items]
    total = sum(values)

    fig, ax = plt.subplots(figsize=(14, 7))
    ax.set_position([0.1, 0.1, 0.5, 0.8])


    colors = [ETIOLOGY_COLORS.get(e, '#CCCCCC') for e, _ in sorted_items]

    wedges, texts, autotexts = ax.pie(values, colors=colors, startangle=90, wedgeprops=dict(width=0.4), autopct=autopct_visible, pctdistance=0.8)
    style_pie_percentages(autotexts)

    legend_labels = [f"{etiology} ({value:.0f} activity, {value / total * 100:.1f}%)" for etiology, value in sorted_items]

    legend = ax.legend(wedges, legend_labels, title="SIGNATURES ETIOLOGY\n", fontsize=LEGEND_FONTSIZE, loc="center left", bbox_to_anchor=(1, 0, 0.5, 1), frameon=False)

    legend.get_title().set_fontweight('bold')
    legend.get_title().set_fontsize(LEGEND_TITLE_FONTSIZE)

    ax.set_title(title_text, fontsize=PIE_TITLE_FONTSIZE, fontweight='bold')

    if show_only:
        plt.show()
    else:
        plt.tight_layout()
        output_path = os.path.join(output_dir, filename)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Etiology activity pie chart saved: {output_path}")

    plt.close()

def create_etiology_combined_piecharts(data, output_dir, dataset_name, sample_id=None, show_only=False):
    """
    Create one figure with two etiology donut charts side by side:
    1. Prevalence: number of samples/signatures per etiology.
    2. Activity: summed mutational signature activity per etiology.
    """
    if sample_id and sample_id not in data.index:
        print(f"Error: Sample {sample_id} not found in data.")
        return

    prevalence = calculate_etiology_prevalence(data, sample_id=sample_id)
    activity = calculate_etiology_activity(data, sample_id=sample_id)

    if not prevalence or not activity:
        print("No active signatures found for combined etiology pie charts.")
        return

    all_etiologies = sorted(
        set(prevalence.keys()) | set(activity.keys()),
        key=lambda e: max(prevalence.get(e, 0), activity.get(e, 0)),
        reverse=True)

    prevalence_values = [prevalence.get(e, 0) for e in all_etiologies]
    activity_values = [activity.get(e, 0) for e in all_etiologies]
    colors = [ETIOLOGY_COLORS.get(e, '#CCCCCC') for e in all_etiologies]

    fig, axes = plt.subplots(1, 2, figsize=(14, 9))


    # Left pie: prevalence
    wedges_prev, _, autotexts_prev = axes[0].pie(
        prevalence_values,
        startangle=90,
        colors=colors,
        radius=1.35,
        wedgeprops=dict(width=0.65, linewidth=0.5, edgecolor="#777777"),
        autopct=autopct_visible,
        pctdistance=0.78)
    style_pie_percentages(autotexts_prev)

    if sample_id:
        prevalence_title = 'Prevalence by etiology\nactive signatures in this sample'
        figure_title = f'Etiology of Active Mutational Signatures - Sample {sample_id}'
        filename = f"{dataset_name}_etiology_combined_piecharts_{sample_id}.png"
        prev_singular, prev_plural = "signature", "signatures"
    else:
        prevalence_title = 'Etiology Prevalence'
        figure_title = f'Etiology of Active Mutational Signatures (n = {data.shape[0]} samples)'
        filename = f"{dataset_name}_etiology_combined_piecharts_overall.png"
        prev_singular, prev_plural = "sample", "samples"

    axes[0].set_title(prevalence_title, fontsize=PIE_TITLE_FONTSIZE, fontweight='bold', pad=35)
    #axes[0].set_position([0.00, 0.54, 0.58, 0.40])
    axes[0].set_position([0.05, 0.30, 0.40, 0.52])

    # Right pie: summed activity
    wedges_act, _, autotexts_act = axes[1].pie(
        activity_values,
        startangle=90,
        colors=colors,
        radius=1.35,
        wedgeprops=dict(width=0.65, linewidth=0.5, edgecolor="#777777"),
        autopct=autopct_visible,
        pctdistance=0.78
    )
    style_pie_percentages(autotexts_act)
    axes[1].set_title('Etiology Activity', fontsize=PIE_TITLE_FONTSIZE, fontweight='bold', pad=35)
    #axes[1].set_position([0.00, 0.08, 0.58, 0.40])
    axes[1].set_position([0.55, 0.30, 0.40, 0.52])

    # Shared legend with both values in one place
    total_prevalence = sum(prevalence_values)
    total_activity = sum(activity_values)
    legend_handles = [Patch(facecolor=ETIOLOGY_COLORS.get(e, '#CCCCCC'), edgecolor='none') for e in all_etiologies]

    legend_labels = []
    for etiology in all_etiologies:
        prev_count = prevalence.get(etiology, 0)
        prev_unit = format_count_unit(prev_count, prev_singular, prev_plural)
        act_value = activity.get(etiology, 0)

        prev_pct = (prev_count / total_prevalence * 100) if total_prevalence else 0
        act_pct = (act_value / total_activity * 100) if total_activity else 0

        legend_labels.append(
            f"{etiology} | {prev_count} {prev_unit} ({prev_pct:.1f}%) | activity: {act_value:.0f} ({act_pct:.1f}%)")

    legend = fig.legend(
        legend_handles,
        legend_labels,
        title="SIGNATURES ETIOLOGY",
        fontsize=LEGEND_FONTSIZE,
        loc='lower center',
        bbox_to_anchor=(0.5, -0.0),
        frameon=False,
        ncol=1,
        handlelength=1.2,
        handletextpad=0.7,
    )
    legend._legend_box.sep = 13
    legend.get_title().set_fontweight('bold')
    legend.get_title().set_fontsize(LEGEND_TITLE_FONTSIZE)

    # plt.subplots_adjust(left=0.01, right=0.66, top=0.88, bottom=0.05, hspace=0.35)
    #plt.subplots_adjust(left=0.04, right=0.96, top=0.90, bottom=0.32, wspace=0.2)
    plt.subplots_adjust(left=0.04, right=0.96, top=0.90, bottom=0.42, wspace=0.2)

    if show_only:
        plt.show()
    else:
        output_path = os.path.join(output_dir, filename)
        plt.savefig(output_path, dpi=300, bbox_inches='tight', pad_inches=0.2)
        print(f"Combined etiology pie charts saved: {output_path}")

    plt.close()


def cluster_signatures(data, output_dir=None, dataset_name=None, method='complete', show_only=False):
    """Cluster signatures hierarchically and visualize as a heatmap with dendrograms."""

    active_data = data.loc[:, (data != 0).any(axis=0)]
    if active_data.shape[1] == 0:
        print("No active signatures found for clustering heatmap.")
        return

    clean_name = dataset_name.replace('.txt', '') if dataset_name else 'heatmap'

    cg = sns.clustermap(
        active_data,
        metric='cosine',
        method=method,
        cmap='viridis',
        figsize=(max(12, active_data.shape[1] * 0.55), 10),
        z_score=1,
        linewidths=0.1,
        cbar_kws={'label': 'Z-score activity'},
        cbar_pos=(1.02, 0.2, 0.01, 0.6)
    )

    cg.ax_heatmap.set_xlabel('Mutational signatures', fontsize=AXIS_LABEL_FONTSIZE, fontweight='bold')
    cg.ax_heatmap.set_ylabel('Samples', fontsize=AXIS_LABEL_FONTSIZE, fontweight='bold')
    cg.ax_heatmap.set_title('Hierarchical Clustering of Mutational Signatures', fontsize=TITLE_FONTSIZE, fontweight='bold')
    cg.ax_heatmap.tick_params(axis='x', labelsize=TICK_LABEL_FONTSIZE, rotation=90)
    cg.ax_heatmap.tick_params(axis='y', labelsize=10)

    if show_only:
        plt.subplots_adjust(right=0.90)
        plt.show()
    else:
        if output_dir is None or dataset_name is None:
            print("Output directory and dataset name required to save clustering heatmap.")
            return
        output_path = os.path.join(output_dir, f"{clean_name}_signature_clustermap.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Signature clustering heatmap saved: {output_path}")
    plt.close()


def generate_summary_report(data, output_dir, dataset_name, sample_sizes=None):
    """ Generate a summary report of the analysis """
    report_path = os.path.join(output_dir, f"{dataset_name}_summary_report.txt")

    with open(report_path, 'w') as f:
        f.write("Mutational Signatures Analysis Summary Report\n")
        f.write("=" * 57 + "\n\n")

        f.write(f"Dataset: {dataset_name}\n")
        f.write(f"Dataset Overview:\n")
        if sample_sizes:
            f.write(f"- Unique samples across all signature types: {data.shape[0]}\n")
            for signature_type, n in sample_sizes.items():
                f.write(f"- {signature_type} samples: {n}\n")
        else:
            f.write(f"- Number of samples: {data.shape[0]}\n")
        f.write(f"- Number of signatures: {data.shape[1]}\n\n")

        # Number of patients with each signature
        prevalence = (data > 0).sum().sort_values(ascending=False)
        f.write("Top 10 Most Prevalent Signatures (by number of samples):\n")
        for i, (sig, count) in enumerate(prevalence.head(10).items(), 1):
            etiology = SIGNATURE_ETIOLOGY.get(sig, 'Unknown')
            if sample_sizes:
                denominator = sample_sizes.get(get_signature_type(sig), data.shape[0])
            else:
                denominator = data.shape[0]
            percentage = (count / denominator) * 100 if denominator else 0
            unit = format_count_unit(int(count), "sample", "samples")
            f.write(f"{i:2d}. {sig}: {count} {unit} ({percentage:.1f}%) - {etiology}\n")
        f.write("\n")

        # Etiology distribution
        etiology_counts = {}
        for sig, activity in data.sum().items():
            if activity > 0:
                etiology = SIGNATURE_ETIOLOGY.get(sig, 'Unknown')
                etiology_counts[etiology] = etiology_counts.get(etiology, 0) + 1

        f.write("Etiology Distribution:\n")
        f.write("-" * 40 + "\n")
        f.write("By number of different signatures:\n")
        for etiology, count in sorted_metric_dict(etiology_counts):
            unit = format_count_unit(count, "signature", "signatures")
            f.write(f"- {etiology}: {count} {unit}\n")

        f.write("\nBy patient prevalence (how many samples have each etiology):\n")
        etiology_prevalence = calculate_etiology_prevalence(data)
        if sample_sizes:
            f.write("Note: percentages are omitted because SBS, DBS and ID have different sample availability.\n")
        for etiology, count in sorted_metric_dict(etiology_prevalence):
            unit = format_count_unit(count, "sample", "samples")
            if sample_sizes:
                f.write(f"- {etiology}: {count} {unit}\n")
            else:
                percentage = (count / data.shape[0]) * 100
                f.write(f"- {etiology}: {count} {unit} ({percentage:.1f}%)\n")

    print(f"Summary report saved: {report_path}")


def main():
    """
    MutSigMA Visualization Tool

    Arguments:
      -i, --input                            Path to one Assignment_Solution folder
      --combined SBS DBS ID                  Paths to SBS, DBS and ID Assignment outputs
      -o, --output                           Output directory. Default: plots
      -b, --boxplot                          Generate boxplot of signature activities
      -y, --y_max                            Optional upper limit for y-axis (applied to boxplot with outliers)
      -x, --outlier_threshold                Filter outliers above this threshold before plotting
      -n, --no_outliers                      Hide outliers in boxplot
      -r, --barplot                          Generate barplot of active signatures count
      -p, --piechart                         Generate pie chart of signature etiologies
      -pa --piechart_activity                Generate pie chart of etiologies based on total signature activity
      -pc --piechart_combined                Generate prevalence and activity pie charts side by side
      -d, --id                               Specific sample ID for individual pie chart
      -f, --first_n                          Generate pie charts for first N patients (default: 0)
      -a, --all                              Generate all visualizations
      -t, --report                           Generate summary report
      -s, --show                             Display plots only (no saving)
      -c, --cluster_signatures               Generate clustering heatmap of signatures
    """

    parser = argparse.ArgumentParser(description='MutSigMA Visualization Tool')

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-i', '--input', help='Path to one Assignment_Solution folder')
    input_group.add_argument('--combined', nargs=3, metavar=('SBS', 'DBS', 'ID'),
                             help='Paths to SBS, DBS and ID Assignment outputs for one combined visualization')
    parser.add_argument('-o', '--output', default='plots', help='Output directory (default: plots)')
    parser.add_argument('-b', '--boxplot', action='store_true', help='Generate boxplot of signature activities')
    parser.add_argument('-y', '--y_max', type=float, default=None, help='Optional upper limit for y-axis (applied to boxplot with outliers)')
    parser.add_argument('-x', '--outlier_threshold', type=float, default=None,
                        help='Filter outliers above this threshold before plotting')
    parser.add_argument('-n', '--no_outliers', action='store_true', help='Hide outliers in boxplot')
    parser.add_argument('-r', '--barplot', action='store_true', help='Generate barplot of active signatures count')
    parser.add_argument('-p', '--piechart', action='store_true', help='Generate pie chart of signature etiologies')
    parser.add_argument('-pa', '--piechart_activity', action='store_true', help='Generate pie chart of signatures etiology based on total signature activity')
    parser.add_argument('-pc', '--piechart_combined', action='store_true', help='Generate prevalence and activity pie charts side by side')
    parser.add_argument('-d', '--id', help='Specific sample ID for individual pie chart')
    parser.add_argument('-f', '--first_n', type=int, default=0,
                        help='Generate pie charts for first N patients (default: 0)')
    parser.add_argument('-a', '--all', action='store_true', help='Generate all visualizations')
    parser.add_argument('-t', '--report', action='store_true', help='Generate summary report')
    parser.add_argument('-s', '--show', action='store_true', help='Display plots only (no saving)')
    parser.add_argument('-c', '--cluster_signatures', action='store_true',
                        help='Generate clustering heatmap of signatures')
    args = parser.parse_args()

    combined_mode = args.combined is not None
    sample_sizes = None

    if combined_mode:
        data, dataset_name, sample_sizes = load_combined_data(args.combined)
        if data is None:
            sys.exit(1)
    else:
        # Original single-input mode
        input_file, dataset_name = find_input_file(args.input)
        if input_file is None:
            sys.exit(1)

        print(f"Using input file: {input_file}")
        print(f"Dataset name: {dataset_name}")

        data = load_data(input_file)
        if data is None:
            sys.exit(1)

    if combined_mode and (args.piechart or args.piechart_combined or args.all):
        print("\nNote: in combined prevalence plots, SBS/DBS/ID may have different sample availability.")

    # Create output directory if not --show arg
    if not args.show:
        os.makedirs(args.output, exist_ok=True)
        print(f"Output directory: {args.output}")

    # If sample ID is specified, generate individual sample plots
    if args.id:
        if args.all or args.boxplot:
            print("\nGenerating activity plot for sample...")
            create_sample_activity_plot(data, args.output, dataset_name, args.id, show_only=args.show)

        if args.all or args.piechart:
            print("\nGenerating pie chart...")
            create_etiology_piechart(data, args.output, dataset_name, sample_id=args.id, show_only=args.show)

        if args.all or args.piechart_activity:
            print("\nGenerating etiology activity pie chart...")
            create_etiology_activity_piechart(data, args.output, dataset_name, sample_id=args.id, show_only=args.show)

        if args.all or args.piechart_combined:
            print("\nGenerating combined etiology pie charts for sample...")
            create_etiology_combined_piecharts(data, args.output, dataset_name, sample_id=args.id, show_only=args.show)

    else:
        if args.all or args.boxplot:
            print("\nGenerating boxplot...")
            if args.all:
                # If -all, generate both
                create_boxplot(data, args.output, dataset_name, show_outliers=True, show_only=args.show, y_max=args.y_max,
                               outlier_threshold=args.outlier_threshold, combined=combined_mode, sample_sizes=sample_sizes)
                create_boxplot(data, args.output, dataset_name, show_outliers=False, show_only=args.show,
                               combined=combined_mode, sample_sizes=sample_sizes)
            else:
                # If --boxplot: only with outliers, if --boxplot --no_outliers: only no outliers
                create_boxplot(data, args.output, dataset_name, show_outliers=not args.no_outliers, show_only=args.show, y_max=args.y_max,
                               outlier_threshold=args.outlier_threshold, combined=combined_mode, sample_sizes=sample_sizes)

        if args.all or args.barplot:
            print("\nGenerating barplot...")
            create_barplot_active_signatures(data, args.output, dataset_name, show_only=args.show,
                                             combined=combined_mode, sample_sizes=sample_sizes)

        if args.first_n > 0:
            sample_list = data.index[:args.first_n]

            if args.all or args.piechart:
                for sample in sample_list:
                    create_etiology_piechart(data, args.output, dataset_name, sample_id=sample, show_only=args.show)

            if args.all or args.piechart_activity:
                print("\nGenerating activity-based etiology pie charts for first samples...")
                for sample in sample_list:
                    create_etiology_activity_piechart(data, args.output, dataset_name, sample_id=sample, show_only=args.show)

            if args.all or args.piechart_combined:
                print("\nGenerating combined etiology pie charts for first samples...")
                for sample in sample_list:
                    create_etiology_combined_piecharts(data, args.output, dataset_name, sample_id=sample, show_only=args.show)

        else:
            if args.all or args.piechart:
                # If --piechart: generate chart for all samples (default)
                print("\nGenerating etiology prevalence pie chart...")
                create_etiology_piechart(data, args.output, dataset_name, sample_id=None, show_only=args.show, sample_sizes=sample_sizes)

            if args.all or args.piechart_activity:
                # Cohort-level activity pie chart
                print("\nGenerating etiology activity pie chart...")
                create_etiology_activity_piechart(data, args.output, dataset_name, sample_id=None, show_only=args.show, sample_sizes=sample_sizes)

            if args.all or args.piechart_combined:
                print("\nGenerating combined etiology pie charts...")
                create_etiology_combined_piecharts(data, args.output, dataset_name, show_only=args.show)

        if args.cluster_signatures: # removing args.all
            print("\nGenerating cluster heatmap...")
            cluster_signatures(data, args.output, dataset_name, show_only=args.show)

        if (args.all or args.report) and not args.show:
            print("\nGenerating summary report...")
            generate_summary_report(data, args.output, dataset_name, sample_sizes=sample_sizes)

    if not any([args.boxplot, args.barplot, args.piechart, args.piechart_activity, args.piechart_combined, args.all, args.report, args.cluster_signatures]):
        print("No visualization option selected. Use --help for available options.")
        print("Quick start: python visualizer.py -i SBS_output --all")
        print("Combined: python visualizer.py --combined SBS_output DBS_output ID_output --all")


if __name__ == "__main__":
    main()
