import pandas as pd

cohorts = {
    "immu": {
        "expression": "R_output/immu_expression_matrix.tsv",
        "signatures": "../SBS_assignment/IMMUcan_after_filtering/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
    },
    "synergy": {
        "expression": "R_output/synergy_expression_matrix.tsv",
        "signatures": "../SBS_assignment/SYNERGY_after_filtering/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
    }
}

def core_id(name):
    name = str(name).split("_vs_")[0]
    return name.rsplit("-", 2)[0]

def filter_cohort(name, expression_file, signatures_file):
    expr = pd.read_csv(expression_file, sep="\t", index_col=0)
    sigs = pd.read_csv(signatures_file, sep="\t", index_col=0).T

    expr_map = {core_id(col): col for col in expr.columns}

    sigs.columns = [col.split("_vs_")[0] for col in sigs.columns]
    sigs = sigs.loc[:, ~pd.Index([core_id(col) for col in sigs.columns]).duplicated()]
    sigs_map = {core_id(col): col for col in sigs.columns}

    common_ids = [x for x in expr_map if x in sigs_map]

    expr = expr[[expr_map[x] for x in common_ids]]
    sigs = sigs[[sigs_map[x] for x in common_ids]]

    expr.columns = common_ids
    sigs.columns = common_ids

    expr.to_csv(f"{name}_expression_filtered.tsv", sep="\t")
    sigs.to_csv(f"{name}_signatures_filtered.tsv", sep="\t")

    print(f"{name}: {len(common_ids)} matched samples")

for name, files in cohorts.items():
    filter_cohort(name, files["expression"], files["signatures"])