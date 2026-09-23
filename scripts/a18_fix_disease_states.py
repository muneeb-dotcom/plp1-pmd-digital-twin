import pandas as pd
import numpy as np
import json

files = {
    "WT10": "../plp1-pmd-project/scripts/GSE111605/extracted/GSM3034722_WT10-D1-T3-2.fpkm_tracking.gz",
    "WT13": "../plp1-pmd-project/scripts/GSE111605/extracted/GSM3034723_WT13-D1-T3-2.fpkm_tracking.gz",
    "WT14": "../plp1-pmd-project/scripts/GSE111605/extracted/GSM3034724_WT14-D1-T3-2.fpkm_tracking.gz",
    "jp14": "../plp1-pmd-project/scripts/GSE111605/extracted/GSM3034725_jp14-D1-T3-1.fpkm_tracking.gz",
    "jp15": "../plp1-pmd-project/scripts/GSE111605/extracted/GSM3034726_jp15-D1-T3-1.fpkm_tracking.gz",
    "jp16": "../plp1-pmd-project/scripts/GSE111605/extracted/GSM3034727_jp-16-D1-T3-2.fpkm_tracking.gz",
}

def read_fpkm(path, sample_name):
    df = pd.read_csv(path, sep="\t")[["gene_short_name", "FPKM"]]
    df.columns = ["gene", sample_name]
    return df.groupby("gene", as_index=False).mean()

dfs = [read_fpkm(p, n) for n, p in files.items()]
expr = dfs[0]
for d in dfs[1:]:
    expr = expr.merge(d, on="gene", how="outer")
expr = expr.set_index("gene").fillna(0)

log_expr = np.log2(expr + 1)

with open("results/state_modules.json") as f:
    ref = json.load(f)
MODULE_GENES = ref["genes"]
baseline = ref["healthy_baseline"]

def score_genes_simple(log_expr_df, genes):
    present = [g for g in genes if g in log_expr_df.index]
    if not present:
        return None
    return log_expr_df.loc[present].mean(axis=0)

wt_cols = ["WT10", "WT13", "WT14"]
jp_cols = ["jp14", "jp15", "jp16"]

rows = []
for module, genes in MODULE_GENES.items():
    module_scores = score_genes_simple(log_expr, genes)
    if module_scores is None:
        continue
    wt_mean = module_scores[wt_cols].mean()
    jp_mean = module_scores[jp_cols].mean()
    rows.append({
        "dataset": "jimpy_vs_WT",
        "mechanism": "misfolding_confirmed",
        "module": module,
        "wt_mean_log2fpkm": wt_mean,
        "jimpy_mean_log2fpkm": jp_mean,
        "delta_log2fpkm": jp_mean - wt_mean,
    })

df = pd.DataFrame(rows)
df.to_csv("results/pmd_state_vectors.csv", index=False)
print(df.to_string(index=False))