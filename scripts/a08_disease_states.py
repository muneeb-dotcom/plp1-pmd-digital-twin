import pandas as pd
import numpy as np
import json

with open("results/state_modules.json") as f:
    ref = json.load(f)

MODULE_GENES = ref["genes"]
BASELINE = ref["healthy_baseline"]

jimpy = pd.read_csv("../plp1-therapeutics/results/plp1_jimpy_deseq_results.csv", index_col=0)

def sample_module_score(de_df, genes):
    present = [g for g in genes if g in de_df.index]
    return de_df.loc[present, "t"].mean() if present else np.nan

rows = []
for name, genes in MODULE_GENES.items():
    t_mean = sample_module_score(jimpy, genes)
    rows.append({
        "dataset": "jimpy_vs_WT",
        "mechanism": "misfolding_confirmed",
        "module": name,
        "mean_t_statistic": t_mean,
    })

df = pd.DataFrame(rows)
df.to_csv("results/pmd_state_vectors.csv", index=False)
print(df.to_string(index=False))