import scanpy as sc
import numpy as np
import pandas as pd
import json

adata = sc.read_h5ad("results/ol_reference_atlas_scored.h5ad")

mature_mask = adata.obs["cell_type"].isin(["MOL1","MOL2","MOL3","MOL4","MOL5","MOL6"])

with open("results/state_modules.json") as f:
    ref = json.load(f)
baseline = ref["healthy_baseline"]

# extend baseline with maturation itself (not in the original 6 modules)
baseline["maturation"] = {
    "mean": float(adata.obs.loc[mature_mask, "maturation"].mean()),
    "std": float(adata.obs.loc[mature_mask, "maturation"].std()),
}
with open("results/state_modules.json", "w") as f:
    json.dump({"genes": ref["genes"], "healthy_baseline": baseline}, f, indent=2)

def zscore(value, module):
    b = baseline[module]
    return (value - b["mean"]) / b["std"]

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

def naive_myelination(row):
    """The assumption most of the field implicitly makes:
    mature + non-apoptotic -> myelinating. Ignores the myelin
    program and ER-stress state entirely."""
    mat_z = zscore(row["maturation"], "maturation")
    apop_z = zscore(row["score_apoptosis"], "apoptosis")
    return sigmoid(mat_z) * sigmoid(-apop_z)

def gated_myelination(row):
    """Survival is necessary but NOT sufficient (Elitt et al. 2018).
    Requires the myelin program to also be active and ER stress low."""
    survival = naive_myelination(row)
    myelin_z = zscore(row["score_myelin_output"], "myelin_output")
    stress_z = zscore(row["score_er_stress"], "er_stress")
    program = sigmoid(myelin_z) * sigmoid(-stress_z)
    return survival * program

adata.obs["naive_myelination"] = adata.obs.apply(naive_myelination, axis=1)
adata.obs["gated_myelination"] = adata.obs.apply(gated_myelination, axis=1)

order = ["OPC","COP","NFOL1","NFOL2","MFOL1","MFOL2","MOL1","MOL2","MOL3","MOL4","MOL5","MOL6"]
summary = adata.obs.groupby("cell_type", observed=True)[["naive_myelination","gated_myelination"]].mean()
summary = summary.reindex(order)

print("Naive vs gated myelination score across the healthy maturation trajectory:")
print(summary)
print("\n(Expected: both rise OPC->MOL and stay close together in HEALTHY cells --")
print(" they only diverge under a pathological state with survival but no myelin")
print(" program, e.g. Elitt et al. 2018's Ro 25-6981 rescue arm. That specific test")
print(" is reserved for Module C, which holds that dataset out for validation.)")

summary.to_csv("results/myelination_gate_healthy_check.csv")
adata.write("results/ol_reference_atlas_scored.h5ad")
print("\nSCRIPT COMPLETE")