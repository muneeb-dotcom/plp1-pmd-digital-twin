import scanpy as sc
import pandas as pd
import json

with open("results/state_modules.json") as f:
    ref = json.load(f)
MODULE_GENES = ref["genes"]

human = sc.read_h5ad("results/human_ingested.h5ad")  # already log-normalized, uppercase genes

order = ["OPCs", "COPs", "ImOlGs", "Oligo1", "Oligo2", "Oligo3", "Oligo4", "Oligo5", "Oligo6"]
human.obs["celltype_human"] = pd.Categorical(human.obs["celltype_human"], categories=order, ordered=True)

for name, genes in MODULE_GENES.items():
    genes_upper = [g.upper() for g in genes]
    present = [g for g in genes_upper if g in human.var_names]
    if present:
        sc.tl.score_genes(human, present, score_name=f"score_{name}")
        print(f"{name}: {len(present)}/{len(genes)} genes found in human data")

print("\nHuman module scores by the dataset's OWN cell type labels (OPC -> mature order):")
score_cols = [f"score_{m}" for m in MODULE_GENES if f"score_{m}" in human.obs.columns]
print(human.obs.groupby("celltype_human", observed=True)[score_cols].mean())

human.write("results/human_scored.h5ad")
print("\nSAVED. SCRIPT COMPLETE")