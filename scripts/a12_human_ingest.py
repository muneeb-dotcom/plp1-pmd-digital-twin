import matplotlib
matplotlib.use("Agg")
import scanpy as sc
import pandas as pd
import numpy as np

OL_TYPES = ["OPCs", "COPs", "ImOlGs", "Oligo1", "Oligo2", "Oligo3", "Oligo4", "Oligo5", "Oligo6"]

print("Loading human annotation...", flush=True)
ann = pd.read_csv("data/reference/human/GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz",
                  sep="\t", compression="gzip")
keep = ann[(ann["Celltypes"].isin(OL_TYPES)) & (ann["Condition"] == "Ctrl")]
print(f"Keeping {len(keep)} of {len(ann)} cells (OL-lineage, Ctrl only)")

print("Loading human expression matrix (this may take a minute)...", flush=True)
expr = pd.read_csv("data/reference/human/GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt.gz",
                   sep="\t", compression="gzip", index_col=0)
expr = expr[keep["Detected"].tolist()]  # subset to kept cells
print("Human matrix subset:", expr.shape)

human = sc.AnnData(expr.T)
human.var_names_make_unique()
human.obs["celltype_human"] = keep.set_index("Detected").loc[human.obs_names, "Celltypes"].values
human.var_names = human.var_names.str.upper()

print("Loading mouse reference...", flush=True)
mouse = sc.read_h5ad("results/ol_reference_atlas_scored.h5ad")
mouse.var_names = mouse.var_names.str.upper()

shared = mouse.var_names.intersection(human.var_names)
print(f"{len(shared)} shared genes between mouse and human")

mouse_shared = mouse[:, shared].copy()
human_shared = human[:, shared].copy()

sc.pp.normalize_total(human_shared, target_sum=1e4)
sc.pp.log1p(human_shared)

sc.pp.scale(mouse_shared, max_value=10)
sc.tl.pca(mouse_shared, n_comps=50)
sc.pp.neighbors(mouse_shared, n_neighbors=15, n_pcs=30)
sc.tl.umap(mouse_shared)

sc.tl.ingest(human_shared, mouse_shared, obs="cell_type")

print("\nConfusion: human's own labels vs mouse-trajectory-predicted labels")
crosstab = pd.crosstab(human_shared.obs["celltype_human"], human_shared.obs["cell_type"])
print(crosstab)

crosstab.to_csv("results/human_mouse_concordance.csv")
human_shared.write("results/human_ingested.h5ad")
print("\nSAVED. SCRIPT COMPLETE", flush=True)