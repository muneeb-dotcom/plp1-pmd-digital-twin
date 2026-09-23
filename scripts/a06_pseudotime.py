import matplotlib
matplotlib.use("Agg")

import scanpy as sc
import numpy as np

sc.settings.figdir = "figures/"

print("Loading clustered atlas...", flush=True)
adata = sc.read_h5ad("data/reference/marques_clustered.h5ad")

# root the trajectory at an OPC cell (the maturation starting point)
opc_idx = np.flatnonzero(adata.obs["cell_type"] == "OPC")
adata.uns["iroot"] = opc_idx[0]

sc.tl.diffmap(adata)
sc.tl.dpt(adata)
adata.obs["maturation"] = adata.obs["dpt_pseudotime"]

print("Writing h5ad...", flush=True)
adata.write("results/ol_reference_atlas.h5ad")
print("SAVED to results/ol_reference_atlas.h5ad", flush=True)

print("\nMean maturation pseudotime by cell type (should rise OPC -> MOL):")
print(adata.obs.groupby("cell_type", observed=True)["maturation"].mean().sort_values())

try:
    sc.pl.umap(adata, color=["maturation", "Pdgfra", "Mbp"], save="_maturation_axis.png")
    print("Figure saved.")
except Exception as e:
    print(f"Plotting failed but data is safe: {e}")

print("SCRIPT COMPLETE", flush=True)