import matplotlib
matplotlib.use("Agg")

import scanpy as sc
import sys

sc.settings.verbosity = 2
sc.settings.figdir = "figures/"

print("Loading QC'd data...", flush=True)
adata = sc.read_h5ad("data/reference/marques_qc.h5ad")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, n_top_genes=2000, layer="counts", flavor="seurat_v3")
adata = adata[:, adata.var.highly_variable].copy()

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.8, key_added="leiden")

print("Writing h5ad BEFORE plotting...", flush=True)
adata.write("data/reference/marques_clustered.h5ad")
print("SAVED SUCCESSFULLY to data/reference/marques_clustered.h5ad", flush=True)

print("\nLeiden cluster counts:")
print(adata.obs["leiden"].value_counts())

print("\nNow attempting to save the figure (data is already safe either way)...", flush=True)
try:
    sc.pl.umap(adata, color=["cell_type", "leiden"], save="_atlas_overview.png")
    print("Figure saved.", flush=True)
except Exception as e:
    print(f"Plotting failed but data is safe: {e}", flush=True)

print("SCRIPT COMPLETE", flush=True)