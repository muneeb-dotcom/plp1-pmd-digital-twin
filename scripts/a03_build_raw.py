import scanpy as sc
import anndata as ad
import pandas as pd

sc.settings.verbosity = 2

counts = pd.read_csv("data/reference/marques_counts.tab.gz",
                     sep="\t", index_col=0, compression="gzip")
labels = pd.read_csv("data/reference/cell_labels.csv")

# genes x cells -> cells x genes
adata = ad.AnnData(counts.T)
adata.var_names_make_unique()

# attach labels, keeping only cells present in both
labels = labels.set_index("sample_title")
common = adata.obs_names.intersection(labels.index)
print(f"Keeping {len(common)} of {adata.n_obs} cells with matched labels")

adata = adata[common].copy()
adata.obs["cell_type"] = labels.loc[common, "cell_type"].values
adata.layers["counts"] = adata.X.copy()

adata.write("data/reference/marques_raw.h5ad")
print(adata)
print("\nCell type counts:")
print(adata.obs["cell_type"].value_counts())