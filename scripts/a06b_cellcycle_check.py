import scanpy as sc

adata = sc.read_h5ad("results/ol_reference_atlas.h5ad")

cc_genes = ["Mki67", "Top2a", "Ccnb1", "Pcna", "Mcm2"]
present = [g in adata.raw.var_names for g in cc_genes]
print(dict(zip(cc_genes, present)))

adata_full = adata.raw.to_adata()
adata_full.obs["cell_type"] = adata.obs["cell_type"]

found = [g for g in cc_genes if g in adata_full.var_names]
sc.tl.score_genes(adata_full, found, score_name="cc_score")

print(adata_full.obs.groupby("cell_type", observed=True)["cc_score"].mean().sort_values())