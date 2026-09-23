import scanpy as sc
import numpy as np
from scipy.stats import median_abs_deviation

adata = sc.read_h5ad("data/reference/marques_raw.h5ad")

# --- calculate QC metrics (mouse mito prefix is lowercase 'mt-') ---
adata.var["mt"] = adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

print("Before filtering:")
print(f"  Cells: {adata.n_obs}, Genes: {adata.n_vars}")
print(f"  Median counts/cell: {adata.obs['total_counts'].median():.0f}")
print(f"  Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
print(f"  Mean mito%: {adata.obs['pct_counts_mt'].mean():.2f}%")

def mad_outlier(values, n_mads):
    median = np.median(values)
    mad = median_abs_deviation(values)
    lower, upper = median - n_mads * mad, median + n_mads * mad
    return (values < lower) | (values > upper), lower, upper

# scverse-standard: 5 MADs on counts/genes (permissive), 3 MADs on mito% (stricter)
out_counts, lo_c, hi_c = mad_outlier(adata.obs["total_counts"], 5)
out_genes, lo_g, hi_g = mad_outlier(adata.obs["n_genes_by_counts"], 5)
out_mito, lo_m, hi_m = mad_outlier(adata.obs["pct_counts_mt"], 3)

print(f"\ntotal_counts outliers: {out_counts.sum()} (bounds {lo_c:.0f}-{hi_c:.0f})")
print(f"n_genes outliers: {out_genes.sum()} (bounds {lo_g:.0f}-{hi_g:.0f})")
print(f"pct_mt outliers: {out_mito.sum()} (bounds {lo_m:.1f}-{hi_m:.1f})")

is_outlier = out_counts | out_genes | out_mito
print(f"\nTotal outliers (any criterion): {is_outlier.sum()} of {adata.n_obs} ({is_outlier.sum()/adata.n_obs*100:.1f}%)")

adata_filtered = adata[~is_outlier].copy()
sc.pp.filter_genes(adata_filtered, min_cells=20)

print(f"\nAfter filtering: {adata_filtered.n_obs} cells, {adata_filtered.n_vars} genes")
print("\nCell type counts after filtering:")
print(adata_filtered.obs["cell_type"].value_counts())

adata_filtered.write("data/reference/marques_qc.h5ad")