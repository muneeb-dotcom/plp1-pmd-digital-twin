import matplotlib
matplotlib.use("Agg")

import scanpy as sc

adata = sc.read_h5ad("data/reference/marques_clustered.h5ad")
adata_full = adata.raw.to_adata()  # full gene set, log-normalized, pre-HVG-subset
adata_full.obs["cell_type"] = adata.obs["cell_type"]

markers = {
    "OPC":       ["Pdgfra", "Cspg4", "Ptprz1"],
    "COP":       ["Bmp4", "Sox6", "Neu4"],
    "NFOL":      ["Tcf7l2", "Itpr2", "Cnksr3"],
    "MFOL":      ["Ctps", "Mal", "Opalin"],
    "MOL":       ["Mog", "Mbp", "Plp1", "Klk6", "Apod"],
}

present = {}
for group, genes in markers.items():
    found = [g for g in genes if g in adata_full.var_names]
    missing = [g for g in genes if g not in adata_full.var_names]
    present[group] = found
    if missing:
        print(f"{group}: missing {missing}")

order = ["OPC", "PPR", "COP", "NFOL1", "NFOL2", "MFOL1", "MFOL2",
         "MOL1", "MOL2", "MOL3", "MOL4", "MOL5", "MOL6"]
adata_full.obs["cell_type"] = adata_full.obs["cell_type"].astype("category")
adata_full.obs["cell_type"] = adata_full.obs["cell_type"].cat.set_categories(order)

sc.settings.figdir = "figures/"
sc.pl.dotplot(adata_full, present, groupby="cell_type", save="_marker_validation.png")
print("\nSaved figures/dotplot_marker_validation.png")

print("\nMean expression check (should rise OPC -> MOL for Mbp/Mog/Plp1, fall for Pdgfra):")
for gene in ["Pdgfra", "Mbp", "Mog", "Plp1"]:
    means = adata_full.obs.groupby("cell_type", observed=True).apply(
        lambda x: adata_full[x.index, gene].X.mean() if gene in adata_full.var_names else None
    )
    print(f"\n{gene}:")
    print(means)