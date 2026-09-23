import scanpy as sc
import pandas as pd

human = sc.read_h5ad("results/human_ingested.h5ad")

markers = ["PDGFRA", "CSPG4", "BMP4", "MOG", "MBP", "PLP1"]
present = [g for g in markers if g in human.var_names]
print("Markers found:", present)

human_full = human.raw.to_adata() if human.raw is not None else human

for gene in present:
    means = human.obs.groupby("celltype_human", observed=True).apply(
        lambda x, g=gene: human[x.index, g].X.mean()
    )
    print(f"\n{gene} by human's own labels:")
    print(means.sort_values())