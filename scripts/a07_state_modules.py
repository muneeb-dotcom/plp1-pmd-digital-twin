import scanpy as sc
import json

adata = sc.read_h5ad("results/ol_reference_atlas.h5ad")
adata = adata[adata.obs["cell_type"] != "PPR"].copy()  # exclude non-trajectory outlier

adata_full = adata.raw.to_adata()
adata_full.obs["cell_type"] = adata.obs["cell_type"]
adata_full.obs["maturation"] = adata.obs["maturation"]

MODULES = {
    "myelin_output": ["Mbp","Mog","Mag","Cnp","Mal","Plp1","Ugt8a","Mobp"],
    "er_stress":    ["Hspa5","Ddit3","Atf4","Atf6","Xbp1","Herpud1","Eif2ak3",
                     "Dnajb9","Edem1","Pdia4","Pdia6","Hyou1","Manf","Sel1l","Derl1"],
    "opc_identity": ["Pdgfra","Cspg4","Ptprz1","Sox10","Olig1","Olig2"],
    "apoptosis":    ["Casp3","Bax","Bcl2l11","Trp53","Cdkn1a"],
    "lipid_synth":  ["Fdft1","Hmgcr","Sqle","Idi1","Srebf2","Cyp51"],
}

module_genes_found = {}
for name, genes in MODULES.items():
    present = [g for g in genes if g in adata_full.var_names]
    module_genes_found[name] = present
    sc.tl.score_genes(adata_full, present, score_name=f"score_{name}")
    print(f"{name}: {len(present)}/{len(genes)} -> {present}")

mature_mask = adata_full.obs["cell_type"].isin(["MOL1","MOL2","MOL3","MOL4","MOL5","MOL6"])
baseline = {}
for name in MODULES:
    col = f"score_{name}"
    baseline[name] = {
        "mean": float(adata_full.obs.loc[mature_mask, col].mean()),
        "std": float(adata_full.obs.loc[mature_mask, col].std()),
    }

with open("results/state_modules.json", "w") as f:
    json.dump({"genes": module_genes_found, "healthy_baseline": baseline}, f, indent=2)

adata_full.write("results/ol_reference_atlas_scored.h5ad")

print("\nHealthy mature-OL baseline (mean/std):")
for k, v in baseline.items():
    print(f"  {k}: mean={v['mean']:.3f}, std={v['std']:.3f}")

print("\nModule scores by cell type:")
print(adata_full.obs.groupby("cell_type", observed=True)[[f"score_{m}" for m in MODULES]].mean())