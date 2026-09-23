import pandas as pd
import re
from scipy.stats import spearmanr

traits = pd.read_csv("data/pathogenic_traits.csv")

def severity_score(trait_str):
    t = str(trait_str)
    if "connatal" in t.lower():
        return 2
    if "pelizaeus-merzbacher disease" in t.lower():
        return 1
    if "hereditary spastic paraplegia 2" in t.lower():
        return 0
    return None

traits["severity"] = traits["traits"].apply(severity_score)
traits = traits.dropna(subset=["severity"])
print(f"{len(traits)} variants with a usable severity label")
print(traits["severity"].value_counts())

# extract HGVS protein change from title to join against Part 1 data
traits["Name"] = traits["title"]

classified = pd.read_csv("../plp1-pmd-project/data/clinvar/plp1_classified.csv")
merged = traits.merge(classified, on="Name", how="inner")
print(f"\n{len(merged)} matched to Part 1 classification")
print(merged["mechanism"].value_counts())

# --- Test 1: within misfolding class, does ddG correlate with severity? ---
ddg = pd.read_csv("../plp1-therapeutics/results/../../plp1-pmd-project/data/structures/plp1_structural_disruption_scores.csv") \
    if False else pd.read_csv("../plp1-pmd-project/data/structures/plp1_structural_disruption_scores.csv")
mis = merged[merged["mechanism"] == "misfolding_candidate"].merge(
    ddg[["Name", "ddG_kcal_mol"]], on="Name", how="inner")
print(f"\nMisfolding variants with both severity AND ddG: {len(mis)}")
if len(mis) >= 3:
    rho, p = spearmanr(mis["ddG_kcal_mol"], mis["severity"])
    print(f"Spearman rho={rho:.3f}, p={p:.4f}")
    print(mis[["Name", "ddG_kcal_mol", "severity"]].to_string(index=False))

# --- Test 2: within duplication class, does copy number correlate with severity? ---
dup = merged[merged["mechanism"] == "duplication"].copy()
dup["copy_number"] = dup["Name"].str.extract(r"x(\d+)\b").astype(float)
dup = dup.dropna(subset=["copy_number"])
print(f"\nDuplication variants with both severity AND copy number: {len(dup)}")
if len(dup) >= 3:
    rho, p = spearmanr(dup["copy_number"], dup["severity"])
    print(f"Spearman rho={rho:.3f}, p={p:.4f}")
    print(dup[["Name", "copy_number", "severity"]].to_string(index=False))

merged.to_csv("results/variant_severity_merged.csv", index=False)