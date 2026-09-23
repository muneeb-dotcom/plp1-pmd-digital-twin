import pandas as pd

path = "data/reference/series_matrix.txt"

titles = None
cell_types = None

with open(path, encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!Sample_title"):
            titles = [x.strip('"') for x in line.strip().split("\t")[1:]]
        if line.startswith("!Sample_characteristics_ch1") and "inferred cell type" in line:
            cell_types = [x.strip('"').replace("inferred cell type: ", "") for x in line.strip().split("\t")[1:]]

assert titles is not None and cell_types is not None, "Failed to find one of the fields"
assert len(titles) == len(cell_types), "Mismatched lengths"

labels = pd.DataFrame({"sample_title": titles, "cell_type": cell_types})
labels.to_csv("data/reference/cell_labels.csv", index=False)

print(f"{len(labels)} labeled samples")
print("\nCell type counts:")
print(labels["cell_type"].value_counts())

# sanity check against the counts matrix columns
counts_cols = pd.read_csv("data/reference/marques_counts.tab.gz", sep="\t",
                          index_col=0, compression="gzip", nrows=0).columns
matched = labels["sample_title"].isin(counts_cols).sum()
print(f"\n{matched} of {len(labels)} labels match a column in the counts matrix")