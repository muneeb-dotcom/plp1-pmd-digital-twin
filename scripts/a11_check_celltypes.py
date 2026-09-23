import pandas as pd

ann = pd.read_csv("data/reference/human/GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz",
                  sep="\t", compression="gzip")

print("Unique Celltypes:")
print(ann["Celltypes"].value_counts())
print("\nUnique Condition values:")
print(ann["Condition"].value_counts())