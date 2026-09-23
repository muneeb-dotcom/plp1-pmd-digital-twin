import pandas as pd

print("=== Annotation table ===")
ann = pd.read_csv("data/reference/human/GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz",
                  sep="\t", compression="gzip")
print(ann.shape)
print(ann.head())
print("\nColumns:", ann.columns.tolist())

print("\n=== Expression matrix (header + few rows only) ===")
expr_head = pd.read_csv("data/reference/human/GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt.gz",
                        sep="\t", compression="gzip", nrows=5)
print(expr_head.shape)
print(expr_head.iloc[:, :5])