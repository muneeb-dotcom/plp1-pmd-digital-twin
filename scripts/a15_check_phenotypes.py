import pandas as pd
import re

df = pd.read_csv("../plp1-pmd-project/data/clinvar/plp1_classified.csv")

print("PhenotypeList value counts (top 20):")
print(df["PhenotypeList"].value_counts().head(20))

print(f"\nNon-null phenotypes: {df['PhenotypeList'].notna().sum()} of {len(df)}")

# check copy-number extractability from duplication variant names
dup = df[df["mechanism"] == "duplication"]
dup_copy = dup["Name"].str.extract(r"x(\d+)\b")
print(f"\nDuplication variants with extractable copy number: {dup_copy[0].notna().sum()} of {len(dup)}")
print(dup_copy[0].value_counts())

# check truncation position extractability for LOF
lof = df[df["mechanism"] == "loss_of_function"]
print(f"\nLOF variants with a protein_pos value: {lof['protein_pos'].notna().sum()} of {len(lof)}")