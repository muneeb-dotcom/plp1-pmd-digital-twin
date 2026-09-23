import pandas as pd
raw = pd.read_csv("../plp1-pmd-project/data/clinvar/plp1_raw.csv")
print(f"Non-null PhenotypeList in RAW pull: {raw['PhenotypeList'].notna().sum()} of {len(raw)}")
print(raw["PhenotypeList"].dropna().unique()[:10])