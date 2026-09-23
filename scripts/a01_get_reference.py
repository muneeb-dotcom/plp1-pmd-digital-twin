import requests
import os
import pandas as pd

os.makedirs("data/reference", exist_ok=True)

url = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75330/suppl/"
       "GSE75330_Marques_et_al_mol_counts2.tab.gz")
gz = "data/reference/marques_counts.tab.gz"

if not os.path.exists(gz):
    print("Downloading Marques 2016 OL atlas...")
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(gz, "wb") as f:
            for chunk in r.iter_content(1024 * 1024):
                f.write(chunk)
    print("Done.")
else:
    print("Already downloaded.")

df = pd.read_csv(gz, sep="\t", index_col=0, compression="gzip")
print("Matrix shape (genes x cells):", df.shape)
print(df.iloc[:5, :3])