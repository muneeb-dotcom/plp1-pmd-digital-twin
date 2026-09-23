import requests
import gzip
import shutil
import os

url = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75330/matrix/"
       "GSE75330_series_matrix.txt.gz")
gz_path = "data/reference/series_matrix.txt.gz"
out_path = "data/reference/series_matrix.txt"

if not os.path.exists(gz_path):
    print("Downloading series matrix...")
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    with open(gz_path, "wb") as f:
        f.write(r.content)

if not os.path.exists(out_path):
    with gzip.open(gz_path, "rb") as fin, open(out_path, "wb") as fout:
        shutil.copyfileobj(fin, fout)

print("Inspecting sample metadata lines...\n")
with open(out_path, encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!Sample_title") or line.startswith("!Sample_characteristics") or line.startswith("!Sample_geo_accession"):
            print(line[:250])
            print("---")