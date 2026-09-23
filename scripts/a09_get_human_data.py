import requests
import os

url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118257/suppl/"
headers = {"User-Agent": "Mozilla/5.0"}

files = [
    "GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt.gz",
    "GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz",
]

os.makedirs("data/reference/human", exist_ok=True)

def get_remote_size(file_url):
    r = requests.head(file_url, headers=headers, timeout=30)
    return int(r.headers.get("content-length", 0))

for f in files:
    out_path = f"data/reference/human/{f}"
    file_url = url + f
    total_size = get_remote_size(file_url)

    for attempt in range(1, 11):
        existing = os.path.getsize(out_path) if os.path.exists(out_path) else 0
        if existing >= total_size:
            break
        print(f"{f}: attempt {attempt}, resuming from {existing/(1024*1024):.0f} MB / {total_size/(1024*1024):.0f} MB")
        range_headers = {**headers, "Range": f"bytes={existing}-"}
        try:
            with requests.get(file_url, headers=range_headers, stream=True, timeout=60) as r:
                mode = "ab" if existing > 0 else "wb"
                with open(out_path, mode) as fout:
                    for chunk in r.iter_content(256 * 1024):
                        fout.write(chunk)
        except Exception as e:
            print(f"Dropped: {e}, retrying...")
            continue

    final = os.path.getsize(out_path)
    print(f"{f}: final size {final/(1024*1024):.1f} MB / {total_size/(1024*1024):.1f} MB {'OK' if final >= total_size else 'INCOMPLETE'}")

print("Done.")