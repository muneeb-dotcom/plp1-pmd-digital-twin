from Bio import Entrez
import json
import time
import pandas as pd
from collections import Counter

Entrez.email = "shahzadmuneeb15@gmail.com"

search = Entrez.esearch(db="clinvar", term="PLP1[gene]", retmax=5000)
record = Entrez.read(search)
ids = record["IdList"]

pathogenic_labels = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
trait_counter = Counter()
rows = []

batch_size = 50
for i in range(0, len(ids), batch_size):
    batch = ids[i:i + batch_size]
    handle = Entrez.esummary(db="clinvar", id=",".join(batch), retmode="json")
    data = json.load(handle)
    for uid in data["result"]["uids"]:
        rec = data["result"][uid]
        gc = rec.get("germline_classification") or {}
        desc = gc.get("description")
        if desc not in pathogenic_labels:
            continue
        traits = gc.get("trait_set") or []
        names = [t.get("trait_name") for t in traits]
        for n in names:
            trait_counter[n] += 1
        rows.append({"uid": uid, "title": rec.get("title"), "traits": "; ".join(names)})
    time.sleep(0.3)

print(f"{len(rows)} pathogenic records checked")
print("\nAll distinct trait_name values and counts:")
for name, count in trait_counter.most_common():
    print(f"  {count}: {name}")

pd.DataFrame(rows).to_csv("data/pathogenic_traits.csv", index=False)