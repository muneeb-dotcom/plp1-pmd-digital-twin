from Bio import Entrez
import json

Entrez.email = "shahzadmuneeb15@gmail.com"

search = Entrez.esearch(db="clinvar", term="PLP1[gene]", retmax=5)
record = Entrez.read(search)
ids = record["IdList"]

handle = Entrez.esummary(db="clinvar", id=",".join(ids), retmode="json")
data = json.load(handle)

for uid in ids:
    rec = data["result"][uid]
    print(f"UID {uid}:")
    print(json.dumps(rec.get("trait_set"), indent=2)[:500])
    print(json.dumps(rec.get("germline_classification"), indent=2)[:500])
    print("---")