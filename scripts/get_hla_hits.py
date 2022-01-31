# From https://github.com/NCBI-Hackathons/TheHumanPangenome/blob/5055b4c01af69483709883b2f82fdf208e75d0ec/MHC/e2e_notebooks/X01_HLA-Typing-from-Assmbly.ipynb

import sys

path_hla_alignment_paf = sys.argv[1]

def get_allele_to_hits(fn):
    allele_to_hits = {}
    with open(fn)  as f:
        for row in f:
            if row[0] in ["#","@"]:
                continue
            row = row.strip().split("\t")
            HLA_allele = row[0].split("_")[1].split("*")[0]

            data = [c.split(":") for c in row[13:]]
            data = dict( ((c[0],c[2]) for c in data) )

            allele_to_hits.setdefault(HLA_allele, [])
            allele_to_hits[HLA_allele].append( (float(data["de"]), row[0], data["cg"], row[5], int(row[7]), row[4] ) )
    return allele_to_hits

for allele, hits in get_allele_to_hits(path_hla_alignment_paf).items():
    hits.sort(key=lambda x: x[0])
    for h in hits[:5]:
        print(allele, *[str(c) for c in h], sep="\t")
