import sys

path_hla_hits_tsv = sys.argv[1]

# Example line:
# A     0.0     HLA:HLA00001_A*01:01:01:01_3503_bp      MaSuRCA_Combo.phap#chr6  29784906  +  3503=
with open(path_hla_hits_tsv) as f:
    for line in f:
        gene, _, allele, contig, _, _, cigar = line.strip().split('\t')

        allele_len = allele.split('_bp')[0].split('_')[-1]
        match_len = cigar[:-1]

        if allele_len == match_len:
            print(line.strip())
