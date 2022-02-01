import sys

path_hla_perfect_hits_tsv = sys.argv[1]
path_haplotypes = sys.argv[2]
gene_as_rows = sys.argv[3]

# Kepp all haplotypes and preserve their order in the list
haplotype_list = list()

with open(path_haplotypes) as f:
    for line in f:
        haplotype_list.append(line.strip())

# Read all possible genes
gene_set = set()
with open(path_hla_perfect_hits_tsv) as f:
    f.readline() # Skip header

    for line in f:
        Gene = line.strip().split('\t')[1]
        gene_set.add(Gene)


haplotype_2_gene_2_allele_contig_dict = {}
for haplotype in haplotype_list:
    haplotype_2_gene_2_allele_contig_dict[haplotype] = {x: [] for x in sorted(gene_set)}

gene_2_haplotype_2_allele_contig_dict = {}

# Parse file
with open(path_hla_perfect_hits_tsv) as f:
    f.readline() # Skip header

    for line in f:
        AbbreviatedName, Gene, Allele, Contig = line.strip().split('\t')

        if Gene not in gene_2_haplotype_2_allele_contig_dict:
            gene_2_haplotype_2_allele_contig_dict[Gene] = {x: [] for x in haplotype_list}
        gene_2_haplotype_2_allele_contig_dict[Gene][AbbreviatedName].append(f'{Allele}({Contig})')

        haplotype_2_gene_2_allele_contig_dict[AbbreviatedName][Gene].append(f'{Allele}({Contig})')


if gene_as_rows == 'Y':
    print('Gene', '\t'.join([x for x in haplotype_list]), sep="\t")
    for gene, haplotype_2_allele_contig_dict in gene_2_haplotype_2_allele_contig_dict.items():
        print(gene, '\t'.join([';'.join(x) for x in haplotype_2_allele_contig_dict.values()]), sep="\t")
else:
    print('AbbreviatedName', '\t'.join([x for x in sorted(gene_set)]), sep="\t")
    for haplotype, gene_2_allele_contig_dict in haplotype_2_gene_2_allele_contig_dict.items():
        print(haplotype, '\t'.join([';'.join(x) for x in gene_2_allele_contig_dict.values()]), sep="\t")
