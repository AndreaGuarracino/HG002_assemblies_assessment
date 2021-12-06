# HG002_assemblies_assessment

Pangenomic assessment of HG002 assemblies

## Versions

# TODO###

## Preparation

Clone the repository:

```shell
cd /lizardfs/guarracino/
git clone --recursive https://github.com/AndreaGuarracino/HG002_assemblies_assessment.git
```

## Download data

Create the `assemblies` folder:

```shell
mkdir -p assemblies/
cd assemblies/
```

Download HG002 assemblies:

```shell
#Url: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=publications/HG002_BAKEOFF_2021/HG002_Assessment/assemblies/

grep s3 data/HGRC_bakeoff_HG002_assemblies_identifiers.tsv | 
  cut -f 6 | 
  sed 's,s3://human-pangenomics/HPRC/HG002_Assessment/assemblies/,https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/HG002_BAKEOFF_2021/HG002_Assessment/assemblies/,g' | 
  cut -f 1 -d ' ' |
  while read f; do
    wget $f;
  done

# Missing link
grep s3 -v data/HGRC_bakeoff_HG002_assemblies_identifiers.tsv
```

Download the new HG002 reference assemblies:

```shell
# Url: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/21edcb42-02c4-4e9f-b226-6773e62484a4--RU-HG002-commons/assembly/curated_round2/
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/21edcb42-02c4-4e9f-b226-6773e62484a4--RU-HG002-commons/assembly/curated_round2/HG002.mat.cur.20211005.fasta.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/21edcb42-02c4-4e9f-b226-6773e62484a4--RU-HG002-commons/assembly/curated_round2/HG002.pat.cur.20211005.fasta.gz
```

