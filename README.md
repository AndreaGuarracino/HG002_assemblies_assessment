# HG002_assemblies_assessment

Pangenomic assessment of HG002 assemblies

## Tools

```shell
# guix install parallel

mkdir -p ~/tools && cd ~/tools

git clone --recursive https://github.com/ekg/fastix.git
cd fastix
git checkout 331c1159ea16625ee79d1a82522e800c99206834
cargo build --release
mv target/release/fastix target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834
cd ..
```

## Preparation

Clone the repository:

```shell
cd /lizardfs/guarracino/
git clone --recursive https://github.com/AndreaGuarracino/HG002_assemblies_assessment.git
```

## Download and prepare the data

Create the `assemblies` folder:

```shell
mkdir -p /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/
cd /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/
```

Download the HG002 assemblies:

```shell
# Url: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=publications/HG002_BAKEOFF_2021/HG002_Assessment/assemblies/

grep s3 ../data/HGRC_bakeoff_HG002_assemblies_identifiers.tsv | 
  cut -f 6 | 
  sed 's,s3://human-pangenomics/HPRC/HG002_Assessment/assemblies/,https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/HG002_BAKEOFF_2021/HG002_Assessment/assemblies/,g' | 
  cut -f 1 -d ' ' |
  while read f; do
    wget $f;
  done

# Fix collisions
mv asm.fa.gz.1 asm.v2.fa.gz
mv flye.scaffolds.fasta.gz paternal.ONT.std.flye.scaffolds.fasta.gz
mv flye.scaffolds.fasta.gz.1 maternal.ONT.std.flye.scaffolds.fasta.gz
mv flye.scaffolds.fasta.gz.2 paternal.ONT.UL.flye.scaffolds.fasta.gz
mv flye.scaffolds.fasta.gz.3 maternal.ONT.UL.flye.scaffolds.fasta.gz
mv hg002_crossstitch_upload.tar.gz hg002_crossstitch_upload.hap1.tar.gz
mv hg002_crossstitch_upload.tar.gz.1 hg002_crossstitch_upload.hap2.tar.gz
mv canu.contigs.fasta.gz maternal.canu.contigs.fasta.gz
mv canu.contigs.fasta.gz.1 paternal.canu.contigs.fasta.gz
mv peregrine.contigs.fasta.gz maternal.peregrine.contigs.fasta.gz
mv peregrine.contigs.fasta.gz.1 paternal.peregrine.contigs.fasta.gz

# Missing link in the file
grep s3 -v ../data/HGRC_bakeoff_HG002_assemblies_identifiers.tsv
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/HG002_BAKEOFF_2021/HG002_Assessment/assemblies/Dovetail_Genomics/new/Dovetail_HG002_phase1_scaffolds_with_X_Y.fa.gz
```

Download the new HG002 reference assemblies:

```shell
# Url: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/21edcb42-02c4-4e9f-b226-6773e62484a4--RU-HG002-commons/assembly/curated_round2/

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/21edcb42-02c4-4e9f-b226-6773e62484a4--RU-HG002-commons/assembly/curated_round2/HG002.mat.cur.20211005.fasta.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/21edcb42-02c4-4e9f-b226-6773e62484a4--RU-HG002-commons/assembly/curated_round2/HG002.pat.cur.20211005.fasta.gz
```

Decompress:

```shell
ls *a.gz | while read f; do gunzip $f; done
ls *tar.gz | while read f; do tar -xvzf $f; done && rm *tar.gz
mv hg002_crossstitch_upload/hap1.fa hg002_crossstitch_upload_hap1.fa
mv hg002_crossstitch_upload/hap2.fa hg002_crossstitch_upload_hap2.fa
rm -rf hg002_crossstitch_upload/
```

Add prefixes:

```shell
# Check existence
#cat ../data/HGRC_bakeoff_HG002_assemblies_v3_renaming.tsv | sed 1,1d | sed 's/"//g' | while read -r Id Filename AbbreviatedName; do ls -l $(echo $Filename); done

cat ../data/HGRC_bakeoff_HG002_assemblies_v3_renaming.tsv | sed 1,1d | sed 's/"//g' | while read -r Id Filename AbbreviatedName; do
  AbbreviatedName2=$(echo $AbbreviatedName | sed 's/ /_/g');
  echo "$Filename -> ${AbbreviatedName2}.fa.gz"
  ~/tools/fastix/target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834 -p "${AbbreviatedName2}#" $Filename | bgzip -@ 48 -c > ${AbbreviatedName2}.fa.gz;
done
```

## Pangenomic assessment of HG002 assemblies

[Link to the `HG002_pangenomic_assessment.md` workflow](workflows/HG002_pangenomic_assessment.md).

## HG002 pangenome graphs of selected _loci_

[Link to the `HG002_pangenome_graphs.md` workflow](workflows/HG002_pangenome_graphs.md).
