# HG002_assemblies_assessment

Pangenomic assessment of HG002 assemblies

## Tools

# TODO###

```shell
mkdir -p ~/tools $$ cd ~/tools

git clone --recursive https://github.com/ekg/fastix.git
cd fastix
git checkout 331c1159ea16625ee79d1a82522e800c99206834
cargo build --release

git clone --recursive https://github.com/ekg/wfmash.git
cd wfmash
git checkout 0a6c3b4ed296c8ef48a7fbd9d8f8af7de00ec0fb
cmake -H. -Bbuild && cmake --build build -- -j 48

git clone --recursive https://github.com/ekg/seqwish.git
cd seqwish
git checkout beb7a3805c5af599b8425db94b9524eecbfc73e0 # speedwish branch
cmake -H. -Bbuild && cmake --build build -- -j 48
```

## Preparation

Clone the repository:

```shell
cd /lizardfs/guarracino/
git clone --recursive https://github.com/AndreaGuarracino/HG002_assemblies_assessment.git
```

## Download and prepare data

Create the `assemblies` folder:

```shell
mkdir -p /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/
cd /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/
```

Download HG002 assemblies:

```shell
#Url: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=publications/HG002_BAKEOFF_2021/HG002_Assessment/assemblies/

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
# Check existance
#cat ../data/HGRC_bakeoff_HG002_assemblies_v3_renaming.tsv | sed 1,1d | sed 's/"//g' | while read -r Id Filename AbbreviatedName; do ls -l $(echo $Filename); done

cat ../data/HGRC_bakeoff_HG002_assemblies_v3_renaming.tsv | sed 1,1d | sed 's/"//g' | while read -r Id Filename AbbreviatedName; do
  AbbreviatedName2=$(echo $AbbreviatedName | sed 's/ /_/g');
  echo "$Filename -> ${AbbreviatedName2}.fa.gz"
  ~/tools/fastix/target/release/fastix -p "${AbbreviatedName2}#" $Filename | bgzip -@ 48 -c > ${AbbreviatedName2}.fa.gz;
done
```

Merge all assemblies:

```shell
# To merge assemblies sorted by ID
FA_GZS=$(cat ../data/HGRC_bakeoff_HG002_assemblies_v3_renaming.tsv | sed 1,1d | sed 's/"//g' | while read -r Id Filename AbbreviatedName; do AbbreviatedName2=$(echo $AbbreviatedName | sed 's/ /_/g'); echo ${AbbreviatedName2}.fa.gz; done | tr '\n' ' ')

zcat $FA_GZS | bgzip -@ 48 -c > HG002_all.fa.gz
samtools faidx HG002_all.fa.gz
```

## Mapping and Alignment

Generate all-vs-all mappings:

```shell
mkdir -p /lizardfs/guarracino/HG002_assemblies_assessment/alignment/

sbatch -p lowmem -c 48 --wrap 'cd /scratch && \time -v ~/tools/wfmash/build/bin/wfmash -X -s 20k -l 60k -p 95 -n 45 -k 16 -t 48 /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/HG002_all.fa.gz /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/HG002_all.fa.gz -m > HG002_all.s20k.l60k.p95.n45.k16.approx.paf && mv HG002_all.s20k.l60k.p95.n45.k16.approx.paf /lizardfs/guarracino/HG002_assemblies_assessment/alignment/'
```

Split the mappings in chunks:

```shell
cd /lizardfs/guarracino/alignment/
python3 ~/tools/wfmash/scripts/split_approx_mappings_in_chunks.py /lizardfs/guarracino/alignment/HG002_all.s20k.l60k.p95.n45.k16.approx.paf 5
```

Run the alignments on multiple nodes:

```shell
seq 0 4 | while read i; do sbatch -p lowmem -c 48 --wrap 'cd /scratch && ~/tools/wfmash/build/bin/wfmash -X -s 20k -l 60k -p 95 -n 45 -k 16 -t 48 /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/HG002_all.fa.gz /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/HG002_all.fa.gz -i /lizardfs/guarracino/HG002_assemblies_assessment/alignment/HG002_all.s20k.l60k.p95.n45.k16.approx.paf.chunk_'$i'.paf | pigz -c > HG002_all.s20k.l60k.p95.n45.k16.approx.paf.chunk_'$i'.paf.gz && mv HG002_all.s20k.l60k.p95.n45.k16.approx.paf.chunk_'$i'.paf.gz /lizardfs/guarracino/HG002_assemblies_assessment/alignment/'; done
```

## Inducing the graph

```shell
mkdir -p /lizardfs/guarracino/graphs/

# list all base-level alignment PAFs
PAFS=$(ls HG002_all.s20k.l60k.p95.n45.k16.approx.paf.chunk_*.paf.gz | tr '\n' ',')
PAFS=${input::-1}

seqwish -s reference.fa -p PAFS -k 47 -g seqwish.gfa -B 10000000 -P

sbatch -p lowmem -c 48 --wrap 'cd /scratch && \time -v ~/tools/seqwish/bin/seqwish -t 48 -s /lizardfs/guarracino/assemblies/HG002_all.fa.gz -p '$PAFS' -k 79 -B 50M -g HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.gfa -P && mv HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.gfa /lizardfs/guarracino/graphs/'
```
