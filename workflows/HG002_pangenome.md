# HG002 pangenome

## Tools

```shell
guix install parallel

mkdir -p ~/tools $$ cd ~/tools

git clone --recursive https://github.com/ekg/fastix.git
cd fastix
git checkout 331c1159ea16625ee79d1a82522e800c99206834
cargo build --release
mv target/release/fastix target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834
cd ..

git clone --recursive https://github.com/ekg/wfmash.git
cd wfmash
git checkout 09e73eb3fcf24b8b7312b8890dd0741933f0d1cd
cmake -H. -Bbuild && cmake --build build -- -j 48
mv build/bin/wfmash build/bin/wfmash-09e73eb3fcf24b8b7312b8890dd0741933f0d1cd
cd ..

git clone --recursive https://github.com/ekg/seqwish.git
cd seqwish
git checkout ccfefb016fcfc9937817ce61dc06bbcf382be75e
cmake -H. -Bbuild && cmake --build build -- -j 48
mv bin/seqwish bin/seqwish-ccfefb016fcfc9937817ce61dc06bbcf382be75e
cd ..

git clone --recursive https://github.com/pangenome/smoothxg.git
cd smoothxg
git checkout d8c0fb8e5a5527936a4e550d6d1cd473a872313a
cmake -H. -Bbuild && cmake --build build -- -j 48
mv bin/smoothxg bin/smoothxg-d8c0fb8e5a5527936a4e550d6d1cd473a872313a
cd ..

git clone --recursive https://github.com/pangenome/odgi.git
cd odgi
git checkout 67a7e5bb2f328888e194845a362cef9c8ccc488f
mv bin/odgi bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f
cmake -H. -Bbuild && cmake --build build -- -j 48
cd ..

git clone --recursive https://github.com/pangenome/pggb.git
cd pggb
git checkout 5d2601127e8d08b39c8b05906240e5b50e46baf3
sed 's,"$fmt" wfmash,"$fmt" ~/tools/wfmash/build/bin/wfmash-09e73eb3fcf24b8b7312b8890dd0741933f0d1cd,g' pggb -i
sed 's,"$fmt" seqwish,"$fmt" ~/tools/seqwish/bin/seqwish-ccfefb016fcfc9937817ce61dc06bbcf382be75e,g' pggb -i
sed 's,"$fmt" smoothxg,"$fmt" ~/tools/smoothxg/bin/smoothxg-d8c0fb8e5a5527936a4e550d6d1cd473a872313a,g' pggb -i
sed 's,"$fmt" odgi,"$fmt" ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f,g' pggb -i
mv pggb pggb-5d2601127e8d08b39c8b05906240e5b50e46baf3
cd ..
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

Download the new HG002 reference assemblies:

```shell
# Url: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/21edcb42-02c4-4e9f-b226-6773e62484a4--RU-HG002-commons/assembly/curated_round2/
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/21edcb42-02c4-4e9f-b226-6773e62484a4--RU-HG002-commons/assembly/curated_round2/HG002.mat.cur.20211005.fasta.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/21edcb42-02c4-4e9f-b226-6773e62484a4--RU-HG002-commons/assembly/curated_round2/HG002.pat.cur.20211005.fasta.gz
```

Decompress:

```shell
ls HG002.*at.cur.20211005.fasta.gz | while read f; do gunzip $f; done
```

Rename:

```shell
( ~/tools/fastix/target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834 -p 'HG002_20211005.mat#' HG002.mat.cur.20211005.fasta | bgzip -@ 48 -c >HG002_mat.20211005.fa.gz ) &
( ~/tools/fastix/target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834 -p 'HG002_20211005.pat#' HG002.pat.cur.20211005.fasta | bgzip -@ 48 -c >HG002_pat.20211005.fa.gz ) &
```

Download and prepare the references:

```shell
mkdir -p /lizardfs/guarracino/HG002_assemblies_assessment/references/
cd /lizardfs/guarracino/HG002_assemblies_assessment/references/

wget -O- https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index | grep 'chm13\|h38' | awk '{ print $2 }' | sed 's%s3://human-pangenomics/working/%https://s3-us-west-2.amazonaws.com/human-pangenomics/working/%g' >refs.urls
cat refs.urls | parallel -j 2 'wget -q {} && echo got {}'
```

Add a prefix to the reference sequences:

```shell
( ~/tools/fastix/target/release/fastix -p 'grch38#' <(zcat GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz) >grch38_full.fa && samtools faidx grch38_full.fa ) &
( ~/tools/fastix/target/release/fastix -p 'chm13#' <(zcat chm13.draft_v1.1.fasta.gz) >chm13.fa && samtools faidx chm13.fa ) &
wait
```

Remove unplaced contigs from grch38 that are (hopefully) represented in chm13:

```shell
samtools faidx grch38_full.fa $(cat grch38_full.fa.fai | cut -f 1 | grep -v _ ) > grch38.fa && samtools faidx grch38.fa
```

Put all together:

```shell
cd /lizardfs/guarracino/HG002_assemblies_assessment/

cat references/grch38.fa references/chm13.fa <(zcat assemblies/HG002_20211005.mat.fa.gz) <(zcat assemblies/HG002_20211005.pat.fa.gz) | bgzip -@ 48 -c > HG002_20211005+refs.fa.gz
samtools faidx HG002_20211005+refs.fa.gz
```

## Build the graph with PGGB

```shell
sbatch -p lowmem -c 48 --wrap 'hostname; cd /scratch && ~/tools/pggb/pggb-5d2601127e8d08b39c8b05906240e5b50e46baf3 -i /lizardfs/guarracino/HG002_assemblies_assessment/HG002_20211005+refs.fa.gz -o HG002_20211005+refs -t 48 -p 98 -s 100000 -n 4 -k 311 -O 0.03 -T 48 -U -v -L -V chm13:#,grch38:# ; mv /scratch/HG002_20211005+refs /lizardfs/guarracino/HG002_assemblies_assessment/'
```

