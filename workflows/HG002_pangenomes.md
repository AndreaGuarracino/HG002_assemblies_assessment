# HG002 pangenomes

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

## Download and prepare the data

Create the `assemblies` folder:

```shell
mkdir -p /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/
cd /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/
```


Download some HG002 assemblies:

```shell
# Url: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=publications/HG002_BAKEOFF_2021/HG002_Assessment/assemblies/

grep s3 ../data/HGRC_bakeoff_HG002_assemblies_identifiers.tsv | grep 'asm6\|asm9\|asm23' | 
  cut -f 6 | 
  sed 's,s3://human-pangenomics/HPRC/HG002_Assessment/assemblies/,https://s3-us-west-2.amazonaws.com/human-pangenomics/publications/HG002_BAKEOFF_2021/HG002_Assessment/assemblies/,g' | 
  cut -f 1 -d ' ' |
  while read f; do
    wget $f;
  done

# Fix collisions
mv flye.scaffolds.fasta.gz paternal.ONT.std.flye.scaffolds.fasta.gz
mv flye.scaffolds.fasta.gz.1 maternal.ONT.std.flye.scaffolds.fasta.gz
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
```

Add prefixes:

```shell
# Check existence
#cat ../data/HGRC_bakeoff_HG002_assemblies_v3_renaming.tsv | sed 1,1d | sed 's/"//g' | grep 'asm6\|asm9\|asm23\|?' | while read -r Id Filename AbbreviatedName; do ls -l $(echo $Filename); done

cat ../data/HGRC_bakeoff_HG002_assemblies_v3_renaming.tsv | sed 1,1d | sed 's/"//g' | grep 'asm6\|asm9\|asm23\|?' | while read -r Id Filename AbbreviatedName; do
  AbbreviatedName2=$(echo $AbbreviatedName | sed 's/ /_/g');
  echo "$Filename -> ${AbbreviatedName2}.fa.gz"
  ~/tools/fastix/target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834 -p "${AbbreviatedName2}#" $Filename | bgzip -@ 48 -c > ${AbbreviatedName2}.fa.gz;
done
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
# GRCH38 + CHM13 + HG002
cat ../references/grch38.fa ../references/chm13.fa <(zcat HG002_20211005.mat.fa.gz) <(zcat HG002_20211005.pat.fa.gz) | bgzip -@ 48 -c > HG002_20211005+refs.fa.gz
samtools faidx HG002_20211005+refs.fa.gz
```


## Partitioning

Map `asm6 + asm9 + asm23 + HG002` assembly against the scaffolded references:

```shell
zcat Trio_Flye_ONT_std.pat.fa.gz Trio_Flye_ONT_std.mat.fa.gz \
  Trio_HiFiasm.mat.fa.gz Trio_HiFiasm.pat.fa.gz \
  Trio_VGP.mat.fa.gz Trio_VGP.pat.fa.gz \
  HG002_20211005.mat.fa.gz HG002_20211005.pat.fa.gz | bgzip -@ 48 -c > asm6+asm9+asm23+HG002_20211005.fa.gz
samtools faidx asm6+asm9+asm23+HG002_20211005.fa.gz

sbatch -p headnode -c 48 --wrap 'cd /scratch && \time -v ~/tools/wfmash/build/bin/wfmash-09e73eb3fcf24b8b7312b8890dd0741933f0d1cd -s 50k -l 150k -p 90 -k 16 -N -t 48 /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13+grch38.fa /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/asm6+asm9+asm23+HG002_20211005.fa.gz -m > /lizardfs/guarracino/HG002_assemblies_assessment/alignment/asm6+asm9+asm23+HG002_20211005.vs.chm13+grch38.s50k.l150k.p90.k16.N.approx.paf'
```

Collect unmapped contigs and remap them in split mode:

```shell
out=/lizardfs/guarracino/HG002_assemblies_assessment/alignment/asm6+asm9+asm23+HG002_20211005.vs.chm13+grch38.s50k.l150k.p90.k16.N.unmapped.names
comm -23 \
  <(cut -f 1 asm6+asm9+asm23+HG002_20211005.fa.gz.fai | sort) \
  <(cut -f 1 /lizardfs/guarracino/HG002_assemblies_assessment/alignment/asm6+asm9+asm23+HG002_20211005.vs.chm13+grch38.s50k.l150k.p90.k16.N.approx.paf | sort) > $out.txt

samtools faidx asm6+asm9+asm23+HG002_20211005.fa.gz $(tr '\n' ' ' <$out.txt) | bgzip -@ 48 -c >$out.fa.gz
samtools faidx $out.fa.gz

~/tools/wfmash/build/bin/wfmash-09e73eb3fcf24b8b7312b8890dd0741933f0d1cd -s 50k -l 150k -p 80 -k 16 -t 48 /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13+grch38.fa $out.fa.gz -m > $out.fa.gz.split.p80.approx.paf
```

Collect our best mapping for our attempted split rescues:

```shell

```

Take chr6's contigs:

```shell
cat \
  <(awk -v chm13="chm13#chr6" -v grch38="grch38#chr6" '($6 == chm13 || $6 == grch38) && $4 >= 100000' /lizardfs/guarracino/HG002_assemblies_assessment/alignment/asm6+asm9+asm23+HG002_20211005.vs.chm13+grch38.s50k.l150k.p90.k16.N.approx.paf | cut -f 1) \
   <(grep chr6 $out.fa.gz.split.p80.approx.paf | awk '{ print $1,$11,$0 }' | tr ' ' '\t' |  sort -n -r -k 1,2 | awk '$1 != last { print; last = $1; }' | awk '{ if($2 >= 100000) print $1}') \
   > /lizardfs/guarracino/HG002_assemblies_assessment/alignment/asm6+asm9+asm23+HG002_20211005.vs.chm13+grch38.s50k.l150k.p90.k16.N.chr6.100kbps.names.txt

cat \
  <(samtools faidx ../references/grch38.fa $(echo grch38#chr6)) \
  <(samtools faidx ../references/chm13.fa $(echo chm13#chr6)) \
  <(samtools faidx asm6+asm9+asm23+HG002_20211005.fa.gz $(cat /lizardfs/guarracino/HG002_assemblies_assessment/alignment/asm6+asm9+asm23+HG002_20211005.vs.chm13+grch38.s50k.l150k.p90.k16.N.chr6.100kbps.names.txt | sort)) | \
  bgzip -@ 48 -c > refs+asm6+asm9+asm23+HG002_20211005.chr6.100kbps.fa.gz
samtools faidx refs+asm6+asm9+asm23+HG002_20211005.chr6.100kbps.fa.gz

```


## Build the PGGB graphs

```shell
# GRCH38 + CHM13 + HG002
sbatch -p lowmem -c 48 --wrap 'hostname; cd /scratch && ~/tools/pggb/pggb-5d2601127e8d08b39c8b05906240e5b50e46baf3 -i /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/HG002_20211005+refs.fa.gz -o HG002_20211005+refs -t 48 -p 98 -s 100000 -n 4 -k 311 -O 0.03 -T 48 -U -v -L -V chm13:#,grch38:# ; mv /scratch/HG002_20211005+refs /lizardfs/guarracino/HG002_assemblies_assessment/'


# GRCH38 + CHM13 + asm6 + asm9 + asm23 + HG002

## contigs >= 1Mbps
sbatch -p lowmem -c 48 --wrap 'hostname; cd /scratch && ~/tools/pggb/pggb-5d2601127e8d08b39c8b05906240e5b50e46baf3 -i /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/refs+asm6+asm9+asm23+HG002_20211005.chr6.1Mbps.fa.gz -o refs+asm6+asm9+asm23+HG002_20211005.chr6.1Mbps -t 48 -p 98 -s 100000 -n 10 -k 311 -O 0.03 -T 48 -U -v -L -V chm13:#,grch38:# ; mv /scratch/refs+asm6+asm9+asm23+HG002_20211005.chr6.1Mbps /lizardfs/guarracino/HG002_assemblies_assessment/'

## contigs >= 500kbps
sbatch -p headnode -c 48 --wrap 'hostname; cd /scratch && ~/tools/pggb/pggb-5d2601127e8d08b39c8b05906240e5b50e46baf3 -i /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/refs+asm6+asm9+asm23+HG002_20211005.chr6.500kbps.fa.gz -o refs+asm6+asm9+asm23+HG002_20211005.chr6.500kbps -t 48 -p 98 -s 100000 -n 10 -k 311 -O 0.03 -T 48 -U -v -L -V chm13:#,grch38:# ; mv /scratch/refs+asm6+asm9+asm23+HG002_20211005.chr6.500kbps /lizardfs/guarracino/HG002_assemblies_assessment/'

## contigs >= 100kbps
sbatch -p headnode -c 48 --wrap 'hostname; cd /scratch && ~/tools/pggb/pggb-5d2601127e8d08b39c8b05906240e5b50e46baf3 -i /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/refs+asm6+asm9+asm23+HG002_20211005.chr6.100kbps.fa.gz -o refs+asm6+asm9+asm23+HG002_20211005.chr6.100kbps -t 48 -p 98 -s 100000 -n 10 -k 311 -O 0.03 -T 48 -U -v -L -V chm13:#,grch38:# ; mv /scratch/refs+asm6+asm9+asm23+HG002_20211005.chr6.100kbps /lizardfs/guarracino/HG002_assemblies_assessment/'

## all contigs
sbatch -p headnode -c 48 --wrap 'hostname; cd /scratch && ~/tools/pggb/pggb-5d2601127e8d08b39c8b05906240e5b50e46baf3 -i /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/refs+asm6+asm9+asm23+HG002_20211005.chr6.fa.gz -o refs+asm6+asm9+asm23+HG002_20211005.chr6 -t 48 -p 98 -s 100000 -n 10 -k 311 -O 0.03 -T 48 -U -v -L -V chm13:#,grch38:# ; mv /scratch/refs+asm6+asm9+asm23+HG002_20211005.chr6 /lizardfs/guarracino/HG002_assemblies_assessment/'
```

Compacted visualization:

```shell
cd /lizardfs/guarracino/HG002_assemblies_assessment/HG002_20211005+refs
~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f paths -i HG002_20211005+refs.fa.gz.3525971.4030258.adf7ed8.smooth.og -L | cut -f 1 -d '#' | uniq > HG002_20211005+refs.fa.gz.3525971.4030258.adf7ed8.smooth.prefix.txt

~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i HG002_20211005+refs.fa.gz.3525971.4030258.adf7ed8.smooth.og -o HG002_20211005+refs.fa.gz.3525971.4030258.adf7ed8.smooth.og.viz_inv.M.png -x 1500 -y 500 -a 10 -z -I Consensus_ -M HG002_20211005+refs/HG002_20211005+refs.fa.gz.3525971.4030258.adf7ed8.smooth.prefix.txt
```


## C4 _locus_

Find C4 coordinates:

```shell
mkdir -p /lizardfs/guarracino/HG002_assemblies_assessment/C4_locus
cd /lizardfs/guarracino/HG002_assemblies_assessment/C4_locus


wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
zgrep 'gene_id "C4A"\|gene_id "C4B"' hg38.ncbiRefSeq.gtf.gz |
  awk '$1 == "chr6"' | cut -f 1,4,5 |
  bedtools sort | bedtools merge -d 15000 | bedtools slop -l 10000 -r 20000 -g hg38.chrom.sizes |
  sed 's/chr6/grch38#chr6/g' > hg38.ncbiRefSeq.C4.coordinates.bed
```

Extraction, explosion, optimization, and sorting:

```shell
path_graph=/lizardfs/guarracino/HG002_assemblies_assessment/HG002_20211005+refs/HG002_20211005+refs.fa.gz.3525971.4030258.adf7ed8.smooth.og
prefix=$(basename ${path_graph} .og)
odgi extract -i ${path_graph} -b hg38.ncbiRefSeq.C4.coordinates.bed -o - --full-range -t 48 -P |
  odgi explode -i - --biggest 1 --sorting-criteria P --optimize -p $prefix.C4
odgi sort -i $prefix.C4.0.og -o $prefix.C4.gYs.x100.og -p gYs -x 100 -t 48 -P

```

# odgi viz: default (binned) mode
~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i $prefix.C4.gYs.x100.og -o "$(echo $prefix | tr '.' '_' )"_C4_sorted.png -c 40 -w 100 -y 50

# odgi viz: color by strand
~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i $prefix.C4.gYs.x100.og -o "$(echo $prefix | tr '.' '_' )"_C4_sorted_z.png -c 40 -w 100 -y 50 -z

# odgi viz: color by position
~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i $prefix.C4.gYs.x100.og -o "$(echo $prefix | tr '.' '_' )"_C4_sorted_du.png -c 40 -w 100 -y 50 -du

# odgi viz: color by depth
~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i $prefix.C4.gYs.x100.og -o "$(echo $prefix| tr '.' '_' )"_C4_sorted_m.png -c 40 -w 100 -y 50 -m -B Spectral:4


