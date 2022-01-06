# HG002 pangenome graphs

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

## Graph with GRCH38 + CHM13 + HG002 assemblies

Put all together:

```shell
cd /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/

cat ../references/grch38.fa ../references/chm13.fa <(zcat HG002_20211005.mat.fa.gz) <(zcat HG002_20211005.pat.fa.gz) | bgzip -@ 48 -c > HG002_20211005+refs.fa.gz
samtools faidx HG002_20211005+refs.fa.gz
```

### Pangenome graph building

Build the PGGB graph:

```shell
mkdir -p /lizardfs/guarracino/HG002_assemblies_assessment/GRCH38_CHM13_HG002/

sbatch -p lowmem -c 48 --wrap 'hostname; cd /scratch && ~/tools/pggb/pggb-5d2601127e8d08b39c8b05906240e5b50e46baf3 -i /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/HG002_20211005+refs.fa.gz -o HG002_20211005+refs -t 48 -p 98 -s 100000 -n 4 -k 311 -O 0.03 -T 48 -U -v -L -V chm13:#,grch38:# ; mv /scratch/HG002_20211005+refs /lizardfs/guarracino/HG002_assemblies_assessment/GRCH38_CHM13_HG002/'
```

Compacted visualization:

```shell
cd /lizardfs/guarracino/HG002_assemblies_assessment/HG002_20211005+refs

~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f paths -i HG002_20211005+refs.fa.gz.3525971.4030258.adf7ed8.smooth.og -L | cut -f 1 -d '#' | uniq > HG002_20211005+refs.fa.gz.3525971.4030258.adf7ed8.smooth.prefix.txt
~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i HG002_20211005+refs.fa.gz.3525971.4030258.adf7ed8.smooth.og -o HG002_20211005+refs.fa.gz.3525971.4030258.adf7ed8.smooth.og.viz_inv.M.png -x 1500 -y 500 -a 10 -z -I Consensus_ -M HG002_20211005+refs/HG002_20211005+refs.fa.gz.3525971.4030258.adf7ed8.smooth.prefix.txt
```


### C4 _locus_

Find C4 coordinates:

```shell
mkdir -p /lizardfs/guarracino/HG002_assemblies_assessment/GRCH38_CHM13_HG002/C4
cd /lizardfs/guarracino/HG002_assemblies_assessment/GRCH38_CHM13_HG002/C4

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
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

Visualization:

```shell
# odgi viz: default (binned) mode
~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i $prefix.C4.gYs.x100.og -o "$(echo $prefix | tr '.' '_' )"_C4_sorted.png -c 40 -w 100 -y 50

# odgi viz: color by strand
~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i $prefix.C4.gYs.x100.og -o "$(echo $prefix | tr '.' '_' )"_C4_sorted_z.png -c 40 -w 100 -y 50 -z

# odgi viz: color by position
~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i $prefix.C4.gYs.x100.og -o "$(echo $prefix | tr '.' '_' )"_C4_sorted_du.png -c 40 -w 100 -y 50 -du

# odgi viz: color by depth
~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i $prefix.C4.gYs.x100.og -o "$(echo $prefix| tr '.' '_' )"_C4_sorted_m.png -c 40 -w 100 -y 50 -m -B Spectral:4
```

## Graph with GRCH38 + CHM13 + asm6 + asm9 + asm23 + HG002 assemblies

### Partitioning

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

Take chr6's contigs:

```shell
cat \
  <(awk -v chm13="chm13#chr6" -v grch38="grch38#chr6" '($6 == chm13 || $6 == grch38) && $4 >= 0' /lizardfs/guarracino/HG002_assemblies_assessment/alignment/asm6+asm9+asm23+HG002_20211005.vs.chm13+grch38.s50k.l150k.p90.k16.N.approx.paf | cut -f 1) \
   <(grep chr6 $out.fa.gz.split.p80.approx.paf | awk '{ print $1,$11,$0 }' | tr ' ' '\t' |  sort -n -r -k 1,2 | awk '$1 != last { print; last = $1; }' | awk '{ if($2 >= 0) print $1}') \
   > /lizardfs/guarracino/HG002_assemblies_assessment/alignment/asm6+asm9+asm23+HG002_20211005.vs.chm13+grch38.s50k.l150k.p90.k16.N.chr6.names.txt

cat \
  <(samtools faidx ../references/grch38.fa $(echo grch38#chr6)) \
  <(samtools faidx ../references/chm13.fa $(echo chm13#chr6)) \
  <(samtools faidx asm6+asm9+asm23+HG002_20211005.fa.gz $(cat /lizardfs/guarracino/HG002_assemblies_assessment/alignment/asm6+asm9+asm23+HG002_20211005.vs.chm13+grch38.s50k.l150k.p90.k16.N.chr6.names.txt | sort)) | \
  bgzip -@ 48 -c > refs+asm6+asm9+asm23+HG002_20211005.chr6.fa.gz
samtools faidx refs+asm6+asm9+asm23+HG002_20211005.chr6.fa.gz
```

### Pangenome graph building

Build the PGGB graphs:

```shell
mkdir -p /lizardfs/guarracino/HG002_assemblies_assessment/GRCH38_CHM13_asm6_asm9_asm23_HG002/


sbatch -p headnode -c 24 --wrap 'hostname; cd /scratch && ~/tools/pggb/pggb-5d2601127e8d08b39c8b05906240e5b50e46baf3 -i /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/refs+asm6+asm9+asm23+HG002_20211005.chr6.fa.gz -o refs+asm6+asm9+asm23+HG002_20211005.chr6.s100k.l300k.p98 -t 24 -p 98 -s 100000 -n 10 -k 311 -O 0.03 -T 24 -U -v -L -V chm13:#,grch38:# ; mv /scratch/refs+asm6+asm9+asm23+HG002_20211005.chr6.s100k.l300k.p98 /lizardfs/guarracino/HG002_assemblies_assessment/GRCH38_CHM13_asm6_asm9_asm23_HG002'
sbatch -p headnode -c 24 --wrap 'hostname; cd /scratch && ~/tools/pggb/pggb-5d2601127e8d08b39c8b05906240e5b50e46baf3 -i /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/refs+asm6+asm9+asm23+HG002_20211005.chr6.fa.gz -o refs+asm6+asm9+asm23+HG002_20211005.chr6.s150k.l450k.p98 -t 24 -p 98 -s 150000 -n 10 -k 311 -O 0.03 -T 24 -U -v -L -V chm13:#,grch38:# ; mv /scratch/refs+asm6+asm9+asm23+HG002_20211005.chr6.s150k.l450k.p98 /lizardfs/guarracino/HG002_assemblies_assessment/GRCH38_CHM13_asm6_asm9_asm23_HG002'

sbatch -p headnode -c 24 --wrap 'hostname; cd /scratch && ~/tools/pggb/pggb-5d2601127e8d08b39c8b05906240e5b50e46baf3 -i /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/refs+asm6+asm9+asm23+HG002_20211005.chr6.fa.gz -o refs+asm6+asm9+asm23+HG002_20211005.chr6.s100k.l300k.p95 -t 24 -p 95 -s 100000 -n 10 -k 311 -O 0.03 -T 24 -U -v -L -V chm13:#,grch38:# ; mv /scratch/refs+asm6+asm9+asm23+HG002_20211005.chr6.s100k.l300k.p95 /lizardfs/guarracino/HG002_assemblies_assessment/GRCH38_CHM13_asm6_asm9_asm23_HG002'
sbatch -p headnode -c 24 --wrap 'hostname; cd /scratch && ~/tools/pggb/pggb-5d2601127e8d08b39c8b05906240e5b50e46baf3 -i /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/refs+asm6+asm9+asm23+HG002_20211005.chr6.fa.gz -o refs+asm6+asm9+asm23+HG002_20211005.chr6.s150k.l450k.p95 -t 24 -p 95 -s 150000 -n 10 -k 311 -O 0.03 -T 24 -U -v -L -V chm13:#,grch38:# ; mv /scratch/refs+asm6+asm9+asm23+HG002_20211005.chr6.s150k.l450k.p95 /lizardfs/guarracino/HG002_assemblies_assessment/GRCH38_CHM13_asm6_asm9_asm23_HG002'

sbatch -p headnode -c 24 --wrap 'hostname; cd /scratch && ~/tools/pggb/pggb-5d2601127e8d08b39c8b05906240e5b50e46baf3 -i /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/refs+asm6+asm9+asm23+HG002_20211005.chr6.fa.gz -o refs+asm6+asm9+asm23+HG002_20211005.chr6.s100k.l300k.p90 -t 24 -p 90 -s 100000 -n 10 -k 311 -O 0.03 -T 24 -U -v -L -V chm13:#,grch38:# ; mv /scratch/refs+asm6+asm9+asm23+HG002_20211005.chr6.s100k.l300k.p90 /lizardfs/guarracino/HG002_assemblies_assessment/GRCH38_CHM13_asm6_asm9_asm23_HG002'
sbatch -p headnode -c 24 --wrap 'hostname; cd /scratch && ~/tools/pggb/pggb-5d2601127e8d08b39c8b05906240e5b50e46baf3 -i /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/refs+asm6+asm9+asm23+HG002_20211005.chr6.fa.gz -o refs+asm6+asm9+asm23+HG002_20211005.chr6.s150k.l450k.p90 -t 24 -p 90 -s 150000 -n 10 -k 311 -O 0.03 -T 24 -U -v -L -V chm13:#,grch38:# ; mv /scratch/refs+asm6+asm9+asm23+HG002_20211005.chr6.s150k.l450k.p90 /lizardfs/guarracino/HG002_assemblies_assessment/GRCH38_CHM13_asm6_asm9_asm23_HG002'


## contigs >= 100kbps
sbatch -p headnode -c 48 --wrap 'hostname; cd /scratch && ~/tools/pggb/pggb-5d2601127e8d08b39c8b05906240e5b50e46baf3 -i /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/refs+asm6+asm9+asm23+HG002_20211005.chr6.100kbps.fa.gz -o refs+asm6+asm9+asm23+HG002_20211005.chr6.100kbps -t 48 -p 98 -s 100000 -n 10 -k 311 -O 0.03 -T 48 -U -v -L -V chm13:#,grch38:# ; mv /scratch/refs+asm6+asm9+asm23+HG002_20211005.chr6.100kbps /lizardfs/guarracino/HG002_assemblies_assessment/'

```


### MHC _locus

```shell
mkdir -p /lizardfs/guarracino/HG002_assemblies_assessment/MHC_locus
cd /lizardfs/guarracino/HG002_assemblies_assessment/MHC_locus

path_graph=/lizardfs/guarracino/HG002_assemblies_assessment/refs+asm6+asm9+asm23+HG002_20211005.chr6.1Mbps/refs+asm6+asm9+asm23+HG002_20211005.chr6.1Mbps.fa.gz.338ca16.4030258.976bd92.smooth.og
path_graph=/lizardfs/guarracino/HG002_assemblies_assessment/refs+asm6+asm9+asm23+HG002_20211005.chr6.500kbps/refs+asm6+asm9+asm23+HG002_20211005.chr6.500kbps.fa.gz.338ca16.4030258.976bd92.smooth.og
path_graph=/lizardfs/guarracino/HG002_assemblies_assessment/refs+asm6+asm9+asm23+HG002_20211005.chr6.100kbps/refs+asm6+asm9+asm23+HG002_20211005.chr6.100kbps.fa.gz.338ca16.4030258.976bd92.smooth.og
path_graph=/lizardfs/guarracino/HG002_assemblies_assessment/refs+asm6+asm9+asm23+HG002_20211005.chr6/refs+asm6+asm9+asm23+HG002_20211005.chr6.fa.gz.338ca16.4030258.976bd92.smooth.og

path_graph=/lizardfs/guarracino/HG002_assemblies_assessment/refs+asm6+asm9+asm23+HG002_20211005.chr6.s150k/refs+asm6+asm9+asm23+HG002_20211005.chr6.fa.gz.181272d.4030258.976bd92.smooth.og
prefix=$(basename ${path_graph} .og)

path_graph=/scratch/refs+asm6+asm9+asm23+HG002_20211005.chr6.s150k.l450k.p95/refs+asm6+asm9+asm23+HG002_20211005.chr6.s150k.l450k.p95.fa.gz.181272d.4030258.seqwish.gfa.prep.gfa
prefix=$(basename ${path_graph} .gfa)

/home/guarracino/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f extract -i ${path_graph}  -r grch38#chr6:29000000-34000000 -o - --full-range -t 48 -P | 
  /home/guarracino/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f sort -i - -o - --optimize |
  /home/guarracino/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f sort -i - -o ${prefix}.mhc.og -p gYs -t 48 -P

# Visualization
/home/guarracino/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i ${prefix}.mhc.og -o ${prefix}.mhc.prefix.png -s '#' -x 1500 -c 44
/home/guarracino/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i ${prefix}.mhc.og -o ${prefix}.mhc.N.png -N -c 44 -x 1500 
/home/guarracino/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i ${prefix}.mhc.og -o ${prefix}.mhc.z.png -z -c 44 -x 1500 
/home/guarracino/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i ${prefix}.mhc.og -o ${prefix}.mhc.m.png -m -c 44 -x 1500 
/home/guarracino/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f paths -i  ${prefix}.mhc.og -L | cut -f 1 -d '#' | uniq > $prefix.prefix.txt
/home/guarracino/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i ${prefix}.mhc.og -o ${prefix}.mhc.M.png -m -c 44 -x 1500 -M $prefix.prefix.txt
```


[comment]: <> (wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz)
[comment]: <> (zgrep '_alt\|_fix\|_random\|chrUn_' hg38.ncbiRefSeq.gtf.gz -v | sed -e 's/^/grch38#/' > hg38.ncbiRefSeq.grch38.coordinates.gtf)

[comment]: <> (## MHC _locus_)
[comment]: <> (```shell)
[comment]: <> (mkdir -p /lizardfs/guarracino/HG002_assemblies_assessment/MHC_locus)
[comment]: <> (cd /lizardfs/guarracino/HG002_assemblies_assessment/MHC_locus)
[comment]: <> (zgrep chr6 hg38.ncbiRefSeq.gtf.gz | awk '{print $1"\t"$4"\t"$5"\t"$10"\t"$6"\t"$7}' | sed 's/"//g' | sed 's/;//g' | sed 's/chr6/grch38#chr6/g' > hg38.chr6.bed)
[comment]: <> (cd /lizardfs/guarracino/HG002_assemblies_assessment)
[comment]: <> (cut -f 1 HG002_20211005+refs.fa.gz.fai | grep S6 > HG002_20211005.S6.txt)
[comment]: <> (/home/guarracino/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f position -i HG002_20211005+refs/HG002_20211005+refs.fa.gz.3525971.4030258.adf7ed8.smooth.og -p grch38#chr6,29000000 -R HG002_20211005.S6.txt)
[comment]: <> (##source.path.pos	target.path.pos	dist.to.ref	strand.vs.ref)
[comment]: <> (grch38#chr6,29000000,+	HG002_20211005.mat#S6,28946544,+	0	+)
[comment]: <> (/home/guarracino/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f position -i HG002_20211005+refs/HG002_20211005+refs.fa.gz.3525971.4030258.adf7ed8.smooth.og -p grch38#chr6,34000000 -R HG002_20211005.S6.txt)
[comment]: <> (#source.path.pos	target.path.pos	dist.to.ref	strand.vs.ref)
[comment]: <> (grch38#chr6,34000000,+	HG002_20211005.mat#S6,33920246,+	0	+)
[comment]: <> (sbatch -p headnode -w octopus01 -c 48 --wrap '\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f extract -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.chr6.og -o /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.chr6.MHC.og -r HG002_20211005.mat#S6:0-171862164:28946544-33920246 -E -t 48 -P')
[comment]: <> (```)
