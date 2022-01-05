# Pangenomic assessment of HG002 assemblies

## Tools

```shell
mkdir -p ~/tools && cd ~/tools

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

git clone --recursive https://github.com/pangenome/odgi.git
cd odgi
git checkout 67a7e5bb2f328888e194845a362cef9c8ccc488f
mv bin/odgi bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f
cmake -H. -Bbuild && cmake --build build -- -j 48
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

Download HG002 assemblies:

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

Merge all assemblies:

```shell
# To merge assemblies sorted by ID
FA_GZS=$(cat ../data/HGRC_bakeoff_HG002_assemblies_v3_renaming.tsv | sed 1,1d | sed 's/"//g' | while read -r Id Filename AbbreviatedName; do AbbreviatedName2=$(echo $AbbreviatedName | sed 's/ /_/g'); echo ${AbbreviatedName2}.fa.gz; done | tr '\n' ' ')

zcat $FA_GZS | bgzip -@ 48 -c > HG002_all.fa.gz
samtools faidx HG002_all.fa.gz
```

## Mapping and alignment

Generate all-vs-all mappings:

```shell
mkdir -p /lizardfs/guarracino/HG002_assemblies_assessment/alignment/

# -s 20k -l 60k -p 95
sbatch -p lowmem -c 48 --wrap 'cd /scratch && \time -v ~/tools/wfmash/build/bin/wfmash-09e73eb3fcf24b8b7312b8890dd0741933f0d1cd -X -s 20k -l 60k -p 95 -n 45 -k 16 -t 48 /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/HG002_all.fa.gz /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/HG002_all.fa.gz -m > HG002_all.s20k.l60k.p95.n45.k16.approx.paf && mv HG002_all.s20k.l60k.p95.n45.k16.approx.paf /lizardfs/guarracino/HG002_assemblies_assessment/alignment/'

# -s 100k -l 300k -p 98
sbatch -p lowmem -c 48 --wrap 'cd /scratch && \time -v ~/tools/wfmash/build/bin/wfmash-09e73eb3fcf24b8b7312b8890dd0741933f0d1cd -X -s 100k -l 300k -p 98 -n 45 -k 16 -t 48 /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/HG002_all.fa.gz /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/HG002_all.fa.gz -m > HG002_all.s100k.l300k.p98.n45.k16.approx.paf && mv HG002_all.s100k.l300k.p98.n45.k16.approx.paf /lizardfs/guarracino/HG002_assemblies_assessment/alignment/'
```

Split the mappings in chunks:

```shell
cd /lizardfs/guarracino/HG002_assemblies_assessment/alignment/

# -s 20k -l 60k -p 95
python3 ~/tools/wfmash/scripts/split_approx_mappings_in_chunks.py /lizardfs/guarracino/HG002_assemblies_assessment/alignment/HG002_all.s20k.l60k.p95.n45.k16.approx.paf 5

# -s 100k -l 300k -p 98
python3 ~/tools/wfmash/scripts/split_approx_mappings_in_chunks.py /lizardfs/guarracino/HG002_assemblies_assessment/alignment/HG002_all.s100k.l300k.p98.n45.k16.approx.paf 5
```

Run the alignments on multiple nodes:

```shell
# -s 20k -l 60k -p 95
seq 0 4 | while read i; do sbatch -p lowmem -c 48 --wrap 'cd /scratch && ~/tools/wfmash/build/bin/wfmash-09e73eb3fcf24b8b7312b8890dd0741933f0d1cd -X -s 20k -l 60k -p 95 -n 45 -k 16 -t 48 /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/HG002_all.fa.gz /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/HG002_all.fa.gz -i /lizardfs/guarracino/HG002_assemblies_assessment/alignment/HG002_all.s20k.l60k.p95.n45.k16.approx.paf.chunk_'$i'.paf | pigz -c > HG002_all.s20k.l60k.p95.n45.k16.approx.paf.chunk_'$i'.paf.gz && mv HG002_all.s20k.l60k.p95.n45.k16.approx.paf.chunk_'$i'.paf.gz /lizardfs/guarracino/HG002_assemblies_assessment/alignment/;'; done

# -s 100k -l 300k -p 98
seq 0 4 | while read i; do sbatch -p lowmem -c 48 --wrap 'cd /scratch && ~/tools/wfmash/build/bin/wfmash-09e73eb3fcf24b8b7312b8890dd0741933f0d1cd -X -s 100k -l 300k -p 98 -n 45 -k 16 -t 48 /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/HG002_all.fa.gz /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/HG002_all.fa.gz -i /lizardfs/guarracino/HG002_assemblies_assessment/alignment/HG002_all.s100k.l300k.p98.n45.k16.approx.paf.chunk_'$i'.paf | pigz -c > HG002_all.s100k.l300k.p98.n45.k16.approx.paf.chunk_'$i'.paf.gz && mv HG002_all.s100k.l300k.p98.n45.k16.approx.paf.chunk_'$i'.paf.gz /lizardfs/guarracino/HG002_assemblies_assessment/alignment/'; done
```

## Induce the SEQWISH graph

```shell
mkdir -p /lizardfs/guarracino/HG002_assemblies_assessment/graphs/

# -s 20k -l 60k -p 95
PAFS=$(ls /lizardfs/guarracino/HG002_assemblies_assessment/alignment/HG002_all.s20k.l60k.p95.n45.k16.approx.paf.chunk_*.paf.gz | tr '\n' ',')
PAFS=${PAFS::-1}
sbatch -p highmem -w octopus02 -c 48 --wrap 'cd /scratch && \time -v ~/tools/seqwish/bin/seqwish-ccfefb016fcfc9937817ce61dc06bbcf382be75e -t 48 -s /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/HG002_all.fa.gz -p '$PAFS' -k 79 -B 50M -g HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.gfa -P && mv HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.gfa /lizardfs/guarracino/HG002_assemblies_assessment/graphs/'

# -s 100k -l 300k -p 98
PAFS=$(ls /lizardfs/guarracino/HG002_assemblies_assessment/alignment/HG002_all.s100k.l300k.p98.n45.k16.approx.paf.chunk_*.paf.gz | tr '\n' ',')
PAFS=${PAFS::-1}
sbatch -p highmem -w octopus11 -c 48 --wrap 'cd /scratch && \time -v ~/tools/seqwish/bin/seqwish-ccfefb016fcfc9937817ce61dc06bbcf382be75e -t 48 -s /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/HG002_all.fa.gz -p '$PAFS' -k 79 -B 50M -g HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.gfa -P && mv HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.gfa /lizardfs/guarracino/HG002_assemblies_assessment/graphs/'
```


## Build the ODGI graph and sort it

Build:

```shell
# -s 20k -l 60k -p 95
\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f build -g /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.gfa -o /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.og -t 48 -P

# -s 100k -l 300k -p 98
\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f build -g /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.gfa -o /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.og -t 48 -P
```

Statistics:

```shell
# -s 20k -l 60k -p 95
\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f stats -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.og -bS > /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.stats.txt

# -s 100k -l 300k -p 98
\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f stats -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.og -bS > /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.stats.txt
```

Sort:

```shell
# -s 20k -l 60k -p 95
# -p Y
sbatch -p highmem -w octopus11 -c 48 --wrap 'hostname; cd /scratch && ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f sort -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.og -o /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Y.x100.og -t 48 -Y -x 100 -P'
# -p Ys
sbatch -p highmem -w octopus11 -c 48 --wrap 'hostname; cd /scratch && ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f sort -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Y.x100.og -o /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Ys.x100.og -t 48 -p s -P'
# -p Ygs
sbatch -p highmem -w octopus02 -c 48 --wrap 'hostname; cd /scratch && ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f sort -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.og -o /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Ygs.x100.og -t 48 -p Ygs -x 100 -P'

# -s 100k -l 300k -p 98
# -p Y
sbatch -p highmem -w octopus02 -c 48 --wrap 'hostname; cd /scratch && ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f sort -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.og -o /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.og -t 48 -Y -x 100 -P'
# -p Ys
sbatch -p highmem -w octopus02 -c 48 --wrap 'hostname; cd /scratch && ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f sort -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.og  -o /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Ys.x100.og -t 48 -p s -P'
```

Plots:

```shell
# -s 20k -l 60k -p 95
sbatch -p headnode -w octopus01 -c 48 --wrap "\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f paths -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.og -L | cut -f 1 -d '#' | sed 's/$/#/g' | uniq > HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.prefix.txt"
# -p Y
sbatch -p headnode -w octopus01 -c 48 --wrap "\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Y.x100.og -o HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Y.x100.z.M.png -z -M HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.prefix.txt -P"
sbatch -p headnode -w octopus01 -c 48 --wrap "\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Y.x100.og -o HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Y.x100.m.Spectral4.M.png -m -B Spectral:4 -M HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.prefix.txt -P"
# -p Ys
sbatch -p headnode -w octopus01 -c 48 --wrap "\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Ys.x100.og -o HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Ys.x100.z.M.png -z -M HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.prefix.txt -P"
sbatch -p headnode -w octopus01 -c 48 --wrap "\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Ys.x100.og -o HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Ys.x100.m.Spectral4.M.png -m -B Spectral:4 -M HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.prefix.txt -P"
# -p Ygs
sbatch -p headnode -w octopus01 -c 48 --wrap "\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Ygs.x100.og -o HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Ygs.x100.z.M.png -z -M HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.prefix.txt -P"
sbatch -p headnode -w octopus01 -c 48 --wrap "\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Ygs.x100.og -o HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Ygs.x100.m.Spectral4.M.png -m -B Spectral:4 -M HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.prefix.txt -P"

# -s 100k -l 300k -p 98
sbatch -p headnode -w octopus01 -c 48 --wrap "\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f paths -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.og -L | cut -f 1 -d '#' | sed 's/$/#/g' | uniq > HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.prefix.txt"
# -p Y
sbatch -p headnode -w octopus01 -c 48 --wrap '\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.og -o /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.z.M.png -z -M HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.prefix.txt -P'
sbatch -p headnode -w octopus01 -c 48 --wrap '\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.og -o /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.m.Spectral4.M.png -m -B Spectral:4 -M HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.prefix.txt -P'
# -p Ys
sbatch -p headnode -w octopus01 -c 48 --wrap '\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Ys.x100.og -o /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Ys.x100.z.M.png -z -M HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.prefix.txt -P'
sbatch -p headnode -w octopus01 -c 48 --wrap '\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f viz -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Ys.x100.og -o /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Ys.x100.m.Spectral4.M.png -m -B Spectral:4 -M HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.prefix.txt -P'
```

## Partitioning

Map each assembly against the scaffolded references:

```shell
sbatch -p lowmem -c 48 --wrap 'cd /scratch && \time -v ~/tools/wfmash/build/bin/wfmash-09e73eb3fcf24b8b7312b8890dd0741933f0d1cd -s 50k -l 150k -p 90 -k 16 -N -t 48 /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13+grch38.fa /lizardfs/guarracino/HG002_assemblies_assessment/assemblies/HG002_all.fa.gz -m > HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.approx.paf && mv HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.approx.paf /lizardfs/guarracino/HG002_assemblies_assessment/alignment/'
```

Partitioning by chromosome, ignoring the unmapped contigs:

```shell
cd /lizardfs/guarracino/HG002_assemblies_assessment/alignment

( seq 22; echo M ) | while read i; do
  echo chr$i
  awk -v chm13="chm13#chr$i" -v grch38="grch38#chr$i" '$6 == chm13 || $6 == grch38' HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.approx.paf | cut -f 1,3,4 > HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chr$i.bed;
  cut -f 1 HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chr$i.bed > HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chr$i.names.txt;
done
```

Partitioning by chromosome type (sex and autosome), ignoring the unmapped contigs:

```shell
grep 'chm13#chrX\|grch38#chrX\|grch38#chrY' HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.approx.paf | cut -f 1,3,4 > HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chrXY.bed
cut -f 1 HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chrXY.bed > HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chrXY.names.txt

grep 'chm13#chrX\|grch38#chrX\|grch38#chrY' HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.approx.paf -v | cut -f 1,3,4 > HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chr1to22.bed
cut -f 1 HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chr1to22.bed > HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chr1to22.names.txt
```

Create partitioned graphs:

```shell
# -s 20k -l 60k -p 95
#
#sbatch -p headnode -w octopus01 -c 48 --wrap '\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f extract -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Y.x100.og -o /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Y.x100.chrXY.og -b /lizardfs/guarracino/HG002_assemblies_assessment/alignmentHG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chrXY.bed -p /lizardfs/guarracino/HG002_assemblies_assessment/alignment/HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chrXY.names.txt -t 48 -P'
#sbatch -p headnode -w octopus01 -c 48 --wrap '\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f extract -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Y.x100.og -o /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.Y.x100.chr1to22.og -b /lizardfs/guarracino/HG002_assemblies_assessment/alignment/HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chr1to22.bed -p /lizardfs/guarracino/HG002_assemblies_assessment/alignment/HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chr1to22.names.txt -t 48 -P'

# -s 100k -l 300k -p 98
( seq 22; echo M ) | while read i; do
  sbatch -p headnode -w octopus01 -c 24 --wrap '\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f extract -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.og -o /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.chr'$i'.og -b /lizardfs/guarracino/HG002_assemblies_assessment/alignment/HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chr'$i'.bed -p /lizardfs/guarracino/HG002_assemblies_assessment/alignment/HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chr'$i'.names.txt -t 24 -P'
done
sbatch -p headnode -w octopus01 -c 48 --wrap '\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f extract -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.og -o /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.chrXY.og -b /lizardfs/guarracino/HG002_assemblies_assessment/alignment/HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chrXY.bed -p /lizardfs/guarracino/HG002_assemblies_assessment/alignment/HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chrXY.names.txt -t 48 -P'
sbatch -p headnode -w octopus01 -c 48 --wrap '\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f extract -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.og -o /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.chr1to22.og -b /lizardfs/guarracino/HG002_assemblies_assessment/alignment/HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chr1to22.bed -p /lizardfs/guarracino/HG002_assemblies_assessment/alignment/HG002_all.vs.chm13+grch38.s50k.l150k.p90.k16.N.chr1to22.names.txt -t 48 -P'
```


##  Multidimensional Scaling analysis

Get matrix of distances on the full graph:

```shell
mkdir -p /lizardfs/guarracino/HG002_assemblies_assessment/matrix/

# -s 20k -l 60k -p 95
sbatch -p lowmem -c 48 --wrap 'hostname; cd /scratch && ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f paths -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.og -d -D "#" -t 48 > HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.chrAll.dist.tsv && mv HG002_all.s20k.l60k.p95.n45.k16.seqwish.k79.B50M.chrAll.dist.tsv /lizardfs/guarracino/HG002_assemblies_assessment/matrix/'

# -s 100k -l 300k -p 98
sbatch -p lowmem -c 48 --wrap 'hostname; cd /scratch && ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f paths -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.og -d -D "#" -t 48 > HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.chrAll.dist.tsv && mv HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.chrAll.dist.tsv /lizardfs/guarracino/HG002_assemblies_assessment/matrix/'
```

Get matrix of distances on the partitioned graphs:

```shell
## -s 20k -l 60k -p 95
#

# -s 100k -l 300k -p 98
( seq 22; echo M ) | while read i; do
  sbatch -p headnode -c 24 --wrap '\time -v ~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f paths -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.chr'$i'.og -d -D "#" -t 24 > /lizardfs/guarracino/HG002_assemblies_assessment/matrix/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.chr'$i'.dist.tsv'
done
sbatch -p headnode -c 24 --wrap '~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f paths -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.chrXY.og -d -D "#" -t 24 > /lizardfs/guarracino/HG002_assemblies_assessment/matrix/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.chrXY.dist.tsv'
sbatch -p headnode -c 24 --wrap '~/tools/odgi/bin/odgi-67a7e5bb2f328888e194845a362cef9c8ccc488f paths -i /lizardfs/guarracino/HG002_assemblies_assessment/graphs/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.chr1to22.og -d -D "#" -t 24 > /lizardfs/guarracino/HG002_assemblies_assessment/matrix/HG002_all.s100k.l300k.p98.n45.k16.seqwish.k79.B50M.Y.x100.chr1to22.dist.tsv'
```

Go to `haccard_mds.R` script.
