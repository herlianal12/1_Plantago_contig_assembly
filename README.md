# Plantago contig assembly
Main steps in generating contigs
1. Genomic DNA extraction
2. PacBio CLR sequencing
3. Converting PacBio unaligned bam files into fastq files
4. Removing contaminant
5. Contig assembly
6. Polishing
7. Purging and clipping haplotig


bamtofastq.sh
```
bam2fastq -c 9 m54078_170831_060817.subreads.bam -o m54078_170831_060817.subreads
bam2fastq -c 9 m54078_170831_160707.subreads.bam -o m54078_170831_160707.subreads
bam2fastq -c 9 m54078_170901_080645.subreads.bam -o m54078_170901_080645.subreads
bam2fastq -c 9 m54078_170901_180552.subreads.bam -o m54078_170901_180552.subreads
bam2fastq -c 9 m54078_170902_041524.subreads.bam -o m54078_170902_041524.subreads
bam2fastq -c 9 m54078_170902_142504.subreads.bam -o m54078_170902_142504.subreads
bam2fastq -c 9 m54078_170903_003441.subreads.bam -o m54078_170903_003441.subreads 

cat *.subreads.fastq > Plantago_pacbio.fastq
bgzip -c -l 9 Plantago_pacbio.fastq > Plantago_pacbio.fastq.gz
```


minimap_index.sh
```
minimap2 -d references/plantago_chloroplast.fasta.gz.mmi references/plantago_chloroplast.fasta.gz
minimap2 -d references/plantago_mitocondria.fasta.gz.mmi references/plantago_mitocondria.fasta.gz
minimap2 -d database_rrna/SILVA_128_LSURef_tax_silva.fasta.gz.mmi database_rrna/SILVA_128_LSURef_tax_silva.fasta.gz
minimap2 -d database_rrna/SILVA_138_SSURef_NR99_tax_silva.fasta.gz.mmi database_rrna/SILVA_138_SSURef_NR99_tax_silva.fasta.gz
minimap2 -d database_rrna/SILVA_138_SSURef_tax_silva.fasta.gz.mmi database_rrna/SILVA_138_SSURef_tax_silva.fasta.gz
```

minimap_samtools.sh
```
data="assembly/raw_reads/Plantago_pacbio.fastq.gz"
index="references/plantago_chloroplast.fasta.gz.mmi"
output_1="assembly/raw_reads/Plantago_pacbio_no_chloro.bam"
output_2="assembly/raw_reads/Plantago_pacbio_chloroplast.bam"

time minimap2 \
-ax map-pb \
-t 2 $index $data \
| samtools view -f 0x04 -u \
| samtools sort \
--threads 2 -l 7 \
-o $output_1

time minimap2 \
-ax map-pb \
-t 2 $index $data \
| samtools view -F 0x04 -u \
| samtools sort \
--threads 2 -l 7 \
-o $output_2


bamToFastq -i assembly/raw_reads/Plantago_pacbio_no_chloro.bam -fq assembly/raw_reads/Plantago_pacbio_no_chloro.fastq
bgzip -c -l 9 Plantago_pacbio_no_chloro.fastq > Plantago_pacbio_no_chloro.fastq.gz


data="assembly/raw_reads/Plantago_pacbio_no_chloro.fastq.gz"
index="references/plantago_mitocondria.fasta.gz.mmi"
output="assembly/raw_reads/Plantago_pacbio_no_mito_chloro.bam"

time minimap2 \
-ax map-pb \
-t 2 $index $data \
| samtools view -f 0x04 -u \
| samtools sort \
--threads 2 -l 7 \
-o $output

bamToFastq -i assembly/raw_reads/Plantago_pacbio_no_mito_chloro.bam -fq assembly/raw_reads/Plantago_pacbio_no_mito_chloro.fastq
bgzip -c -l 9 Plantago_pacbio_no_mito_chloro.fastq > Plantago_pacbio_no_mito_chloro.fastq.gz


data="assembly/raw_reads/Plantago_pacbio_no_mito_chloro.bam.gz"
index="database_rrna/SILVA_128_LSURef_tax_silva.fasta.gz.mmi"
output="assembly/raw_reads/Plantago_pacbio_no_mito_chloro_rrna1.bam"

time minimap2 \
-ax map-pb \
-t 2 $index $data \
| samtools view -f 0x04 -u \
| samtools sort \
--threads 2 -l 7 \
-o $output

bamToFastq -i assembly/raw_reads/Plantago_pacbio_no_mito_chloro_rrna1_2_3.bam -fq assembly/raw_reads/Plantago_pacbio_no_mito_chloro_rrna1.fastq
bgzip -c -l 9 Plantago_pacbio_no_mito_chloro_rrna1_2_3.fastq > Plantago_pacbio_no_mito_chloro_rrna1.fastq.gz

data="assembly/raw_reads/Plantago_pacbio_no_mito_chloro_rrna1.bam.gz"
index="database_rrna/SILVA_138_SSURef_NR99_tax_silva.fasta.gz.mmi"
output="assembly/raw_reads/Plantago_pacbio_no_mito_chloro_rrna1_2.bam"

time minimap2 \
-ax map-pb \
-t 2 $index $data \
| samtools view -f 0x04 -u \
| samtools sort \
--threads 2 -l 7 \
-o $output

bamToFastq -i assembly/raw_reads/Plantago_pacbio_no_mito_chloro_rrna1_2_3.bam -fq assembly/raw_reads/Plantago_pacbio_no_mito_chloro_rrna1_2.fastq
bgzip -c -l 9 Plantago_pacbio_no_mito_chloro_rrna1_2_3.fastq > Plantago_pacbio_no_mito_chloro_rrna1_2.fastq.gz


data="assembly/raw_reads/Plantago_pacbio_no_mito_chloro_rrna1_2.bam.gz"
index="database_rrna/SILVA_138_SSURef_tax_silva.fasta.gz.mmi"
output="assembly/raw_reads/Plantago_pacbio_no_mito_chloro_rrna1_2_3.bam"

time minimap2 \
-ax map-pb \
-t 2 $index $data \
| samtools view -f 0x04 -u \
| samtools sort \
--threads 2 -l 7 \
-o $output

bamToFastq -i assembly/raw_reads/Plantago_pacbio_no_mito_chloro_rrna1_2_3.bam -fq assembly/raw_reads/Plantago_pacbio_no_mito_chloro_rrna1_2_3.fastq
bgzip -c -l 9 Plantago_pacbio_no_mito_chloro_rrna1_2_3.fastq > Plantago_pacbio_no_mito_chloro_rrna1_2_3.fastq.gz
```

canu.sh
```
canu -p "Po_new" -d assembly/canu_po genomeSize=522m \
-correct -pacbio assembly/raw_reads/Plantago_pacbio_no_mito_chloro_rrna1_2_3.fastq.gz \
corMhapSensitivity=high corMinCoverage=0 corOutCoverage=200 correctedErrorRate=0.105

canu -p "Po_new" -d assembly/canu_po genomeSize=522m \
-trim -pacbio-corrected assembly/canu_po/Po_new.correctedReads.fasta.gz \
correctedErrorRate=0.105

canu -p "Po_new" -d assembly/canu_po genomeSize=522m \
-assemble -pacbio-corrected assembly/canu_po/Po_new.trimmedReads.fasta.gz \
correctedErrorRate=0.105 batMemory=9 ovbMemory=16 ovsMemory=16 executiveMemory=1 cnsMemory=20 cnsThreads=8 \
gridOptions="--partition=batch,highmem --time=72:00:00" "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50" utgovlMemory=30
```

polish.sh

```
for i in *.subreads.bam
do
pbindex $i
done

ls *.subreads.bam > mymovies.fofn
dataset create --type SubreadSet --name Plantago PlantagoGenomeSet.subreadset.xml mymovies.fofn
grep ’@’ ?.fastq > PlantagoGenome.txt
sed 's|[@,]||g' PlantagoGenome.txt > PlantagoGenome_final.txt
dataset filter PlantagoGenomeSet.subreadset.xml Plantago_filter.subreadset.xml 'qname=PlantagoGenome_final.txt'
pbmm2 align --log-level INFO --log-file pbmm2_log --sample Plantago /fast/users/a1697274/genome/assembly/canu_po/additional_files/Po_new_new.contigs.fasta.mmi Plantago_filter.subreadset.xml Plantago.aligned.bam

samtools sort -m 10G -o Plantago_sorted_aligned.bam -T tmp.ali Plantago.aligned.bam
samtools index -b Plantago_sorted_aligned.bam Plantago_sorted_aligned.bam.bai


awk '{ print $1 }' Po_new_new.contigs.fasta.gz.fai > contigs.txt

for line in `cat contigs.txt`
do
/hpcfs/users/a1697274/.conda/envs/pacbio/bin/samtools view -bh Plantago_sorted_aligned.bam ${line} > split/${line}.bam
done

for i in split/*.bam
do
samtools index -b $i $i.bai
done

tig00000006-tig00002996
tig00003003-tig00005999
tig00006000-tig00210499
tig00210500-tig00211778
```
```
#!/usr/bin/env python
import sys
"""
usage: python Fasta_splitter.py genome.fasta 2> /dev/null
"""
f=open(sys.argv[1],"r");
opened = False
for line in f :
  if(line[0] == ">") :
    if(opened) :
      of.close()
    opened = True
    of=open("%s.fasta" % (line[1:12].rstrip()), "w")
    print(line[1:12].rstrip())
  of.write(line)
of.close()

####modifying from https://raw.githubusercontent.com/harish0201/General_Scripts/master/Fasta_splitter.py
tig00000006-tig00002996
tig00003003-tig00005999
tig00006000-tig00210499
tig00210500-tig00211778


for i in *.fasta.fai
do
OUTPUT=$(basename "$i" .fasta.fai).bed
awk -F'\t' '{ print $1,$2}' $i > $OUTPUT
done

N=100
(
for file in *.fasta
do
((i=i%N)); ((i++==0)) && wait
OUTPUT=$(basename "$file").mmi
pbmm2 index $file $OUTPUT &
done
)

```
```
awk '{ print $1":0-"$2 }' Po_new_new.contigs.fasta.gz.fai > window.txt
tig00000006-tig00002996 --> window_1.txt
tig00003003-tig00005999 --> window_2.txt
tig00006000-tig00210499 --> window_3.txt
tig00210500-tig00211778 --> window_4.txt

Head –n 1100
Head –n 2194 file.txt | tail –n 1094
Head –n 3114 file.txt | tail –n 920
Head –n 4365 file.txt | tail –n 1251


for window in `cat window_1.txt`
do
line=$(echo ${window} | cut -c 1-11)
echo "gcpp --max-iterations 4 --log-level INFO --log-file polished_seqs/file_1/${line}.log \
-w ${window} -r fasta_split/fasta_1/${line}.fasta \
-o polished_seqs/file_1/${line}.polished.fastq,polished_seqs/file_1/${line}.polished.fasta,polished_seqs/file_1/${line}.polished.gff \   
split/bam_1/${line}.bam"
done | parallel -j7 --tmpdir TMPDIR_1 --compress

for window in `cat window_2.txt`
do
line=$(echo ${window} | cut -c 1-11)
echo "gcpp --max-iterations 4 --log-level INFO --log-file polished_seqs/file_2/${line}.log \
-w ${window} -r fasta_split/fasta_2/${line}.fasta \
-o polished_seqs/file_2/${line}.polished.fastq,polished_seqs/file_2/${line}.polished.fasta,polished_seqs/file_2/${line}.polished.gff \   
split/bam_2/${line}.bam"
done | parallel -j7 --tmpdir TMPDIR_2 --compress

for window in `cat window_3.txt`
do
line=$(echo ${window} | cut -c 1-11)
echo "gcpp --max-iterations 4 --log-level INFO --log-file polished_seqs/file_3/${line}.log \
-w ${window} -r fasta_split/fasta_3/${line}.fasta \
-o polished_seqs/file_3/${line}.polished.fastq,polished_seqs/file_3/${line}.polished.fasta,polished_seqs/file_3/${line}.polished.gff \   
split/bam_3/${line}.bam"
done | parallel -j7 --tmpdir TMPDIR_3 --compress

for window in `cat window_4.txt`
do
line=$(echo ${window} | cut -c 1-11)
echo "gcpp --max-iterations 4 --log-level INFO --log-file polished_seqs/file_4/${line}.log \
-w ${window} -r fasta_split/fasta_4/${line}.fasta \
-o polished_seqs/file_4/${line}.polished.fastq,polished_seqs/file_4/${line}.polished.fasta,polished_seqs/file_4/${line}.polished.gff \   
split/bam_4/${line}.bam"
done | parallel -j7 --tmpdir TMPDIR_4 --compress

cat file_1/*.fasta > polished_1.fasta
cat file_2/*.fasta > polished_2.fasta
cat file_3/*.fasta > polished_3.fasta
cat file_4/*.fasta > polished_4.fasta

cat polished_*.fasta > polished.fasta
bgzip -c -l 9 polished.fasta > polished.fasta.gz
```

purging.sh
```
minimap2 -d /hpcfs/users/a1697274/genome/assembly/plantago_genome_sequences/sequel/polished_seqs/polished.fasta.mmi /hpcfs/users/a1697274/genome/assembly/plantago_genome_sequences/sequel/polished_seqs/polished.fasta

minimap2 -t 4 -ax map-pb /hpcfs/users/a1697274/genome/assembly/plantago_genome_sequences/sequel/polished_seqs/polished.fasta.mmi \
/hpcfs/users/a1697274/genome/assembly/raw_reads/Plantago_pacbio_no_mito_chloro_rrna1_2_3.fasta --secondary=no \
| samtools sort -m 1G -o after_polishing_aligned.bam -T after_polishing_tmp.ali

purge_haplotigs  hist  -b after_polishing_aligned.bam  \
-g /hpcfs/users/a1697274/genome/assembly/plantago_genome_sequences/sequel/polished_seqs/polished.fasta  -t 4 -d 200

purge_haplotigs cov -i ./after_polishing_aligned.bam.gencov -l 5 -m 70 -h 190 -o try1_coverage_stats.csv -j 80 -s 80

purge_haplotigs purge -g /hpcfs/users/a1697274/genome/assembly/plantago_genome_sequences/sequel/polished_seqs/polished.fasta \
-c try1_coverage_stats.csv -t 4 -o try1_final

purge_haplotigs  clip  -p try1_final.fasta -h try1_final.haplotigs.fasta
```
