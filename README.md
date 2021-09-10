# **Plantago contig assembly**

This respository aims to document every step in generating contigs from PacBio CLR reads.

Main steps in generating contigs:
1. Genomic DNA extraction and PacBio CLR sequencing
2. Installing bioinformatic tools
3. Converting PacBio unaligned bam files into fastq files
4. Removing contaminant
5. Genome size prediction and contig assembly
6. Polishing
7. Purging

**Step 1. Genomic DNA extraction and sequencing process were explained in this publication.**

**Step 2. Installing softwares**

I used conda to install all tools and I created several environments due to incompatibility of softwares.
An example how to create conda environment:
```
conda create -n pacbio
conda activate pacbio
conda install -c bioconda bam2fastx  
```
List of main tools for contig assembly:
- Pacbio tools (bam2fastx 1.3.1, pbbam 1.6.0, pbcommand 2.1.1, pbcopper 1.9.1, pbcore 2.1.2, pbcoretools 0.8.1, pbgcpp 1.0.0, pbmm2 1.4.0, pbzip2 1.1.13) :   https://github.com/PacificBiosciences/pbbioconda
- tabix 0.2.6 : https://github.com/samtools/tabix
- canu 2.1.1 : https://github.com/marbl/canu 
- minimap2 2.17 : https://github.com/lh3/minimap2
- samtools 1.11 :https://github.com/samtools/samtools
- bedtools 2.29.2 : https://github.com/arq5x/bedtools2 or https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html


**Step 3. Converting PacBio unaligned bam files into fastq files**

Raw data from PacBio CLR are in bam format. For downstream analysis, we need to convert them into fastq files. 

Here are examples of our data:
```
bam2fastq -c 9 m54078_170831_060817.subreads.bam -o m54078_170831_060817.subreads
bam2fastq -c 9 m54078_170831_160707.subreads.bam -o m54078_170831_160707.subreads
bam2fastq -c 9 m54078_170901_080645.subreads.bam -o m54078_170901_080645.subreads
bam2fastq -c 9 m54078_170901_180552.subreads.bam -o m54078_170901_180552.subreads
bam2fastq -c 9 m54078_170902_041524.subreads.bam -o m54078_170902_041524.subreads
bam2fastq -c 9 m54078_170902_142504.subreads.bam -o m54078_170902_142504.subreads
bam2fastq -c 9 m54078_170903_003441.subreads.bam -o m54078_170903_003441.subreads 

### if you have a lot of files, you can do this:

for i in *.bam
do
outfile=$(basename ${i} .bam)
bam2fastq -c 9 ${i} -o ${outfile}
done


### concatenating and compressing all fastq files
cat *.subreads.fastq > Plantago_pacbio.fastq
bgzip -c -l 9 Plantago_pacbio.fastq > Plantago_pacbio.fastq.gz
```

**Step 4. Removing contaminant**

I found removing contaminants from PacBio raw reads helped me to solve my problem in contig assembly (Canu). We interested in nuclear genome, so chloroplast and mitochondrial reads are considered as contaminants. Plantago chloroplast genome can be found at https://www.ncbi.nlm.nih.gov/nuccore/MH205737.1/) and a mithochondrial gene is in here https://www.ncbi.nlm.nih.gov/nuccore/EU069524.1/). Only one mitochondrial gene was found in NCBI database (mitochondrial genome is still not available in May 2021). 

Creating index file
```
minimap2 -d references/plantago_chloroplast.fasta.gz.mmi references/plantago_chloroplast.fasta.gz
minimap2 -d references/plantago_mitocondria.fasta.gz.mmi references/plantago_mitocondria.fasta.gz
```
Removing chroloplast reads
```
data="assembly/raw_reads/Plantago_pacbio.fastq.gz"
index="references/plantago_chloroplast.fasta.gz.mmi"
output="assembly/raw_reads/Plantago_pacbio_no_chloro.bam"

time minimap2 \
-ax map-pb \
-t 2 $index $data \
| samtools view -f 0x04 -u \
| samtools sort \
--threads 2 -l 7 \
-o $output

### Converting a bam file to fastq file then compressing it
bamToFastq -i assembly/raw_reads/Plantago_pacbio_no_chloro.bam -fq assembly/raw_reads/Plantago_pacbio_no_chloro.fastq
bgzip -c -l 9 Plantago_pacbio_no_chloro.fastq > Plantago_pacbio_no_chloro.fastq.gz
```

Removing mithochondrial reads
```
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

### Converting a bam file to fastq file then compressing it
bamToFastq -i assembly/raw_reads/Plantago_pacbio_no_mito_chloro.bam -fq assembly/raw_reads/Plantago_pacbio_no_mito_chloro.fastq
bgzip -c -l 9 Plantago_pacbio_no_mito_chloro.fastq > Plantago_pacbio_no_mito_chloro.fastq.gz
```

**Step 5. Genome size prediction and contig assembly**

To predict Plantago ovata genome size, I utilized publicly short read genomic Illumina data (SRR10076762) using genomescope2 (https://github.com/tbenavi1/genomescope2.0).
I run three steps of Canu on clean PacBio reads to generate contig assembly.

```
canu -p Po_2021 -d canu_2021 genomeSize=584m \
-correct -pacbio Plantago_pacbio_no_mito_chloro.fastq.gz \
corMhapSensitivity=high corMinCoverage=0 corOutCoverage=200 correctedErrorRate=0.105 \
gridEngineArrayOption="-a ARRAY_JOBS%20" gridOptions="--partition=batch --nodes=1 --time=24:00:00"

canu -p Po_2021 -d canu_2021 genomeSize=584m \
-trim -corrected -pacbio "canu_2021/Po_2021.correctedReads.fasta.gz" \
correctedErrorRate=0.105 gridOptions="--partition=batch --nodes=1 --time=72:00:00"

canu -p Po_2021 -d canu_2021 genomeSize=584m \
-assemble -trimmed -corrected -pacbio "canu_2021/Po_2021.trimmedReads.fasta.gz" \
correctedErrorRate=0.105 batMemory=9 ovbMemory=16 ovsMemory=16 executiveMemory=1 cnsMemory=20 cnsThreads=8 \
gridOptions="--partition=batch --nodes=1 --time=24:00:00" "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50" utgovlMemory=30
```

**Step 6. Polishing**

After contig assembly, I polished the genome with clean raw data. As far as I am aware PacBio tools accept only files generated from their sequencer or processed using their tools. This means I cannot use clean PacBio reads in fastq format. I do not want chloroplast and mithochondrial reads that were still present in PacBio raw reads polished the assembled contigs. To prevent this, we needed to filter original reads (native bam files).

This is how I did it:

Creating index files
```
for i in *.subreads.bam
do
pbindex $i
done
```

Filtering reads
```
### creating xml file from all PacBio raw read
ls *.subreads.bam > mymovies.fofn
dataset create --type SubreadSet --name Plantago PlantagoGenomeSet.subreadset.xml mymovies.fofn

### listing clean reads from fastq file
grep ’@’ Plantago_pacbio_no_mito_chloro.fastq > PlantagoGenome.txt
sed 's|[@,]||g' PlantagoGenome.txt > PlantagoGenome_final.txt

### filtering PacBio reads using list of clean reads
dataset filter PlantagoGenomeSet.subreadset.xml Plantago_filter.subreadset.xml 'qname=PlantagoGenome_final.txt'
```


```
###
pbmm2 align --log-level INFO --log-file pbmm2_log --sample Plantago /hpcfs/users/a1697274/canu_2021/Po_2021.assembled.unassembled.fasta.mmi Plantago_filter.subreadset.xml Plantago.aligned.bam

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

mkdir bam_0 bam_1 bam_2 bam_3 bam_4 bam_5 bam_6 bam_7
mv tig00007* bam_7
mv tig00006* bam_6
mv tig00005* bam_5
mv tig00004* bam_4
mv tig00003* bam_3
mv tig00002* bam_2
mv tig00001* bam_1
mv tig0000* bam_0
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
mkdir fasta_0 fasta_1 fasta_2 fasta_3 fasta_4 fasta_5 fasta_6 fasta_7
mv tig00007* fasta_7
mv tig00006* fasta_6
mv tig00005* fasta_5
mv tig00004* fasta_4
mv tig00003* fasta_3
mv tig00002* fasta_2
mv tig00001* fasta_1
mv tig0000* fasta_0

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

mkdir fasta_0 fasta_1 fasta_2 fasta_3 fasta_4 fasta_5 fasta_6 fasta_7
mv tig00007* fasta_7
mv tig00006* fasta_6
mv tig00005* fasta_5
mv tig00004* fasta_4
mv tig00003* fasta_3
mv tig00002* fasta_2
mv tig00001* fasta_1
mv tig0000* fasta_0

```
```
awk '{ print $1":0-"$2 }' Po_new_new.contigs.fasta.gz.fai > window.txt
window_7.txt
window_6.txt
window_5.txt
window_4.txt
window_3.txt
window_2.txt
window_1.txt
window_0.txt

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



cat file_0/*.fasta > polished_0.fasta
cat file_1/*.fasta > polished_1.fasta
cat file_2/*.fasta > polished_2.fasta
cat file_3/*.fasta > polished_3.fasta
cat file_4/*.fasta > polished_4.fasta
cat file_5/*.fasta > polished_5.fasta
cat file_6/*.fasta > polished_6.fasta
cat file_7/*.fasta > polished_7.fasta

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
