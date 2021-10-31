# ***Plantago ovata* contig assembly**

Continuous long reads (CLRs) from the PacBio platform (~76X coverage) were used to assemble the *P. ovata* genome.

Main steps in generating contigs:
1. Genomic DNA extraction and PacBio CLR sequencing
2. Installing bioinformatic tools
3. Converting PacBio unaligned bam files into fastq files
4. Removing unwanted reads
5. Genome size prediction
6. Contig assembly
7. Polishing
8. Purging and clipping

**Step 1. Genomic DNA extraction and sequencing processes were explained in the publication (link).**

**Step 2. Installing software**

I used conda to install all tools and created several environments due to software incompatibility.
An example of how to create a conda environment:
```
conda create -n pacbio
conda activate pacbio
conda install -c bioconda bam2fastx  
```
List of main tools for contig assembly:
- Pacbio tools (bam2fastx v1.3.1, pbbam v1.6.0, pbcommand v2.1.1, pbcopper v1.9.1, pbcore v2.1.2, pbcoretools v0.8.1, pbgcpp v1.0.0, pbmm2 v1.4.0, pbzip2 v1.1.13) :   https://github.com/PacificBiosciences/pbbioconda
- tabix v0.2.6 : https://github.com/samtools/tabix
- canu v2.1.1 : https://github.com/marbl/canu 
- minimap2 v2.17 : https://github.com/lh3/minimap2
- samtools v1.11 :https://github.com/samtools/samtools
- bedtools v2.29.2 : https://github.com/arq5x/bedtools2 or https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html
- fastqc v0.11.9 : https://github.com/s-andrews/FastQC

More tools can be found in Supplementary File Table 6 (link)

**Step 3. Converting PacBio unaligned bam files into fastq files**


Raw data has been deposited in SRA NCBI : SRR14643405

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

### if you have many files, you can do this:

for i in *.bam
do
outfile=$(basename ${i} .bam)
bam2fastq -c 9 ${i} -o ${outfile}
done


### concatenating and compressing all fastq files
cat *.subreads.fastq > Plantago_pacbio.fastq
bgzip -c -l 9 Plantago_pacbio.fastq > Plantago_pacbio.fastq.gz
```

**Step 4. Removing unwanted reads**

I found that removing contaminants from PacBio raw reads helped me solve my contig assembly problem using Canu. We are interested in nuclear genome, so chloroplast and mitochondrial reads are considered unwanted reads. *Plantago* chloroplast genome can be found at https://www.ncbi.nlm.nih.gov/nuccore/MH205737.1/ and a mitochondrial gene is in here https://www.ncbi.nlm.nih.gov/nuccore/EU069524.1/. Only one mitochondrial gene was found in the NCBI database. The mitochondrial genome is still not available (October 2021). 


Creating index file
```
minimap2 -d plantago_chloroplast.fasta.gz.mmi plantago_chloroplast.fasta.gz
minimap2 -d plantago_mitocondria.fasta.gz.mmi plantago_mitocondria.fasta.gz
```
Removing chroloplast reads
```
data="Plantago_pacbio.fastq.gz"
index="plantago_chloroplast.fasta.gz.mmi"
output="Plantago_pacbio_no_chloro.bam"

time minimap2 \
-ax map-pb \
-t 2 $index $data \
| samtools view -f 0x04 -u \
| samtools sort \
--threads 2 -l 7 \
-o $output

### Converting a bam file to fastq file then compressing it
bamToFastq -i Plantago_pacbio_no_chloro.bam -fq Plantago_pacbio_no_chloro.fastq
bgzip -c -l 9 Plantago_pacbio_no_chloro.fastq > Plantago_pacbio_no_chloro.fastq.gz
```

Removing mitochondrial reads
```
data="Plantago_pacbio_no_chloro.fastq.gz"
index="plantago_mitocondria.fasta.gz.mmi"
output="Plantago_pacbio_no_mito_chloro.bam"

time minimap2 \
-ax map-pb \
-t 2 $index $data \
| samtools view -f 0x04 -u \
| samtools sort \
--threads 2 -l 7 \
-o $output

### Converting a bam file to fastq file then compressing it
bamToFastq -i Plantago_pacbio_no_mito_chloro.bam -fq Plantago_pacbio_no_mito_chloro.fastq
bgzip -c -l 9 Plantago_pacbio_no_mito_chloro.fastq > Plantago_pacbio_no_mito_chloro.fastq.gz
```

**Step 5. Genome size prediction**

To predict Plantago ovata genome size, I utilized publicly short read genomic Illumina data (SRR10076762) using genomescope2 (https://github.com/tbenavi1/genomescope2.0).

```
jellyfish count -C -m 21 -s 1000000000 -t 10 *.fastq -o reads_21.jf
jellyfish count -C -m 31 -s 3000000000 -t 20 *.fastq -o reads_31.jf
jellyfish count -C -m 39 -s 3000000000 -t 10 *.fastq -o reads_39.jf
jellyfish count -C -m 45 -s 3000000000 -t 20 *.fastq -o reads_45.jf
jellyfish count -C -m 50 -s 3000000000 -t 20 *.fastq -o reads_50.jf
jellyfish count -C -m 70 -s 3000000000 -t 30 *.fastq -o reads_70.jf
jellyfish count -C -m 100 -s 3000000000 -t 50 *.fastq -o reads_100.jf

jellyfish histo -t 10 reads_21.jf > reads_21.histo
jellyfish histo -t 10 reads_31.jf > reads_31.histo
jellyfish histo -t 10 reads_39.jf > reads_39.histo
jellyfish histo -t 10 reads_45.jf > reads_45.histo
jellyfish histo -t 10 reads_50.jf > reads_50.histo
jellyfish histo -t 10 reads_70.jf > reads_70.histo
jellyfish histo -t 10 reads_100.jf > reads_100.histo

Rscript genomescope.R reads_21.histo 21 150 genomescope_plantago_21
Rscript genomescope.R reads_31.histo 31 150 genomescope_plantago_31
Rscript genomescope.R reads_39.histo 39 150 genomescope_plantago_39
Rscript genomescope.R reads_45.histo 45 150 genomescope_plantago_45
Rscript genomescope.R reads_50.histo 50 150 genomescope_plantago_50
Rscript genomescope.R reads_70.histo 70 150 genomescope_plantago_70
Rscript genomescope.R reads_100.histo 100 150 genomescope_plantago_100
```

**Step 6. Contig assembly**

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

**Step 7. Polishing**

After contig assembly, I polished the genome with clean raw data. As far as I am aware PacBio tools accept only files generated from their sequencer or processed using their tools. This means I cannot use clean PacBio reads in fastq format. I do not want chloroplast and mithochondrial reads still present in PacBio raw reads to polish the assembled contigs. To prevent this, we needed to filter original reads (native bam files).

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

Generating a bam file for each contig
```
### aligning or mapping filtered reads to draft contig assembly
pbmm2 align --log-level INFO --log-file pbmm2_log -j 5 --sample Plantago canu_2021/Po_2021.contigs.fasta.mmi \
Plantago_filter.subreadset.xml Plantago.aligned.bam

### sorting and indexing
samtools sort -m 10G --threads 5 -l 7 -o Plantago_sorted_aligned.bam -T tmp.ali Plantago.aligned.bam
samtools index -b Plantago_sorted_aligned.bam Plantago_sorted_aligned.bam.bai

### getting conting names
samtools faidx Po_2021.contigs.fasta > Po_2021.contigs.fasta.fai
awk '{ print $1 }' Po_2021.contigs.fasta.fai > contigs.txt


### splitting a bam file by contigs
for line in `cat contigs.txt`
do
samtools view -bh Plantago_sorted_aligned.bam ${line} > split/${line}.bam
done

### indexing each bam file
for i in split/*.bam
do
samtools index -b $i $i.bai
done

### creating new directories and distributing bam files to the new folders
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

Generating a fasta file for each contig
```
### splitting unpolished assembled genome by contigs

### modifying from https://raw.githubusercontent.com/harish0201/General_Scripts/master/Fasta_splitter.py

#!/usr/bin/env python
import sys
"""
usage: python Fasta_splitter.py Po_2021.contigs.fasta 2> /dev/null
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

### indexing each fasta file using samtools
for i in *.fasta
do
samtools faidx $i
done

### creating bed files
for i in *.fasta.fai
do
OUTPUT=$(basename "$i" .fasta.fai).bed
awk -F'\t' '{ print $1,$2}' $i > $OUTPUT
done


### indexing fasta files using a PacBio tool
N=100
(
for file in *.fasta
do
((i=i%N)); ((i++==0)) && wait
OUTPUT=$(basename "$file").mmi
pbmm2 index $file $OUTPUT &
done
)

### creating new directories and distributing fasta files to the new folders

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

Parallelising polishing step
```

### getting contig names and their starts and ends (e.g tig00000002:0-23132)
awk '{ print $1":0-"$2 }' Po_2021.contigs.fasta.fai > window.txt
grep "tig00007" window.txt > window_7.txt
grep "tig00006" window.txt > window_6.txt
grep "tig00005" window.txt > window_5.txt
grep "tig00004" window.txt > window_4.txt
grep "tig00003" window.txt > window_3.txt
grep "tig00002" window.txt > window_2.txt
grep "tig00001" window.txt > window_1.txt
grep "tig00000" window.txt > window_0.txt


### polishing step by running individual script (8 scripts) below:
### script 1
for window in `cat window_0.txt`
do
line=$(echo ${window} | cut -c 1-11)
echo "gcpp --max-iterations 4 --log-level INFO --log-file polished_seqs/file_0/${line}.log \
-w ${window} -r fasta_split/fasta_0/${line}.fasta \
-o polished_seqs/file_0/${line}.polished.fastq,polished_seqs/file_0/${line}.polished.fasta,polished_seqs/file_0/${line}.polished.gff \   
split/bam_0/${line}.bam"
done | parallel -j7 --tmpdir TMPDIR_1 --compress

### script 2
for window in `cat window_1.txt`
do
line=$(echo ${window} | cut -c 1-11)
echo "gcpp --max-iterations 4 --log-level INFO --log-file polished_seqs/file_1/${line}.log \
-w ${window} -r fasta_split/fasta_1/${line}.fasta \
-o polished_seqs/file_1/${line}.polished.fastq,polished_seqs/file_1/${line}.polished.fasta,polished_seqs/file_1/${line}.polished.gff \   
split/bam_1/${line}.bam"
done | parallel -j7 --tmpdir TMPDIR_1 --compress

### script 3
for window in `cat window_2.txt`
do
line=$(echo ${window} | cut -c 1-11)
echo "gcpp --max-iterations 4 --log-level INFO --log-file polished_seqs/file_2/${line}.log \
-w ${window} -r fasta_split/fasta_2/${line}.fasta \
-o polished_seqs/file_2/${line}.polished.fastq,polished_seqs/file_2/${line}.polished.fasta,polished_seqs/file_2/${line}.polished.gff \   
split/bam_2/${line}.bam"
done | parallel -j7 --tmpdir TMPDIR_2 --compress

### script 4
for window in `cat window_3.txt`
do
line=$(echo ${window} | cut -c 1-11)
echo "gcpp --max-iterations 4 --log-level INFO --log-file polished_seqs/file_3/${line}.log \
-w ${window} -r fasta_split/fasta_3/${line}.fasta \
-o polished_seqs/file_3/${line}.polished.fastq,polished_seqs/file_3/${line}.polished.fasta,polished_seqs/file_3/${line}.polished.gff \   
split/bam_3/${line}.bam"
done | parallel -j7 --tmpdir TMPDIR_3 --compress

### script 5
for window in `cat window_4.txt`
do
line=$(echo ${window} | cut -c 1-11)
echo "gcpp --max-iterations 4 --log-level INFO --log-file polished_seqs/file_4/${line}.log \
-w ${window} -r fasta_split/fasta_4/${line}.fasta \
-o polished_seqs/file_4/${line}.polished.fastq,polished_seqs/file_4/${line}.polished.fasta,polished_seqs/file_4/${line}.polished.gff \   
split/bam_4/${line}.bam"
done | parallel -j7 --tmpdir TMPDIR_4 --compress

### script 6
for window in `cat window_5.txt`
do
line=$(echo ${window} | cut -c 1-11)
echo "gcpp --max-iterations 4 --log-level INFO --log-file polished_seqs/file_5/${line}.log \
-w ${window} -r fasta_split/fasta_5/${line}.fasta \
-o polished_seqs/file_5/${line}.polished.fastq,polished_seqs/file_5/${line}.polished.fasta,polished_seqs/file_5/${line}.polished.gff \   
split/bam_5/${line}.bam"
done | parallel -j7 --tmpdir TMPDIR_4 --compress

### script 7
for window in `cat window_6.txt`
do
line=$(echo ${window} | cut -c 1-11)
echo "gcpp --max-iterations 4 --log-level INFO --log-file polished_seqs/file_6/${line}.log \
-w ${window} -r fasta_split/fasta_6/${line}.fasta \
-o polished_seqs/file_6/${line}.polished.fastq,polished_seqs/file_6/${line}.polished.fasta,polished_seqs/file_6/${line}.polished.gff \   
split/bam_6/${line}.bam"
done | parallel -j7 --tmpdir TMPDIR_4 --compress

### script 8
for window in `cat window_7.txt`
do
line=$(echo ${window} | cut -c 1-11)
echo "gcpp --max-iterations 4 --log-level INFO --log-file polished_seqs/file_7/${line}.log \
-w ${window} -r fasta_split/fasta_7/${line}.fasta \
-o polished_seqs/file_7/${line}.polished.fastq,polished_seqs/file_7/${line}.polished.fasta,polished_seqs/file_7/${line}.polished.gff \   
split/bam_7/${line}.bam"
done | parallel -j7 --tmpdir TMPDIR_4 --compress


### concatenating and compressing
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


**Step 8. Purging haplotigs (alternative contigs)**
```
### indexing polished assembly
minimap2 -d polished.fasta.mmi polished.fasta

### mapping clean read to polished assembly
minimap2 -t 4 -ax map-pb polished.fasta.mmi Plantago_pacbio_no_mito_chloro.fasta --secondary=no \
| samtools sort -m 1G -o polished_aligned.bam -T polished_tmp.ali

### purging and clipping
purge_haplotigs  hist  -b polished_aligned.bam -g polished.fasta  -t 4 -d 200

purge_haplotigs cov -i ./polished_aligned.bam.gencov -l 5 -m 70 -h 190 -o polished_coverage_stats.csv -j 80 -s 80

purge_haplotigs purge -g polished.fasta -c polishde_coverage_stats.csv -t 4 -r repeat.bed -o purge_polished

purge_haplotigs  clip  -p purge_polished.fasta -h purge_polished.haplotigs.fasta

### renaming and compressing file
mv clip.fasta Plantago.fasta
bgzip -c -l 9 Plantago.fasta > Plantago.fasta.gz
```
