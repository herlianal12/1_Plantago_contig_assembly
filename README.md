# Plantago contig assembly
Main steps in generating contigs
1. Genomic DNA extraction
2. PacBio CLR sequencing
3. Removing contaminant
3. Contig assembly
4. Polishing
5. Purging and clipping haplotig





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
