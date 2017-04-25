#!/bin/bash -l

## Adapted from Robinson Research Institute RNAseq pipeline by Jimmy Breen (jimmybreen@gmail.com)
## https://github.com/jimmybgammyknee/gammytools/blob/master/RRI_RNAseq_pipeline.sh


## Requires:
## AdapterRemoval for trimming raw reads
## Salmon for transcript-level quantification of trimmed reads
## STAR for aligning to genome
## StringTie for transcript assembly
## featureCounts for gene-level quantification of trimmed reads

# Setting Up 
## Module loading
. /opt/shared/Modules/3.2.7/init/bash
module load gnu/4.9.2
module load java/java-1.7.09
module load fastQC/0.11.2
module load parallel

## Directory setup
base=/data/biohub/2017_Lardelli_Q96K97del
RAWDATA=$base/rawData
TRIMDATA=$base/trimmedData
QUANTDATA_TRANSCRIPT=$base/quantifiedData/salmon
ALIGNDATA=$base/alignedData
QUANTDATA_GENE=$base/quantifiedData/featureCounts
ASSEMDATA=$base/assembledData
QC_RAW=$RAWDATA/fastqc_reports
QC_TRIMMED=$TRIMDATA/fastqc_reports

SCRIPTS=$base/scripts
ALIGN_SLURM_HEADER=$base/scripts/starAlignment.sh

## Programs (location on Phoenix)
AdapterRemoval=/data/biohub/local/adapterremoval/build/AdapterRemoval
salmon=/data/biohub/local/Salmon-0.8.2_linux_x86_64/bin/salmon
featureCounts=/data/biohub/local/subread-1.5.2-Linux-x86_64/bin/featureCounts
stringtie=/data/biohub/local/stringtie-1.3.3b.Linux_x86_64/stringtie
module load STAR/2.5.1a-foss-2015b
module load SAMtools/1.3.1-GCC-5.3.0-binutils-2.25 
module load fastqc/0.11.4

## Threads (usually use one thread per sample)
threads=32

## Transcript references for zebrafish 
#### The following downloaded from ENSEMBL: http://www.ensembl.org/info/data/ftp/index.html
z10_dna=/data/biohub/Refs/zebrafish/Danio_rerio.GRCz10.dna.toplevel.fa.gz
z10_gtf=/data/biohub/Refs/zebrafish/Danio_rerio.GRCz10.88.gtf
z10_cdna=/data/biohub/Refs/zebrafish/Danio_rerio.GRCz10.cdna.all.fa.gz

## Indexes for zebrafish 
#### STAR index was built using the following command:
#### STAR --runThreadN 32 --runMode genomeGenerate --genomeDir /data/biohub/Refs/zebrafish/STAR_index --genomeFastaFiles /data/biohub/Refs/zebrafish/Danio_rerio.GRCz10.dna.toplevel.fa --sjdbGTFfile /data/biohub/Refs/zebrafish/Danio_rerio.GRCz10.88.gtf --sjdbOverhang 100
z10_star=/data/biohub/Refs/zebrafish/STAR_index
#### salmon index was built using the following command:
#### /data/biohub/local/Salmon-0.8.2_linux_x86_64/bin/salmon index -t /data/biohub/Refs/zebrafish/Danio_rerio.GRCz10.cdna.all.fa.gz -i /data/biohub/Refs/zebrafish/salmon_index
z10_salmon=/data/biohub/Refs/zebrafish/salmon_index




# 0. FastQC for raw data
mkdir -p $QC_RAW
cd $RAWDATA
fastqc *.fastq.gz -o $QC_RAW/

# 1.1. Adapter trimming with AdapterRemoval
## Note: paired-end files should be named like *_R1.fastq.gz and *_R2.fastq.gz
mkdir -p $TRIMDATA
cd $RAWDATA
for firstpair in *R1.fastq.gz
do
	pairID="${firstpair%%_R1*}"
	$AdapterRemoval --file1 ${pairID}_R1.fastq.gz \
	--file2 ${pairID}_R2.fastq.gz \
	--output1 $TRIMDATA/${pairID}_trimmed_1.fastq.gz \
	--output2 $TRIMDATA/${pairID}_trimmed_2.fastq.gz \
	--threads $threads --gzip \
	--trimqualities --trimns --minquality 10 --minlength 25
done

## 1.2. FastQC for trimmed data
mkdir -p $QC_TRIMMED
cd $TRIMDATA
fastqc *.fastq.gz -o $QC_TRIMMED

# 2.1. Transcript level quantification from trimmed fastq with salmon using their quasi-alignment mode
mkdir $QUANTDATA_TRANSCRIPT
cd $TRIMDATA
for firstpair in *_trimmed_1.fastq.gz
do
  pairID="${firstpair%%_trimmed_1*}"
  $salmon quant -i $z10_salmon \
  -l A \ 
  -1 ${pairID}_trimmed_1.fastq.gz \
  -2 ${pairID}_trimmed_2.fastq.gz \
  -o $QUANTDATA_TRANSCRIPT/${pairID}_transcripts_quant \
  --numBootstraps 100
  --numThreads $threads
done

# 2.2. Aligning trimmed reads to reference genome using STAR
mkdir -p $ALIGNDATA
cd $TRIMDATA
for firstpair in *trimmed_1.fastq.gz
do
	pairID="${firstpair%%_trimmed_1*}"
	STAR --runThreadN $threads \
	--genomeDir $z10_star \
	--readFilesIn ${pairID}_trimmed_1.fastq.gz ${pairID}_trimmed_2.fastq.gz \
	--readFilesCommand gunzip -c \
	--outSAMtype BAM Unsorted SortedByCoordinate \
	--outBAMcompression 5 \
	--outFileNamePrefix $ALIGNDATA/${pairID}_
done

## Index sorted .bam files
cd $ALIGNDATA
for bam in *Aligned.sortedByCoord.out.bam
do
#	samtools sort ${file} ${file}.sorted   #Note: no longer needed since STAR also produces .bam files that are sorted by coordinate, apparently the same as samtools-sorted bam files. These just need to be indexed.
	samtools index $bam
done

# 3.1. Gene level quantification of sorted and indexed bam files using featureCounts
cd $ALIGNDATA

sampleList=`find $ALIGNDATA -name "*Aligned.sortedByCoord.out.bam" | tr '\n' ' '`
$featureCounts -Q 10 -s 2 -T $threads -p -a $z10_gtf -o $QUANTDATA_GENE/project_genes1.out $sampleList
cut -f1,7- $QUANTDATA_GENE/project_genes1.out | sed 1d > $QUANTDATA_GENE/project_genes1.txt


# 3.2. Transcript assembly of each sorted bam file with stringtie
mkdir -p $ASSEMDATA
cd $ALIGNDATA
for bam in *Aligned.sortedByCoord.out.bam
do 
	$stringtie -p $threads -m 30 -G $z10_gtf \
	-o $ASSEMDATA/${bam%%Aligned.sortedByCoord.out.bam}assembly.gtf \
	-A $ASSEMDATA/${bam%%Aligned.sortedByCoord.out.bam}gene_abund.tab \
	-C $ASSEMDATA/${bam%%Aligned.sortedByCoord.out.bam}cov_refs.gtf $bam
done

## Merge the assemblies using stringtie's merge mode
assemblies_to_merge=$(ls $ASSEMDATA/*assembly.gtf)
$stringtie --merge -G $z10_gtf -o $ASSEMDATA/merged.gtf $assemblies_to_merge
awk '{if($3=="transcript")print}' $ASSEMDATA/merged.gtf > $ASSEMDATA/merged.transcript.gtf




