#########################################################################################
# Paired-End RNA-seq Pipeline for Phoenix HPC 
#########################################################################################

threads=16

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Setup 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
DATA=/data/biohub/2017_Lardelli_K97Gfs #All data here
LOCAL=/data/biohub/local #Local programs
REFS=/data/biohub/Refs/zebrafish #Change if not using zebrafish.

#---------------------------------------------------------------------------------------#
# Data Directories 
#---------------------------------------------------------------------------------------#
# Note: All raw paired-end fastq files should be in $RAWDATA/fastq and named like 
# `*_R1.fastq.gz` and `*_R2.fastq.gz`.

# Specify data directory names.
RAWDATA=$DATA/0_rawData
TRIMDATA=$DATA/1_trimmedData
RRRDATA=$DATA/2_rRNAremovedData
ALIGNDATA=$DATA/3_alignedData
DEDUPDATA=$DATA/4_dedupData
CLEANDATA=$DATA/5_cleanRData
QUANTDATA=$DATA/6_quantifiedData
QUANTDATA_GENE=$QUANTDATA/geneLevel
QUANTDATA_TRANSCRIPT=$QUANTDATA/transcriptLevel
ASSEMDATA=$DATA/7_assembledData

# Make directories if they do not already exist. 
DIRS=($TRIMDATA $RRRDATA $ALIGNDATA $DEDUPDATA $CLEANDATA \
	$QUANTDATA $ASSEMDATA)
for dir in ${DIRS[*]}; do mkdir -pv $dir; done

#---------------------------------------------------------------------------------------#
# Rename raw data 									#
#---------------------------------------------------------------------------------------#
# Rename paired-end raw fastq files to have an ID number in front. This is required to 
# use the SLURM Job Arrays to run jobs in parallel. This following bit is only run once.
# More info on Job Arrays: https://slurm.schedmd.com/job_array.html

################## NOTE: To list the files in numeric order, use ls -v ##################

cd ${RAWDATA}/fastq
if [ ! -d "$RAWDATA/fastqc" ]; then
	ID=1
	for firstread in *_R1.fastq.gz
	do
		mv firstread ${ID}_${firstread}
		ID=$(expr $ID + 1)
	done

	ID=1
	for secondread in *_R2.fastq.gz
	do
		mv $secondread ${ID}_${secondread}
		ID=$(expr $ID + 1)
	done
fi

#---------------------------------------------------------------------------------------#
# Create additional directories to organise data 
#---------------------------------------------------------------------------------------#
# 0. Quality checking of raw data
mkdir -pv $RAWDATA/fastqc

# 1. Trim and filter raw data, then check quality of trimmed reads. 
mkdir -pv $TRIMDATA/fastq
mkdir -pv $TRIMDATA/fastqc

# 2. Align trimmed/filtered reads to known zebrafish rRNA sequences and write all 
# unaligned reads (what we want) to /fastqc. Check quality of reads.
# /sams contains rRNA sequences and can be deleted to free space. 
mkdir -pv $RRRDATA/fastq
mkdir -pv $RRRDATA/fastqc
mkdir -pv $RRRDATA/sams
mkdir -pv $RRRDATA/info

# 3. Align rRNA-free reads to reference zebrafish genome, convert the .sam files to 
# .bam files, then check quality of aligned reads. /sams can be deleted to free space.
mkdir -pv $ALIGNDATA/sams
mkdir -pv $ALIGNDATA/bams
mkdir -pv $ALIGNDATA/info
mkdir -pv $ALIGNDATA/fastqc

# 4. Deduplicate aligned reads, then check quality of aligned reads.
mkdir -pv $DEDUPDATA/bams
mkdir -pv $DEDUPDATA/fastqc
mkdir -pv $DEDUPDATA/info

# 5. Remove aligned reads corresponding to DNA contamination, and then check quality 
# of aligned reads. 
mkdir -pv $CLEANDATA/bams
mkdir -pv $CLEANDATA/fastqc
mkdir -pv $CLEANDATA/info

# 6. Quantification of gene expression at gene and transcript levels.
mkdir -pv $QUANTDATA_GENE
mkdir -pv $QUANTDATA_TRANSCRIPT

#---------------------------------------------------------------------------------------#
# Zebrafish References and Indexes 		
#---------------------------------------------------------------------------------------#
# Downloaded from: http://www.ensembl.org/info/data/ftp/index.html
z10_gtf=$REFS/Danio_rerio.GRCz10.88.gtf
z10_gff=$REFS/Danio_rerio.GRCz10.88.gff3
z10_cdna=$REFS/Danio_rerio.GRCz10.cdna.all.fa
z10_dna=$REFS/Danio_rerio.GRCz10.dna.toplevel.fa

# Splice sites file for HISAT2, generated using:
# cd /data/biohub/Refs/zebrafish;
# python hisat2_extract_splice_sites.py \
# Danio_rerio.GRCz10.88.gtf > Danio_rerio.GRCz10.88_splice_sites.txt
z10_splicesites=$REFS/Danio_rerio.GRCz10.88_splice_sites.txt 

# HISAT2 index, generated using:
# cd /data/biohub/Refs/zebrafish; mkdir hisat2_index;
# hisat2-build Danio_rerio.GRCz10.dna.toplevel.fa \
# hisat2_index/Danio_rerio.GRCz10.dna_run.toplevel
z10_hisat2index=$REFS/hisat2_index/Danio_rerio.GRCz10.dna_run.toplevel

# HISAT2 index using the zebrafish rRNA sequences, generated using:
# cd /data/biohub/Refs/zebrafish; mkdir rRNA_hisat2_index;
# hisat2-build Danio_rerio_rRNA_sequences_SILVA_20170524.fa \
# rRNA_hisat2_index/Danio_rerio_rRNA_sequences_SILVA_20170524_hisat2index
z10_rRNAindex=$REFS/rRNA_hisat2_index/Danio_rerio_rRNA_sequences_SILVA_20170524_hisat2index

# Salmon index is required to use its quasi-mapping mode. Index was generated using:
# /data/biohub/local/Salmon-0.8.2_linux_x86_64/bin/salmon index \
# -t /data/biohub/Refs/zebrafish/Danio_rerio.GRCz10.cdna.all.fa \
# -i /data/biohub/Refs/zebrafish/salmon_index --type quasi -k 31
z10_salmonindex=$REFS/salmon_index

#---------------------------------------------------------------------------------------#
# Programs 								
#---------------------------------------------------------------------------------------#
module load Java/1.8.0_121
module load fastqc/0.11.4
module load GCC/6.3.0-2.27
module load HISAT2/2.0.5-foss-2017a
module load SAMtools/1.3.1-GCC-5.3.0-binutils-2.25
module load R/3.3.0-foss-2016uofa
AdapterRemoval=$LOCAL/adapterremoval/build/AdapterRemoval
featureCounts=$LOCAL/subread-1.5.2-Linux-x86_64/bin/featureCounts
salmon=$LOCAL/Salmon-0.8.2_linux_x86_64/bin/salmon



