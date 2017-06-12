#!/bin/bash
#SBATCH -p batch	            	                            
#SBATCH -N 1               	                                
#SBATCH -n 8              	                                
#SBATCH --time=2:00:00    	                                
#SBATCH --mem=32GB 
#SBATCH --array=1-12

# Use Hisat2 to align the rRNA-removed fastq files to the reference
# zebrafish genome, writing out the aligned .sam file and log file.

source RNASEQ_CONFIG.sh

cd ${RRRDATA}/fastq
FIRSTREAD=`ls ${SLURM_ARRAY_TASK_ID}_*_RR1.fastq.gz`
SECONDREAD=`ls ${SLURM_ARRAY_TASK_ID}_*_RR2.fastq.gz`
ALIGNEDSAM=${ALIGNDATA}/sams/${FIRSTREAD%%_RR1.fastq.gz}.sam
ALIGNEDINFO=${ALIGNDATA}/info/${FIRSTREAD%%_RR1.fastq.gz}.sam.info

# Run HISAT2
hisat2 --known-splicesite-infile ${z10_splicesites} \
-x ${z10_hisat2index} \
-p ${threads} \
-1 ${FIRSTREAD} -2 ${SECONDREAD} \
-S ${ALIGNEDSAM} &> ${ALIGNEDINFO}


# Convert the aligned sam files to bam and then sort by coordinate.
SORTEDBAM=${ALIGNDATA}/bams/${FIRSTREAD%%_RR1.fastq.gz}.bam

java -jar /data/biohub/local/picard-2.9.0.jar SortSam \
      I=${ALIGNEDSAM} \
      O=${SORTEDBAM} \
      SORT_ORDER=coordinate

# Index the sorted bam files, should only take ~3 mins.
java -jar /data/biohub/local/picard-2.9.0.jar BuildBamIndex \
      I=${SORTEDBAM}
      O=${SORTEDBAM}.bam.bai

# Run FastQC on the sorted .bam files, only on the mapped reads.
fastqc -t ${threads} -f bam_mapped -o ${ALIGNDATA}/fastqc ${SORTEDBAM}
