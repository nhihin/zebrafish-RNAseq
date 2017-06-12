#!/bin/bash
#SBATCH -p batch            	                                
#SBATCH -N 1               	                                
#SBATCH -n 8              	                                
#SBATCH --time=2:00:00    	                                
#SBATCH --mem=32GB         	                                
#SBATCH --array=1-12

source RNASEQ_CONFIG.sh

cd ${ALIGNDATA}/bams

SORTEDBAM=`ls ${SLURM_ARRAY_TASK_ID}_*_sorted.bam`
DEDUPBAM=${DEDUPDATA}/bams/${SORTEDBAM%%.bam}_dedup.bam
DEDUPINFO=${DEDUPDATA}/info/${SORTEDBAM%%.bam}_dedup.txt

# Run picard to detect and remove PCR + optical duplicates and write 
# the non-duplicated reads to a bam file.
# Note: If only PCR duplicates should be removed, replace the 
# `REMOVE_DUPLICATES=TRUE` with `REMOVE_SEQUENCING_DUPLICATES=TRUE`.
java -jar /data/biohub/local/picard-2.9.0.jar MarkDuplicates \
      I=${SORTEDBAM} \
      O=${DEDUPBAM} \
      M=${DEDUPINFO} \
      REMOVE_DUPLICATES=TRUE

# Run Samtools to filter the deduplicated bams for unique reads.
# The -1 is used for fast bam compression with 16 threads.
# The -F 1036 flag is used to filter out unmapped reads/mates along
# with any reads marked as duplicates but not removed. 
# The -q 10 option filters out reads with MAPQ<10 indicating they
# are likely to be not unique. 
UNIQUEBAM=${DEDUPBAM%%.bam}_unique.bam

samtools view -b -1 -@ 16 -F 1036 -q 10 ${DEDUPBAM} > ${UNIQUEBAM}

# Build index for the deduplicated + filtered bam files
java -jar /data/biohub/local/picard-2.9.0.jar BuildBamIndex \
      I=${UNIQUEBAM}
      O=${UNIQUEBAM}.bam.bai

# Run FastQC on the deduplicated + filtered bam files, only on the mapped reads.
fastqc -t ${threads} -f bam_mapped -o ${DEDUPDATA}/fastqc ${UNIQUEBAM}