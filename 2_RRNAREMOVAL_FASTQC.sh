#!/bin/bash
#SBATCH -p cpu	            	                            
#SBATCH -N 1               	                                
#SBATCH -n 8              	                                
#SBATCH --time=3:00:00    	                                
#SBATCH --mem=32GB         	                                
#SBATCH --array=1-12

source RNASEQ_CONFIG.sh

## This step is useful if there are significant rRNA sequences in some 
## libraries but not others, probably due to bad library prep. This 
## might result in false-positive DE genes later. To remove rRNA 
## sequences from the trimmed reads, we align the reads to the known 
## rRNA sequences in zebrafish and write unaligned reads (the non-rRNA 
## ones) to paired-end .fastq format. 

cd ${TRIMDATA}/fastq
FIRSTREAD=`ls ${SLURM_ARRAY_TASK_ID}_*_R1_T1.fastq.gz`
SECONDREAD=`ls ${SLURM_ARRAY_TASK_ID}_*_R2_T2.fastq.gz`
OUTPUTSAM=${RRRDATA}/sams/${FIRSTREAD%%_R1_T1.fastq.gz}.sam
OUTPUTINFO=${RRRDATA}/info/${FIRSTREAD%%_R1_T1.fastq.gz}.sam.info
OUTPUTFASTQ=${RRRDATA}/fastq/${FIRSTREAD%%_R1_T1.fastq.gz}_RR%.fastq.gz

hisat2 --known-splicesite-infile ${z10_splicesites} \
-x ${z10_rRNAindex} \
-p 16 \
-1 ${FIRSTREAD} -2 ${SECONDREAD} \
-S ${OUTPUTSAM} &> ${OUTPUTINFO} \
--un-conc-gz ${OUTPUTFASTQ}


# Check to see if GC content has improved via Fastqc

cd ${RRRDATA}/fastq
FIRSTREAD=`ls ${SLURM_ARRAY_TASK_ID}_*_RR1.fastq.gz`
SECONDREAD=`ls ${SLURM_ARRAY_TASK_ID}_*_RR2.fastq.gz`

fastqc -t 16 -f fastq -o ${RRRDATA}/fastqc ${FIRSTREAD}
fastqc -t 16 -f fastq -o ${RRRDATA}/fastqc ${SECONDREAD}
