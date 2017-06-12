#!/bin/bash
#SBATCH -p cpu	            	                            
#SBATCH -N 1               	                                
#SBATCH -n 8              	                                
#SBATCH --time=2:00:00    	                                
#SBATCH --mem=32GB         	                             
#SBATCH --array=1-12  		

source RNASEQ_CONFIG.sh

cd ${RAWDATA}/fastq					

FIRSTREAD=`ls ${SLURM_ARRAY_TASK_ID}_*_R1.fastq.gz`
SECONDREAD=`ls ${SLURM_ARRAY_TASK_ID}_*_R2.fastq.gz`
FIRSTTRIM=${TRIMDATA}/fastq/${FIRSTREAD%%.fastq.gz}_T1.fastq.gz
SECONDTRIM=${TRIMDATA}/fastq/${SECONDREAD%%.fastq.gz}_T2.fastq.gz

$AdapterRemoval --file1 ${FIRSTREAD} \
--file2 ${SECONDREAD} \
--output1 ${FIRSTTRIM} \
--output2 ${SECONDTRIM} \
--threads 16 --gzip \
--trimqualities --trimns --minquality 20 --minlength 35

fastqc -t 8 -f fastq -o ${TRIMDATA}/fastqc ${FIRSTTRIM}
fastqc -t 8 -f fastq -o ${TRIMDATA}/fastqc ${SECONDTRIM}