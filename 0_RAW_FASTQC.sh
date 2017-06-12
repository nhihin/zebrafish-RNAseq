#!/bin/bash
#SBATCH -p cpu	            	                            
#SBATCH -N 1               	                                
#SBATCH -n 8              	                                
#SBATCH --time=0:40:00    	                                
#SBATCH --mem=32GB         	                                
#SBATCH --array=1-12

source RNASEQ_CONFIG.sh

cd ${RAWDATA}/fastq
FIRSTREAD=`ls ${SLURM_ARRAY_TASK_ID}_*_R1.fastq.gz`
SECONDREAD=`ls ${SLURM_ARRAY_TASK_ID}_*_R2.fastq.gz`

fastqc -t ${threads} -f fastq -o ${RAWDATA}/fastqc ${FIRSTREAD}
fastqc -t ${threads} -f fastq -o ${RAWDATA}/fastqc ${SECONDREAD}