#!/bin/bash
#SBATCH -p batch            	                            
#SBATCH -N 1               	                                
#SBATCH -n 8              	                               
#SBATCH --time=3:00:00    	                                
#SBATCH --mem=32GB       
#SBATCH --array=1-12	  

# Run Salmon to quantify gene expression at the transcript level.                               

source RNASEQ_CONFIG.sh

cd ${RRRDATA}/fastq

FIRSTREAD=`ls ${SLURM_ARRAY_TASK_ID}_*_RR1.fastq.gz`
SECONDREAD=`ls ${SLURM_ARRAY_TASK_ID}_*_RR2.fastq.gz`
OUTPUTCOUNTS=$QUANTDATA_TRANSCRIPT/${FIRSTREAD%%_RR1.fastq.gz}_transcripts_quant

$salmon quant -i ${z10_salmonindex} \
-p ${threads} \
-l A \
-1 ${FIRSTREAD} \
-2 ${SECONDREAD} \
-o ${OUTPUTCOUNTS} \
--numBootstraps 100
