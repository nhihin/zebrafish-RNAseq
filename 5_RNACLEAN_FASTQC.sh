#!/bin/bash
#SBATCH -p batch            	                                
#SBATCH -N 1               	                                
#SBATCH -n 8              	                                
#SBATCH --time=2:00:00    	                                
#SBATCH --mem=64GB         	                                
#SBATCH --array=1-12

source RNASEQ_CONFIG.sh

# The following script must be run ONCE in console to generate one R script 
# for each .bam file. Scripts were saved in a folder in my home directory 
# on Phoenix. 
# mkdir -pv ~/RNASEQ/5_RNACLEAN_FASTQC_SLURM/scripts/
# SCRIPTDIR=~/RNASEQ/5_RNACLEAN_FASTQC_SLURM/scripts/
# source ~/RNASEQ/RNASEQ_CONFIG.sh
# cd $DEDUPDATA/bams
# for uniquebam in *_unique.bam
# do 
# 	echo "library(rnaCleanR); library(Rsamtools); filterOnePairs(bamfilein=\"$DEDUPDATA/bams/$uniquebam\", bamfileout=\"${CLEANDATA}/bams/${uniquebam%%.bam}_cleanR.bam\", statfile=\"${CLEANDATA}/info/${uniquebam%%.bam}_cleanR.txt\", readLength=150)" > $SCRIPTDIR/${uniquebam}_cleanRScript.R
# done

# Run Hien's R package rnaCleanR for each .bam file. 
# This removes sequences corresponding to DNA contamination. 
# Because the $SLURM_ARRAY_TASK_ID environmental variable
# cannot be used in an R script easily , I prepared separate 
# R scripts for each .bam file. See the comments above.

cd ~/RNASEQ/5_RNACLEAN_FASTQC_SLURM/scripts 

CLEANRSCRIPT=`ls ${SLURM_ARRAY_TASK_ID}_*.R`

# Pass the script to R, which has already been loaded
# from the RNASEQ_CONFIG.sh file.
R < ${CLEANRSCRIPT} --no-save

# Get the name from the unique, deduplicated bams.
cd ${DEDUPDATA}/bams
UNIQUEBAM=`ls ${SLURM_ARRAY_TASK_ID}_*_unique.bam`
# Specify the name of the clean bam file. 
cd ${CLEANDATA}/bams
CLEANBAM=${UNIQUEBAM%%.bam}_cleanR.bam
# Run FastQC on the DNA-contamination removed .bam files. 
fastqc -t ${threads} -f bam_mapped -o ${CLEANDATA}/fastqc ${CLEANBAM}

# Index the cleaned bam
java -jar /data/biohub/local/picard-2.9.0.jar BuildBamIndex \
      I=${CLEANBAM}
      O=${CLEANBAM}.bam.bai

# R
# library(rnaStrand)
# bamfilein = "your bam file"


# plotFileWin = win.pdf
# plotFileHist = hist.pdf

# win <- getWinPairs(bamfilein)
# plotWin(win,save=TRUE,file=plotFileWin)
# plotHist(win,save=TRUE,file=plotFileHist)





