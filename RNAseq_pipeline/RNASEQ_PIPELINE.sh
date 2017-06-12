#!/bin/bash

#########################################################################################
# Paired-End RNA-seq Pipeline for Phoenix HPC 
#########################################################################################

## IMPORTANT:
## This script (RNASEQ_PIPELINE.sh) should be run in console and NOT submitted to SLURM.

## DESCRIPTION:
## The pipeline consists of several steps (see below). Each step is contained in a 
## separate bash script which will be run for each paired-end RNA-seq library in
## parallel using Phoenix's SLURM Job Arrays. When all RNA-seq libraries finish that 
## particular step, the next step will begin. If a job fails (e.g. not enough time, 
## missing file, etc.) the job will not have an exit status of COMPLETED on SLURM, 
## and so the subsequent jobs will not run. When this happens, simply use sbatch in 
## console to submit the jobs manually. 

## Run RNASEQ_CONFIG.sh to (1) import the variables for data directories, zebrafish 
## references, indexes and programs, (2) rename the raw paired-end .fastq data for
## compatibility with Phoenix's SLURM Job Arrays, which allows this script to be
## automatically run once per paired-end data in parallel, and (3) specify the
## number of threads to be used for each job. 
source RNASEQ_CONFIG.sh

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. RAW DATA QUALITY CHECK								
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
RAW=`sbatch 0_RAW_FASTQC.sh`
RAW_SLURM_ID=$(echo "$RAW" | sed 's/Submitted batch job //')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. TRIMMING + QUALITY CHECK
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
TRIM=`sbatch --dependency=afterok:${RAW_SLURM_ID} 1_TRIM_FASTQC.sh`
TRIM=`sbatch 1_TRIM_FASTQC.sh`
TRIM_SLURM_ID=$(echo "$TRIM" | sed 's/Submitted batch job //')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. rRNA REMOVAL 									
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
RRNA=`sbatch --dependency=afterok:${TRIM_SLURM_ID} 2_RRNAREMOVAL_FASTQC.sh`
RRNA_SLURM_ID=$(echo "$RRNA" | sed 's/Submitted batch job //')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. ALIGN rRNA-REMOVED READS TO REFERENCE GENOME
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ALIGN=`sbatch --dependency=afterok:${RRNA_SLURM_ID} 3_ALIGN_SAMTOBAM_INDEX.sh`
ALIGN_SLURM_ID=$(echo "$ALIGN" | sed 's/Submitted batch job //')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. GENE EXPRESSION QUANTIFICATION AT TRANSCRIPT LEVEL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
QUANT_TR=`sbatch --dependency=afterok:${RRNA_SLURM_ID} 3_SALMON.sh`
QUANT_TR_SLURM_ID=$(echo "$QUANT_TR" | sed 's/Submitted batch job //')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. DEDUPLICATION OF ALIGNED BAM FILES 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
DEDUP=`sbatch --dependency=afterok:${ALIGN_SLURM_ID} 4_DEDUPLICATION_UNIQUE_INDEX.sh`
DEDUP_SLURM_ID=$(echo "$DEDUP" | sed 's/Submitted batch job //')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. REMOVING DNA CONTAMINATION FROM DEDUPLICATED BAM FILES	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
CLEAN=`sbatch --dependency=afterok:${DEDUP_SLURM_ID} 5_RNACLEAN_FASTQC.sh`
CLEAN_SLURM_ID=$(echo "$CLEAN" | sed 's/Submitted batch job //')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. GENE EXPRESSION QUANTIFICATION AT GENE LEVEL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
QUANT_GENE=`sbatch --dependency=afterok:${CLEAN_SLURM_ID} 6_FEATURECOUNTS.sh`
QUANT_GENE_SLURM_ID=$(echo "$CLEAN" | sed 's/Submitted batch job //')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SLURM JOB IDs	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SLURM Job IDs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "0 - Quality check raw data: $RAW_SLURM_ID"
echo "1 - Trim and filter raw data: $TRIM_SLURM_ID"
echo "2 - rRNA removal: $RRNA_SLURM_ID"
echo "3 - Align to reference genome: $ALIGN_SLURM_ID"
echo "3 - Transcript quantification: $QUANT_TR_SLURM_ID"
echo "4 - Deduplication of aligned reads: $DEDUP_SLURM_ID"
echo "5 - Remove DNA contamination from aligned reads: $CLEAN_SLURM_ID"
echo "6 - Gene quantification: $QUANT_GENE_SLURM_ID"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

