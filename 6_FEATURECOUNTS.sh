#!/bin/bash
#SBATCH -p batch            	                                
#SBATCH -N 1               	                                
#SBATCH -n 8              	                                
#SBATCH --time=3:00:00    	                                
#SBATCH --mem=32GB         	                                

source RNASEQ_CONFIG.sh

# Run FeatureCounts to quantify gene-counts in the bam files.
# Only needs to be run once.
cd ${CLEANDATA}/bams

sampleList=`find $CLEANDATA/bams -name "*_cleanR.bam" | tr '\n' ' '`
$featureCounts -Q 10 -s 2 -T ${threads} -p -a $z10_gtf -o $QUANTDATA_GENE/project_genes.out ${sampleList}
cut -f1,7- $QUANTDATA_GENE/project_genes.out | sed 1d > $QUANTDATA_GENE/project_genes.txt





#samp_list=`find $ALIGNDATA -name "*.sorted.bam" | tr '\n' ' '`

#featureCounts -Q 10 -s 1 -T $threads -a $gff -o $QUANT1/project_genes.out $samp_list
#cut -f1,7- $QUANT1/project_genes.out | sed 1d > $QUANT1/project_genes.txt

# cd $CLEANDATA
# CLEANBAM=`ls ${SLURM_ARRAY_TASK_ID}_*_cleanR.bam`
# GENECOUNTS=$QUANTDATA_GENE/${CLEANBAM%%_sorted_dedup_unique_cleanR.bam}.out $CLEANBAM


# $featureCounts -Q 10 -s 2 -T 16 -p -a $z10_gtf -o $QUANTDATA_GENE/

# cd $CLEANDATA


# threads=16
# z10_gtf=/data/biohub/Refs/zebrafish/Danio_rerio.GRCz10.88.gtf
# featureCounts=/data/biohub/local/subread-1.5.2-Linux-x86_64/bin/featureCounts