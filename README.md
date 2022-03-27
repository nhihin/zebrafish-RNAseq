# RNA-seq pipeline for processing of zebrafish RNA-seq data

- This is a collection of bash scripts for processing paired-end RNA-seq data on a SLURM HPC environment. The scripts use SLURM array jobs to parallelise running multiple samples through at the same time. 

## SLURM array jobs

- The raw .fastq.gz files need to be named such that there is a number in front of the filename. e.g. `1_SAMPLE_r1.fastq.gz`, `1_SAMPLE_r2.fastq.gz`, `2_SAMPLE_r1.fastq.gz`, `2_SAMPLE_r2.fastq.gz` ... `12_SAMPLE_r1.fastq.gz`, `12_SAMPLE_r2.fastq.gz`. 

- In the SLURM header, the `#SBATCH --array=1-12` controls how many samples to run the script with. 

## Pipeline steps

- Define relevant parameters and directories in `RNASEQ_CONFIG.sh`
- Run pipeline using `RNASEQ_PIPELINE.sh`.

1. Quality check on raw fastq.gz files `0_RAW_FASTQC.sh`
2. Adapter and quality trimming + quality check `1_TRIM_FASTQC.sh`
3. Removal of incompletely removed rRNA `2_RRNAREMOVAL_FASTQC.sh` (The rRNA sequences used for rRNA depletion for these total RNA libraries was not perfect, resulting in rRNA sequences dominating the downstream differential expression analysis). This step aligns the fastq's against known zebrafish rRNA sequences to try to remove these. 
4. Alignment to reference zebrafish genome assembly (`3_ALIGN_SAMTOBAM_INDEX.sh`) and also running transcript-level quantification (`3_SALMON.sh). 
5. Deduplication based on unique indexes `4_DEDUPLICATION_UNIQUE_INDEX.sh`
6. Run R package `rnaCleanR` to remove small amount of DNA contaminants `5_RNACLEAN_FASTQC.sh`
7. Gene expression quantification `6_FEATURECOUNTS.sh`


