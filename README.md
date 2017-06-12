# RNA-seq analysis of zebrafish models of aging and Alzheimer's disease.

## Data Description

Total RNA-seq was done on whole zebrafish brains. The biological factors being varied include `sibship`, `presenilin1 genotype` and `age`. The experimental design is summarised below:

![experimental design for zebrafish experiment](http://i.imgur.com/lJFc7fU.png)

## Completed so far

- [x] Create an RNA-seq pipeline to clean and process raw paired-end reads, align them to the reference zebrafish genome, and quantify gene expression at the gene and transcript levels. 
- [x] Apply generalised linear models (using the R/Bioconductor package **limma**) to perform differential gene expression analysis on the dataset and identify differentially expressed genes for various comparisons, to determine the effect of `presenilin1 genotype` and `age` on gene expression in the zebrafish brain. 
- [ ] Network analysis
- [ ] Comparison with MCI and Alzheimer's disease in RNA-seq and microarray human data
