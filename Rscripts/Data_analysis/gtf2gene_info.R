rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("rtracklayer")
#BiocManager::install("GenomicFeatures")
#BiocManager::install("AcidGenomes")

## Setup the environment

# Import packages
library(rstudioapi)
library(rtracklayer)
library(GenomicFeatures)
library(AcidGenomes)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# GTF directory
gtf_dir=file.path(work_dir,"Data/refs/gencode.v45.primary_assembly.ERCC.M81_DFLR.chrEBV.inverted.gtf.gz")

# Import GTF file using rtracklayer
gtf <- import(gtf_dir)

# Get ENSEMBL, SYMBOL, gene_type, strand, and chromosome information
ENSEMBL = data.frame(gtf$gene_id)
SYMBOLS = data.frame(gtf$gene_name)
gene_type = data.frame(gtf$gene_type)
strand = data.frame(strand(gtf))
chrom = data.frame(chrom(gtf))

# Combine information into data frame
gtf_df = cbind(ENSEMBL,SYMBOLS,gene_type,strand,chrom)
colnames(gtf_df) = c("ENSEMBL","SYMBOL","gene_type","strand","chrom")

# Strip ENSEMBL gene versions and remove duplicates
gtf_df$ENSEMBL = stripGeneVersions(gtf_df$ENSEMBL)
gtf_df = gtf_df[!duplicated(gtf_df$ENSEMBL),]

# Export
write.csv(gtf_df,file.path(dirname(gtf_dir),"gene_info.csv"), row.names = FALSE)

## print session info ##
print("Session Info below: ")
sessionInfo()
