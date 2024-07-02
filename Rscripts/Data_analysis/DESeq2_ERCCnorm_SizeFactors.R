rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("tidyverse")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("tximport")
#BiocManager::install("DESeq2")
#BiocManager::install("AcidGenomes")

## Setup the environment

# Import packages
library(rstudioapi)
library(tximport) # imports data from RSEM
library(DESeq2) # DGE analysis
library(tidyverse) # formatting data
library(AcidGenomes) # strip gene versions

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Set directory containing RSEM counts files
counts_dir=file.path(work_dir,"Data/RSEM")

# Set directory for sampleinfo matrix file
sampleinfo_dir=file.path(work_dir,"Data/refs/sampleinfo.txt") # 3-column file: id, condition, replicate

# Set directory for gene info file
geneinfo_dir=file.path(work_dir,"Data/refs/gene_info.csv") # see gtf2gene_info.R script
if (!file.exists(geneinfo_dir)) {stop(" gene_info.csv file does not exist.\n\tRun gtf2gene_info.R script to generate.", call. = FALSE)} # check file exists

# Set desired comparisons
contrasts_list=list(c("WT_Bulk","WT_Latent"),
                    c("WT_Early","WT_Latent"),c("WT_Late","WT_Latent"), c("WT_Late","WT_Early"),
                    c("dBGLF5_Early","dBGLF5_Latent"),c("dBGLF5_Late","dBGLF5_Latent"), c("dBGLF5_Late","dBGLF5_Early"),
                    c("WT_Early","dBGLF5_Early"),c("WT_Late","dBGLF5_Late"))

## Differential Gene Expression Analysis

# Create output directory
out_dir=file.path(work_dir, "Data/DESeq2")
dir.create(out_dir, showWarnings=FALSE)

# Import sampleInfo.txt file
sampleinfo=read.delim(sampleinfo_dir)
sampleinfo=sampleinfo[grep(paste(unique(unlist(contrasts_list)), collapse="|"), sampleinfo[[2]]),] # isolates samples included in contrasts
sampleinfo # print

# Create data frame containing all samples and respective factors
study=sampleinfo[2]
rownames(study)=sampleinfo[[1]]

if (dim(study)[1] >= 2){
  group=apply(study,1,paste,collapse=" & ") # concatenate multiple factors into one condition per sample
} else{
  group=study[,1]
}

group_names=paste0("(",group,")",sep="") # human readable group names
group=make.names(group) # group naming compatible with R models
names(group)=group_names
rm(group_names)

# Format contrasts table, defining pairwise comparisons for all groups
contrasts.tmp=contrasts_list
for (i in 1:length(contrasts.tmp)) {
  if (i == 1) {
    contrasts=combn(((unlist(contrasts.tmp[i]))),2)
    contrast.names=paste("(", paste(contrasts[,i], collapse=")v("), ")", sep="")
    
  } else {
    contrasts=cbind(contrasts, combn(((unlist(contrasts.tmp[i]))),2))
    contrast.names=cbind(contrast.names, paste("(", paste(contrasts[,i], collapse=")v("), ")", sep=""))
  }
}
contrast.names=as.character(contrast.names)
colnames(contrasts)=contrast.names
rm(contrasts.tmp, contrast.names)
contrasts # print comparisons

# Import RSEM raw gene counts data
files=list.files(file.path(counts_dir),pattern=".genes.results", full.names=TRUE)

# Reorder list according to sampleinfo
temp_list=as.list(sampleinfo[[1]])
temp_files=c()
for (i in as.character(temp_list)) {
  temp_files=append(temp_files, c(print(files[grep(i, files)])))
}
files=temp_files
rm(temp_list, temp_files)

# Name Samples
names(files)=paste0(sampleinfo[[2]], "_", sampleinfo[[3]])

# Filter files for desired contrasts
files=files[grep(paste(unique(as.vector(contrasts)), collapse="|"), names(files))]

# Import files with tximport as counts
txi.rsem=tximport(files, type="rsem", txIn=FALSE, txOut=FALSE)

# Import files with tximport as FPKM for filtering
txi.rsem_2 <- tximport(files, type = "none", txIn = FALSE, txOut = FALSE, geneIdCol = "gene_id", abundanceCol = "FPKM",
                       countsCol = "expected_count", lengthCol = "effective_length", importer = function(x) read_tsv(x))
rm(files)

# Add 1 to genes with lengths of zero - needed to make DESeqDataSet object 
txi.rsem$length[txi.rsem$length == 0]=1

## Make DESeqDataSet object

# Create data frame defining which group each sample belongs to
sampleTable=data.frame(condition=factor(group))
rownames(sampleTable)=colnames(txi.rsem$counts)

# Create DESeqDataSet Object (dds)
dds=DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)
rm(sampleTable)

## Filter out genes with FPKM < 1 or < 5 reads across more than the smallest group's sample size
# This makes the pipeline fast by reducing the size of the object in memory.
# This also increases the power of multiple testing correciton.
# rowSums(FPKM) >= x) >= y), where x is 1 FPKM and 5 counts and y is the smallest group's sample size. Similar to recommendation from edgeR and limma.

# Also, filter by gene_type
gene_info = read.csv(geneinfo_dir)

FPKM <- txi.rsem_2$abundance # FPKM
rownames(FPKM) = stripGeneVersions(rownames(FPKM))
keep <- rowSums(FPKM >= 1) >= min(table(study)) &  rowSums(counts(dds) >= 5) >= min(table(study)) & (!grepl("ENSG", rownames(FPKM)) | rownames(FPKM) %in% gene_info$ENSEMBL[gene_info$gene_type == "protein_coding"])
dds=dds[keep,]
rm(txi.rsem,txi.rsem_2,FPKM,study)

## Perform DESeq analysis

# Create list of rows containing ERCC group B genes to use for ERCC-normalization
# Note: ERCC genes should be at the same concentration across comparisons
ercc_rows <- grep("ERCC-",rownames(dds))
ercc_dds <- dds[ercc_rows,]

# Run DESeq analysis with ERCC-normalization by replacing size factor object with ERCC size factors for rescaling
# This includes normalization, dispersion estimation, linear model fitting.
dds_ERCC <- estimateSizeFactors(dds, controlGenes=ercc_rows)
dds_ERCC <- estimateDispersions(dds_ERCC)
dds_ERCC <- nbinomWaldTest(dds_ERCC)

# Calculate size factors for IGV and export as sample matrix file
ERCC_size_factors=1/colMeans(normalizationFactors(dds_ERCC))
sampleinfo$ERCC_norm_factors = ERCC_size_factors

# Export samplinfo file with ERCC normalization factors
write.table(sampleinfo,file.path(out_dir, 'sampleinfo_IGV.txt'), row.names = F, quote = F, sep = '\t')

## print session info ##
print("Session Info below: ")
sessionInfo()
