rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("tidyverse")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("tximport")
#BiocManager::install("DESeq2")
#BiocManager::install("AcidGenomes")
#BiocManager::install("apeglm")

## Setup the environment

# Import packages
library(rstudioapi)
library(tximport) # imports data from RSEM
library(DESeq2) # DGE analysis
library(tidyverse) # formatting data
library(AcidGenomes) # strip gene versions
library(apeglm) # shrinkage

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# ARTDeco readthrough counts
counts_dir=file.path(work_dir,"Data/ARTDeco/dogs/all_dogs.raw.txt")
counts_fpkm_dir=file.path(work_dir,"Data/ARTDeco/dogs/all_dogs.fpkm.txt")

# DESeq2 differential expressed gene (DEG) files
deg_dir=file.path(work_dir,"Data/DESeq2")

# Set directory for sampleinfo matrix file
sampleinfo_dir=file.path(work_dir,"Data/refs/sampleinfo.txt") # 3-column file: id, condition, replicate

# Set directory for gene info file
geneinfo_dir=file.path(work_dir,"Data/refs/gene_info.csv") # see gtf2gene_info.R script
if (!file.exists(geneinfo_dir)) {stop(" gene_info.csv file does not exist.\n\tRun gtf2gene_info.R script to generate.", call. = FALSE)} # check file exists

# RSEM counts
rsem_dir=file.path(work_dir,"Data/RSEM")

# Comparisons to keep
deg_files = list.files(deg_dir, full.names=TRUE)
deg_files = deg_files[grepl("WT_Early.vs.WT_Latent|WT_Late.vs.WT_Latent|dBGLF5_Early.vs.dBGLF5_Latent|dBGLF5_Late.vs.dBGLF5_Latent",deg_files)]

# Prepare contrasts
contrast_cols = sub("\\.csv","",sub("ERCCnorm_DE\\.","",basename(deg_files)))
contrast_1 = sub("\\.vs\\..*","",contrast_cols)
contrast_2 = sub(".*\\.vs\\.","",contrast_cols)
contrasts = data.frame(matrix(ncol = length(contrast_cols), nrow = 0))
contrasts = rbind(contrasts,contrast_1)
contrasts = rbind(contrasts,contrast_2)
colnames(contrasts) = contrast_cols
rownames(contrasts)= c("Treatment", "Control")
contrasts

## Differential Gene Expression Analysis

# Create working directory
out_dir=file.path(work_dir, "ARTDeco_DESeq2/readthrough")
dir.create(out_dir, showWarnings=FALSE, recursive = TRUE)
for (i in 1:length(contrasts)) {
  contrast = contrasts[i]
  
  # Import DEG file
  deg = read.csv(deg_files[grepl(colnames(contrasts)[i], deg_files)])
  
  # Import sampleInfo.txt file
  sampleinfo=read.delim(sampleinfo_dir)
  sampleinfo=sampleinfo[grepl(paste(paste0("\\b",unique(unlist(contrast)),"\\b"), collapse="|"), sampleinfo[[2]]),] # isolates samples included in contrast
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
  
  # Format contrast table, defining pairwise comparisons for all groups
  contrast.tmp=contrast
  for (i in 1:length(contrast.tmp)) {
    if (i == 1) {
      contrast=combn(((unlist(contrast.tmp[i]))),2)
      contrast.names=paste("(", paste(contrast[,i], collapse=")v("), ")", sep="")
      
    } else {
      contrast=cbind(contrast, combn(((unlist(contrast.tmp[i]))),2))
      contrast.names=cbind(contrast.names, paste("(", paste(contrast[,i], collapse=")v("), ")", sep=""))
    }
  }
  contrast.names=as.character(contrast.names)
  colnames(contrast)=contrast.names
  rm(contrast.tmp, contrast.names)
  contrast # print comparisons
  
  # Import RSEM raw gene counts data
  files=list.files(file.path(rsem_dir),pattern=".genes.results", full.names=TRUE)
  
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
  
  # Filter files for desired contrast
  files=files[grep(paste(unique(as.vector(contrast)), collapse="|"), names(files))]
  
  # Import files with tximport as counts
  txi.rsem=tximport(files, type="rsem", txIn=FALSE, txOut=FALSE)
  ercc_counts = txi.rsem$abundance
  ercc_counts = ercc_counts[grepl("ERCC-",rownames(ercc_counts)),]
  
  # Import readthrough matrix
  countData <- read.delim(counts_dir)
  rownames(countData) = countData[,1]
  countData = countData[,-1:-2]
  colnames(countData) = sub("\\..*","",colnames(countData))
  rownames(countData) = stripGeneVersions(rownames(countData))
  countData = countData[,colnames(countData) %in% sampleinfo[,1]]
  colnames(countData)
  
  # Import FPKM matrix
  FPKM <- read.delim(counts_fpkm_dir)
  rownames(FPKM) = FPKM[,1]
  FPKM = FPKM[,-1:-2]
  colnames(FPKM) = sub("\\..*","",colnames(FPKM))
  rownames(FPKM) = stripGeneVersions(rownames(FPKM))
  FPKM = FPKM[,colnames(FPKM) %in% sampleinfo[,1]]
  colnames(FPKM)
  
  # Order
  FPKM = FPKM[match(rownames(countData), rownames(FPKM)),]
  
  # Filter
  keep = rowSums(FPKM >= 0.2) >= min(table(study)) & rowSums(countData >= 5) >= min(table(study))
  countData = countData[keep,]
  
  # Reorder columns according to sampleinfo
  countData = countData[,sampleinfo[[1]]]
  colnames(countData) = paste(sampleinfo[[2]], sampleinfo[[3]], sep = "_")
  
  # Combine readthrough/readin counts and ERCC gene counts
  ercc_counts = round(ercc_counts,0)
  ercc_counts = ercc_counts[rownames(ercc_counts) %in% deg$SYMBOL,]
  countData = rbind(countData, ercc_counts)
  
  # Get only human and ERCC genes
  countData = countData[grepl("ENSG|ERCC-",rownames(countData)),]
  
  ## Make DESeqDataSet object
  
  # Create data frame defining which group each sample belongs to
  sampleTable=data.frame(condition=factor(group))
  
  ## Round data
  countData = round(countData)
  
  # Filter based on DESeq2
  deg = deg[,1:5]
  deg = deg[!is.na(deg$ENSEMBL) | grepl("ERCC-",deg$SYMBOL),]
  deg$ENSEMBL = ifelse(is.na(deg$ENSEMBL), deg$SYMBOL, deg$ENSEMBL)
  #countData = countData[rownames(countData) %in% deg$ENSEMBL,]
  
  # Create DESeqDataSet Object (dds)
  dds <- DESeqDataSetFromMatrix(countData=countData, 
                                colData=sampleTable, 
                                design=~condition, tidy = FALSE)
  
  
  ## Make a DESeqDataSet object using only filtered ERCC genes, which will be used to generate ERCC counts table
  ## Create list of rows containing ERCC group B genes to use for ERCC-normalization
  ## Note: ERCC genes should be the same concentration in all samples
  ercc_rows <- grep("ERCC-",rownames(dds))
  
  ## Perform DESeq analysis
  
  ## Run DESeq analysis with ERCC-normalization by replacing size factor object with ERCC size factors for rescaling
  ## This includes normalization, dispersion estimation, linear model fitting.
  dds_ERCC <- estimateSizeFactors(dds, controlGenes=ercc_rows)
  dds_ERCC <- estimateDispersions(dds_ERCC)
  dds_ERCC <- nbinomWaldTest(dds_ERCC)
  
  ## Run DESeq analysis without considering ERCC genes
  ## This includes normalization, dispersion estimation, linear model fitting.
  dds <- DESeq(dds)
  
  # Get normalized counts data
  ERCC_normCounts=as.data.frame(counts(dds_ERCC, normalized=TRUE))
  normCounts=as.data.frame(counts(dds, normalized=TRUE))
  
  # Strip gene versions
  rownames(ERCC_normCounts) = stripGeneVersions(rownames(ERCC_normCounts))
  rownames(normCounts) = stripGeneVersions(rownames(normCounts))
  
  ## Create output table with normalized sample expression values used to generate DGE table
  ERCC_output_table=ERCC_normCounts
  output_table=normCounts
  
  ## Iterate through Wald Tests to generate pairwise comparisons of all groups with apeglm shrinkage
  # In the Wald test, the estimated standard error of a log2 fold change is used
  #    to compare the differences between two groups. This generates the adjusted p values
  LFC_coef <- c("condition",contrast[1,1],contrast[2,1])
  dds_ERCC$condition <- relevel(dds_ERCC$condition, ref = LFC_coef[[3]])
  dds_ERCC <- nbinomWaldTest(dds_ERCC)
  res <- results(dds_ERCC, contrast=LFC_coef)
  res <- lfcShrink(dds_ERCC, coef=which(resultsNames(dds_ERCC) == paste("condition_", LFC_coef[[2]], "_vs_", LFC_coef[[3]], sep ="")), 
                   res=res, type = "apeglm")
  res <- as.data.frame(res)[,c(2,4,5)]
  rownames(res) <- 1:NROW(res)
  ERCC_output_table<-cbind(ERCC_output_table,res)
  
  rm(res)
  
  # Add gene symbol annotations from GTF file (generated separately)
  ERCC_output_table = merge(x = deg, y = ERCC_output_table, by.x = "ENSEMBL", by.y = 0, all.y = TRUE)
  ERCC_output_table$SYMBOL = ifelse(!grepl("ENSG",ERCC_output_table$ENSEMBL), ERCC_output_table$ENSEMBL, ERCC_output_table$SYMBOL)
  ERCC_output_table$ENSEMBL = ifelse(!grepl("ENSG",ERCC_output_table$ENSEMBL), NA, ERCC_output_table$ENSEMBL)
  
  ERCC_output_table[[paste0(contrast[2,1],"_groupMean")]] = rowMeans(ERCC_output_table[,grepl(paste0(contrast[2,1],"_c"),colnames(ERCC_output_table))])
  ERCC_output_table[[paste0(contrast[1,1],"_groupMean")]] = rowMeans(ERCC_output_table[,grepl(paste0(contrast[1,1],"_c"),colnames(ERCC_output_table))])
  ERCC_output_table$baseMean = rowMeans(ERCC_output_table[,grepl("_c",colnames(ERCC_output_table))])
  
  # Export DGE table
  write.csv(ERCC_output_table,file.path(out_dir, paste0("ERCCnorm_DE.",LFC_coef[[2]], ".vs.", LFC_coef[[3]],".csv")), row.names=FALSE)
  
}

## print session info ##
print("Session Info below: ")
sessionInfo()

