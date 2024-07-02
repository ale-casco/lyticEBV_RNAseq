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

# Set directory containing RSEM counts files
counts_dir=file.path(work_dir,"Data/RSEM")

# Set directory for sampleinfo matrix file
sampleinfo_dir=file.path(work_dir,"Data/refs/sampleinfo.txt") # 3-column file: id, condition, replicate

# Set directory for gene info file
geneinfo_dir = file.path(work_dir, "Data/refs/gene_info.csv") # see gtf2gene_info.R script
if (!file.exists(geneinfo_dir)) {stop(" gene_info.csv file does not exist.\n\tRun gtf2gene_info.R script to generate.", call. = FALSE)} # check file exists

# Set desired comparisons
contrasts_list=list(c("WT_Latent","dBGLF5_Latent"))

## Differential Gene Expression Analysis

# Create working directory
out_dir=file.path(work_dir, "Data/DESeq2")
dir.create(out_dir, showWarnings=FALSE)

for (i in 1:length(contrasts_list)) {
  contrasts = contrasts_list[i]
  
  # Import sampleInfo.txt file
  sampleinfo=read.delim(sampleinfo_dir)
  sampleinfo=sampleinfo[grepl(paste(paste0("\\b",unique(unlist(contrasts)),"\\b"), collapse="|"), sampleinfo[[2]]),] # isolates samples included in contrasts
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
  contrasts.tmp=contrasts
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
 
  # Get normalized counts data
  ERCC_normCounts=as.data.frame(counts(dds_ERCC, normalized=TRUE))
  
  # Strip gene versions
  rownames(ERCC_normCounts) = stripGeneVersions(rownames(ERCC_normCounts))
  
  ## Create output table with normalized sample expression values used to generate DGE table
  ERCC_output_table=ERCC_normCounts
  
  ## Iterate through Wald Tests to generate pairwise comparisons of all groups with apeglm shrinkage
  # Test alterantive hypothesis
  LFC_coef <- c("condition",contrasts[1,1],contrasts[2,1])
  dds_ERCC$condition <- relevel(dds_ERCC$condition, ref = LFC_coef[[3]])
  dds_ERCC <- nbinomWaldTest(dds_ERCC, betaPrior = FALSE)
  res <- results(dds_ERCC, contrast=LFC_coef, lfcThreshold = 1, altHypothesis = "lessAbs")
  res <- lfcShrink(dds_ERCC, coef=which(resultsNames(dds_ERCC) == paste("condition_", LFC_coef[[2]], "_vs_", LFC_coef[[3]], sep ="")), 
                   res=res, type = "apeglm")
  res <- as.data.frame(res)[,c(2,4,5)]
  rownames(res) <- 1:NROW(res)
  ERCC_output_table<-cbind(ERCC_output_table,res)
  rm(res)

  # Add gene symbol annotations from GTF file (generated separately)
  gene_info$ENSEMBL = ifelse(is.na(gene_info$ENSEMBL), gene_info$SYMBOL, gene_info$ENSEMBL)
  ERCC_output_table = merge(x = gene_info, y = ERCC_output_table, by.x = "ENSEMBL", by.y = 0, all.y = TRUE)
  ERCC_output_table$SYMBOL = ifelse(!grepl("ENSG",ERCC_output_table$ENSEMBL), ERCC_output_table$ENSEMBL, ERCC_output_table$SYMBOL)
  ERCC_output_table$ENSEMBL = ifelse(!grepl("ENSG",ERCC_output_table$ENSEMBL), NA, ERCC_output_table$ENSEMBL)
  
  ERCC_output_table[[paste0(contrasts[2,1],"_groupMean")]] = rowMeans(ERCC_output_table[,grepl(paste0(contrasts[2,1],"_c"),colnames(ERCC_output_table))])
  ERCC_output_table[[paste0(contrasts[1,1],"_groupMean")]] = rowMeans(ERCC_output_table[,grepl(paste0(contrasts[1,1],"_c"),colnames(ERCC_output_table))])
  ERCC_output_table$baseMean = rowMeans(ERCC_output_table[,grepl("_c",colnames(ERCC_output_table))])
  
  # Export DGE table
  write.csv(ERCC_output_table,file.path(out_dir, paste0("AH_ERCCnorm_DE.",LFC_coef[[2]], ".vs.", LFC_coef[[3]],".csv")), row.names=FALSE)

}

## print session info ##
print("Session Info below: ")
sessionInfo()

