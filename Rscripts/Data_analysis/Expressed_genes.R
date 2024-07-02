rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("openxlsx")

## Setup the environment

# Import packages
library(rstudioapi)
library(tximport) # imports data from RSEM
library(tidyverse) # formatting data
library(AcidGenomes) # strip gene versions

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Set directory containing RSEM counts files
counts_dir=file.path(work_dir, "Data/RSEM")

# Set directory for sampleinfo matrix file
sampleinfo_dir=file.path(work_dir,"Data/refs/sampleinfo.txt") # 3-column file: id, condition, replicate

# Set directory for gene info file
geneinfo_dir=file.path(work_dir,"Data/refs/gene_info.csv") # see gtf2gene_info.R script
if (!file.exists(geneinfo_dir)) {stop(" gene_info.csv file does not exist.\n\tRun gtf2gene_info.R script to generate.", call. = FALSE)} # check file exists

## Differential Gene Expression Analysis

# Create working directory
out_dir=file.path(work_dir, "Expressed_Genes")
dir.create(out_dir, showWarnings=FALSE)

# Import sampleInfo.txt file
sampleinfo=read.delim(sampleinfo_dir)

expressed_genes_list = list()
for (i in 1:length(unique(sampleinfo[[2]]))) {
  # Sample to subset
  sample = unique(sampleinfo[[2]])[i]
  
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
  
  # Filter files for desired samples
  
  files=files[grep(paste0(sample,"_c"), names(files))]
  
  # Import files with tximport as counts
  txi.rsem=tximport(files, type="rsem", txIn=FALSE, txOut=FALSE)
  
  # Import files with tximport as FPKM for filtering
  txi.rsem_2 <- tximport(files, type = "none", txIn = FALSE, txOut = FALSE, geneIdCol = "gene_id", abundanceCol = "FPKM",
                         countsCol = "expected_count", lengthCol = "effective_length", importer = function(x) read_tsv(x))
  rm(files)
  
  # Get read counts and FPKM
  counts = txi.rsem$abundance # read counts
  FPKM <- txi.rsem_2$abundance # FPKM
  
  ## Filter out genes with FPKM < 1 or < 5 reads across more than the smallest group's sample size
  # This makes the pipeline fast by reducing the size of the object in memory.
  # This also increases the power of multiple testing correciton.
  # rowSums(FPKM) >= x) >= y), where x is 1 FPKM and 5 counts and y is the smallest group's sample size. Similar to recommendation from edgeR and limma.
  
  # Also, filter by gene_type
  gene_info = read.csv(geneinfo_dir)
  
  # Strip gene versions
  rownames(FPKM) = stripGeneVersions(rownames(FPKM))
  rownames(counts) = stripGeneVersions(rownames(counts))
  
  # Filter
  keep <- rowSums(FPKM >= 1) >= NCOL(FPKM) &  rowSums(counts >= 5) >= NCOL(FPKM) & (!grepl("ENSG", rownames(FPKM)) | rownames(FPKM) %in% gene_info$ENSEMBL[gene_info$gene_type == "protein_coding"])
  counts = as.data.frame(counts[keep,])
  
  # Add filtered genes to list
  expressed_genes_list = append(expressed_genes_list, list(rownames(counts)))
  names(expressed_genes_list)[length(expressed_genes_list)] = sample
  
}

# Convert list consisting of vector of different lengths to a data frame
expressed_gene_df = t(plyr::ldply(expressed_genes_list, rbind))
colnames(expressed_gene_df) = expressed_gene_df[1,]
expressed_gene_df = expressed_gene_df[-1,]

# Export Expressed Genes data frame
write.csv(expressed_gene_df, file.path(out_dir, "expressed_genes.csv"), row.names=FALSE)

## print session info ##
print("Session Info below: ")
sessionInfo()
