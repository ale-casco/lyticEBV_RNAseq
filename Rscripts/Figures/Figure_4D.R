rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("pheatmap")
#install.packages("devEMF")
#install.packages("DT")
#install.packages("svglite")
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("SpliceWiz")
#BiocManager::install("fgsea")

## Setup the environment

# Import packages
library(rstudioapi)
library(SpliceWiz)
library(pheatmap)
library(devEMF)
library(DT)
library(svglite)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Output directory
out_dir=file.path(work_dir, "Figures/Figure_4")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Sample information directory
sampleinfo_dir=file.path(work_dir,"Data/refs/sampleinfo.txt")

# WT SpliceWiz experiment directory
WT_nxtse_dir=file.path(work_dir,"Data/SpliceWiz/NxtSE/WT_NxtSE")
if (!file.exists(list.files(WT_nxtse_dir, full.names=TRUE)[1])) {stop(" NxtSE objects do not exist.\n\tRun SpliceWiz.R script to generate.", call. = FALSE)} # check file exists

# DESeq2 differential gene expression directory
deg_dir=file.path(work_dir,"Data/DESeq2/ERCCnorm_DE.WT_Early.vs.WT_Latent.csv")
if (!file.exists(deg_dir)) {stop(" No DESeq2 file exists.\n\tRun DESeq2.R script to generate.", call. = FALSE)} # check file exists

# DESeq2 diffential alternative splicing events directory
de_ase_dir=file.path(work_dir,"Data/SpliceWiz/DESeq2/DE.WT_Early.vs.WT_Latent.csv")
if (!file.exists(de_ase_dir)) {stop(" No SpliceWiz DESeq2 file exists.\n\tRun SpliceWiz.R script to generate.", call. = FALSE)} # check file exists

# Import the experiments
WT_se = makeSE(WT_nxtse_dir)

# Load sampleinfo file
sampleinfo <- read.delim(file.path(sampleinfo_dir))

# Assign annotations to samples
colData(WT_se)$condition = sampleinfo[which(sampleinfo$UNIQUE_ID %in% colnames(WT_se)),2]

# Filter for high-confidence events and remove WT_Bulk and WT_Late condition
WT_se.filtered = WT_se[applyFilters(WT_se),which(colnames(WT_se) %in% colnames(WT_se)[WT_se$condition != "WT_Bulk" & WT_se$condition != "WT_Late"])]

# Import differential expression/splicing data
deg = read.csv(deg_dir)
de_ase = read.csv(de_ase_dir)

# Filter genes based on differential expression data
de_ase = de_ase[de_ase$Exc_gene_id %in% deg$ENSEMBL,]

# Get only skipped exon events (lfc <= -1 & padj <= 0.05 & abs_deltaPSI >= 0.05)
se = de_ase[!is.na(de_ase$abs_deltaPSI),]
se = se[se$log2FoldChange <= -1 & se$padj <= 0.05 & se$abs_deltaPSI >= 0.05 & se$EventType == "SE",]
se = se[order(se$deltaPSI, decreasing = TRUE),]

# Perform over-representation analysis of exon skipping events using background genes from analyzed ASEs using nominal p-value
colnames(WT_se.filtered@metadata$ref$ontology)[4] = "gene_id"
go_byASE = goASE(
  enrichedEventNames = se$EventName,
  universeEventNames = de_ase$EventName,
  se = WT_se.filtered,
  pAdjustMethod = "none"
)

# Filter GO terms by nominal pvalue <= 0.05 and >= 10 genes
# GO terms with fewer than 10 genes may reflect mere chance rather than evidence for enrichment
go_byASE = go_byASE[go_byASE$pval <= 0.05 & go_byASE$size >= 10,]

# Table of GO terms
go_byASE$pval = formatC(go_byASE$pval, format = "e", digits = 2)
datatable(
  data = go_byASE[,c("go_term","pval","overlap","size","foldEnrichment")]
)

# Obtain matrix of PSI values from the top differentiall expressed events
go_genes = unlist(go_byASE$overlapGenes[go_byASE$go_term == "defense response to virus"])
se_go_genes = se[se$Exc_gene_id %in% go_genes,]

# Order by PSI
se_go_genes = se_go_genes[order(se_go_genes$deltaPSI, decreasing = FALSE),]

mat = makeMatrix(
  WT_se.filtered,
  event_list = se_go_genes$EventName,
  method = "PSI"
)

# Plot this matrix of values in a heatmap

anno_col_df <- as.data.frame(colData(WT_se.filtered))
anno_col_df <- anno_col_df[, 1, drop=FALSE]
ann_colors = list(condition = c(WT_Latent="gray", WT_Early="green"))

heatmap = pheatmap(mat[1:10,],
         annotation_col = anno_col_df,
         annotation_colors = ann_colors)


emf(file=file.path(out_dir,"Figure_4D.emf"), width = 7, height = 3.5)
heatmap
dev.off()
svg(file=file.path(out_dir,"Figure_4D.svg"), width = 8, height = 3.5)
heatmap
dev.off()
png(file=file.path(out_dir,"Figure_4D.png"), width = 7, height = 3.5, units = "in", res = 1200)
heatmap
dev.off()

## print session info ##
print("Session Info below: ")
sessionInfo()

