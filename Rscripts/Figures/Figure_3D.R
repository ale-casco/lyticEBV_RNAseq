rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("ggplot2")
#install.packages("readr")
#install.packages("svglite")

## Setup the environment

# Import packages
library(rstudioapi)
library(ggplot2)
library(readr)
library(svglite)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Output directory
out_dir=file.path(work_dir, "Figures/Figure_3")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# DAVID results
david_dir=file.path(work_dir, "Data/DAVID")

# DESeq2 differential expressed gene (DEG) files
deg_dir=file.path(work_dir, "Data/DESeq2")
file = file.path(deg_dir,"ERCCnorm_DE.WT_Early.vs.WT_Latent.csv")
if (!file.exists(file)) {stop(" No DESeq2 file exist.\n\tRun DESeq2.R script to generate.", call. = FALSE)} # check file exists

# Import DESeq2 DEG file into data frame
df = read.csv(file)

# keep only human genes that are protein coding
df = df[grepl("ENSG",df$ENSEMBL) & df$gene_type == "protein_coding",]

# Get genes that are not significantly downregulated
df$padj = ifelse(is.na(df$padj), 1, df$padj)
df$Sig = ifelse(df$log2FoldChange < 0 & df$padj <= 0.05, "Downregulated", "Non-downregulated")
#df = df[df$Sig == "Non-downregulated",]

## Run DAVID online (https://david.ncifcrf.gov/summary.jsp)
# 1. Copy ENSEMBL genes into Step 1: Gene List
writeClipboard(df$ENSEMBL[df$Sig == "Non-downregulated"])
writeClipboard(df$ENSEMBL)
# 2. For Step 2: Select Identifier, selected ENSEMBL_GENE_ID
# 3. For Step 3: List Type, selected Gene List
# 4. Select only GOTERM_BP_DIRECT and KEGG_PATHWAY
# 5. Perform Functional Annotation Clustering  & Download file to directory specified by david_dir

# Import DAVID results
david_results_path = file.path(work_dir, "Data/DAVID/Early_clusters.txt")
if (!file.exists(david_results_path)) {stop(" No DAVID_results file exists.\n\tRun DAVID online according to instructions above.", call. = FALSE)} # check file exists
df_cluster = data.frame(read_delim(david_results_path, 
                        delim = "\t", escape_double = FALSE, 
                        col_names = FALSE, trim_ws = TRUE))

# Format DAVID results
cluster_rows = which(grepl("Annotation Cluster", df_cluster[,1]))[-1]
for (i in 1:length(cluster_rows)) {
  if (i == 1) {
    df_tmp = df_cluster[2:(cluster_rows[i]-1),]
    df_tmp = data.frame(t(data.frame(strsplit(df_tmp[,2], '\t'))))
    colnames(df_tmp) = df_tmp[1,]
    df_tmp = df_tmp[-1,]
    rownames(df_tmp) = 1:NROW(df_tmp)
    df_tmp$Cluster = i
    df = df_tmp
  } else {
    df_tmp = df_cluster[(cluster_rows[i-1]+1):(cluster_rows[i]-1),]
    df_tmp = data.frame(t(data.frame(strsplit(df_tmp[,2], '\t'))))
    colnames(df_tmp) = df_tmp[1,]
    df_tmp = df_tmp[-1,]
    rownames(df_tmp) = 1:NROW(df_tmp)
    df_tmp$Cluster = i
    df = rbind(df, df_tmp)
  }
}

# Calculate GeneRatio
colnames(df) = sub("%","GeneRatio",colnames(df))
df$GeneRatio = as.numeric(df$GeneRatio)/100

# Make GeneRatio and FDR columns numeric
df$GeneRatio = as.numeric(df$GeneRatio)
df$FDR = signif(as.numeric(df$FDR),digits = 1)

# Keep only significant pathways
df = df[df$FDR <= 0.05,]
df$FDR_transformed = -log10(df$FDR)

# Determine clusters based on Functional ANnotation Clustering file
df$Term_fixed = sub(".*\\:","",sub(".*\\~","",df$Term))

reorder(df$Term_fixed, df$FDR_transformed)

# Plot
ggplot(df, aes(x = GeneRatio, y = reorder(Term_fixed, FDR_transformed), fill = FDR_transformed)) +
  geom_bar(stat="identity",alpha = 0.8) +
  #coord_flip() +
  theme_classic() +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_size(range = c(2,8)) +
  theme(text = element_text(size = 15), axis.title = element_blank(), axis.text = element_text(color = "black")) +
  facet_grid(Cluster~., scales = "free", space='free_y')

# Plot as SVG file
ggsave(file.path(out_dir,"Figure_3D.svg"), dpi = 300, width = 12, height = 6, units = "in")

## print session info ##
print("Session Info below: ")
sessionInfo()
