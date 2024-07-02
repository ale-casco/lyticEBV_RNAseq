rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("svglite")
#install.packages("dplyr")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")

## Setup the environment

# Import packages
library(rstudioapi)
library(svglite)
library(dplyr)
library(ComplexHeatmap)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))


# Output directory
out_dir=file.path(work_dir, "Figures/Supplementary_Figure_2")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Import DGE file
deg_dir = file.path(work_dir, "Data/DESeq2")
if (!file.exists(list.files(deg_dir,full.names = TRUE)[1])) {stop(" No DESeq2 files exist.\n\tRun DESeq2.R script to generate.", call. = FALSE)} # check file exists
files = list.files(deg_dir, full.names=TRUE)

# Comparisons to keep
files = files[grepl("WT_Bulk.vs.WT_Latent|WT_Early.vs.WT_Latent|WT_Late.vs.WT_Latent",files)]

# Import DESeq2 DEG files into list
df_list = list()
for (i in 1:length(files)) {
  df_tmp = read.csv(files[i]) # read deg file
  df_tmp = df_tmp[!grepl("ENSG", df_tmp$ENSEMBL) & !grepl("ERCC-", df_tmp$SYMBOL),] # keep only ebv genes
  df_tmp = df_tmp[,!grepl("ENSEMBL|gene_type|strand|chrom|log2FoldChange|pvalue|padj|Mean",colnames(df_tmp))] # keep SYMBOL and read counts
  df_list = append(df_list, list(df_tmp)) # append to list
}

# Merge list into data frame
df <- Reduce(function(...) merge(..., by="SYMBOL", all = TRUE, no.dups = FALSE), df_list)

# Remove duplicated columns
df_dup = df[,grepl("\\.x$",colnames(df))]
colnames(df_dup) = sub("\\.x","",colnames(df_dup))
df = cbind(df[,!grepl(paste(unique(colnames(df_dup)),collapse="|"),colnames(df))], df_dup)

# Replace NA values with 0
df[is.na(df)] <- 0

# Convert to matrix style
rownames(df) = df[,1]
df = df[,-1]

# Rename some EBV genes
rownames(df) = sub("BILF2-tPT2A-mScarlet-I", "BILF2-mSI", rownames(df))
rownames(df) = sub("mGreenLantern", "mGL", rownames(df))
rownames(df) = sub("BSLF2/BMLF1", "SM", rownames(df))
rownames(df) = sub("BBLF2/BBLF3", "BBLF2/3", rownames(df))
rownames(df) = sub("BGRF1/BDRF1", "BG/DRF1", rownames(df))

# Remove genes
df = df[!grepl("EBER1|EBER2|BWRF1", rownames(df)),]

##### get lists of EBV genes by their expression kinetics
Early = c("BZLF1","BRLF1","mGL", "BHLF1","BFLF2","BORF2","BaRF1","BMRF1","BMRF2","SM","BSLF1","BSRF1","BLLF3","BLLF2","BRRF1","BKRF3",
          "BBLF4","BBLF2/3","BGLF5","BGLF4","BGLF3.5","BGLF3","BDLF3.5","BVRF1","BVLF1","BALF2","BALF1","BARF1",
          "BCRF1","BFRF2","BG/DRF1","BTRF1","LF3","LF2","LF1","BILF1","BALF5","BNLF2a/b")
Leaky = c("BFLF1","BORF1","BLRF1","BLRF2","BRRF2","BKRF2","BBRF1","BBRF2","BBRF3","BBLF1","BDLF3","BcRF1","BXLF2","BdRF1",
          "BALF4","BALF3","BFRF0.5","BKRF4")
Late = c("BNRF1","BFRF1","BFRF3","BPLF1","BOLF1","BLLF1","BZLF2","BGLF2","BGLF1","BDLF2","BDLF1","BcLF1","BXRF1",
         "BVRF2","BILF2-mSI","BDLF4","LMP1_lytic")
Latent = c("LMP2", "LMP2A", "LMP2B","EBER1","EBER2","EBNALP","EBNA2","EBNA3A","EBNA3B","EBNA3C","EBNA1","BART","LMP1","BHRF1")

# get vector of EBV genes by their expression kinetics
gene = ifelse(rownames(df) %in% Early, "Early",
                     ifelse(rownames(df) %in% Leaky, "Leaky",
                            ifelse(rownames(df) %in% Late, "Late",
                                   ifelse(rownames(df) %in% Latent, "Latent", "Unclassified"))))

rm(list=setdiff(ls(), c("df","gene")))

# Rename columns
colnames(df) = sub("_"," ", sub("WT_","",colnames(df)))

# Get df of conditions
conditions = data.frame(gsub('.{5}$', '', colnames(df)))
colnames(conditions) = "condition"

# Color conditions
conditions$color = ifelse(conditions$condition == "Bulk", "yellow",
                                     ifelse(conditions$condition == "Latent", "gray",
                                            ifelse(conditions$condition == "Early", "green",
                                                   ifelse(conditions$condition == "Late","red",NA))))

# Factor condition
conditions$condition = factor(conditions$condition, levels = c("Bulk","Latent","Early","Late"))

# Z transform by gene
df = scale(t(df), scale=TRUE, center=TRUE)

# Get annotation colors
row_ha <- HeatmapAnnotation(conditions = conditions %>% pull(condition), 
                            col = list(conditions = c("Bulk"="yellow","Latent"="gray","Early"="green","Late"="red")),
                            which = "row")

heatmap = Heatmap(df, name = "z-score",
                  left_annotation = row_ha,
                  column_split = factor(c(gene), levels = rev(c("Late", "Leaky", "Early","Latent", "Unclassified"))),
                  #row_split = factor(colsplit,levels = c("Bulk","Latent","Early","Late")),
                  cluster_row_slices = FALSE,
                  cluster_column_slices = FALSE,
                  show_row_dend = TRUE,
                  show_column_dend = TRUE,
                  show_row_names = TRUE,
                  cluster_columns = TRUE,
                  column_names_side = "bottom", 
                  column_dend_side = "top",
                  row_gap = unit(1, "mm"), column_gap = unit(1, "mm"), border = TRUE,
                  show_parent_dend_line = TRUE,
                  column_names_rot = 90,
                  row_title_rot = 0)
heatmap

# Save as SVG

## print session info ##
print("Session Info below: ")
sessionInfo()

