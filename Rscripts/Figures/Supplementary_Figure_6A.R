rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("svglite")
#install.packages("ggplot2")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexUpset")

## Setup the environment

# Import packages
library(rstudioapi)
library(svglite)
library(ggplot2)
library(ComplexUpset)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Output directory
out_dir=file.path(work_dir, "Figures/Supplementary_Figure_6")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Import files
LCL_early_path = file.path(work_dir, "Data/DESeq2/ERCCnorm_DE.WT_Early.vs.WT_Latent.csv")
if (!file.exists(LCL_early_path)) {stop(" No LCL DESeq2 file exists.\n\tRun DESeq2.R script to generate.", call. = FALSE)} # check file exists
LCL_early = read.csv(LCL_early_path,check.names = TRUE)
Raji_early_path = file.path(work_dir, "Data_Raji/DESeq2/ERCCnorm_DE.Early.vs.Latent.csv")
if (!file.exists(Raji_early_path)) {stop(" No Raji DESeq2 file exists.\n\tRun Raji_DESeq2.R script to generate.", call. = FALSE)} # check file exists
Raji_early = read.csv(Raji_early_path,check.names = TRUE)

# Get list
early_list = list(LCL_early,Raji_early)
names(early_list) = c("LCL","Raji")

# Keep only human genes
for (i in 1:length(early_list)) {  early_list[[i]] = early_list[[i]][grepl("ENSG",early_list[[i]]$ENSEMBL),] }

# Filter for similar genes
filter="yes" # <no|yes>
if (filter == "yes") {
  early_list[[1]] = early_list[[1]][early_list[[1]]$ENSEMBL %in% early_list[[2]]$ENSEMBL,]
  early_list[[2]] = early_list[[2]][early_list[[2]]$ENSEMBL %in% early_list[[1]]$ENSEMBL,]
}

# Get sig
early_list[[1]]$Sig = ifelse(early_list[[1]]$log2FoldChange < 0 & early_list[[1]]$padj <= 0.05, "Down", "Escapee")
early_list[[2]]$Sig = ifelse(early_list[[2]]$log2FoldChange < 0 & early_list[[2]]$padj <= 0.05, "Down", "Escapee")

# Get overlapping genes list
degs = list(
  "LCL down" = early_list[[1]]$ENSEMBL[early_list[[1]]$Sig == "Down"],
  "LCL Escapee" = early_list[[1]]$ENSEMBL[early_list[[1]]$Sig == "Escapee"],
  "Raji down" = early_list[[2]]$ENSEMBL[early_list[[2]]$Sig == "Down"],
  "Raji Escapee" = early_list[[2]]$ENSEMBL[early_list[[2]]$Sig == "Escapee"]
)

# Get intersections
fromList_function = function (input) {
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  return(data)
}

# upset function
upset_function = function(input_list,ylabel="% Gene Intersections") {
  upset(
    fromList_function(input_list), rev(names(input_list)), name='Cell Type', width_ratio=0.1, min_size=10,wrap=TRUE,
    base_annotations=list(
      'Intersection size'=intersection_size(
        text_mapping=aes(label=paste0(
          paste0(round(!!get_size_mode("distinct")/nrow(fromList_function(input_list)) * 100,2), "%"),
          '\n',
          '(',
          !!get_size_mode('exclusive_intersection'),
          ')'
        )),
        bar_number_threshold=1,
        text=list(vjust=0.55)
      )
    )
  ) +
    ggtitle(ylabel)
}

# Early
out_plot = upset_function(degs, "Overlap between host shutoff genes across cell types during early lytic phase")

# Export as SVG
ggsave(file.path(out_dir, "Supplementary_Figure_6A.svg"), plot = out_plot, bg="transparent")

## print session info ##
print("Session Info below: ")
sessionInfo()
