rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("svglite")
#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("reshape2")
#install.packages("eulerr")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("AcidGenomes")

## Setup the environment

# Import packages
library(rstudioapi)
library(svglite)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(eulerr)
library(AcidGenomes)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Output directory
out_dir=file.path(work_dir, "Figures/Supplementary_Figure_8")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Directory containing DoG BED files
dogs_dir = file.path(work_dir, "Data/ARTDeco/dogs")

# Sampleinfo file path
sampleinfo_path = file.path(work_dir, "Data/refs/sampleinfo.txt")

# Expressed genes file path
expressed_genes_path = file.path(work_dir, "Data/Expressed_genes/expressed_genes.csv")
if (!file.exists(expressed_genes_path)) {stop(" expressed_genes.csv file does not exist.\n\tRun Expressed_genes.R script to generate.", call. = FALSE)} # check file exists

# Import DoG files
files = list.files(dogs_dir, pattern = "*out.dogs.bed", full.names=TRUE)

# Import sampleinfo file
sampleinfo = read.delim(sampleinfo_path)

# Reorder DoG list according to sampleinfo
names(files) = sub("\\.Aligned.*","",basename(files))
sampleinfo = sampleinfo[sampleinfo[,1] %in% names(files),]
temp_list=as.list(sampleinfo[[1]])
temp_files=c()
for (i in as.character(temp_list)) {
  temp_files=append(temp_files, c(print(files[grep(i, files)])))
}
files=temp_files
rm(temp_list, temp_files)

# Name Samples
names(files)=paste0(sampleinfo[[2]], "_", sampleinfo[[3]])

# Import expressed genes
expressed_genes = read.csv(expressed_genes_path)

# Get expressed genes
latent_expressed = na.omit(unique(expressed_genes$dBGLF5_Latent))
early_expressed = na.omit(unique(expressed_genes$dBGLF5_Early))
late_expressed = na.omit(unique(expressed_genes$dBGLF5_Late))

# Import files into list, filtering for expressed genes
latent_list = list()
for (i in which(grepl("dBGLF5_Latent_c", names(files)))) {
  tmp = stripGeneVersions(read.delim(files[i])[[4]])
  tmp = tmp[tmp %in% latent_expressed]
  latent_list = append(latent_list, list(tmp))
  names(latent_list)[length(latent_list)] = names(files)[i]
}
early_list = list()
for (i in which(grepl("dBGLF5_Early_c", names(files)))) {
  tmp = stripGeneVersions(read.delim(files[i])[[4]])
  tmp = tmp[tmp %in% early_expressed]
  early_list = append(early_list, list(tmp))
  names(early_list)[length(early_list)] = names(files)[i]
}
late_list = list()
for (i in which(grepl("dBGLF5_Late_c", names(files)))) {
  tmp = stripGeneVersions(read.delim(files[i])[[4]])
  tmp = tmp[tmp %in% late_expressed]
  late_list = append(late_list, list(tmp))
  names(late_list)[length(late_list)] = names(files)[i]
}

# Create list
latent_list = list(
  "c2-1" = latent_list[[1]],
  "c2-2" = latent_list[[2]],
  "c2-3" = latent_list[[3]]
)

early_list = list(
  "c2-1" = early_list[[1]],
  "c2-2" = early_list[[2]],
  "c2-3" = early_list[[3]]
)

late_list = list(
  "c2-1" = late_list[[1]],
  "c2-2" = late_list[[2]],
  "c2-3" = late_list[[3]]
)

# Convert lists to data frames
latent_df = t(plyr::ldply(latent_list, rbind))
colnames(latent_df) = latent_df[1,]
latent_df = latent_df[-1,]
for (i in 2:NCOL(latent_df)) {
  latent_df[,i] = latent_df[match(latent_df[,1], latent_df[,i]),i]
}

early_df = t(plyr::ldply(early_list, rbind))
colnames(early_df) = early_df[1,]
early_df = early_df[-1,]
for (i in 2:NCOL(early_df)) {
  early_df[,i] = early_df[match(early_df[,1], early_df[,i]),i]
}

late_df = t(plyr::ldply(late_list, rbind))
colnames(late_df) = late_df[1,]
late_df = late_df[-1,]
for (i in 2:NCOL(late_df)) {
  late_df[,i] = late_df[match(late_df[,1], late_df[,i]),i]
}

# Get data frame containing total and DoGs across sample replicates
df = as.data.frame(matrix(ncol = 4, nrow = 2))
colnames(df) = c("Gene_Summary", "Latent", "Early", "Late")
df$Gene_Summary = c("Variable", "Common")
df[,2] = c(length(unique(unlist(latent_list)))-NROW(as.data.frame(na.omit(latent_df))), NROW(as.data.frame(na.omit(latent_df))))
df[,3] = c(length(unique(unlist(early_list)))-NROW(as.data.frame(na.omit(early_df))), NROW(as.data.frame(na.omit(early_df))))
df[,4] = c(length(unique(unlist(late_list)))-NROW(as.data.frame(na.omit(late_df))), NROW(as.data.frame(na.omit(late_df))))

# Melt
df = melt(df)

# Factor
df$Gene_Summary = factor(df$Gene_Summary, levels = c("Variable", "Common"))
df$variable = factor(df$variable, levels = c("Latent", "Early", "Late"))

##### Plot Figure #####

# Plot
ggbarplot(position=position_stack(),label = TRUE, label.pos = "out",
          df, x = "variable", y = "value",
          color = "Gene_Summary", fill = "Gene_Summary", palette = c("gray","blue"),
          ylab = "DoGs", xlab = FALSE, add.params = list(color = "black"))

# Export as SVG
ggsave(file.path(out_dir, "Supplementary_Figure_8F.svg"), width = 5, height = 5)

## print session info ##
print("Session Info below: ")
sessionInfo()
