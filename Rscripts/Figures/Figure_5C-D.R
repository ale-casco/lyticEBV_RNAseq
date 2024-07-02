rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("reshape2")
#install.packages("eulerr")
#install.packages("svglite")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("AcidGenomes")

## Setup the environment

# Import packages
library(rstudioapi)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(eulerr)
library(AcidGenomes)
library(svglite)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Output directory
out_dir=file.path(work_dir, "Figures/Figure_5")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Directory containing DoG BED files
dogs_dir=file.path(work_dir, "Data/ARTDeco/dogs")

# Sampleinfo file path
sampleinfo_path=file.path(work_dir, "Data/refs/sampleinfo.txt")

# Expressed genes file path
expressed_genes_path=file.path(work_dir, "Data/Expressed_genes/expressed_genes.csv")
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
latent_expressed = na.omit(unique(expressed_genes$WT_Latent))
early_expressed = na.omit(unique(expressed_genes$WT_Early))
late_expressed = na.omit(unique(expressed_genes$WT_Late))

# Import files into list, filtering for expressed genes
latent_list = list()
for (i in which(grepl("WT_Latent_c", names(files)))) {
  tmp = stripGeneVersions(read.delim(files[i])[[4]])
  tmp = tmp[tmp %in% latent_expressed]
  latent_list = append(latent_list, list(tmp))
  names(latent_list)[length(latent_list)] = names(files)[i]
}
early_list = list()
for (i in which(grepl("WT_Early_c", names(files)))) {
  tmp = stripGeneVersions(read.delim(files[i])[[4]])
  tmp = tmp[tmp %in% early_expressed]
  early_list = append(early_list, list(tmp))
  names(early_list)[length(early_list)] = names(files)[i]
}
late_list = list()
for (i in which(grepl("WT_Late_c", names(files)))) {
  tmp = stripGeneVersions(read.delim(files[i])[[4]])
  tmp = tmp[tmp %in% late_expressed]
  late_list = append(late_list, list(tmp))
  names(late_list)[length(late_list)] = names(files)[i]
}

# Create list
latent_list = list(
  "c3-1" = latent_list[[1]],
  "c3-2" = latent_list[[2]],
  "c3-3" = latent_list[[3]],
  "c4-1" = latent_list[[4]],
  "c4-2" = latent_list[[5]]
)

early_list = list(
  "c3-1" = early_list[[1]],
  "c3-2" = early_list[[2]],
  "c3-3" = early_list[[3]],
  "c4-1" = early_list[[4]],
  "c4-2" = early_list[[5]]
)

late_list = list(
  "c3-1" = late_list[[1]],
  "c3-2" = late_list[[2]],
  "c3-3" = late_list[[3]],
  "c4-1" = late_list[[4]],
  "c4-2" = late_list[[5]]
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

##### Barplot Figure 5C #####

# Plot
ggbarplot(position=position_stack(),label = TRUE, label.pos = "out",
          df, x = "variable", y = "value",
          color = "Gene_Summary", fill = "Gene_Summary", palette = c("gray","blue"),
          ylab = "DoGs", xlab = FALSE, add.params = list(color = "black"))

# Export as SVG
ggsave(file.path(out_dir, "Figure_5C.svg"), width = 5, height = 5)

##### Venn Diagram Figure 5D #####

vennfun <- function(x) { 
  x$id <- seq(1, nrow(x))  #add a column of numbers (required for melt)
  xm <- melt(x, id.vars="id", na.rm=TRUE)  #melt table into two columns (value & variable)
  xc <- dcast(xm, value~variable, fun.aggregate=length)  #remove NA's, list presence/absence of each value for each variable (1 or 0)
  rownames(xc) <- xc$value  #value column = rownames (required for Venneuler)
  xc$value <- NULL  #remove redundent value column
  xc  #output the new dataframe
}
venn_dat = function(x) {
  x = plyr::ldply(x, rbind)
  x = t(x)
  colnames(x) = x[1,]
  x = as.data.frame(x[-1,])
  x = vennfun(x)
  x = euler(x)
}

# venn diagram list
vd_list = list(
  "Latent" = unique(as.vector(na.omit(latent_df))),
  "Early" = unique(as.vector(na.omit(early_df))),
  "Late" = unique(as.vector(na.omit(late_df)))
)

# Plot venn diagram
venn = venn_dat(vd_list)
vd_plot = plot(venn, quantities = TRUE,
     fills = c("gray", "green", "red"))

# Export as SVG
ggsave(file.path(out_dir, "Figure_5D.svg"), plot = vd_plot, width = 5, height = 5)

## print session info ##
print("Session Info below: ")
sessionInfo()
