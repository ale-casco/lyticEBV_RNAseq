rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("stringr")
#install.packages("dplyr")
#install.packages("reshape2")
#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("svglite")

## Setup the environment

# Import packages
library(rstudioapi)
library(stringr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(svglite)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Output directory
out_dir=file.path(work_dir, "Figures/Figure_3")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Load molecules per cell (MPC) data
MPC_path = file.path(work_dir, "Data/MPC/MPC.csv")
if (!file.exists(MPC_path)) {stop(" MPC.csv file does not exist.\n\tRun MPC.R script to generate.", call. = FALSE)} # check file exists
MPC=read.csv(MPC_path,check.names=FALSE)

# keep only human protein coding genes
MPC = MPC[grepl("ENSG",MPC$ENSEMBL) & MPC$gene_type == "protein_coding" | is.na(MPC$ENSEMBL),]

# Convert to matrix style
MPC$ENSEMBL = ifelse(is.na(MPC$ENSEMBL), MPC$SYMBOL, MPC$ENSEMBL) # if ENSEMBL ID does not exist, use SYMBOL
rownames(MPC) = MPC$ENSEMBL
MPC = MPC[,-1:-4]

# Subset host and ebv genes
MPC_host = MPC[grepl("ENSG",rownames(MPC)),]
MPC_ebv = MPC[!grepl("ENSG|ERCC",rownames(MPC)),]

# Create data frame of summed valeus
MPC_list = list(MPC_host, MPC_ebv)
organism = c("Host","EBV")
for (i in 1:length(MPC_list)) {
  MPC_list[[i]]=colSums(MPC_list[[i]]) # sum the MPC per condition
  MPC_list[[i]]=melt(MPC_list[[i]]) # melt the data
  MPC_list[[i]]$Condition=sub(".*_","",gsub('.{5}$', '', rownames(MPC_list[[i]]))) # Condition column
  MPC_list[[i]]$Mutant=sub("_.*","",gsub('.{5}$', '', rownames(MPC_list[[i]]))) # Mutant (WT or dBGLF5) column
  MPC_list[[i]]$Organism=organism[i] # Organism (host or EBV) column
  rownames(MPC_list[[i]])=NULL
}
MPC=bind_rows(MPC_list) # combine host and ebv data frames into one
rm(list=setdiff(ls(), c("MPC","out_dir")))

# Factor data
MPC$Condition=factor(MPC$Condition, levels=c("Bulk","Latent","Refractory","Early","Late"))
MPC$Mutant=factor(MPC$Mutant, levels=c("WT","dBGLF5"))

################ FIGURE ################

# Host + EBV WT (350 x 350)
MPC = MPC[!grepl("Bulk",MPC$Condition),] # Remove Bulk condition
ggbarplot(position=position_stack(),label = TRUE, label.pos = "out",
          MPC[grepl("WT",MPC$Mutant),], x = "Condition", y = "value", add = "mean_se",
          color = "Organism", fill = "Organism", palette = c("tomato","skyblue"),
          ylab = "Mean RNA Molecules Per Cell", xlab = FALSE, add.params = list(color = "black"))

# Export plot
ggsave(file.path(out_dir,"Figure_3A.svg"), dpi = 300, width = 6.2, height = 6.2, units = "in")

## print session info ##
print("Session Info below: ")
sessionInfo()
