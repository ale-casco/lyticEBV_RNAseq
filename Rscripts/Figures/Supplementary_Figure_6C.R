rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("svglite")
#install.packages("reshape2")
#install.packages("ggplot2")
#install.packages("ggstatsplot")

## Setup the environment

# Import packages
library(rstudioapi)
library(svglite)
library(reshape2)
library(ggplot2)
library(ggstatsplot)

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

# Keep only human genes, subset LFC column, and rename column
LCL_early = LCL_early[grepl("ENSG",LCL_early$ENSEMBL),c("ENSEMBL","log2FoldChange")]
Raji_early = Raji_early[grepl("ENSG",Raji_early$ENSEMBL),c("ENSEMBL","log2FoldChange")]
colnames(LCL_early)[2] = c("LCL")
colnames(Raji_early)[2] = c("Raji")

# Get list
early_list = list(LCL_early,Raji_early)
names(early_list) = c("LCL","Raji")

# Merge list into data frame, keeping all genes
df <- Reduce(function(...) merge(..., by="ENSEMBL", all = TRUE), early_list)

# Melt data
df = melt(df)
df = df[!is.na(df$value),]

# Plot
plot_out = ggbetweenstats(data = df, x = variable, y = value, type = "p", var.equal = FALSE, p.adjust.method = "fdr",point.args = list(alpha = 0)) +
  coord_cartesian(ylim = c(-6,3)) +
  geom_hline(yintercept=1,linetype="dashed", color = "green", linewidth = 1) +
  geom_hline(yintercept=-1,linetype="dashed", color = "green", linewidth = 1) +
  theme(text = element_text(size = 12))

# Welch’s t-test
t_test.result = t.test(value ~ variable, data = df, var.equal = FALSE)
t_test.result

# Export
plot_out$layers[[1]] <- NULL
plot_out
ggsave(file.path(out_dir, "Supplementary_Figure_6C.svg"),
       dpi = 10, width = 7.64, height = 4.91, units = "in")

## print session info ##
print("Session Info below: ")
sessionInfo()
