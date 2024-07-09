rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("svglite")
#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("scales")

## Setup the environment

# Import packages
library(rstudioapi)
library(svglite)
library(ggplot2)
library(ggpubr)
library(scales)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Output directory
out_dir=file.path(work_dir, "Figures/Supplementary_Figure_9")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# MPC path
MPC_path = file.path(work_dir, "Data/MPC/MPC.csv")
if (!file.exists(MPC_path)) {stop(" MPC.csv file does not exist.\n\tRun MPC.R script to generate.", call. = FALSE)} # check file exists

# Path to LCL half-life data from PMID: 23422947 
df_hl_path = file.path(work_dir, "Data/LCL_mRNA_Half-Life/LCL_mRNA_half-lives.csv")

# DESeq2 differential expressed gene (DEG) files
deg_dir = file.path(work_dir, "Data/DESeq2")
if (!file.exists(list.files(deg_dir,full.names = TRUE)[1])) {stop(" No DESeq2 files exist.\n\tRun DESeq2.R script to generate.", call. = FALSE)} # check file exists

# Import differential expressed genes
deg_early = read.csv(file.path(deg_dir, "ERCCnorm_DE.dBGLF5_Early.vs.dBGLF5_Latent.csv"),check.names = TRUE)
deg_late = read.csv(file.path(deg_dir, "ERCCnorm_DE.dBGLF5_Late.vs.dBGLF5_Latent.csv"),check.names = TRUE)

# Import LCL half-life data from PMID: 23422947 
df_hl = read.csv(df_hl_path)

# Take average of samples with biological and tehcnical replicates (GM07029, GM10835, GM12813)
df_hl$HL = rowMeans(df_hl[,grepl("GM07029|GM10835|GM12813",colnames(df_hl))])

# Keep only human genes
deg_early = deg_early[grepl("ENSG",deg_early$ENSEMBL),]
deg_late = deg_late[grepl("ENSG",deg_late$ENSEMBL),]

# subset columns
deg_early = deg_early[,grepl("ENSEMBL|SYMBOL|log2FoldChange|padj|Latent_groupMean",colnames(deg_early))]
deg_late = deg_late[,grepl("ENSEMBL|SYMBOL|log2FoldChange|padj|Latent_groupMean",colnames(deg_late))]

# Merge deg files with half-life data
df_early = merge(x = deg_early, y = df_hl[,c("SYMBOL","HL")], by = "SYMBOL")
df_late = merge(x = deg_late, y = df_hl[,c("SYMBOL","HL")], by = "SYMBOL")

# Remove duplicates
df_early = df_early[!duplicated(df_early$SYMBOL),]
df_late = df_late[!duplicated(df_late$SYMBOL),]

# Plot Latent expression comparison
expression_plot_function = function(tmp_df) {
  tmp_df = tmp_df[tmp_df$log2FoldChange < 0 & tmp_df$padj <= 0.05,]
  colnames(tmp_df)[5] = "Latent_groupMean"
  tmp_df$Shape = ifelse(abs(tmp_df$log2FoldChange) > 8 | log2(tmp_df$Latent_groupMean) > 18, 24, 21)
  tmp_df$Size = ifelse(tmp_df$Shape == 21, 2, 1.5)
  ggplot(tmp_df, aes(x = log2(Latent_groupMean), y = log2FoldChange)) +
    geom_point(color = "black",fill="#1465AC", shape = tmp_df$Shape, size=tmp_df$Size) +
    theme_classic() +
    theme(text = element_text(size = 18)) +
    labs(x = "Log2 Latent Normalized Counts", y = "Log2 Fold Change (Relative Latent)") +
    scale_y_continuous(limits = c(-8,0), breaks = seq(-8,0,2), expand = c(0,0.15), oob = squish) +
    scale_x_continuous(limits = c(0,18), breaks = seq(0,18,3), expand = c(0,0.15), oob = squish)
}

# Plot half-life comparison function
hl_plot_function = function(tmp_df) {
  tmp_df = tmp_df[tmp_df$log2FoldChange < 0 & tmp_df$padj <= 0.05,]
  tmp_df$Shape = ifelse(abs(tmp_df$log2FoldChange) > 8 | log2(tmp_df$HL) > 5, 24, 21)
  tmp_df$Size = ifelse(tmp_df$Shape == 21, 2, 1.5)
  ggplot(tmp_df, aes(x = log2(HL), y = log2FoldChange)) +
    geom_point(color = "black",fill="#1465AC", shape = tmp_df$Shape, size=tmp_df$Size) +
    theme_classic() +
    theme(text = element_text(size = 18)) +
    labs(x = "Log2 mRNA Half-Life", y = "Log2 Fold Change (Relative Latent)") +
    scale_y_continuous(limits = c(-8,0), breaks = seq(-8,0,2), expand = c(0,0.15), oob = squish) +
    scale_x_continuous(limits = c(-1,5), breaks = seq(-1,5,2), expand = c(0,0.15), oob = squish)
}

expression_plot_early = expression_plot_function(df_early)
expression_plot_late = expression_plot_function(df_late)

hl_plot_early = hl_plot_function(df_early)
hl_plot_late = hl_plot_function(df_late)

# Get theme function
theme_white = theme(text = element_text(color = "white"), axis.text = element_text(color = "white"),
                    axis.line = element_line(color = "white"), axis.ticks = element_line(color = "white"))

# expression early
expression_plot_early + theme_white
ggsave(file.path(out_dir, "Supplementary_Figure_9A.Early.png"), dpi = 300, width = 5.21, height = 5.21, units = "in")

# expression late
expression_plot_late + theme_white
ggsave(file.path(out_dir, "Supplementary_Figure_9A.Late.png"), dpi = 300, width = 5.21, height = 5.21, units = "in")

# Half-life early
hl_plot_early + theme_white
ggsave(file.path(out_dir, "Supplementary_Figure_9B.Early.png"), dpi = 300, width = 5.21, height = 5.21, units = "in")

# Half-life late
hl_plot_late + theme_white
ggsave(file.path(out_dir, "Supplementary_Figure_9B.Late.png"), dpi = 300, width = 5.21, height = 5.21, units = "in")

# Get background
expression_bg_plot_early = expression_plot_early
expression_bg_plot_early$layers <- NULL
expression_bg_plot_early + geom_smooth(method = "lm", se = FALSE, color = "black") + stat_cor(method = "pearson", aes(label = after_stat(r.label)), label.y = 7, label.x = -1)
ggsave(file.path(out_dir, "Supplementary_Figure_9A.Early.svg"), dpi = 300, width = 5.21, height = 5.21, units = "in")

expression_bg_plot_late = expression_plot_late
expression_bg_plot_late$layers <- NULL
expression_bg_plot_late + geom_smooth(method = "lm", se = FALSE, color = "black") + stat_cor(method = "pearson", aes(label = after_stat(r.label)), label.y = 7, label.x = -1)
ggsave(file.path(out_dir, "Supplementary_Figure_9A.Late.svg"), dpi = 300, width = 5.21, height = 5.21, units = "in")

hl_bg_plot_early = hl_plot_early
hl_bg_plot_early$layers <- NULL
hl_bg_plot_early + geom_smooth(method = "lm", se = FALSE, color = "black") + stat_cor(method = "pearson", aes(label = after_stat(r.label)), label.y = 7, label.x = -1)
ggsave(file.path(out_dir, "Supplementary_Figure_9B.Early.svg"), dpi = 300, width = 5.21, height = 5.21, units = "in")

hl_bg_plot_late = hl_plot_late
hl_bg_plot_late$layers <- NULL
hl_bg_plot_late + geom_smooth(method = "lm", se = FALSE, color = "black") + stat_cor(method = "pearson", aes(label = after_stat(r.label)), label.y = 7, label.x = -1)
ggsave(file.path(out_dir, "Supplementary_Figure_9B.Late.svg"), dpi = 300, width = 5.21, height = 5.21, units = "in")

## print session info ##
print("Session Info below: ")
sessionInfo()
