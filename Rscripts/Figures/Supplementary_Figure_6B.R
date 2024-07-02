rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("svglite")
#install.packages("ggpubr")
#install.packages("scales")

## Setup the environment

# Import packages
library(rstudioapi)
library(svglite)
library(ggpubr)
library(scales)

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
colnames(LCL_early)[2] = c("LCL_LFC")
colnames(Raji_early)[2] = c("Raji_LFC")

# Get list
early_list = list(LCL_early,Raji_early)
names(early_list) = c("LCL","Raji")

# Merge list into data frame, keeping only genes present in all comparisons
df <- Reduce(function(...) merge(..., by="ENSEMBL"), early_list)

df$Shape = ifelse(abs(df$LCL_LFC) > 6 | abs(df$Raji_LFC) > 6, 17, 19)
int = lm(LCL_LFC ~ Raji_LFC, df)[[1]][[1]]

plot_out = ggplot(df, aes(x = Raji_LFC, y = LCL_LFC)) +
  #scale_color_manual(values = c("blue","green4","orange","gray")) +
  geom_point(shape = df$Shape) +
  theme_classic() +
  theme(legend.position = "none", text = element_text(size = 12)) +
  scale_y_continuous(limits = c(-6,6), breaks = seq(-6,6,2), expand = c(0,0.15), oob = squish) +
  scale_x_continuous(limits = c(-6,6), breaks = seq(-6,6,2), expand = c(0,0.15), oob = squish) +
  labs(x = "Raji LFC", y = "LCL LFC")

# Get theme function
theme_white = theme(text = element_text(color = "white"), axis.text = element_text(color = "white"),
                    axis.line = element_line(color = "white"), axis.ticks = element_line(color = "white"),
                    panel.background = element_rect(fill='transparent'),
                    plot.background = element_rect(fill='transparent', color=NA))

# Get point data
plot_out + theme_white
ggsave(file.path(out_dir, "Supplementary_Figure_6B.png"), dpi = 300, width = 5.21, height = 5.21, units = "in", bg = NULL)

# Get background
background_plot = plot_out
background_plot$layers <- NULL
background_plot +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), label.y = 5, label.x = -5) +
  geom_smooth(method = "lm", color = "purple", linewidth = 1, se = FALSE) +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  geom_abline(slope = 1, intercept = int, color = "black", linetype = "dashed")

ggsave(file.path(out_dir, "Supplementary_Figure_6B.Background.svg"), dpi = 300, width = 5.21, height = 5.21, units = "in", bg = NULL)

## print session info ##
print("Session Info below: ")
sessionInfo()
