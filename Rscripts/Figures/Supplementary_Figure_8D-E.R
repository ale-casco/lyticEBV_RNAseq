rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("svglite")
#install.packages("ggplot2")
#install.packages("scales")
#install.packages("reshape2")
#install.packages("ggpubr")

## Setup the environment

# Import packages
library(rstudioapi)
library(svglite)
library(ggplot2)
library(scales)
library(reshape2)
library(ggpubr)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Output directory
out_dir=file.path(work_dir, "Figures/Supplementary_Figure_8")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# DESeq2 differential alternative splicing events files
de_ase_dir = file.path(work_dir, "Data/SpliceWiz/DESeq2")
if (!file.exists(list.files(de_ase_dir,full.names = TRUE)[1])) {stop(" No ASE DESeq2 files exist.\n\tRun SpliceWiz.R script to generate.", call. = FALSE)} # check file exists

# DESeq2 differential expressed gene (DEG) files
deg_dir= file.path(work_dir, "Data/DESeq2")
if (!file.exists(list.files(deg_dir,full.names = TRUE)[1])) {stop(" No DESeq2 files exist.\n\tRun DESeq2.R script to generate.", call. = FALSE)} # check file exists

# Load differential ASEs
df_WT_early = read.csv(file.path(de_ase_dir, "DE.WT_Early.vs.WT_Latent.csv"))
df_WT_late = read.csv(file.path(de_ase_dir, "DE.WT_Late.vs.WT_Latent.csv"))
df_dBGLF5_early = read.csv(file.path(de_ase_dir, "DE.dBGLF5_Early.vs.dBGLF5_Latent.csv"))
df_dBGLF5_late = read.csv(file.path(de_ase_dir, "DE.dBGLF5_Late.vs.dBGLF5_Latent.csv"))

# Load differential expressed genes
deg_AH = read.csv(file.path(deg_dir, "AH_ERCCnorm_DE.WT_Latent.vs.dBGLF5_Latent.csv"))
deg_WT_early = read.csv(file.path(deg_dir, "ERCCnorm_DE.WT_Early.vs.WT_Latent.csv"))
deg_WT_late = read.csv(file.path(deg_dir, "ERCCnorm_DE.WT_Late.vs.WT_Latent.csv"))
deg_dBGLF5_early = read.csv(file.path(deg_dir, "ERCCnorm_DE.dBGLF5_Early.vs.dBGLF5_Latent.csv"))
deg_dBGLF5_late = read.csv(file.path(deg_dir, "ERCCnorm_DE.dBGLF5_Late.vs.dBGLF5_Latent.csv"))

# Keep only human genes and non-differential expressed genes between latent WT and dBGLF5
deg_AH = deg_AH[grepl("ENSG",deg_AH$ENSEMBL),]
deg_WT_early = deg_WT_early[grepl("ENSG",deg_WT_early$ENSEMBL) & deg_WT_early$ENSEMBL %in% deg_AH$ENSEMBL,]
deg_WT_late = deg_WT_late[grepl("ENSG",deg_WT_late$ENSEMBL) & deg_WT_late$ENSEMBL %in% deg_AH$ENSEMBL,]
deg_dBGLF5_early = deg_dBGLF5_early[grepl("ENSG",deg_dBGLF5_early$ENSEMBL) & deg_dBGLF5_early$ENSEMBL %in% deg_AH$ENSEMBL,]
deg_dBGLF5_late = deg_dBGLF5_late[grepl("ENSG",deg_dBGLF5_late$ENSEMBL) & deg_dBGLF5_late$ENSEMBL %in% deg_AH$ENSEMBL,]

# Filter genes in DEGs
df_WT_early = df_WT_early[,c("EventName","EventType","log2FoldChange","padj","abs_deltaPSI","Inc_gene_id", "Exc_gene_id")]
df_WT_late = df_WT_late[,c("EventName","EventType","log2FoldChange","padj","abs_deltaPSI", "Inc_gene_id", "Exc_gene_id")]
df_dBGLF5_early = df_dBGLF5_early[,c("EventName","EventType","log2FoldChange","padj","abs_deltaPSI","Inc_gene_id", "Exc_gene_id")]
df_dBGLF5_late = df_dBGLF5_late[,c("EventName","EventType","log2FoldChange","padj","abs_deltaPSI", "Inc_gene_id", "Exc_gene_id")]

df_WT_early = df_WT_early[(df_WT_early$Inc_gene_id %in% deg_WT_early$ENSEMBL) | (df_WT_early$Exc_gene_id %in% deg_WT_early$ENSEMBL),]
df_WT_late = df_WT_late[(df_WT_late$Inc_gene_id %in% deg_WT_late$ENSEMBL) | (df_WT_late$Exc_gene_id %in% deg_WT_late$ENSEMBL),]
df_dBGLF5_early = df_dBGLF5_early[(df_dBGLF5_early$Inc_gene_id %in% deg_dBGLF5_early$ENSEMBL) | (df_dBGLF5_early$Exc_gene_id %in% deg_dBGLF5_early$ENSEMBL),]
df_dBGLF5_late = df_dBGLF5_late[(df_dBGLF5_late$Inc_gene_id %in% deg_dBGLF5_late$ENSEMBL) | (df_dBGLF5_late$Exc_gene_id %in% deg_dBGLF5_late$ENSEMBL),]

# Combine
colnames(df_WT_early)[3:5] = paste0("WT_Early_", colnames(df_WT_early)[3:5])
colnames(df_WT_late)[3:5] = paste0("WT_Late_", colnames(df_WT_late)[3:5])
colnames(df_dBGLF5_early)[3:5] = paste0("dBGLF5_Early_", colnames(df_dBGLF5_early)[3:5])
colnames(df_dBGLF5_late)[3:5] = paste0("dBGLF5_Late_", colnames(df_dBGLF5_late)[3:5])
df_early = merge(df_WT_early, df_dBGLF5_early[,c(1,3:5)], by = "EventName")
df_late = merge(df_WT_late, df_dBGLF5_late[,c(1,3:5)], by = "EventName")

# Remove values with LFC 0 in both conditions
df_early$Remove = ifelse(df_early$WT_Early_log2FoldChange == 0 & df_early$dBGLF5_Early_log2FoldChange == 0, "Remove", "Keep")
df_early = df_early[df_early$Remove != "Remove",]
df_late$Remove = ifelse(df_late$WT_Late_log2FoldChange == 0 & df_late$dBGLF5_Late_log2FoldChange == 0, "Remove", "Keep")
df_late = df_late[df_late$Remove != "Remove",]

# Subset intron retention and exon skipping events only
df_early = df_early[df_early$EventType == "IR" | df_early$EventType == "SE",]
df_late = df_late[df_late$EventType == "IR" | df_late$EventType == "SE",]

# Combine
colnames(df_early) = sub("Early_","", colnames(df_early))
colnames(df_late) = sub("Late_","", colnames(df_late))
df_early$Condition = "Early"
df_late$Condition = "Late"
df = rbind(df_early, df_late)

# If values are beyond the limits, change shape to triangle
df$Shape = ifelse(abs(df$WT_log2FoldChange) > 10 | abs(df$dBGLF5_log2FoldChange) > 10, 24, 21)
df$Size = ifelse(df$Shape == 21, 2, 1.5)

# Split
df_IR = df[df$EventType == "IR",]
df_SE = df[df$EventType == "SE",]

plot_IR = ggplot(df_IR, aes(x = dBGLF5_log2FoldChange, y = WT_log2FoldChange)) +
  geom_point(shape = df_IR$Shape, size=df_IR$Size,  color = "black", fill = "black", alpha = 0.8) +
  #geom_smooth(method = "lm", se = FALSE) +
  #geom_abline(slope = 1, color = "red", linetype = "dashed") +
  facet_grid(Condition ~ .) +
  scale_y_continuous(limits = c(-10,10), breaks = seq(-10,10,5), expand = c(0,0.15), oob = squish) +
  scale_x_continuous(limits = c(-10,10), breaks = seq(-10,10,5), expand = c(0,0.15), oob = squish) +
  theme(legend.position="none",axis.line=element_line(colour="black"),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),panel.border=element_blank(),panel.background=element_blank(),text=element_text(size=18),
        strip.background = element_rect(color="black", fill="white"))

plot_SE = ggplot(df_SE, aes(x = dBGLF5_log2FoldChange, y = WT_log2FoldChange)) +
  geom_point(shape = df_SE$Shape, size=df_SE$Size,  color = "black", fill = "black", alpha = 0.8) +
  #geom_smooth(method = "lm", se = FALSE) +
  #geom_abline(slope = 1, color = "red", linetype = "dashed") +
  facet_grid(Condition ~ .) +
  scale_y_continuous(limits = c(-10,10), breaks = seq(-10,10,5), expand = c(0,0.15), oob = squish) +
  scale_x_continuous(limits = c(-10,10), breaks = seq(-10,10,5), expand = c(0,0.15), oob = squish) +
  theme(legend.position="none",axis.line=element_line(colour="black"),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),panel.border=element_blank(),panel.background=element_blank(),text=element_text(size=18),
        strip.background = element_rect(color="black", fill="white"))

# Get theme function
theme_white = theme(text = element_text(color = "white"), axis.text = element_text(color = "white"),
                    axis.line = element_line(color = "white"), axis.ticks = element_line(color = "white"),
                    strip.background = element_rect(color="white", fill="white"), strip.text = element_text(color = "white"))

# Intron Retention
plot_IR + theme_white
ggsave(file.path(out_dir,"Supplementary_Figure_8D.png"), dpi = 300, width = 4.34, height = 7.27, units = "in")

# Exon Skipping
plot_SE + theme_white
ggsave(file.path(out_dir,"Supplementary_Figure_8E.png"), dpi = 300, width = 4.34, height = 7.27, units = "in")

# Background
background_plot = plot_IR
background_plot$layers <- NULL
background_plot + geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), label.y = 7, label.x = -1) +
  geom_abline(slope = 1, color = "red", linetype = "dashed")
ggsave(file.path(out_dir,"Supplementary_Figure_8D.Background.svg"), dpi = 300, width = 4.34, height = 7.27, units = "in")

background_plot = plot_SE
background_plot$layers <- NULL
background_plot + geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), label.y = 7, label.x = -1) +
  geom_abline(slope = 1, color = "red", linetype = "dashed")
ggsave(file.path(out_dir,"Supplementary_Figure_8E.Background.svg"), dpi = 300, width = 4.34, height = 7.27, units = "in")

## print session info ##
print("Session Info below: ")
sessionInfo()
