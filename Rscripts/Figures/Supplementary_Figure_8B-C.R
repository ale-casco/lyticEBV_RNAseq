rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("svglite")
#install.packages("ggplot2")
#install.packages("scales")
#install.packages("reshape2")

## Setup the environment

# Import packages
library(rstudioapi)
library(svglite)
library(ggplot2)
library(scales)
library(reshape2)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Output directory
out_dir=file.path(work_dir, "Figures/Supplementary_Figure_8")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# DESeq2 differential alternative splicing events files
de_ase_dir=file.path(work_dir, "Data/SpliceWiz/DESeq2")
if (!file.exists(list.files(de_ase_dir,full.names = TRUE)[1])) {stop(" No ASE DESeq2 files exist.\n\tRun SpliceWiz.R script to generate.", call. = FALSE)} # check file exists

# DESeq2 differential expressed gene (DEG) files
deg_dir=file.path(work_dir, "Data/DESeq2")
if (!file.exists(list.files(deg_dir,full.names = TRUE)[1])) {stop(" No DESeq2 files exist.\n\tRun DESeq2.R script to generate.", call. = FALSE)} # check file exists

# Load differential ASEs
df_early=read.csv(file.path(de_ase_dir, "DE.dBGLF5_Early.vs.dBGLF5_Latent.csv"))
df_late=read.csv(file.path(de_ase_dir, "DE.dBGLF5_Late.vs.dBGLF5_Latent.csv"))

# Load differential expressed genes
deg_early = read.csv(file.path(deg_dir, "ERCCnorm_DE.dBGLF5_Early.vs.dBGLF5_Latent.csv"))
deg_late = read.csv(file.path(deg_dir, "ERCCnorm_DE.dBGLF5_Late.vs.dBGLF5_Latent.csv"))

# Keep only human genes
deg_early = deg_early[grepl("ENSG",deg_early$ENSEMBL),]
deg_late = deg_late[grepl("ENSG",deg_late$ENSEMBL),]

# Filter genes in DEGs
df_early = df_early[,c("EventName","EventType","log2FoldChange","padj","abs_deltaPSI","Inc_gene_id", "Exc_gene_id")]
df_late = df_late[,c("EventName","EventType","log2FoldChange","padj","abs_deltaPSI", "Inc_gene_id", "Exc_gene_id")]
df_early = df_early[(df_early$Inc_gene_id %in% deg_early$ENSEMBL) | (df_early$Exc_gene_id %in% deg_early$ENSEMBL),]
df_late = df_late[(df_late$Inc_gene_id %in% deg_early$ENSEMBL) | (df_late$Exc_gene_id %in% deg_early$ENSEMBL),]

# Combine early and late
df_early$Condition = "dBGLF5_Early"
df_late$Condition = "dBGLF5_Late"
df = rbind(df_early, df_late)

# Volcano plot function
ggvaplot_fun = function(tmp_df,ASE, xlim = 7.5, ylim = 20) {
  
  # Filter only genes with PSI >= 5%
  tmp_df  = tmp_df[!is.na(tmp_df$abs_deltaPSI),]
  tmp_df = tmp_df[tmp_df$abs_deltaPSI >= 0.05,]
  
  # Determine differentially expressed genes
  tmp_df$Sig = ifelse(tmp_df$log2FoldChange >= 1 & tmp_df$padj <= 0.05, "Up",
                      ifelse(tmp_df$log2FoldChange <= -1 & tmp_df$padj <= 0.05, "Down", "NS"))
  
  # If values are beyond the limits, change shape to triangle
  tmp_df$Shape = ifelse(tmp_df$log2FoldChange < -xlim | tmp_df$log2FoldChange > xlim | -log10(tmp_df$padj) > ylim, 24, 21)
  tmp_df$Size = ifelse(tmp_df$Shape == 21, 2, 1.5)
  
  # Filter for ASE
  tmp_df = tmp_df[tmp_df$EventType == ASE,]
  
  # Factor by fraction
  tmp_df$Condition = factor(tmp_df$Condition, levels = c("dBGLF5_Early","dBGLF5_Late"))
  
  
  ggplot(tmp_df, aes(x = log2FoldChange, y = -log10(padj), fill = Sig)) +
    facet_grid(Condition~EventType) +
    #geom_point(alpha = 1, shape = 21, color = "black") +
    geom_point(shape = tmp_df$Shape, size=tmp_df$Size,  color = "black", alpha = 0.8) +
    scale_fill_manual(values = c("#1465AC","darkgray","#B31B21")) +
    scale_y_continuous(limits = c(0,ylim), breaks = seq(0,ylim,5), expand = c(0,0.25), oob = squish) +
    scale_x_continuous(limits = c(-xlim,xlim), breaks = seq(-xlim,xlim,2.5), expand = c(0,0.15), oob = squish) +
    geom_hline(yintercept = -log10(0.05), color = "black") +
    geom_vline(xintercept = 1, color = "black") +
    geom_vline(xintercept = -1, color = "black") +
    theme(legend.position="none",axis.line=element_line(colour="black"),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),panel.border=element_blank(),panel.background=element_blank(),text=element_text(size=18),
          strip.background = element_rect(color="black", fill="white"))
}

# Differential ASE table function
aseTable_fun = function(tmp_df, ASE) {
  tmp_df = ggvaplot_fun(tmp_df, ASE)$data
  tmp_df = as.data.frame(table(tmp_df$Condition,tmp_df$Sig))
  out_table = dcast(tmp_df, Var1 ~ Var2, value.var = "Freq")
  colnames(out_table)[1] = ASE
  return(out_table)
}

# Export ASE table
ir_out_table = aseTable_fun(df, "IR")
se_out_table = aseTable_fun(df, "SE")
write.csv(ir_out_table,file.path(out_dir,"Supplementary_Figure_8B.Table_events.csv"), row.names = FALSE)
write.csv(se_out_table,file.path(out_dir,"Supplementary_Figure_8C.Table_events.csv"), row.names = FALSE)

# Genes with differential ASEs table function
aseGeneTable_fun = function(tmp_df, ASE) {
  tmp_df = ggvaplot_fun(tmp_df, ASE)$data
  tmp_df = tmp_df[tmp_df$Sig != "NS",]
  
  if (ASE == "IR") {tmp_df$SYMBOL = tmp_df$Inc_gene_id}
  if (ASE == "SE") {tmp_df$SYMBOL = tmp_df$Exc_gene_id}
  
  tmp_df1 = tmp_df[tmp_df$Condition == "dBGLF5_Early" & tmp_df$Sig == "Down",]
  tmp_df2 = tmp_df[tmp_df$Condition == "dBGLF5_Early" & tmp_df$Sig == "Up",]
  tmp_df3 = tmp_df[tmp_df$Condition == "dBGLF5_Late" & tmp_df$Sig == "Down",]
  tmp_df4 = tmp_df[tmp_df$Condition == "dBGLF5_Late" & tmp_df$Sig == "Up",]
  
  tmp_df1 = tmp_df1[!duplicated(tmp_df1$SYMBOL),]
  tmp_df2 = tmp_df2[!duplicated(tmp_df2$SYMBOL),]
  tmp_df3 = tmp_df3[!duplicated(tmp_df3$SYMBOL),]
  tmp_df4 = tmp_df4[!duplicated(tmp_df4$SYMBOL),]
  
  tmp_df = rbind(tmp_df1, tmp_df2, tmp_df3, tmp_df4)
  
  tmp_df = as.data.frame(table(tmp_df$Condition,tmp_df$Sig))
  out_table = dcast(tmp_df, Var1 ~ Var2, value.var = "Freq")
  colnames(out_table)[1] = ASE
  return(out_table)
}

# Export ASE Gene table
ir_gene_out_table = aseGeneTable_fun(df, "IR")
se_gene_out_table = aseGeneTable_fun(df, "SE")
write.csv(ir_gene_out_table,file.path(out_dir,"Supplementary_Figure_8B.Table_genes.csv"), row.names = FALSE)
write.csv(se_gene_out_table,file.path(out_dir,"Supplementary_Figure_8C.Table_genes.csv"), row.names = FALSE)

# Get theme function
theme_white = theme(text = element_text(color = "white"), axis.text = element_text(color = "white"),
                    axis.line = element_line(color = "white"), axis.ticks = element_line(color = "white"),
                    strip.background = element_rect(color="white", fill="white"), strip.text = element_text(color = "white"))

# Intron Retention
ggvaplot_fun(df, "IR",10,15) + theme_white
ggsave(file.path(out_dir,"Supplementary_Figure_8B.png"), dpi = 300, width = 4.25, height = 7.52, units = "in")

# Exon Skipping
ggvaplot_fun(df, "SE",10,15) + theme_white
ggsave(file.path(out_dir,"Supplementary_Figure_8C.png"), dpi = 300, width = 4.25, height = 7.52, units = "in")

# Background
background_plot = ggvaplot_fun(df, "IR")
background_plot$layers <- NULL
background_plot
ggsave(file.path(out_dir,"Supplementary_Figure_8B-C.Background.svg"), dpi = 300, width = 4.25, height = 7.52, units = "in")

## print session info ##
print("Session Info below: ")
sessionInfo()
