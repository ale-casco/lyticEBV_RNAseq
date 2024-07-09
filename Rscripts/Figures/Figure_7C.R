rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("ggplot2")
#install.packages("scales")
#install.packages("svglite")

## Setup the environment

# Import packages
library(rstudioapi)
library(ggplot2)
library(scales)
library(svglite)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Output directory
out_dir=file.path(work_dir, "Figures/Figure_7")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Import DGE file
deg_dir = file.path(work_dir, "Data/DESeq2")
if (!file.exists(list.files(deg_dir,full.names = TRUE)[1])) {stop(" No DESeq2 files exist.\n\tRun DESeq2.R script to generate.", call. = FALSE)} # check file exists
files = list.files(deg_dir, full.names=TRUE)
files = files[grepl("WT_Early.vs.dBGLF5_Early|WT_Late.vs.dBGLF5_Late|.vs.dBGLF5_Latent",files) & !grepl("WT_Latent",files)]

# Import alternative hypothesis DGE file
deg_ah_path = file.path(deg_dir, "AH_ERCCnorm_DE.WT_Latent.vs.dBGLF5_Latent.csv")
if (!file.exists(deg_ah_path)) {stop(" No DESeq2 Alt Hypothesis file exist.\n\tRun DESeq2_AltHypothesis.R script to generate.", call. = FALSE)} # check file exists
deg_ah = read.csv(deg_ah_path)
deg_ah = deg_ah[!is.na(deg_ah$ENSEMBL) & deg_ah$padj <= 0.1, "ENSEMBL"]

# Prepare contrasts
contrast_cols = sub("\\.csv","",sub("ERCCnorm_DE\\.","",basename(files)))
contrast_1 = sub("\\.vs\\..*","",contrast_cols)
contrast_2 = sub(".*\\.vs\\.","",contrast_cols)
contrasts = data.frame(matrix(ncol = length(contrast_cols), nrow = 0))
contrasts = rbind(contrasts,contrast_1)
contrasts = rbind(contrasts,contrast_2)
colnames(contrasts) = contrast_cols
rownames(contrasts)= c("Treatment", "Control")
contrasts
rm(contrast_cols,contrast_1,contrast_2)

# load files
df_list = list()
for (i in 1:length(files)) {
  # Import file
  df_tmp = read.csv(files[i])
  
  # Keep only human genes and subset columns
  df_tmp = df_tmp[grepl("ENSG",df_tmp$ENSEMBL),grepl("ENSEMBL|log2FoldChange|padj|baseMean",colnames(df_tmp))]
  colnames(df_tmp)[-1] = paste0(colnames(contrasts[i]),"_",colnames(df_tmp)[-1])
  
  # Filter genes by alternative hypothesis padj <= 0.1
  df_tmp = df_tmp[df_tmp$ENSEMBL %in% deg_ah,]
  
  # List
  df_list = append(df_list, list(df_tmp))
  names(df_list)[i] = colnames(contrasts)[i]
}

# Merge keeping only genes present in each comparison
df_early <- Reduce(function(...) merge(..., by="ENSEMBL"), df_list[grepl("dBGLF5_Early",names(df_list))])
df_late <- Reduce(function(...) merge(..., by="ENSEMBL"), df_list[grepl("dBGLF5_Late\\b",names(df_list))])

# Simplify col names
colnames(df_early) = sub(".*Early","BGLF5",sub(".*Latent", "BGLF5_independent", colnames(df_early)))
colnames(df_late) = sub(".*Late","BGLF5",sub(".*Latent", "BGLF5_independent", colnames(df_late)))

# Log-Log Plot Function
plot_logLog = function(tmp_df, threshold = 4, title = "", center = FALSE) {
  # Categorize based on adjusted p-value 
  tmp_df$BGLF5_Category = ifelse(tmp_df$BGLF5_log2FoldChange < 0 & tmp_df$BGLF5_padj <= 0.05, "Shutoff","Escpaee")
  tmp_df$BGLF5_Independent_Category = ifelse(tmp_df$BGLF5_independent_log2FoldChange < 0 & tmp_df$BGLF5_independent_padj <= 0.05, "Shutoff","Escpaee")
  tmp_df$Category = ifelse(tmp_df$BGLF5_Category == "Shutoff" & tmp_df$BGLF5_Independent_Category != "Shutoff", "BGLF5",
                           ifelse(tmp_df$BGLF5_Category == "Shutoff" & tmp_df$BGLF5_Independent_Category == "Shutoff", "Both",
                                  ifelse(tmp_df$BGLF5_Category != "Shutoff" & tmp_df$BGLF5_Independent_Category == "Shutoff", "BGLF5-Independent", "Escapee")))
  tmp_df$Category = factor(tmp_df$Category, levels = c("BGLF5","BGLF5-Independent","Both","Escapee"))
  
  # Center data
  tmp_df$BGLF5_log2FoldChange = scale(tmp_df$BGLF5_log2FoldChange,center = T, scale = F)
  tmp_df$BGLF5_independent_log2FoldChange = scale(tmp_df$BGLF5_independent_log2FoldChange,center = T, scale = F)
  
  # Shape and size
  tmp_df$Shape = ifelse(tmp_df$BGLF5_log2FoldChange < -threshold | tmp_df$BGLF5_log2FoldChange > threshold | 
                          tmp_df$BGLF5_independent_log2FoldChange < -threshold | tmp_df$BGLF5_independent_log2FoldChange > threshold, 24, 21)
  tmp_df$Size = ifelse(tmp_df$Shape == 21, 2, 1.5)
  
  # Plot
  ggplot(tmp_df, aes(x = BGLF5_log2FoldChange, y = BGLF5_independent_log2FoldChange, fill = Category)) +
    theme_classic() +
    labs(x = "BGLF5-mediated log2FC", y = "BGLF5-independent log2FC") +
    geom_point(shape = tmp_df$Shape, size=tmp_df$Size,  color = "black", alpha = 0.8) +
    scale_fill_manual(values = c("red","blue","purple","orange")) +
    theme(text = element_text(size = 18), legend.position = "none") +
    scale_y_continuous(limits = c(-threshold,threshold), breaks = seq(-threshold,threshold,threshold/2), expand = c(0,0.15), oob = squish) +
    scale_x_continuous(limits = c(-threshold,threshold), breaks = seq(-threshold,threshold,threshold/2), expand = c(0,0.15), oob = squish)
}

# Generate table
early_plot = plot_logLog(df_early,4, NA,TRUE)
late_plot = plot_logLog(df_late,4, NA,TRUE)
out_table = cbind(data.frame(table(early_plot$data$Category)),data.frame(table(late_plot$data$Category))[,2])
colnames(out_table) = c("Category","Early", "Late")

# Export table
write.csv(out_table,file.path(out_dir,"Figure_7C.table.csv"), row.names = FALSE)

# Export tables with ENSEMBL and host shutoff classifications
write.csv(early_plot$data[,c("ENSEMBL","Category")], file.path(work_dir,"Data/DESeq2/hs_classification.early.csv"), row.names = FALSE)
write.csv(late_plot$data[,c("ENSEMBL","Category")], file.path(work_dir,"Data/DESeq2/hs_classification.late.csv"), row.names = FALSE)

# Get theme function
theme_white = theme(text = element_text(color = "white"), axis.text = element_text(color = "white"),
                    axis.line = element_line(color = "white"), axis.ticks = element_line(color = "white"))

## Get point data
# Early
early_plot + theme_white
ggsave(file.path(out_dir, "Figure_7C.Early.png"), dpi = 300, width = 5.21, height = 5.21, units = "in")

# Late
late_plot + theme_white
ggsave(file.path(out_dir, "Figure_7C.Late.png"), dpi = 300, width = 5.21, height = 5.21, units = "in")

# Background plot
background_plot = early_plot
background_plot$layers <- NULL
background_plot
ggsave(file.path(out_dir, "Figure_7C.Background.svg"), dpi = 300, width = 5.21, height = 5.21, units = "in")

## print session info ##
print("Session Info below: ")
sessionInfo()
