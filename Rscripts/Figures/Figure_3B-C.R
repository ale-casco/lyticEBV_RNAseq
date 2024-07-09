rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("ggplot2")
#install.packages("scales")
#install.packages("ggstatsplot")
#install.packages("svglite")

## Setup the environment

# Import packages
library(rstudioapi)
library(ggplot2)
library(scales)
library(ggstatsplot)
library(svglite)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Output directory
out_dir=file.path(work_dir, "Figures/Figure_3")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# DESeq2 differential expressed gene (DEG) files
deg_dir=file.path(work_dir, "Data/DESeq2")
if (!file.exists(list.files(deg_dir,full.names = TRUE)[1])) {stop(" No DESeq2 files exist.\n\tRun DESeq2.R script to generate.", call. = FALSE)} # check file exists
files = list.files(deg_dir, full.names=TRUE)

# Comparisons to keep
files = files[grepl("WT_Early.vs.WT_Latent|WT_Late.vs.WT_Latent",files)]

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

# Import DESeq2 DEG files into list
df_list = list()
for (i in 1:length(files)) {
  df_tmp = read.csv(files[i]) # read deg file
  df_tmp = df_tmp[grepl("ENSG",df_tmp$ENSEMBL),] # keep only human genes
  df_tmp = df_tmp[,c("ENSEMBL","log2FoldChange","padj","baseMean")] # columns to keep
  colnames(df_tmp)[-1] = paste0(colnames(contrasts[i]),"_",colnames(df_tmp)[-1]) # paste comparison to each column name
  df_list = append(df_list, list(df_tmp)) # append to list
}

# Merge list into data frame, keeping all genes
df <- Reduce(function(...) merge(..., by="ENSEMBL", all = TRUE), df_list)

# Prepare list of comparisons
df_list = list()
for (i in 1:length(contrasts)) {
  tmp = df[,grepl(paste0("ENSEMBL|",contrasts[1,i],".vs.",contrasts[2,i]),colnames(df))]  # Subset by comparison
  colnames(tmp)[-1] = c("log2FoldChange", "padj", "baseMean")
  tmp = tmp[!is.na(tmp$baseMean),]
  
  # Append to a list
  df_list=append(df_list,list(tmp))
  names(df_list)[i]=colnames(contrasts)[i]
}

rm(list=setdiff(ls(), c("df_list","contrasts","out_dir")))

################### Figure 2B ###################

# Prepare melted data frame with all comparisons
df = data.frame(matrix(ncol = 5))
colnames(df) = c("ENSEMBL","log2FoldChange","padj","baseMean","Contrast")
for (i in 1:length(contrasts)) {
  tmp_df = df_list[[i]]
  tmp_df$Contrast = paste0(contrasts[1,i], ".vs.", contrasts[2,i]) # Create comparison column
  df = rbind(df,tmp_df)
}
df = df[-1,]

# Factor contrast
df$Contrast = factor(df$Contrast, levels = colnames(contrasts))

# Due to large sample, we will assume normality.
plot_export = ggbetweenstats(data = df, x = Contrast, y = log2FoldChange, type = "p", var.equal = FALSE, p.adjust.method = "fdr", point.args = list(alpha = 0)) +
  scale_y_continuous(limits = c(-8,12.5), breaks = seq(-8,8,2), expand = c(0,0.15), oob = squish)

# Welch’s t-test
t_test.result = t.test(log2FoldChange ~ Contrast, data = df, var.equal = FALSE)
t_test.result

# Export violin plot as SVG
plot_export$labels$title
plot_export$layers[[1]] <- NULL # removes hidden points to decrease file size

plot_export
ggsave(file.path(out_dir,"Figure_3B.svg"), dpi = 10, width = 3.82, height = 4.91, units = "in")

rm(list=setdiff(ls(), c("df_list","contrasts","out_dir")))

################### Figure 2C ###################

# MA Plot function
ggmaplot_fun = function(tmp_df,i) {
  
  # Determine differentially expressed genes
  tmp_df$Sig = ifelse(tmp_df$log2FoldChange >= 0 & tmp_df$padj <= 0.05, "Up",
                       ifelse(tmp_df$log2FoldChange <= 0 & tmp_df$padj <= 0.05, "Down","NS"))
  tmp_df$Sig = ifelse(is.na(tmp_df$Sig), "NS", tmp_df$Sig)
  
  # If values are beyond the limits, change shape to triangle
  tmp_df$Shape = ifelse(tmp_df$log2FoldChange < -8 | tmp_df$log2FoldChange > 8 | log2(tmp_df$baseMean) > 18, 24, 21)
  tmp_df$Size = ifelse(tmp_df$Shape == 21, 2, 1.5)
  
  ggplot(tmp_df, aes(x = log2(baseMean), y = log2FoldChange, fill = Sig)) +
    geom_hline(yintercept=0, color = "black", linewidth = 0.5) +
    geom_point(shape = tmp_df$Shape, size=tmp_df$Size,  color = "black", alpha = 0.8) + scale_fill_manual(values = c("#1465AC","darkgray","#B31B21")) +
    scale_y_continuous(limits = c(-8,8), breaks = seq(-8,8,2), expand = c(0,0.15), oob = squish) +
    scale_x_continuous(limits = c(0,18), breaks = seq(0,18,3), expand = c(0,0.15), oob = squish) +
    labs(x = "Log2 mean normalized counts", y = "Log2 fold change") + theme_bw() +
    theme(legend.position="none",axis.line=element_line(colour="black"),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),panel.border=element_blank(),panel.background=element_blank(),text=element_text(size=18))
}

# Prepare DEG table
if (exists("plot_list")) {  rm(plot_list) }
`%notin%` <- Negate(`%in%`)
deg_table = data.frame(matrix(nrow = 0, ncol = 3))
for (i in 1:length(contrasts)) {
  
  # Get df
  tmp_df <- df_list[[i]]
  
  # run MA plot function
  tmp_plot <- ggmaplot_fun(tmp_df,i)
  
  # Get data from plot (contains DEGs)
  deg_tmp = tmp_plot$data
  
  # Prepare table containing DEG info
  if ("Up" %notin% unique(deg_tmp$Sig)) { Up  = 0 } else { Up = as.numeric(table(deg_tmp$Sig)["Up"]) }
  if ("NS" %notin% unique(deg_tmp$Sig)) { NS = 0 } else { NS = as.numeric(table(deg_tmp$Sig)["NS"]) }
  if ("Down" %notin% unique(deg_tmp$Sig)) { Down = 0 } else { Down = as.numeric(table(deg_tmp$Sig)["Down"]) }
  deg_table = rbind(deg_table,c(Up,NS,Down))
  rownames(deg_table)[i] = colnames(contrasts[i])
  
  if (exists("plot_list") == FALSE) {
    plot_list <- list(tmp_plot)
  } else {
    plot_list[length((plot_list))+1] <- list(tmp_plot)
  }
  
  names(plot_list)[i] = names(df_list)[i]
}
colnames(deg_table) = c("Up","NS","Down")

# Export DEG table
write.csv(deg_table,file.path(out_dir,"Figure_3C.Table.csv"), row.names = TRUE)

## PLOT GENERATION (500 x 500)

# Get theme function
theme_white = theme(text = element_text(color = "white"), axis.text = element_text(color = "white"),
                    axis.line = element_line(color = "white"), axis.ticks = element_line(color = "white"))

## Get point data and save as PNG
# WT Early
plot_list[[1]] + theme_white

ggsave(file.path(out_dir,"Figure_3C.Early.png"), dpi = 300, width = 5.21, height = 5.21, units = "in")

# WT Late
plot_list[[2]] + theme_white
ggsave(file.path(out_dir,"Figure_3C.Late.png"), dpi = 300, width = 5.21, height = 5.21, units = "in")
dev.off()

## Get bounderies+labels and save as SVG
background_plot = plot_list[[1]]
background_plot$layers <- NULL
background_plot
ggsave(file.path(out_dir,"Figure_3C.Background.svg"), dpi = 300, width = 5.21, height = 5.21, units = "in")

# to make final figures, overlay plots

## print session info ##
print("Session Info below: ")
sessionInfo()
