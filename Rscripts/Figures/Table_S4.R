rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("openxlsx")

## Setup the environment

# Import packages
library(rstudioapi)
library(openxlsx)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Import DGE file
deg_dir = file.path(work_dir, "Data/DESeq2")
if (!file.exists(list.files(deg_dir,pattern="^ERCCnorm_DE*",full.name=TRUE)[1])) {
  stop(" No DESeq2 files exist.\n\tRun DESeq2.R script to generate.", call. = FALSE)} # check file exists
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
  df_tmp = df_tmp[grepl("ENSG",df_tmp$ENSEMBL),grepl("ENSEMBL|SYMBOL|log2FoldChange|padj",colnames(df_tmp))]
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
colnames(df_early) = sub(".*Early","BGLF5_dependent",sub(".*Latent", "BGLF5_independent", colnames(df_early)))
colnames(df_late) = sub(".*Late","BGLF5_dependent",sub(".*Latent", "BGLF5_independent", colnames(df_late)))

# Keep one SYMBOL column
symbol_cols = which(grepl("SYMBOL",colnames(df_early)))
df_early = df_early[,-symbol_cols[2]]
df_late = df_late[,-symbol_cols[2]]
colnames(df_early)[symbol_cols[1]] = "SYMBOL"
colnames(df_late)[symbol_cols[1]] = "SYMBOL"

# Host shutoff classification function
hs_class_fun = function(tmp_df) {
  # Categorize based on adjusted p-value 
  tmp_df$BGLF5_dependent_Category = ifelse(tmp_df$BGLF5_dependent_log2FoldChange < 0 & tmp_df$BGLF5_dependent_padj <= 0.05, "Shutoff","Escpaee")
  tmp_df$BGLF5_Independent_Category = ifelse(tmp_df$BGLF5_independent_log2FoldChange < 0 & tmp_df$BGLF5_independent_padj <= 0.05, "Shutoff","Escpaee")
  tmp_df$HS_class = ifelse(tmp_df$BGLF5_dependent_Category == "Shutoff" & tmp_df$BGLF5_Independent_Category != "Shutoff", "BGLF5-dependent",
                           ifelse(tmp_df$BGLF5_dependent_Category == "Shutoff" & tmp_df$BGLF5_Independent_Category == "Shutoff", "BGLF5-dependent + BGLF5-independent",
                                  ifelse(tmp_df$BGLF5_dependent_Category != "Shutoff" & tmp_df$BGLF5_Independent_Category == "Shutoff", "BGLF5-independent", "Escapee")))
  tmp_df = tmp_df[,-7:-8]
  
  for (i in c("BGLF5-dependent", "BGLF5-independent", "BGLF5-dependent + BGLF5-independent", "Escapee")) {
    tmp_class_df = tmp_df[tmp_df$HS_class == i,]
    if (i == "BGLF5-dependent" | i == "BGLF5-independent + BGLF5-Independent") {
      tmp_class_df = tmp_class_df[order(tmp_class_df$BGLF5_dependent_log2FoldChange, decreasing = FALSE), ] # Order by increasing BGLF5 log2FoldChange
    } else if (i == "BGLF5-independent") {
      tmp_class_df = tmp_class_df[order(tmp_class_df$BGLF5_independent_log2FoldChange, decreasing = FALSE), ] # Order by increasing BGLF5-Independent log2FoldChange
    } else if (i == "Escapee") {
      tmp_class_df = tmp_class_df[order(tmp_class_df$BGLF5_dependent_log2FoldChange, decreasing = TRUE), ] # Order by decreasing BGLF5 log2FoldChange
    }
    if (i == "BGLF5-dependent") {tmp_out_df = tmp_class_df} else {tmp_out_df = rbind(tmp_out_df, tmp_class_df)}
  }
  tmp_out_df = tmp_out_df[,c("ENSEMBL","SYMBOL","HS_class","BGLF5_dependent_log2FoldChange","BGLF5_dependent_padj","BGLF5_independent_log2FoldChange","BGLF5_independent_padj")]
  colnames(tmp_out_df) = sub("_"," ",sub("_independent","-independent",sub("_dependent", "-dependent", colnames(tmp_out_df))))
  tmp_out_df$SYMBOL = ifelse(grepl("ENSG",tmp_out_df$SYMBOL), NA, tmp_out_df$SYMBOL)
  return(tmp_out_df)
}

# Generate output data frames
early_out = hs_class_fun(df_early)
late_out = hs_class_fun(df_late)


# Add data into Excel workbook object
wb = createWorkbook() # create empty Excel workbook object
addWorksheet(wb, "Early lytic")
writeData(wb, "Early lytic", early_out, rowNames = FALSE)
addWorksheet(wb, "Late lytic")
writeData(wb, "Late lytic", late_out, rowNames = FALSE)

# Export
dir.create(file.path(work_dir, "Figures"), showWarnings=FALSE)
saveWorkbook(wb, file = file.path(work_dir, "Figures/Table_S4.xlsx"), overwrite = TRUE)

## print session info ##
print("Session Info below: ")
sessionInfo()
