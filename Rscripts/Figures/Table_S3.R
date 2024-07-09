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

# DESeq2 differential expressed gene (DEG) files
deg_dir=file.path(work_dir, "Data/DESeq2")
if (!file.exists(list.files(deg_dir,pattern="^ERCCnorm_DE*",full.name=TRUE)[1])) {
  stop(" No DESeq2 files exist.\n\tRun DESeq2.R script to generate.", call. = FALSE)} # check file exists
files = list.files(deg_dir, full.names=TRUE, pattern = "*.csv")

# Rename files
names(files) = gsub("_"," ",paste0("(",sub(".vs.",")vs(",sub(".csv","",sub("ERCCnorm_DE.","",basename(files)))),")"))
names(files) = ifelse(grepl("AH ",names(files)),paste0(sub("AH ","",names(files))," Ha"),names(files))

# Reorder files
names(files)
files = files[c(5,7,10,9,2,4,3,6,8,1)]

# Import DESeq2 DEG files into list
wb = createWorkbook() # create empty Excel workbook object
for (i in 1:length(files)) {
  df_tmp = read.csv(files[i], check.names=FALSE) # read deg file
  df_tmp = df_tmp[order(df_tmp$log2FoldChange, decreasing = TRUE),] # Order by decreasing log2FoldChange
  df_host = df_tmp[grepl("ENSG", df_tmp$ENSEMBL),] # get human genes
  df_ebv = df_tmp[!grepl("ENSG", df_tmp$ENSEMBL) & !grepl("ERCC-", df_tmp$SYMBOL),] # get EBV genes
  df_ercc = df_tmp[grepl("ERCC-", df_tmp$SYMBOL),] # get ERCC genes
  
  # Replace gene biotypes
  df_ebv$gene_type = ifelse(grepl("BART|EBER", df_ebv$gene_type), "ncRNA", "protein_coding")
  df_ercc$gene_type = "spike_in"

  # Replace ENSEMBLE gene IDs in SYMBOL column with NA
  df_host$SYMBOL = ifelse(grepl("ENSG",df_host$SYMBOL), NA, df_host$SYMBOL)
  
  # rbind
  df_out = rbind(df_host, df_ebv, df_ercc)

  # Add data into Excel workbook object
  addWorksheet(wb, names(files)[i])
  writeData(wb, names(files)[i], df_out, rowNames = FALSE)
}

# Export
dir.create(file.path(work_dir, "Figures"), showWarnings=FALSE)
saveWorkbook(wb, file = file.path(work_dir, "Figures/Table_S3.xlsx"), overwrite = TRUE)

## print session info ##
print("Session Info below: ")
sessionInfo()
