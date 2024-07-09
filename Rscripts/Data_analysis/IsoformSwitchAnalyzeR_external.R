rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("IsoformSwitchAnalyzeR")

## Setup the environment

# Import packages
library(rstudioapi)
library(IsoformSwitchAnalyzeR)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Set directory for sampleinfo matrix file
sampleinfo_dir=file.path(work_dir, "Data/refs/sampleinfo.txt") # 3-column file: id, condition, replicate

# Set directory containing rsem isoforms.results files
counts_dir=file.path(work_dir, "Data/RSEM")

# set path to gtf file
gtf_file=file.path(work_dir, "Data/refs/gencode.v45.primary_assembly.ERCC.M81_DFLR.chrEBV.inverted.gtf.gz")

# set path to rsem_index.transcripts.fa.gz file
fa_file=file.path(work_dir, "Data/refs/rsem_index.transcripts.fa.gz")

# DESeq2 differential expressed gene (DEG) files
deg_dir=file.path(work_dir, "Data/DESeq2")
if (!file.exists(list.files(deg_dir,pattern="^ERCCnorm_DE*",full.name=TRUE)[1])) {
  stop(" No DESeq2 files exist.\n\tRun DESeq2.R script to generate.", call. = FALSE)} # check file exists

# list of genes retained from DESeq2 analysis
deg_files = list.files(deg_dir, full.names=TRUE)

# Comparisons to keep
deg_files = deg_files[grepl("WT_Early.vs.WT_Latent|WT_Late.vs.WT_Latent",deg_files)]

# Directory to output extracted sequences and save RDS objects
out_dir=file.path(work_dir, "Data/IsoformSwitchAnalyzeR")

# Import DESeq2 DEG files into list
SYMBOL_list = list()
for (i in 1:length(deg_files)) {
  df_tmp = read.csv(deg_files[i]) # read deg file
  df_tmp = df_tmp[grepl("ENSG",df_tmp$ENSEMBL),] # keep only human genes
  SYMBOL = df_tmp$SYMBOL # list of SYMBOL IDs
  SYMBOL_list = append(SYMBOL_list, list(SYMBOL)) # append to list
}

# Get vector of SYMBOL genes retained from DESeq2 analysis for downstream analysis
SYMBOL_vec = unique(unlist(SYMBOL_list))
rm(deg_dir,deg_files,SYMBOL_list,df_tmp,i,SYMBOL)

# Import sample treatment matrix
myDesign <- read.delim(file.path(sampleinfo_dir))
colnames(myDesign) <- c("sampleID", "condition", "replicates")
myDesign <- myDesign[,c(-3)]
myDesign <- myDesign[!grepl("Bulk|dBGLF5",myDesign$condition),] # Remove Bulk and dBGLF5 condition

# List RSEM isoform.results files
files = list.files(file.path(counts_dir),pattern = ".isoforms.results", full.names = TRUE)
files = files[grepl(paste(myDesign$sampleID,collapse = "|"),files)]

# Import RSEM data 
rsemQuant <- importIsoformExpression(
  sampleVector = files,
  addIsofomIdAsColumn = TRUE
)

# Create switchAnalyzeRlist
aSwitchList <- importRdata(
  isoformCountMatrix   = rsemQuant$counts,
  isoformRepExpression = rsemQuant$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = gtf_file,
  isoformNtFasta       = fa_file,
  showProgress = FALSE,
  comparisonsToMake = as.data.frame(cbind(condition_1 = c("WT_Latent", "WT_Latent", "WT_Early"),
                                          condition_2 = c("WT_Early", "WT_Late", "WT_Late")),
                                    stringsAsFactors = F)
)

# Subset based on gene list
aSwitchList = subsetSwitchAnalyzeRlist(
  switchAnalyzeRlist = aSwitchList,
  subset = aSwitchList$isoformFeatures$gene_name %in% SYMBOL_vec
)
rm(SYMBOL_vec)

# Filtering - removes things like single-isoform and non-expressed genes that may make downstream processes much longer
aSwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 0,
  isoformExpressionCutoff = 1,
  removeSingleIsoformGenes = TRUE
)

# Test for Isoform Switches via DEXSeq
aSwitchListAnalyzedExt <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = aSwitchListFiltered,
  reduceToSwitchingGenes=TRUE,
  quiet = FALSE,
  showProgress = TRUE
)

# Analyze ORF
aSwitchListAnalyzedExt <- analyzeORF(aSwitchListAnalyzedExt)

# Extract Sequences
dir.create(file.path(out_dir, "extractedSequences"), recursive = TRUE, showWarnings=FALSE)
aSwitchListAnalyzedExt <- extractSequence( aSwitchListAnalyzedExt ,
                                        pathToOutput = file.path(out_dir, "extractedSequences"),
                                        writeToFile=TRUE,
                                        alsoSplitFastaFile=TRUE
)

# Save objects as RDS
dir.create(file.path(out_dir,"RDS_objects"), showWarnings=FALSE)
saveRDS(aSwitchListAnalyzedExt, file.path(out_dir, "RDS_objects/aSwitchListAnalyzedExt.RDS"))

# Run external softwares

## print session info ##
print("Session Info below: ")
sessionInfo()
