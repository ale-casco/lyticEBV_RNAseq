rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("svglite")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("IsoformSwitchAnalyzeR")

## Setup the environment

# Import packages
library(rstudioapi)
library(IsoformSwitchAnalyzeR)
library(svglite)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Output directory
out_dir=file.path(work_dir, "Figures/Figure_4")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Directory containing saved RDS objects
RDS_dir=file.path(work_dir, "Data/IsoformSwitchAnalyzeR/RDS_objects/aSwitchListAnalyzedConsequences.RDS")
if (!file.exists(RDS_dir)) {stop(" RDS object does not exist.\n\tRun IsoformSwitchAnalyzeR_external.R and IsoformSwitchAnalyzeR.R scripts to generate.", call. = FALSE)} # check file exists

# Read in analyzed RDS object from IsoformSwitchAnalyzeR_external.R script
aSwitchListAnalyzedConsequences = readRDS(RDS_dir)

# Figure 4F
switchPlotIsoUsage(aSwitchListAnalyzedConsequences, 
           gene = "ILF3", 
           condition1 = "WT_Latent", 
           condition2 = "WT_Early",
           IFcutoff = 0.075 # default = 0.05
)
ggsave(file.path(out_dir, "Figure_4F.svg"), dpi = 300, width = 5.21, height = 5.21, units = "in")

# Figure 4G
extractConsequenceEnrichment(
  aSwitchListAnalyzedConsequences,
  consequencesToAnalyze=c(
    'coding_potential',
    'intron_retention',
    'ORF_seq_similarity',
    'NMD_status',
    'domains_identified',
    'IDR_identified',
    'signal_peptide_identified'
  ),
  analysisOppositeConsequence = TRUE,
)
ggsave(file.path(out_dir, "Figure_4G.svg"), dpi = 300, width = 10.35, height = 3.73, units = "in")

## print session info ##
print("Session Info below: ")
sessionInfo()
