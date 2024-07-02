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

# CPC2 results path
cpc2_path=file.path(work_dir, "Data/IsoformSwitchAnalyzeR/CPC2/cpc2output.txt")

# PFAM results path
pfam_path=file.path(work_dir, "Data/IsoformSwitchAnalyzeR/PFAM/pfam_output.tabular")

# IUPred2A results path
iurpred2a_path=file.path(work_dir, "Data/IsoformSwitchAnalyzeR/IUPred2A/merged.result")

# SignalP6 results path
signalp6_path=file.path(work_dir, "Data/IsoformSwitchAnalyzeR/SignalP/prediction_results.txt")

# DeepLoc2 results path
deeploc2_path=file.path(work_dir, "Data/IsoformSwitchAnalyzeR/DeepLoc2/deeploc_results.csv")

# DeepTMHMM results path
deeptmhmm_path=file.path(work_dir, "Data/IsoformSwitchAnalyzeR/DeepTMHMM/TMRs.gff3")

# Directory containing saved RDS objects
out_dir=file.path(work_dir, "Data/IsoformSwitchAnalyzeR/RDS_objects")

# Read in analyzed RDS object from IsoformSwitchAnalyzeR_external.R script
RDS_path = file.path(out_dir, "aSwitchListAnalyzedExt.RDS")
if (!file.exists(RDS_path)) {stop(" No aSwitchListAnalyzedExt.RDS file exists.\n\tRun IsoformSwitchAnalyzeR_external.R script to generate.", call. = FALSE)} # check file exists
aSwitchListAnalyzedExt = readRDS(RDS_path)

# Add CPC2 analysis
aSwitchListAnalyzed <- analyzeCPC2(
  switchAnalyzeRlist   = aSwitchListAnalyzedExt,
  pathToCPC2resultFile = cpc2_path,
  removeNoncodinORFs   = FALSE   # because ORF was predicted de novo
)

# Add PFAM analysis
aSwitchListAnalyzed <- analyzePFAM(
  switchAnalyzeRlist   = aSwitchListAnalyzed,
  pathToPFAMresultFile = pfam_path,
  showProgress=FALSE
)

# Add IUPred2A analysis
aSwitchListAnalyzed <- analyzeIUPred2A(
  switchAnalyzeRlist        = aSwitchListAnalyzed,
  pathToIUPred2AresultFile = iurpred2a_path,
  showProgress = FALSE
)

# Add SignalP analysis
aSwitchListAnalyzed <- analyzeSignalP(
  switchAnalyzeRlist       = aSwitchListAnalyzed,
  pathToSignalPresultFile  = signalp6_path
)

# Alternative splicing
aSwitchListAnalyzed <- analyzeAlternativeSplicing( aSwitchListAnalyzed )

# Analyze intron retention
aSwitchListAnalyzed <- analyzeIntronRetention( aSwitchListAnalyzed)

# Add DeepLoc2 analysis
aSwitchListAnalyzed <- analyzeDeepLoc2(
  switchAnalyzeRlist = aSwitchListAnalyzed,
  pathToDeepLoc2resultFile = deeploc2_path
)

# Add DeepTMHMM analysis
aSwitchListAnalyzed <- analyzeDeepTMHMM(
  switchAnalyzeRlist   = aSwitchListAnalyzed,
  pathToDeepTMHMMresultFile = deeptmhmm_path,
  showProgress=FALSE
)

# Analyze consequences
aSwitchListAnalyzedConsequences <- analyzeSwitchConsequences(
  aSwitchListAnalyzed,
  consequencesToAnalyze='all',
  alpha=0.05,
  dIFcutoff=0.1,
  onlySigIsoforms=FALSE,
  ntCutoff=50,
  ntFracCutoff=NULL,
  ntJCsimCutoff=0.8,
  AaCutoff=10,
  AaFracCutoff=0.5,
  AaJCsimCutoff=0.9,
  removeNonConseqSwitches=TRUE,
  showProgress=TRUE,
  quiet=FALSE
)

# Save objects as RDS
saveRDS(aSwitchListAnalyzed, file.path(out_dir, "aSwitchListAnalyzed.RDS"))
saveRDS(aSwitchListAnalyzedConsequences, file.path(out_dir, "aSwitchListAnalyzedConsequences.RDS"))

## print session info ##
print("Session Info below: ")
sessionInfo()
