rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("SpliceWiz")
#BiocManager::install("DESeq2")
#BiocManager::install("Biostrings")

## Setup the environment

# Import packages
library(rstudioapi)
library(SpliceWiz)
library(Biostrings)
library(DESeq2)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Output directory
out_dir=file.path(work_dir,"Data/SpliceWiz")

# BAM files directory
bams_dir=file.path(work_dir,"Data/bams")

# FASTA file
fasta_path=file.path(work_dir,"Data/refs/GRCh38.p14.ERCC.M81_DFLR.chrEBV.inverted.fa.gz")

# GTF file
gtf_dir=file.path(work_dir,"Data/refs/gencode.v45.primary_assembly.ERCC.M81_DFLR.chrEBV.inverted.gtf.gz")

# sample info file
sampleinfo_dir=file.path(work_dir,"Data/refs/sampleinfo.txt")

# gene info dir
geneinfo_dir=file.path(work_dir,"Data/refs/gene_info.csv")
if (!file.exists(geneinfo_dir)) {stop(" gene_info.csv file does not exist.\n\tRun gtf2gene_info.R script to generate.", call. = FALSE)} # check file exists

## Fix abmiguities in fasta file

# Read FASTA file
dna_sequences = readDNAStringSet(fasta_path)

# Replace ambiguities in FASTA file to ensure sequences only contain A, C, G, T, or N
clean_sequences = replaceAmbiguities(dna_sequences) # replaceAmbiguities() is a simple wrapper around chartr() that replaces all IUPAC ambiguities with N.

# Write the cleaned sequences to a new FASTA file
fasta_clean_path = file.path(dirname(fasta_path),sub("\\.fa","\\.no_ambiguities\\.fa",basename(fasta_path)))
writeXStringSet(clean_sequences, fasta_clean_path)
rm(dna_sequences,clean_sequences,fasta_path)

## Set directories within out_dir
ref_dir = file.path(out_dir, "ref")
WT_bams_dir = file.path(bams_dir, "WT")
dBGLF5_bams_dir = file.path(bams_dir, "dBGLF5")
WT_out_dir = file.path(out_dir, "output","WT")
dBGLF5_out_dir = file.path(out_dir, "output","dBGLF5")
NxtSE_dir = file.path(out_dir, "NxtSE")
WT_nxtse_dir = file.path(NxtSE_dir, "WT_NxtSE")
dBGLF5_nxtse_dir = file.path(NxtSE_dir, "dBGLF5_NxtSE")
#WT_nxtse_novel_dir = file.path(NxtSE_dir, "WT_NxtSE_novel")
#dBGLF5_nxtse_novel_dir = file.path(NxtSE_dir, "dBGLF5_NxtSE_novel")
DESeq2_dir = file.path(out_dir, "DESeq2")

## Create directories
dir.create(ref_dir, showWarnings=FALSE)
dir.create(file.path(out_dir, "output"), showWarnings=FALSE)
dir.create(WT_out_dir, showWarnings=FALSE)
dir.create(dBGLF5_out_dir, showWarnings=FALSE)
dir.create(NxtSE_dir,showWarnings=FALSE)
dir.create(WT_nxtse_dir,showWarnings=FALSE)
dir.create(dBGLF5_nxtse_dir,showWarnings=FALSE)
#dir.create(WT_nxtse_novel_dir,showWarnings=FALSE)
#dir.create(dBGLF5_nxtse_novel_dir,showWarnings=FALSE)
dir.create(DESeq2_dir,showWarnings=FALSE)

## Build the SpliceWiz reference
buildRef(
  reference_path = ref_dir,
  fasta = fasta_clean_path,
  gtf = gtf_dir,
  genome_type = "hg38"
)

## get bams
WT_bams = findBAMS(WT_bams_dir, level = 0) 
dBGLF5_bams = findBAMS(dBGLF5_bams_dir, level = 0) 
WT_bams$sample = sub(".Aligned.sortedByCoord.out","",WT_bams$sample)
dBGLF5_bams$sample = sub(".Aligned.sortedByCoord.out","",dBGLF5_bams$sample)

## Process WT BAM files using SpliceWiz
processBAM(
  bamfiles = WT_bams$path,
  sample_names = WT_bams$sample,
  reference_path = ref_dir,
  output_path = WT_out_dir,
  n_threads = 4,
  overwrite = FALSE,
  run_featureCounts = FALSE #  runs featureCounts to obtain gene counts (which outputs results as an RDS file)
)

## Process dBGLF5 BAM files using SpliceWiz
processBAM(
  bamfiles = dBGLF5_bams$path,
  sample_names = dBGLF5_bams$sample,
  reference_path = ref_dir,
  output_path = dBGLF5_out_dir,
  n_threads = 4,
  overwrite = FALSE,
  run_featureCounts = FALSE #  runs featureCounts to obtain gene counts (which outputs results as an RDS file)
)

# Organize output files of SpliceWiz's processBAM() function
WT_expr = findSpliceWizOutput(WT_out_dir)
dBGLF5_expr = findSpliceWizOutput(dBGLF5_out_dir)

## Collate the experiments

# WT
collateData(
  Experiment = WT_expr,
  reference_path = ref_dir,
  output_path = WT_nxtse_dir
)

# dBGLF5
collateData(
  Experiment = dBGLF5_expr,
  reference_path = ref_dir,
  output_path = dBGLF5_nxtse_dir
)

## Collate the experiments with novel ASE discovery

# WT
#collateData(
#  Experiment = WT_expr,
#  reference_path = ref_dir,
#  output_path = WT_nxtse_novel_dir,
#  novelSplicing = TRUE
#)

# dBGLF5
#collateData(
#  Experiment = dBGLF5_expr,
#  reference_path = ref_dir,
#  output_path = dBGLF5_nxtse_novel_dir,
#  novelSplicing = TRUE
#)

## Differential analysis

# Import the experiments
WT_se = makeSE(WT_nxtse_dir)
dBGLF5_se = makeSE(dBGLF5_nxtse_dir)
#WT_novel_se = makeSE(WT_nxtse_novel_dir)
#dBGLF5_novel_se = makeSE(dBGLF5_nxtse_novel_dir)

# Load sampleinfo file
sampleinfo <- read.delim(file.path(sampleinfo_dir))
colnames(WT_se)[match(colnames(WT_se),sampleinfo$UNIQUE_ID)]

# Assign annotations to samples
colData(WT_se)$condition = sampleinfo[which(sampleinfo$UNIQUE_ID %in% colnames(WT_se)),2]
colData(dBGLF5_se)$condition = sampleinfo[which(sampleinfo$UNIQUE_ID %in% colnames(dBGLF5_se)),2]
#colData(WT_novel_se)$condition = sampleinfo[which(sampleinfo$UNIQUE_ID %in% colnames(WT_novel_se)),2]
#colData(dBGLF5_novel_se)$condition = sampleinfo[which(sampleinfo$UNIQUE_ID %in% colnames(dBGLF5_novel_se)),2]

# Filter for high-confidence events
WT_se.filtered = WT_se[applyFilters(WT_se),]
dBGLF5_se.filtered = dBGLF5_se[applyFilters(dBGLF5_se),]
#WT_novel_se.filtered = WT_novel_se[applyFilters(WT_novel_se),]
#dBGLF5_novel_se.filtered = dBGLF5_novel_se[applyFilters(dBGLF5_novel_se),]

# Perform differential analysis with DESeq2 (WT early vs. latent)
WT_EvL_res = ASE_DESeq(
  se = WT_se.filtered,
  test_factor = "condition",
  test_nom = "WT_Early",
  test_denom = "WT_Latent",
  n_threads = 4
)

# Perform differential analysis with DESeq2 (dBGLF5 early vs. latent)
dBGLF5_EvL_res = ASE_DESeq(
  se = dBGLF5_se.filtered,
  test_factor = "condition",
  test_nom = "dBGLF5_Early",
  test_denom = "dBGLF5_Latent",
  n_threads = 4
)

# Perform differential analysis with DESeq2 (WT early vs. latent) including novel splices
#WT_EvL_novel_res = ASE_DESeq(
#  se = WT_novel_se.filtered,
#  test_factor = "condition",
#  test_nom = "WT_Early",
#  test_denom = "WT_Latent",
#  n_threads = 4
#)

# Perform differential analysis with DESeq2 (dBGLF5 early vs. latent) including novel splices
#dBGLF5_EvL_novel_res = ASE_DESeq(
#  se = dBGLF5_novel_se.filtered,
#  test_factor = "condition",
#  test_nom = "dBGLF5_Early",
#  test_denom = "dBGLF5_Latent",
#  n_threads = 4
#)

# Perform differential analysis with DESeq2 (WT Late vs. latent)
WT_LvL_res = ASE_DESeq(
  se = WT_se.filtered,
  test_factor = "condition",
  test_nom = "WT_Late",
  test_denom = "WT_Latent",
  n_threads = 4
)

# Perform differential analysis with DESeq2 (dBGLF5 Late vs. latent)
dBGLF5_LvL_res = ASE_DESeq(
  se = dBGLF5_se.filtered,
  test_factor = "condition",
  test_nom = "dBGLF5_Late",
  test_denom = "dBGLF5_Latent",
  n_threads = 4
)

# Perform differential analysis with DESeq2 (WT Late vs. latent) including novel splices
#WT_LvL_novel_res = ASE_DESeq(
#  se = WT_novel_se.filtered,
#  test_factor = "condition",
#  test_nom = "WT_Late",
#  test_denom = "WT_Latent",
#  n_threads = 4
#)

# Perform differential analysis with DESeq2 (dBGLF5 Late vs. latent) including novel splices
#dBGLF5_LvL_novel_res = ASE_DESeq(
#  se = dBGLF5_novel_se.filtered,
#  test_factor = "condition",
#  test_nom = "dBGLF5_Late",
#  test_denom = "dBGLF5_Latent",
#  n_threads = 4
#)

# Perform differential analysis with DESeq2 (WT late vs. early)
WT_LvE_res = ASE_DESeq(
  se = WT_se.filtered,
  test_factor = "condition",
  test_nom = "WT_Late",
  test_denom = "WT_Early",
  n_threads = 4
)

# Perform differential analysis with DESeq2 (WT late vs. early) including novel splices
#WT_LvE_novel_res = ASE_DESeq(
#  se = WT_novel_se.filtered,
#  test_factor = "condition",
#  test_nom = "WT_Late",
#  test_denom = "WT_Early",
#  n_threads = 4
#)

# Include gene ids in output
df <- viewASE(ref_dir) # save SpliceWIz reference as data frame
df$Inc_gene_id = sub("\\..*","",df$Inc_gene_id)
df$Exc_gene_id = sub("\\..*","",df$Exc_gene_id)

WT_EvL_res = merge(x = WT_EvL_res, y = df[,c("EventName","Inc_gene_id","Exc_gene_id")], by = "EventName", all.x = TRUE)
dBGLF5_EvL_res = merge(x = dBGLF5_EvL_res, y = df[,c("EventName","Inc_gene_id","Exc_gene_id")], by = "EventName", all.x = TRUE)
#WT_EvL_novel_res = merge(x = WT_EvL_novel_res, y = df[,c("EventName","Inc_gene_id","Exc_gene_id")], by = "EventName", all.x = TRUE)
#dBGLF5_EvL_novel_res = merge(x = dBGLF5_EvL_novel_res, y = df[,c("EventName","Inc_gene_id","Exc_gene_id")], by = "EventName", all.x = TRUE)
WT_LvL_res = merge(x = WT_LvL_res, y = df[,c("EventName","Inc_gene_id","Exc_gene_id")], by = "EventName", all.x = TRUE)
dBGLF5_LvL_res = merge(x = dBGLF5_LvL_res, y = df[,c("EventName","Inc_gene_id","Exc_gene_id")], by = "EventName", all.x = TRUE)
#WT_LvL_novel_res = merge(x = WT_LvL_novel_res, y = df[,c("EventName","Inc_gene_id","Exc_gene_id")], by = "EventName", all.x = TRUE)
#dBGLF5_LvL_novel_res = merge(x = dBGLF5_LvL_novel_res, y = df[,c("EventName","Inc_gene_id","Exc_gene_id")], by = "EventName", all.x = TRUE)
WT_LvE_res = merge(x = WT_LvE_res, y = df[,c("EventName","Inc_gene_id","Exc_gene_id")], by = "EventName", all.x = TRUE)
#WT_LvE_novel_res = merge(x = WT_LvE_novel_res, y = df[,c("EventName","Inc_gene_id","Exc_gene_id")], by = "EventName", all.x = TRUE)

## Add ENSEMBL IDs to intron retention events

# Import gene info
gene_info <- read.csv(geneinfo_dir)
gene_info = gene_info[!duplicated(gene_info$SYMBOL),]

# WT EvL
WT_EvL_res$SYMBOL = WT_EvL_res$Inc_gene_id
WT_EvL_res$SYMBOL[WT_EvL_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",WT_EvL_res$EventName[WT_EvL_res$EventType == "IR"])))
WT_EvL_res = merge(x = WT_EvL_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
WT_EvL_res$ENSEMBL = ifelse(is.na(WT_EvL_res$ENSEMBL), WT_EvL_res$Inc_gene_id, WT_EvL_res$ENSEMBL)
WT_EvL_res$Inc_gene_id = WT_EvL_res$ENSEMBL
WT_EvL_res = WT_EvL_res[,-30]
WT_EvL_res$SYMBOL = WT_EvL_res$Exc_gene_id
WT_EvL_res$SYMBOL[WT_EvL_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",WT_EvL_res$EventName[WT_EvL_res$EventType == "IR"])))
WT_EvL_res = merge(x = WT_EvL_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
WT_EvL_res$ENSEMBL = ifelse(is.na(WT_EvL_res$ENSEMBL), WT_EvL_res$Exc_gene_id, WT_EvL_res$ENSEMBL)
WT_EvL_res$Exc_gene_id = WT_EvL_res$ENSEMBL
WT_EvL_res = WT_EvL_res[,-30]

# dBGLF5 EvL
dBGLF5_EvL_res$SYMBOL = dBGLF5_EvL_res$Inc_gene_id
dBGLF5_EvL_res$SYMBOL[dBGLF5_EvL_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",dBGLF5_EvL_res$EventName[dBGLF5_EvL_res$EventType == "IR"])))
dBGLF5_EvL_res = merge(x = dBGLF5_EvL_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
dBGLF5_EvL_res$ENSEMBL = ifelse(is.na(dBGLF5_EvL_res$ENSEMBL), dBGLF5_EvL_res$Inc_gene_id, dBGLF5_EvL_res$ENSEMBL)
dBGLF5_EvL_res$Inc_gene_id = dBGLF5_EvL_res$ENSEMBL
dBGLF5_EvL_res = dBGLF5_EvL_res[,-30]
dBGLF5_EvL_res$SYMBOL = dBGLF5_EvL_res$Exc_gene_id
dBGLF5_EvL_res$SYMBOL[dBGLF5_EvL_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",dBGLF5_EvL_res$EventName[dBGLF5_EvL_res$EventType == "IR"])))
dBGLF5_EvL_res = merge(x = dBGLF5_EvL_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
dBGLF5_EvL_res$ENSEMBL = ifelse(is.na(dBGLF5_EvL_res$ENSEMBL), dBGLF5_EvL_res$Exc_gene_id, dBGLF5_EvL_res$ENSEMBL)
dBGLF5_EvL_res$Exc_gene_id = dBGLF5_EvL_res$ENSEMBL
dBGLF5_EvL_res = dBGLF5_EvL_res[,-30]

# WT EvL novel
#WT_EvL_novel_res$SYMBOL = WT_EvL_novel_res$Inc_gene_id
#WT_EvL_novel_res$SYMBOL[WT_EvL_novel_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",WT_EvL_novel_res$EventName[WT_EvL_novel_res$EventType == "IR"])))
#WT_EvL_novel_res = merge(x = WT_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_EvL_novel_res$ENSEMBL = ifelse(is.na(WT_EvL_novel_res$ENSEMBL), WT_EvL_novel_res$Inc_gene_id, WT_EvL_novel_res$ENSEMBL)
#WT_EvL_novel_res$Inc_gene_id = WT_EvL_novel_res$ENSEMBL
#WT_EvL_novel_res = WT_EvL_novel_res[,-30]
#WT_EvL_novel_res$SYMBOL = WT_EvL_novel_res$Exc_gene_id
#WT_EvL_novel_res$SYMBOL[WT_EvL_novel_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",WT_EvL_novel_res$EventName[WT_EvL_novel_res$EventType == "IR"])))
#WT_EvL_novel_res = merge(x = WT_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_EvL_novel_res$ENSEMBL = ifelse(is.na(WT_EvL_novel_res$ENSEMBL), WT_EvL_novel_res$Exc_gene_id, WT_EvL_novel_res$ENSEMBL)
#WT_EvL_novel_res$Exc_gene_id = WT_EvL_novel_res$ENSEMBL
#WT_EvL_novel_res = WT_EvL_novel_res[,-30]
#WT_EvL_novel_res$SYMBOL[is.na(WT_EvL_novel_res$Inc_gene_id)] = sub("\\..*","",sub("\\-novel.*","",sub("\\-2.*","",sub(".*;","",WT_EvL_novel_res$EventName[is.na(WT_EvL_novel_res$Inc_gene_id)]))))
#WT_EvL_novel_res = merge(x = WT_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_EvL_novel_res$Inc_gene_id  = ifelse(is.na(WT_EvL_novel_res$Inc_gene_id), WT_EvL_novel_res$ENSEMBL, WT_EvL_novel_res$Inc_gene_id)
#WT_EvL_novel_res = WT_EvL_novel_res[,-30]
#WT_EvL_novel_res$SYMBOL[is.na(WT_EvL_novel_res$Exc_gene_id)] = sub("\\..*","",sub("\\-novel.*","",sub("\\-2.*","",sub(".*;","",WT_EvL_novel_res$EventName[is.na(WT_EvL_novel_res$Exc_gene_id)]))))
#WT_EvL_novel_res = merge(x = WT_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_EvL_novel_res$Exc_gene_id  = ifelse(is.na(WT_EvL_novel_res$Exc_gene_id), WT_EvL_novel_res$ENSEMBL, WT_EvL_novel_res$Exc_gene_id)
#WT_EvL_novel_res = WT_EvL_novel_res[,-30]
#WT_EvL_novel_res$SYMBOL[is.na(WT_EvL_novel_res$Inc_gene_id)] = sub("_.*","",sub(".*;","",WT_EvL_novel_res$EventName[is.na(WT_EvL_novel_res$Inc_gene_id)]))
#WT_EvL_novel_res = merge(x = WT_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_EvL_novel_res$Inc_gene_id  = ifelse(is.na(WT_EvL_novel_res$Inc_gene_id), WT_EvL_novel_res$ENSEMBL, WT_EvL_novel_res$Inc_gene_id)
#WT_EvL_novel_res = WT_EvL_novel_res[,-30]
#WT_EvL_novel_res$SYMBOL[is.na(WT_EvL_novel_res$Exc_gene_id)] = sub("_.*","",sub(".*;","",WT_EvL_novel_res$EventName[is.na(WT_EvL_novel_res$Exc_gene_id)]))
#WT_EvL_novel_res = merge(x = WT_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_EvL_novel_res$Exc_gene_id  = ifelse(is.na(WT_EvL_novel_res$Exc_gene_id), WT_EvL_novel_res$ENSEMBL, WT_EvL_novel_res$Exc_gene_id)
#WT_EvL_novel_res = WT_EvL_novel_res[,-30]
#WT_EvL_novel_res$SYMBOL[is.na(WT_EvL_novel_res$Inc_gene_id)] = sub("_.*","",sub(".*:","",WT_EvL_novel_res$EventName[is.na(WT_EvL_novel_res$Inc_gene_id)]))
#WT_EvL_novel_res = merge(x = WT_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_EvL_novel_res$Inc_gene_id  = ifelse(is.na(WT_EvL_novel_res$Inc_gene_id), WT_EvL_novel_res$ENSEMBL, WT_EvL_novel_res$Inc_gene_id)
#WT_EvL_novel_res = WT_EvL_novel_res[,-30]
#WT_EvL_novel_res$SYMBOL[is.na(WT_EvL_novel_res$Exc_gene_id)] = sub("_.*","",sub(".*:","",WT_EvL_novel_res$EventName[is.na(WT_EvL_novel_res$Exc_gene_id)]))
#WT_EvL_novel_res = merge(x = WT_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_EvL_novel_res$Exc_gene_id  = ifelse(is.na(WT_EvL_novel_res$Exc_gene_id), WT_EvL_novel_res$ENSEMBL, WT_EvL_novel_res$Exc_gene_id)
#WT_EvL_novel_res = WT_EvL_novel_res[,-30]
#WT_EvL_novel_res$SYMBOL[is.na(WT_EvL_novel_res$Exc_gene_id)]  =  sub("\\-.*","",WT_EvL_novel_res$SYMBOL[is.na(WT_EvL_novel_res$Exc_gene_id)])
#WT_EvL_novel_res = merge(x = WT_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_EvL_novel_res$Inc_gene_id  = ifelse(is.na(WT_EvL_novel_res$Inc_gene_id), WT_EvL_novel_res$ENSEMBL, WT_EvL_novel_res$Inc_gene_id)
#WT_EvL_novel_res = WT_EvL_novel_res[,-30]
#WT_EvL_novel_res$SYMBOL[is.na(WT_EvL_novel_res$Exc_gene_id)]  =  sub("\\-.*","",WT_EvL_novel_res$SYMBOL[is.na(WT_EvL_novel_res$Exc_gene_id)])
#WT_EvL_novel_res = merge(x = WT_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_EvL_novel_res$Exc_gene_id  = ifelse(is.na(WT_EvL_novel_res$Exc_gene_id), WT_EvL_novel_res$ENSEMBL, WT_EvL_novel_res$Exc_gene_id)
#WT_EvL_novel_res = WT_EvL_novel_res[,-30]

# dBGLF5 EvL novel
#dBGLF5_EvL_novel_res$SYMBOL = dBGLF5_EvL_novel_res$Inc_gene_id
#dBGLF5_EvL_novel_res$SYMBOL[dBGLF5_EvL_novel_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",dBGLF5_EvL_novel_res$EventName#[dBGLF5_EvL_novel_res$EventType == "IR"])))
#dBGLF5_EvL_novel_res = merge(x = dBGLF5_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_EvL_novel_res$ENSEMBL = ifelse(is.na(dBGLF5_EvL_novel_res$ENSEMBL), dBGLF5_EvL_novel_res$Inc_gene_id, dBGLF5_EvL_novel_res$ENSEMBL)
#dBGLF5_EvL_novel_res$Inc_gene_id = dBGLF5_EvL_novel_res$ENSEMBL
#dBGLF5_EvL_novel_res = dBGLF5_EvL_novel_res[,-30]
#dBGLF5_EvL_novel_res$SYMBOL = dBGLF5_EvL_novel_res$Exc_gene_id
#dBGLF5_EvL_novel_res$SYMBOL[dBGLF5_EvL_novel_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",dBGLF5_EvL_novel_res$EventName[dBGLF5_EvL_novel_res$EventType == "IR"])))
#dBGLF5_EvL_novel_res = merge(x = dBGLF5_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_EvL_novel_res$ENSEMBL = ifelse(is.na(dBGLF5_EvL_novel_res$ENSEMBL), dBGLF5_EvL_novel_res$Exc_gene_id, dBGLF5_EvL_novel_res$ENSEMBL)
#dBGLF5_EvL_novel_res$Exc_gene_id = dBGLF5_EvL_novel_res$ENSEMBL
#dBGLF5_EvL_novel_res = dBGLF5_EvL_novel_res[,-30]
#dBGLF5_EvL_novel_res$SYMBOL[is.na(dBGLF5_EvL_novel_res$Inc_gene_id)] = sub("\\..*","",sub("\\-novel.*","",sub("\\-2.*","",sub(".*;","",dBGLF5_EvL_novel_res$EventName[is.na(dBGLF5_EvL_novel_res$Inc_gene_id)]))))
#dBGLF5_EvL_novel_res = merge(x = dBGLF5_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_EvL_novel_res$Inc_gene_id  = ifelse(is.na(dBGLF5_EvL_novel_res$Inc_gene_id), dBGLF5_EvL_novel_res$ENSEMBL, dBGLF5_EvL_novel_res$Inc_gene_id)
#dBGLF5_EvL_novel_res = dBGLF5_EvL_novel_res[,-30]
#dBGLF5_EvL_novel_res$SYMBOL[is.na(dBGLF5_EvL_novel_res$Exc_gene_id)] = sub("\\..*","",sub("\\-novel.*","",sub("\\-2.*","",sub(".*;","",dBGLF5_EvL_novel_res$EventName[is.na(dBGLF5_EvL_novel_res$Exc_gene_id)]))))
#dBGLF5_EvL_novel_res = merge(x = dBGLF5_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_EvL_novel_res$Exc_gene_id  = ifelse(is.na(dBGLF5_EvL_novel_res$Exc_gene_id), dBGLF5_EvL_novel_res$ENSEMBL, dBGLF5_EvL_novel_res$Exc_gene_id)
#dBGLF5_EvL_novel_res = dBGLF5_EvL_novel_res[,-30]
#dBGLF5_EvL_novel_res$SYMBOL[is.na(dBGLF5_EvL_novel_res$Inc_gene_id)] = sub("_.*","",sub(".*;","",dBGLF5_EvL_novel_res$EventName[is.na(dBGLF5_EvL_novel_res$Inc_gene_id)]))
#dBGLF5_EvL_novel_res = merge(x = dBGLF5_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_EvL_novel_res$Inc_gene_id  = ifelse(is.na(dBGLF5_EvL_novel_res$Inc_gene_id), dBGLF5_EvL_novel_res$ENSEMBL, dBGLF5_EvL_novel_res$Inc_gene_id)
#dBGLF5_EvL_novel_res = dBGLF5_EvL_novel_res[,-30]
#dBGLF5_EvL_novel_res$SYMBOL[is.na(dBGLF5_EvL_novel_res$Exc_gene_id)] = sub("_.*","",sub(".*;","",dBGLF5_EvL_novel_res$EventName[is.na(dBGLF5_EvL_novel_res$Exc_gene_id)]))
#dBGLF5_EvL_novel_res = merge(x = dBGLF5_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_EvL_novel_res$Exc_gene_id  = ifelse(is.na(dBGLF5_EvL_novel_res$Exc_gene_id), dBGLF5_EvL_novel_res$ENSEMBL, dBGLF5_EvL_novel_res$Exc_gene_id)
#dBGLF5_EvL_novel_res = dBGLF5_EvL_novel_res[,-30]
#dBGLF5_EvL_novel_res$SYMBOL[is.na(dBGLF5_EvL_novel_res$Inc_gene_id)] = sub("_.*","",sub(".*:","",dBGLF5_EvL_novel_res$EventName[is.na(dBGLF5_EvL_novel_res$Inc_gene_id)]))
#dBGLF5_EvL_novel_res = merge(x = dBGLF5_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_EvL_novel_res$Inc_gene_id  = ifelse(is.na(dBGLF5_EvL_novel_res$Inc_gene_id), dBGLF5_EvL_novel_res$ENSEMBL, dBGLF5_EvL_novel_res$Inc_gene_id)
#dBGLF5_EvL_novel_res = dBGLF5_EvL_novel_res[,-30]
#dBGLF5_EvL_novel_res$SYMBOL[is.na(dBGLF5_EvL_novel_res$Exc_gene_id)] = sub("_.*","",sub(".*:","",dBGLF5_EvL_novel_res$EventName[is.na(dBGLF5_EvL_novel_res$Exc_gene_id)]))
#dBGLF5_EvL_novel_res = merge(x = dBGLF5_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_EvL_novel_res$Exc_gene_id  = ifelse(is.na(dBGLF5_EvL_novel_res$Exc_gene_id), dBGLF5_EvL_novel_res$ENSEMBL, dBGLF5_EvL_novel_res$Exc_gene_id)
#dBGLF5_EvL_novel_res = dBGLF5_EvL_novel_res[,-30]
#dBGLF5_EvL_novel_res$SYMBOL[is.na(dBGLF5_EvL_novel_res$Exc_gene_id)]  =  sub("\\-.*","",dBGLF5_EvL_novel_res$SYMBOL[is.na(dBGLF5_EvL_novel_res$Exc_gene_id)])
#dBGLF5_EvL_novel_res = merge(x = dBGLF5_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_EvL_novel_res$Inc_gene_id  = ifelse(is.na(dBGLF5_EvL_novel_res$Inc_gene_id), dBGLF5_EvL_novel_res$ENSEMBL, dBGLF5_EvL_novel_res$Inc_gene_id)
#dBGLF5_EvL_novel_res = dBGLF5_EvL_novel_res[,-30]
#dBGLF5_EvL_novel_res$SYMBOL[is.na(dBGLF5_EvL_novel_res$Exc_gene_id)]  =  sub("\\-.*","",dBGLF5_EvL_novel_res$SYMBOL[is.na(dBGLF5_EvL_novel_res$Exc_gene_id)])
#dBGLF5_EvL_novel_res = merge(x = dBGLF5_EvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_EvL_novel_res$Exc_gene_id  = ifelse(is.na(dBGLF5_EvL_novel_res$Exc_gene_id), dBGLF5_EvL_novel_res$ENSEMBL, dBGLF5_EvL_novel_res$Exc_gene_id)
#dBGLF5_EvL_novel_res = dBGLF5_EvL_novel_res[,-30]

# WT LvL
WT_LvL_res$SYMBOL = WT_LvL_res$Inc_gene_id
WT_LvL_res$SYMBOL[WT_LvL_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",WT_LvL_res$EventName[WT_LvL_res$EventType == "IR"])))
WT_LvL_res = merge(x = WT_LvL_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
WT_LvL_res$ENSEMBL = ifelse(is.na(WT_LvL_res$ENSEMBL), WT_LvL_res$Inc_gene_id, WT_LvL_res$ENSEMBL)
WT_LvL_res$Inc_gene_id = WT_LvL_res$ENSEMBL
WT_LvL_res = WT_LvL_res[,-30]
WT_LvL_res$SYMBOL = WT_LvL_res$Exc_gene_id
WT_LvL_res$SYMBOL[WT_LvL_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",WT_LvL_res$EventName[WT_LvL_res$EventType == "IR"])))
WT_LvL_res = merge(x = WT_LvL_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
WT_LvL_res$ENSEMBL = ifelse(is.na(WT_LvL_res$ENSEMBL), WT_LvL_res$Exc_gene_id, WT_LvL_res$ENSEMBL)
WT_LvL_res$Exc_gene_id = WT_LvL_res$ENSEMBL
WT_LvL_res = WT_LvL_res[,-30]

# dBGLF5 LvL
dBGLF5_LvL_res$SYMBOL = dBGLF5_LvL_res$Inc_gene_id
dBGLF5_LvL_res$SYMBOL[dBGLF5_LvL_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",dBGLF5_LvL_res$EventName[dBGLF5_LvL_res$EventType == "IR"])))
dBGLF5_LvL_res = merge(x = dBGLF5_LvL_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
dBGLF5_LvL_res$ENSEMBL = ifelse(is.na(dBGLF5_LvL_res$ENSEMBL), dBGLF5_LvL_res$Inc_gene_id, dBGLF5_LvL_res$ENSEMBL)
dBGLF5_LvL_res$Inc_gene_id = dBGLF5_LvL_res$ENSEMBL
dBGLF5_LvL_res = dBGLF5_LvL_res[,-30]
dBGLF5_LvL_res$SYMBOL = dBGLF5_LvL_res$Exc_gene_id
dBGLF5_LvL_res$SYMBOL[dBGLF5_LvL_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",dBGLF5_LvL_res$EventName[dBGLF5_LvL_res$EventType == "IR"])))
dBGLF5_LvL_res = merge(x = dBGLF5_LvL_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
dBGLF5_LvL_res$ENSEMBL = ifelse(is.na(dBGLF5_LvL_res$ENSEMBL), dBGLF5_LvL_res$Exc_gene_id, dBGLF5_LvL_res$ENSEMBL)
dBGLF5_LvL_res$Exc_gene_id = dBGLF5_LvL_res$ENSEMBL
dBGLF5_LvL_res = dBGLF5_LvL_res[,-30]

# WT LvL novel
#WT_LvL_novel_res$SYMBOL = WT_LvL_novel_res$Inc_gene_id
#WT_LvL_novel_res$SYMBOL[WT_LvL_novel_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",WT_LvL_novel_res$EventName[WT_LvL_novel_res$EventType == "IR"])))
#WT_LvL_novel_res = merge(x = WT_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvL_novel_res$ENSEMBL = ifelse(is.na(WT_LvL_novel_res$ENSEMBL), WT_LvL_novel_res$Inc_gene_id, WT_LvL_novel_res$ENSEMBL)
#WT_LvL_novel_res$Inc_gene_id = WT_LvL_novel_res$ENSEMBL
#WT_LvL_novel_res = WT_LvL_novel_res[,-30]
#WT_LvL_novel_res$SYMBOL = WT_LvL_novel_res$Exc_gene_id
#WT_LvL_novel_res$SYMBOL[WT_LvL_novel_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",WT_LvL_novel_res$EventName[WT_LvL_novel_res$EventType == "IR"])))
#WT_LvL_novel_res = merge(x = WT_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvL_novel_res$ENSEMBL = ifelse(is.na(WT_LvL_novel_res$ENSEMBL), WT_LvL_novel_res$Exc_gene_id, WT_LvL_novel_res$ENSEMBL)
#WT_LvL_novel_res$Exc_gene_id = WT_LvL_novel_res$ENSEMBL
#WT_LvL_novel_res = WT_LvL_novel_res[,-30]
#WT_LvL_novel_res$SYMBOL[is.na(WT_LvL_novel_res$Inc_gene_id)] = sub("\\..*","",sub("\\-novel.*","",sub("\\-2.*","",sub(".*;","",WT_LvL_novel_res$EventName[is.na(WT_LvL_novel_res$Inc_gene_id)]))))
#WT_LvL_novel_res = merge(x = WT_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvL_novel_res$Inc_gene_id  = ifelse(is.na(WT_LvL_novel_res$Inc_gene_id), WT_LvL_novel_res$ENSEMBL, WT_LvL_novel_res$Inc_gene_id)
#WT_LvL_novel_res = WT_LvL_novel_res[,-30]
#WT_LvL_novel_res$SYMBOL[is.na(WT_LvL_novel_res$Exc_gene_id)] = sub("\\..*","",sub("\\-novel.*","",sub("\\-2.*","",sub(".*;","",WT_LvL_novel_res$EventName[is.na(WT_LvL_novel_res$Exc_gene_id)]))))
#WT_LvL_novel_res = merge(x = WT_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvL_novel_res$Exc_gene_id  = ifelse(is.na(WT_LvL_novel_res$Exc_gene_id), WT_LvL_novel_res$ENSEMBL, WT_LvL_novel_res$Exc_gene_id)
#WT_LvL_novel_res = WT_LvL_novel_res[,-30]
#WT_LvL_novel_res$SYMBOL[is.na(WT_LvL_novel_res$Inc_gene_id)] = sub("_.*","",sub(".*;","",WT_LvL_novel_res$EventName[is.na(WT_LvL_novel_res$Inc_gene_id)]))
#WT_LvL_novel_res = merge(x = WT_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvL_novel_res$Inc_gene_id  = ifelse(is.na(WT_LvL_novel_res$Inc_gene_id), WT_LvL_novel_res$ENSEMBL, WT_LvL_novel_res$Inc_gene_id)
#WT_LvL_novel_res = WT_LvL_novel_res[,-30]
#WT_LvL_novel_res$SYMBOL[is.na(WT_LvL_novel_res$Exc_gene_id)] = sub("_.*","",sub(".*;","",WT_LvL_novel_res$EventName[is.na(WT_LvL_novel_res$Exc_gene_id)]))
#WT_LvL_novel_res = merge(x = WT_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvL_novel_res$Exc_gene_id  = ifelse(is.na(WT_LvL_novel_res$Exc_gene_id), WT_LvL_novel_res$ENSEMBL, WT_LvL_novel_res$Exc_gene_id)
#WT_LvL_novel_res = WT_LvL_novel_res[,-30]
#WT_LvL_novel_res$SYMBOL[is.na(WT_LvL_novel_res$Inc_gene_id)] = sub("_.*","",sub(".*:","",WT_LvL_novel_res$EventName[is.na(WT_LvL_novel_res$Inc_gene_id)]))
#WT_LvL_novel_res = merge(x = WT_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvL_novel_res$Inc_gene_id  = ifelse(is.na(WT_LvL_novel_res$Inc_gene_id), WT_LvL_novel_res$ENSEMBL, WT_LvL_novel_res$Inc_gene_id)
#WT_LvL_novel_res = WT_LvL_novel_res[,-30]
#WT_LvL_novel_res$SYMBOL[is.na(WT_LvL_novel_res$Exc_gene_id)] = sub("_.*","",sub(".*:","",WT_LvL_novel_res$EventName[is.na(WT_LvL_novel_res$Exc_gene_id)]))
#WT_LvL_novel_res = merge(x = WT_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvL_novel_res$Exc_gene_id  = ifelse(is.na(WT_LvL_novel_res$Exc_gene_id), WT_LvL_novel_res$ENSEMBL, WT_LvL_novel_res$Exc_gene_id)
#WT_LvL_novel_res = WT_LvL_novel_res[,-30]
#WT_LvL_novel_res$SYMBOL[is.na(WT_LvL_novel_res$Exc_gene_id)]  =  sub("\\-.*","",WT_LvL_novel_res$SYMBOL[is.na(WT_LvL_novel_res$Exc_gene_id)])
#WT_LvL_novel_res = merge(x = WT_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvL_novel_res$Inc_gene_id  = ifelse(is.na(WT_LvL_novel_res$Inc_gene_id), WT_LvL_novel_res$ENSEMBL, WT_LvL_novel_res$Inc_gene_id)
#WT_LvL_novel_res = WT_LvL_novel_res[,-30]
#WT_LvL_novel_res$SYMBOL[is.na(WT_LvL_novel_res$Exc_gene_id)]  =  sub("\\-.*","",WT_LvL_novel_res$SYMBOL[is.na(WT_LvL_novel_res$Exc_gene_id)])
#WT_LvL_novel_res = merge(x = WT_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvL_novel_res$Exc_gene_id  = ifelse(is.na(WT_LvL_novel_res$Exc_gene_id), WT_LvL_novel_res$ENSEMBL, WT_LvL_novel_res$Exc_gene_id)
#WT_LvL_novel_res = WT_LvL_novel_res[,-30]

# dBGLF5 LvL novel
#dBGLF5_LvL_novel_res$SYMBOL = dBGLF5_LvL_novel_res$Inc_gene_id
#dBGLF5_LvL_novel_res$SYMBOL[dBGLF5_LvL_novel_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",dBGLF5_LvL_novel_res$EventName[dBGLF5_LvL_novel_res$EventType == "IR"])))
#dBGLF5_LvL_novel_res = merge(x = dBGLF5_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_LvL_novel_res$ENSEMBL = ifelse(is.na(dBGLF5_LvL_novel_res$ENSEMBL), dBGLF5_LvL_novel_res$Inc_gene_id, dBGLF5_LvL_novel_res$ENSEMBL)
#dBGLF5_LvL_novel_res$Inc_gene_id = dBGLF5_LvL_novel_res$ENSEMBL
#dBGLF5_LvL_novel_res = dBGLF5_LvL_novel_res[,-30]
#dBGLF5_LvL_novel_res$SYMBOL = dBGLF5_LvL_novel_res$Exc_gene_id
#dBGLF5_LvL_novel_res$SYMBOL[dBGLF5_LvL_novel_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",dBGLF5_LvL_novel_res$EventName[dBGLF5_LvL_novel_res$EventType == "IR"])))
#dBGLF5_LvL_novel_res = merge(x = dBGLF5_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_LvL_novel_res$ENSEMBL = ifelse(is.na(dBGLF5_LvL_novel_res$ENSEMBL), dBGLF5_LvL_novel_res$Exc_gene_id, dBGLF5_LvL_novel_res$ENSEMBL)
#dBGLF5_LvL_novel_res$Exc_gene_id = dBGLF5_LvL_novel_res$ENSEMBL
#dBGLF5_LvL_novel_res = dBGLF5_LvL_novel_res[,-30]
#dBGLF5_LvL_novel_res$SYMBOL[is.na(dBGLF5_LvL_novel_res$Inc_gene_id)] = sub("\\..*","",sub("\\-novel.*","",sub("\\-2.*","",sub(".*;","",dBGLF5_LvL_novel_res$EventName[is.na(dBGLF5_LvL_novel_res$Inc_gene_id)]))))
#dBGLF5_LvL_novel_res = merge(x = dBGLF5_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_LvL_novel_res$Inc_gene_id  = ifelse(is.na(dBGLF5_LvL_novel_res$Inc_gene_id), dBGLF5_LvL_novel_res$ENSEMBL, dBGLF5_LvL_novel_res$Inc_gene_id)
#dBGLF5_LvL_novel_res = dBGLF5_LvL_novel_res[,-30]
#dBGLF5_LvL_novel_res$SYMBOL[is.na(dBGLF5_LvL_novel_res$Exc_gene_id)] = sub("\\..*","",sub("\\-novel.*","",sub("\\-2.*","",sub(".*;","",dBGLF5_LvL_novel_res$EventName[is.na(dBGLF5_LvL_novel_res$Exc_gene_id)]))))
#dBGLF5_LvL_novel_res = merge(x = dBGLF5_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_LvL_novel_res$Exc_gene_id  = ifelse(is.na(dBGLF5_LvL_novel_res$Exc_gene_id), dBGLF5_LvL_novel_res$ENSEMBL, dBGLF5_LvL_novel_res$Exc_gene_id)
#dBGLF5_LvL_novel_res = dBGLF5_LvL_novel_res[,-30]
#dBGLF5_LvL_novel_res$SYMBOL[is.na(dBGLF5_LvL_novel_res$Inc_gene_id)] = sub("_.*","",sub(".*;","",dBGLF5_LvL_novel_res$EventName[is.na(dBGLF5_LvL_novel_res$Inc_gene_id)]))
#dBGLF5_LvL_novel_res = merge(x = dBGLF5_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_LvL_novel_res$Inc_gene_id  = ifelse(is.na(dBGLF5_LvL_novel_res$Inc_gene_id), dBGLF5_LvL_novel_res$ENSEMBL, dBGLF5_LvL_novel_res$Inc_gene_id)
#dBGLF5_LvL_novel_res = dBGLF5_LvL_novel_res[,-30]
#dBGLF5_LvL_novel_res$SYMBOL[is.na(dBGLF5_LvL_novel_res$Exc_gene_id)] = sub("_.*","",sub(".*;","",dBGLF5_LvL_novel_res$EventName[is.na(dBGLF5_LvL_novel_res$Exc_gene_id)]))
#dBGLF5_LvL_novel_res = merge(x = dBGLF5_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_LvL_novel_res$Exc_gene_id  = ifelse(is.na(dBGLF5_LvL_novel_res$Exc_gene_id), dBGLF5_LvL_novel_res$ENSEMBL, dBGLF5_LvL_novel_res$Exc_gene_id)
#dBGLF5_LvL_novel_res = dBGLF5_LvL_novel_res[,-30]
#dBGLF5_LvL_novel_res$SYMBOL[is.na(dBGLF5_LvL_novel_res$Inc_gene_id)] = sub("_.*","",sub(".*:","",dBGLF5_LvL_novel_res$EventName[is.na(dBGLF5_LvL_novel_res$Inc_gene_id)]))
#dBGLF5_LvL_novel_res = merge(x = dBGLF5_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_LvL_novel_res$Inc_gene_id  = ifelse(is.na(dBGLF5_LvL_novel_res$Inc_gene_id), dBGLF5_LvL_novel_res$ENSEMBL, dBGLF5_LvL_novel_res$Inc_gene_id)
#dBGLF5_LvL_novel_res = dBGLF5_LvL_novel_res[,-30]
#dBGLF5_LvL_novel_res$SYMBOL[is.na(dBGLF5_LvL_novel_res$Exc_gene_id)] = sub("_.*","",sub(".*:","",dBGLF5_LvL_novel_res$EventName[is.na(dBGLF5_LvL_novel_res$Exc_gene_id)]))
#dBGLF5_LvL_novel_res = merge(x = dBGLF5_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_LvL_novel_res$Exc_gene_id  = ifelse(is.na(dBGLF5_LvL_novel_res$Exc_gene_id), dBGLF5_LvL_novel_res$ENSEMBL, dBGLF5_LvL_novel_res$Exc_gene_id)
#dBGLF5_LvL_novel_res = dBGLF5_LvL_novel_res[,-30]
#dBGLF5_LvL_novel_res$SYMBOL[is.na(dBGLF5_LvL_novel_res$Exc_gene_id)]  =  sub("\\-.*","",dBGLF5_LvL_novel_res$SYMBOL[is.na(dBGLF5_LvL_novel_res$Exc_gene_id)])
#dBGLF5_LvL_novel_res = merge(x = dBGLF5_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_LvL_novel_res$Inc_gene_id  = ifelse(is.na(dBGLF5_LvL_novel_res$Inc_gene_id), dBGLF5_LvL_novel_res$ENSEMBL, dBGLF5_LvL_novel_res$Inc_gene_id)
#dBGLF5_LvL_novel_res = dBGLF5_LvL_novel_res[,-30]
#dBGLF5_LvL_novel_res$SYMBOL[is.na(dBGLF5_LvL_novel_res$Exc_gene_id)]  =  sub("\\-.*","",dBGLF5_LvL_novel_res$SYMBOL[is.na(dBGLF5_LvL_novel_res$Exc_gene_id)])
#dBGLF5_LvL_novel_res = merge(x = dBGLF5_LvL_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#dBGLF5_LvL_novel_res$Exc_gene_id  = ifelse(is.na(dBGLF5_LvL_novel_res$Exc_gene_id), dBGLF5_LvL_novel_res$ENSEMBL, dBGLF5_LvL_novel_res$Exc_gene_id)
#dBGLF5_LvL_novel_res = dBGLF5_LvL_novel_res[,-30]

# WT LvE
WT_LvE_res$SYMBOL = WT_LvE_res$Inc_gene_id
WT_LvE_res$SYMBOL[WT_LvE_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",WT_LvE_res$EventName[WT_LvE_res$EventType == "IR"])))
WT_LvE_res = merge(x = WT_LvE_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
WT_LvE_res$ENSEMBL = ifelse(is.na(WT_LvE_res$ENSEMBL), WT_LvE_res$Inc_gene_id, WT_LvE_res$ENSEMBL)
WT_LvE_res$Inc_gene_id = WT_LvE_res$ENSEMBL
WT_LvE_res = WT_LvE_res[,-30]
WT_LvE_res$SYMBOL = WT_LvE_res$Exc_gene_id
WT_LvE_res$SYMBOL[WT_LvE_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",WT_LvE_res$EventName[WT_LvE_res$EventType == "IR"])))
WT_LvE_res = merge(x = WT_LvE_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
WT_LvE_res$ENSEMBL = ifelse(is.na(WT_LvE_res$ENSEMBL), WT_LvE_res$Exc_gene_id, WT_LvE_res$ENSEMBL)
WT_LvE_res$Exc_gene_id = WT_LvE_res$ENSEMBL
WT_LvE_res = WT_LvE_res[,-30]

# WT LvE novel
#WT_LvE_novel_res$SYMBOL = WT_LvE_novel_res$Inc_gene_id
#WT_LvE_novel_res$SYMBOL[WT_LvE_novel_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",WT_LvE_novel_res$EventName[WT_LvE_novel_res$EventType == "IR"])))
#WT_LvE_novel_res = merge(x = WT_LvE_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvE_novel_res$ENSEMBL = ifelse(is.na(WT_LvE_novel_res$ENSEMBL), WT_LvE_novel_res$Inc_gene_id, WT_LvE_novel_res$ENSEMBL)
#WT_LvE_novel_res$Inc_gene_id = WT_LvE_novel_res$ENSEMBL
#WT_LvE_novel_res = WT_LvE_novel_res[,-30]
#WT_LvE_novel_res$SYMBOL = WT_LvE_novel_res$Exc_gene_id
#WT_LvE_novel_res$SYMBOL[WT_LvE_novel_res$EventType == "IR"] = sub("_.*","",sub("\\/.*","",sub("RI\\:","",WT_LvE_novel_res$EventName[WT_LvE_novel_res$EventType == "IR"])))
#WT_LvE_novel_res = merge(x = WT_LvE_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvE_novel_res$ENSEMBL = ifelse(is.na(WT_LvE_novel_res$ENSEMBL), WT_LvE_novel_res$Exc_gene_id, WT_LvE_novel_res$ENSEMBL)
#WT_LvE_novel_res$Exc_gene_id = WT_LvE_novel_res$ENSEMBL
#WT_LvE_novel_res = WT_LvE_novel_res[,-30]
#WT_LvE_novel_res$SYMBOL[is.na(WT_LvE_novel_res$Inc_gene_id)] = sub("\\..*","",sub("\\-novel.*","",sub("\\-2.*","",sub(".*;","",WT_LvE_novel_res$EventName[is.na(WT_LvE_novel_res$Inc_gene_id)]))))
#WT_LvE_novel_res = merge(x = WT_LvE_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvE_novel_res$Inc_gene_id  = ifelse(is.na(WT_LvE_novel_res$Inc_gene_id), WT_LvE_novel_res$ENSEMBL, WT_LvE_novel_res$Inc_gene_id)
#WT_LvE_novel_res = WT_LvE_novel_res[,-30]
#WT_LvE_novel_res$SYMBOL[is.na(WT_LvE_novel_res$Exc_gene_id)] = sub("\\..*","",sub("\\-novel.*","",sub("\\-2.*","",sub(".*;","",WT_LvE_novel_res$EventName[is.na(WT_LvE_novel_res$Exc_gene_id)]))))
#WT_LvE_novel_res = merge(x = WT_LvE_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvE_novel_res$Exc_gene_id  = ifelse(is.na(WT_LvE_novel_res$Exc_gene_id), WT_LvE_novel_res$ENSEMBL, WT_LvE_novel_res$Exc_gene_id)
#WT_LvE_novel_res = WT_LvE_novel_res[,-30]
#WT_LvE_novel_res$SYMBOL[is.na(WT_LvE_novel_res$Inc_gene_id)] = sub("_.*","",sub(".*;","",WT_LvE_novel_res$EventName[is.na(WT_LvE_novel_res$Inc_gene_id)]))
#WT_LvE_novel_res = merge(x = WT_LvE_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvE_novel_res$Inc_gene_id  = ifelse(is.na(WT_LvE_novel_res$Inc_gene_id), WT_LvE_novel_res$ENSEMBL, WT_LvE_novel_res$Inc_gene_id)
#WT_LvE_novel_res = WT_LvE_novel_res[,-30]
#WT_LvE_novel_res$SYMBOL[is.na(WT_LvE_novel_res$Exc_gene_id)] = sub("_.*","",sub(".*;","",WT_LvE_novel_res$EventName[is.na(WT_LvE_novel_res$Exc_gene_id)]))
#WT_LvE_novel_res = merge(x = WT_LvE_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvE_novel_res$Exc_gene_id  = ifelse(is.na(WT_LvE_novel_res$Exc_gene_id), WT_LvE_novel_res$ENSEMBL, WT_LvE_novel_res$Exc_gene_id)
#WT_LvE_novel_res = WT_LvE_novel_res[,-30]
#WT_LvE_novel_res$SYMBOL[is.na(WT_LvE_novel_res$Inc_gene_id)] = sub("_.*","",sub(".*:","",WT_LvE_novel_res$EventName[is.na(WT_LvE_novel_res$Inc_gene_id)]))
#WT_LvE_novel_res = merge(x = WT_LvE_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvE_novel_res$Inc_gene_id  = ifelse(is.na(WT_LvE_novel_res$Inc_gene_id), WT_LvE_novel_res$ENSEMBL, WT_LvE_novel_res$Inc_gene_id)
#WT_LvE_novel_res = WT_LvE_novel_res[,-30]
#WT_LvE_novel_res$SYMBOL[is.na(WT_LvE_novel_res$Exc_gene_id)] = sub("_.*","",sub(".*:","",WT_LvE_novel_res$EventName[is.na(WT_LvE_novel_res$Exc_gene_id)]))
#WT_LvE_novel_res = merge(x = WT_LvE_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvE_novel_res$Exc_gene_id  = ifelse(is.na(WT_LvE_novel_res$Exc_gene_id), WT_LvE_novel_res$ENSEMBL, WT_LvE_novel_res$Exc_gene_id)
#WT_LvE_novel_res = WT_LvE_novel_res[,-30]
#WT_LvE_novel_res$SYMBOL[is.na(WT_LvE_novel_res$Exc_gene_id)]  =  sub("\\-.*","",WT_LvE_novel_res$SYMBOL[is.na(WT_LvE_novel_res$Exc_gene_id)])
#WT_LvE_novel_res = merge(x = WT_LvE_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvE_novel_res$Inc_gene_id  = ifelse(is.na(WT_LvE_novel_res$Inc_gene_id), WT_LvE_novel_res$ENSEMBL, WT_LvE_novel_res$Inc_gene_id)
#WT_LvE_novel_res = WT_LvE_novel_res[,-30]
#WT_LvE_novel_res$SYMBOL[is.na(WT_LvE_novel_res$Exc_gene_id)]  =  sub("\\-.*","",WT_LvE_novel_res$SYMBOL[is.na(WT_LvE_novel_res$Exc_gene_id)])
#WT_LvE_novel_res = merge(x = WT_LvE_novel_res, y = gene_info[,1:2], by = "SYMBOL", all.x = TRUE)
#WT_LvE_novel_res$Exc_gene_id  = ifelse(is.na(WT_LvE_novel_res$Exc_gene_id), WT_LvE_novel_res$ENSEMBL, WT_LvE_novel_res$Exc_gene_id)
#WT_LvE_novel_res = WT_LvE_novel_res[,-30]

# Export DESeq2 ASE files
write.csv(WT_EvL_res, file.path(DESeq2_dir,"DE.WT_Early.vs.WT_Latent.csv"))
write.csv(dBGLF5_EvL_res, file.path(DESeq2_dir,"DE.dBGLF5_Early.vs.dBGLF5_Latent.csv"))
#write.csv(WT_EvL_novel_res, file.path(DESeq2_dir,"DE.novel.WT_Early.vs.WT_Latent.csv"))
#write.csv(dBGLF5_EvL_novel_res, file.path(DESeq2_dir,"DE.novel.dBGLF5_Early.vs.dBGLF5_Latent.csv"))
write.csv(WT_LvL_res, file.path(DESeq2_dir,"DE.WT_Late.vs.WT_Latent.csv"))
write.csv(dBGLF5_LvL_res, file.path(DESeq2_dir,"DE.dBGLF5_Late.vs.dBGLF5_Latent.csv"))
#write.csv(WT_LvL_novel_res, file.path(DESeq2_dir,"DE.novel.WT_Late.vs.WT_Latent.csv"))
#write.csv(dBGLF5_LvL_novel_res, file.path(DESeq2_dir,"DE.novel.dBGLF5_Late.vs.dBGLF5_Latent.csv"))
write.csv(WT_LvE_res, file.path(DESeq2_dir,"DE.WT_Late.vs.WT_Early.csv"))
#write.csv(WT_LvE_novel_res, file.path(DESeq2_dir,"DE.novel.WT_Late.vs.WT_Early.csv"))

## print session info ##
print("Session Info below: ")
sessionInfo()
