rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("tidyverse")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("tximport")
#BiocManager::install("AcidGenomes")

## Setup the environment

# Import packages
library(rstudioapi)
library(tximport) # imports data from RSEM
library(tidyverse) # formatting data
library(AcidGenomes) # strip gene versions

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Set output directory
out_dir=file.path(work_dir,"Data/MPC")

# Set directory containing RSEM counts files
counts_dir=file.path(work_dir,"Data/RSEM")

# Set path to expressed genes file
deg_dir=file.path(work_dir,"Data/DESeq2")
if (!file.exists(list.files(deg_dir,pattern="^ERCCnorm_DE*",full.name=TRUE)[1])) 
  {stop(" No DESeq2 files exist.\n\tRun DESeq2.R script to generate.", call. = FALSE)} # check file exists

# Set directory to sampleinfo matrix file
sampleinfo_dir=file.path(work_dir,"Data/refs/sampleinfo.txt")

# List of DESeq2 outputs to filter final data by
DESeq2_filter=list(c("WT_Early","WT_Latent"),c("WT_Late","WT_Latent"),c("dBGLF5_Early","dBGLF5_Latent"),c("dBGLF5_Late","dBGLF5_Latent"))

# Set ERCC spike-in dilution factor and number of cells
cells=40000
dilution=(195/1.04)/1.5 # 1.04 uL of ERCC RNA spiked into 195 mL of TRIzol LS. 40,000 cells were sorted into 1.5 mL aliquots.

## Molecule per cell (MPC) quantification

# Get list of RSEM gene count files data within counts directory
files <- list.files(file.path(counts_dir),pattern = ".genes.results", full.names = TRUE)

# Reorder list according to sampleinfo
sampleinfo <- read.delim(file.path(sampleinfo_dir))
temp_list <- as.list(sampleinfo[[1]])
temp_files <- c()
for (i in as.character(temp_list)) {
  temp_files <- append(temp_files, c(print(files[grep(i, files)])))
}
files <- temp_files
rm(temp_list, temp_files)

# Name Samples
names(files) <- paste0(sampleinfo[[2]], "_", sampleinfo[[3]])

# Import files with tximport as FPKM for filtering
txi.rsem <- tximport(files, type = "none", txIn = FALSE, txOut = FALSE, geneIdCol = "gene_id",
                     abundanceCol = "FPKM", countsCol = "expected_count", lengthCol = "effective_length",
                     importer = function(x) read_tsv(x))
FPKM = txi.rsem$abundance
rownames(FPKM) = stripGeneVersions(rownames(FPKM))

#  Make matrix of ERCC genes
ERCC_FPKM <- FPKM[grep("ERCC-", rownames(FPKM)),]

# Filter by genes from DESeq2 analysis
for (i in 1:length(DESeq2_filter)) {
  file = paste0("ERCCnorm_DE.",paste(DESeq2_filter[[i]], collapse = ".vs."),".csv")
  dge = read.csv(file.path(deg_dir,file))[,1:4] # import
  dge$ENSEMBL = ifelse(is.na(dge$ENSEMBL), dge$SYMBOL, dge$ENSEMBL) # if ENSEMBL ID does not exist, use SYMBOL
  if (i == 1) {dge_out = dge} else {
    dge_out = rbind(dge_out, dge[!dge$ENSEMBL %in% dge_out$ENSEMBL,]) # this keeps all genes present across comparisons
  }
}
dge_ercc = dge_out[grepl("ERCC-",dge_out$ENSEMBL),] # get ERCC genes present in DESeq2 analysis across comparisons
dge_out = dge_out[!grepl("ERCC-",dge_out$ENSEMBL),] # remove ERCC genes

# Filter by genes present in DESeq2 analysis across comparisons
FPKM = FPKM[rownames(FPKM) %in% dge_out$ENSEMBL,]
ERCC_FPKM = ERCC_FPKM[rownames(ERCC_FPKM) %in% dge_ercc$SYMBOL,]

# Replace ERCC genes containing FPKM 0 with NA. Note: there really shouldn't be any after DESeq2 filtering
ERCC_FPKM[ERCC_FPKM == 0] <- NA

# Calculate the log-base-2 of the ERCC FPKM values from RSEM
ERCC_FPKM <- log2(ERCC_FPKM)

## Generate a standard curve that relates ERCC Spike-In FPKM to the absolute amount of each ERCC transcript added

# Prepare ERCC concentration file
ERCC_Mix <- as.data.frame(1:92)
colnames(ERCC_Mix) <- "ID"
ERCC_Mix[,1] <- c("ERCC-00130", "ERCC-00004", "ERCC-00136", "ERCC-00108", "ERCC-00116", "ERCC-00092", "ERCC-00095", "ERCC-00131", "ERCC-00062", 
                  "ERCC-00019", "ERCC-00144", "ERCC-00170", "ERCC-00154", "ERCC-00085", "ERCC-00028", "ERCC-00033", "ERCC-00134", "ERCC-00147", "ERCC-00097", 
                  "ERCC-00156", "ERCC-00123", "ERCC-00017", "ERCC-00083", "ERCC-00096", "ERCC-00171", "ERCC-00009", "ERCC-00042", "ERCC-00060", "ERCC-00035", 
                  "ERCC-00025", "ERCC-00051", "ERCC-00053", "ERCC-00148", "ERCC-00126", "ERCC-00034", "ERCC-00150", "ERCC-00067", "ERCC-00031", "ERCC-00109", 
                  "ERCC-00073", "ERCC-00158", "ERCC-00104", "ERCC-00142", "ERCC-00138", "ERCC-00117", "ERCC-00075", "ERCC-00074", "ERCC-00113", "ERCC-00145", 
                  "ERCC-00111", "ERCC-00076", "ERCC-00044", "ERCC-00162", "ERCC-00071", "ERCC-00084", "ERCC-00099", "ERCC-00054", "ERCC-00157", "ERCC-00143", 
                  "ERCC-00039", "ERCC-00058", "ERCC-00120", "ERCC-00040", "ERCC-00164", "ERCC-00024", "ERCC-00016", "ERCC-00012", "ERCC-00098", "ERCC-00057", 
                  "ERCC-00002", "ERCC-00046", "ERCC-00003", "ERCC-00043", "ERCC-00022", "ERCC-00112", "ERCC-00165", "ERCC-00079", "ERCC-00078", "ERCC-00163", 
                  "ERCC-00059", "ERCC-00160", "ERCC-00014", "ERCC-00077", "ERCC-00069", "ERCC-00137", "ERCC-00013", "ERCC-00168", "ERCC-00041", "ERCC-00081", 
                  "ERCC-00086", "ERCC-00061", "ERCC-00048")
ERCC_Mix$'attomoles/ul' <- c(30000, 7500, 1875, 937.5, 468.75, 234.375, 117.1875, 117.1875, 58.59375, 29.296875, 29.296875, 
                             14.6484375, 7.32421875, 7.32421875, 3.66210938, 1.83105469, 1.83105469, 0.91552734, 0.45776367, 
                             0.45776367, 0.22888184, 0.11444092, 0.02861023, 15000, 3750, 937.5, 468.75, 234.375, 117.1875, 
                             58.59375, 58.59375, 29.296875, 14.6484375, 14.6484375, 7.32421875, 3.66210938, 3.66210938, 
                             1.83105469, 0.91552734, 0.91552734, 0.45776367, 0.22888184, 0.22888184, 0.11444092, 0.05722046, 
                             0.01430512, 15000, 3750, 937.5, 468.75, 234.375, 117.1875, 58.59375, 58.59375, 29.296875, 
                             14.6484375, 14.6484375, 7.32421875, 3.66210938, 3.66210938, 1.83105469, 0.91552734, 0.91552734, 
                             0.45776367, 0.22888184, 0.22888184, 0.11444092, 0.05722046, 0.01430512, 15000, 3750, 937.5, 468.75, 
                             234.375, 117.1875, 58.59375, 58.59375, 29.296875, 14.6484375, 14.6484375, 7.32421875, 3.66210938, 
                             3.66210938, 1.83105469, 0.91552734, 0.91552734, 0.45776367, 0.22888184, 0.22888184, 0.11444092, 
                             0.05722046, 0.01430512)

# Calculate ERCC attomoles added
ERCC_Mix$attomoles_added <- (ERCC_Mix$`attomoles/ul`/as.numeric(dilution))

# Calculate ERCC Moles
ERCC_Mix$moles <- ERCC_Mix$attomoles_added/1e18

# Calculate ERCC molecules
ERCC_Mix$molecules <- ERCC_Mix$moles*6.0221408e23

# Calculate ERCC molecules per cell
ERCC_Mix$MPC <- ERCC_Mix$molecules/as.numeric(cells)

# Calculate the log-base-2 of the MPC values
ERCC_Mix$ERCC_LogMPC <- log2(ERCC_Mix$MPC)

# Merge ERCC log-base-2 MPC values and RSEM FPKM values for ERCC transcripts
StdCurve <- merge(x = ERCC_Mix[,c(1,7)], y = ERCC_FPKM, by.x = "ID", by.y = 0, all.y = TRUE)
rownames(StdCurve) <- StdCurve[,1]
StdCurve <- StdCurve[,-1]

## Calculate molecules per cell for genes of interest
# Calculate log2(MPC) for each gene from RSEM using x = (y-b)/m
# Calculate MPC using exponential 2x of log2(MPC)

MPC <- as.data.frame((FPKM))
MPC[MPC == 0] <- NA # replace FPKM of genes with NA if FPKM = 0
MPC <- log2(MPC)

for(i in 2:length(StdCurve)) {
  temp <- lm(StdCurve[[i]]~ERCC_LogMPC, data = StdCurve)
  new_var <- (StdCurve[[i]]-temp$coefficients[[1]])/(temp$coefficients[[2]])
  a <- (i-1)
  MPC[,(a)] <- (2^((MPC[[a]]-temp$coefficients[[1]])/(temp$coefficients[[2]])))
}
MPC[is.na(MPC)] <- 0

# Round data
FPKM <- round(FPKM,2)
MPC <- round(MPC,2)

# annotate data using DESeq2 output
FPKM = merge(x = dge_out, y = FPKM, by.x = "ENSEMBL", by.y = 0, all.y = TRUE)
MPC = merge(x = dge_out, y = MPC, by.x = "ENSEMBL", by.y = 0, all.y = TRUE)

# Convert non-Ensembl genes to NA in ENSEMBL column
FPKM$ENSEMBL = ifelse(!grepl("ENSG",FPKM$ENSEMBL), NA, FPKM$ENSEMBL)
MPC$ENSEMBL = ifelse(!grepl("ENSG",MPC$ENSEMBL), NA, MPC$ENSEMBL)

# Create output directory
dir.create(out_dir, showWarnings=FALSE)

## Export MPC and FPKM data
#write.csv(FPKM,file=file.path(out_dir,"FPKM.csv"),row.names = FALSE)
write.csv(MPC,file=file.path(out_dir,"MPC.csv"),row.names = FALSE)

## print session info ##
print("Session Info below: ")
sessionInfo()
