rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("reshape2")
#install.packages("ggpubr")
#install.packages("svglite")

## Setup the environment

# Import packages
library(rstudioapi)
library(reshape2)
library(ggpubr)
library(svglite)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Output directory
out_dir=file.path(work_dir, "Figures/Figure_6")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Load molecules per cell (MPC) data
MPC_path = file.path(work_dir, "Data/MPC/MPC.csv")
if (!file.exists(MPC_path)) {stop(" MPC.csv file does not exist.\n\tRun MPC.R script to generate.", call. = FALSE)} # check file exists
MPC=read.csv(MPC_path,check.names=FALSE)

# Convert to matrix style
MPC$ENSEMBL = ifelse(is.na(MPC$ENSEMBL), MPC$SYMBOL, MPC$ENSEMBL) # if ENSEMBL ID does not exist, use SYMBOL
rownames(MPC) = MPC$ENSEMBL
MPC = MPC[,-1:-4]

# Keep only EBV genes
MPC = MPC[!grepl("ENSG|ERCC",rownames(MPC)),]

# Rename some EBV genes
rownames(MPC) = sub("BILF2-tPT2A-mScarlet-I", "BILF2-mSI", rownames(MPC))
rownames(MPC) = sub("mGreenLantern", "mGL", rownames(MPC))
rownames(MPC) = sub("BSLF2/BMLF1", "SM", rownames(MPC))
rownames(MPC) = sub("BBLF2/BBLF3", "BBLF2/3", rownames(MPC))
rownames(MPC) = sub("BGRF1/BDRF1", "BG/DRF1", rownames(MPC))

# Remove genes, including BGLF5
MPC = MPC[!grepl("EBER1|EBER2|BWRF1|BGLF5", rownames(MPC)),]

# get lists of EBV genes by their expression kinetics
Early = c("BZLF1","BRLF1","mGL", "BHLF1","BFLF2","BORF2","BaRF1","BMRF1","BMRF2","SM","BSLF1","BSRF1","BLLF3","BLLF2","BRRF1","BKRF3",
          "BBLF4","BBLF2/3","BGLF5","BGLF4","BGLF3.5","BGLF3","BDLF3.5","BVRF1","BVLF1","BALF2","BALF1","BARF1",
          "BCRF1","BFRF2","BG/DRF1","BTRF1","LF3","LF2","LF1","BILF1","BALF5","BNLF2a/b")
Leaky = c("BFLF1","BORF1","BLRF1","BLRF2","BRRF2","BKRF2","BBRF1","BBRF2","BBRF3","BBLF1","BDLF3","BcRF1","BXLF2","BdRF1",
          "BALF4","BALF3","BFRF0.5","BKRF4")
Late = c("BNRF1","BFRF1","BFRF3","BPLF1","BOLF1","BLLF1","BZLF2","BGLF2","BGLF1","BDLF2","BDLF1","BcLF1","BXRF1",
         "BVRF2","BILF2-mSI","BDLF4","LMP1_lytic")
Latent = c("LMP2", "LMP2A", "LMP2B","EBER1","EBER2","EBNALP","EBNA2","EBNA3A","EBNA3B","EBNA3C","EBNA1","BART","LMP1","BHRF1")
MPC$Kinetics =  ifelse(rownames(MPC) %in% Early, "Early",
                       ifelse(rownames(MPC) %in% Leaky, "Leaky",
                              ifelse(rownames(MPC) %in% Late, "Late",
                                     ifelse(rownames(MPC) %in% Latent, "Latent", NA))))

rm(list=setdiff(ls(), c("MPC", "out_dir")))

################ FIGURE 6C ################

# Sum by the kinetics
for (i in 1:(length(MPC)-1)) {
  tmp_df = MPC[,c(i,length(MPC))] # subset sample and Kinetics columns (last column)
  tmp_df = aggregate(tmp_df[,1] ~ tmp_df$Kinetics, FUN = sum) # sum the data by kinetics
  colnames(tmp_df) = c("Kinetics",colnames(MPC)[i]) # rename column by kinetics and sample
  
  if ( i == 1) {
    df = tmp_df # create output df
  } else {
    df = merge(x = df, y = tmp_df, by = "Kinetics") # merge temporary df with output df
  }
}
MPC = df
rm(list=setdiff(ls(), c("MPC", "out_dir")))

# Factor by kinetics and melt data. Create new columns
MPC$Kinetics <- factor(MPC$Kinetics, c("Late", "Leaky", "Early", "IE", "Latent"))
MPC=melt(MPC)
MPC$Condition = sub("_c.*", "", MPC$variable)

# Remove Bulk condition
MPC = MPC[!grepl("Bulk", MPC$variable),]

# Factor data
MPC$Condition = factor(MPC$Condition, levels = c("WT_Latent", "dBGLF5_Latent", "WT_Early", "dBGLF5_Early", "WT_Late", "dBGLF5_Late"))

ggbarplot(position = position_stack(),label = TRUE, label.pos = "out",
          MPC, x = "Condition", y = "value", add = "mean_se",
          color = "Kinetics", fill = "Kinetics",
          #palette = c("skyblue","tomato"),
          ylab = "Mean polyA RNA Molecules Per Cell", xlab = FALSE, add.params = list(color = "black")) +
  #theme(legend.position = "none") +
  scale_fill_manual(values = c("red","orange","green","gray")) + 
  scale_color_manual(values = c("red","orange","green","gray")) +
  scale_y_continuous(breaks = seq(0, 1.3e5, by = 3e4))

# Export
ggsave(file.path(out_dir,"Figure_6C.svg"), dpi = 300, width = 6.2, height = 6.2, units = "in")

## print session info ##
print("Session Info below: ")
sessionInfo()
