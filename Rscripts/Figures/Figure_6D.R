rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("ggplot2")
#install.packages("svglite")

## Setup the environment

# Import packages
library(rstudioapi)
library(ggplot2)
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

# Import DESeq2 differential expression files
dge_early_path = file.path(work_dir, "Data/DESeq2/ERCCnorm_DE.WT_Early.vs.dBGLF5_Early.csv")
dge_late_path = file.path(work_dir, "Data/DESeq2/ERCCnorm_DE.WT_Late.vs.dBGLF5_Late.csv")
if (!file.exists(dge_early_path)) {stop(" No DESeq2 files exist.\n\tRun DESeq2.R script to generate.", call. = FALSE)} # check file exists
df_early = read.csv(dge_early_path,check.names = TRUE)
df_late = read.csv(dge_late_path,check.names = TRUE)

# Calculate mean MPC
MPC$Early_baseMean = rowMeans(MPC[,grepl("WT_Early_c|dBGLF5_Early_c",colnames(MPC))])
MPC$Late_baseMean = rowMeans(MPC[,grepl("WT_Late_c|dBGLF5_Late_c",colnames(MPC))])
MPC = MPC[,c("SYMBOL","Early_baseMean","Late_baseMean")]

# Keep only EBV genes and subset columns
df_early = df_early[!grepl("ENSG",df_early$ENSEMBL) & !grepl("ERCC-",df_early$SYMBOL),c("SYMBOL","log2FoldChange","padj")]
df_late = df_late[!grepl("ENSG",df_late$ENSEMBL) & !grepl("ERCC-",df_late$SYMBOL),c("SYMBOL","log2FoldChange","padj")]
colnames(df_early)[-1] = paste0("Early_",colnames(df_early)[-1])
colnames(df_late)[-1] = paste0("Late_",colnames(df_late)[-1])
df = merge(x = df_early, y = df_late, by = "SYMBOL", all = TRUE)
rm(df_early,df_late)

# remove genes, including BGLF5 
df = df[!grepl("EBER1|EBER2|BWRF1|BGLF5", df$SYMBOL),]

# Merge MPC and differential expression data
df = merge(x = df, y = MPC, by = "SYMBOL")

# get lists of EBV genes by their expression kinetics
df$SYMBOL = sub("BILF2-tPT2A-mScarlet-I", "BILF2-mSI", df$SYMBOL)
df$SYMBOL = sub("mGreenLantern", "mGL", df$SYMBOL)
df$SYMBOL = sub("BSLF2/BMLF1", "SM", df$SYMBOL)
df$SYMBOL = sub("BBLF2/BBLF3", "BBLF2/3", df$SYMBOL)
df$SYMBOL = sub("BGRF1/BDRF1", "BG/DRF1", df$SYMBOL)

# Get kinetic class
Early = c("BZLF1","BRLF1","mGL", "BHLF1","BFLF2","BORF2","BaRF1","BMRF1","BMRF2","SM","BSLF1","BSRF1","BLLF3","BLLF2","BRRF1","BKRF3",
          "BBLF4","BBLF2/3","BGLF5","BGLF4","BGLF3.5","BGLF3","BDLF3.5","BVRF1","BVLF1","BALF2","BALF1","BARF1",
          "BCRF1","BFRF2","BG/DRF1","BTRF1","LF3","LF2","LF1","BILF1","BALF5","BNLF2a/b")
Leaky = c("BFLF1","BORF1","BLRF1","BLRF2","BRRF2","BKRF2","BBRF1","BBRF2","BBRF3","BBLF1","BDLF3","BcRF1","BXLF2","BdRF1",
          "BALF4","BALF3","BFRF0.5","BKRF4")
Late = c("BNRF1","BFRF1","BFRF3","BPLF1","BOLF1","BLLF1","BZLF2","BGLF2","BGLF1","BDLF2","BDLF1","BcLF1","BXRF1",
         "BVRF2","BILF2-mSI","BDLF4","LMP1_lytic")
Latent = c("LMP2", "LMP2A", "LMP2B","EBER1","EBER2","EBNALP","EBNA2","EBNA3A","EBNA3B","EBNA3C","EBNA1","BART","LMP1","BHRF1")


# Get MA plot info
df$Kinetics = ifelse(df$SYMBOL %in% Early, "Early",
                     ifelse(df$SYMBOL %in% Leaky, "Leaky",
                            ifelse(df$SYMBOL %in% Late, "Late",
                                   ifelse(df$SYMBOL %in% Latent, "Latent", "Unclassified"))))

df$Early = "Early"
df$Early_log2FoldChange = ifelse(is.na(df$Early_log2FoldChange), 0,df$Early_log2FoldChange)
df$Early_padj = ifelse(is.na(df$Early_padj), 1, df$Early_padj)
df$Early_Sig = ifelse(df$Early_padj <= 0.05, "p < 0.05","NS")
df$Kinetics = factor(df$Kinetics, levels = c("Latent","Early","Leaky","Late"))

df$Late = "Late"
df$Late_log2FoldChange = ifelse(is.na(df$Late_log2FoldChange), 0,df$Late_log2FoldChange)
df$Late_padj = ifelse(is.na(df$Late_padj), 1, df$Late_padj)
df$Late_Sig = ifelse(df$Late_padj <= 0.05, "p < 0.05","NS")

df$SYMBOL_Early = reorder(df$SYMBOL, -df$Early_log2FoldChange)
df$SYMBOL_Late = reorder(df$SYMBOL, -df$Late_log2FoldChange)

# Create legend within graph (otherwise early and late plots are not identical dimensions)
df = rbind(df, c("BGLF1", -2.6, 0.01, -2.6, 0.01, 10, 10, "Late", "Early", "Legend", "Late", "Legend", "BGLF1","BFRF1"))
df = rbind(df, c("BGLF1", -2.2, 0.01, -2.2, 0.01, 1000, 1000, "Late", "Early", "Legend", "Late", "Legend", "BGLF1","BFRF1"))
df = rbind(df, c("BGLF1", -1.8, 0.01, -1.8, 0.01, 3000, 3000, "Late", "Early", "Legend", "Late", "Legend", "BGLF1","BFRF1"))
df = rbind(df, c("BGLF1", -1.4, 0.01, -1.4, 0.01, 6000, 6000, "Late", "Early", "Legend", "Late", "Legend", "BGLF1","BFRF1"))
df = rbind(df, c("BGLF1", -1.0, 0.01, -1.0, 0.01, 11000, 11000, "Late", "Early", "Legend", "Late", "Legend", "BGLF1","BFRF1"))

# Make values numeric
for (i in 2:7) { df[,i] = as.numeric(df[,i])}

# Generate plots
early_plot = ggplot(df, aes(y = -Early_log2FoldChange, x = SYMBOL_Early, color = Early_Sig, size = Early_baseMean)) +
  geom_point(alpha = 0.6) +
  #scale_color_gradient(low = "red", high = "blue") +
  scale_color_manual(values = c("#666666","#bebebe","blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), text = element_text(size = 10), legend.position = "none") +
  ylab("") + 
  xlab("") + 
  ggtitle("Late dBGLF5 versus WT") +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
  geom_hline(yintercept = -1, color = "black", linetype = "dashed") +
  facet_grid(~Kinetics,scales = "free", space='free') +
  scale_size(breaks = c(10,1000,3000,6000,11000), range = c(3,12)) +
  ylim(-2.6,2.6) +
  guides(size=guide_legend(title="baseMean"),color=guide_legend(title="Sig"))

late_plot = ggplot(df, aes(y = -Late_log2FoldChange, x = SYMBOL_Late, color = Late_Sig, size = Late_baseMean)) +
  geom_point(alpha = 0.6) +
  #scale_color_gradient(low = "red", high = "blue") +
  scale_color_manual(values = c("#666666","#bebebe","blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), text = element_text(size = 10), legend.position = "none") +
  ylab("") + 
  xlab("") + 
  ggtitle("Late dBGLF5 versus WT") +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
  geom_hline(yintercept = -1, color = "black", linetype = "dashed") +
  facet_grid(~Kinetics,scales = "free", space='free') +
  scale_size(breaks = c(10,1000,3000,6000,11000), range = c(3,12)) +
  ylim(-2.6,2.6) +
  guides(size=guide_legend(title="baseMean"),color=guide_legend(title="Sig"))

## Get point data
# Early
early_plot
ggsave(file.path(out_dir,"Figure_6D.Early.svg"), dpi = 300, width = 12, height = 3.5, units = "in")

# Late
late_plot
ggsave(file.path(out_dir,"Figure_6D.Late.svg"), dpi = 300, width = 12, height = 3.5, units = "in")

## print session info ##
print("Session Info below: ")
sessionInfo()
