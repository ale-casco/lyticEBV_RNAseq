rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("svglite")
#install.packages("ggplot2")
#install.packages("scales")

## Setup the environment

# Import packages
library(rstudioapi)
library(svglite)
library(ggplot2)
library(scales)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Output directory
out_dir=file.path(work_dir, "Figures/Supplementary_Figure_3")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Load molecules per cell (MPC) data
MPC_path = file.path(work_dir, "Data/MPC/MPC.csv")
if (!file.exists(MPC_path)) {stop(" MPC.csv file does not exist.\n\tRun MPC.R script to generate.", call. = FALSE)} # check file exists
MPC = read.csv(MPC_path,check.names=FALSE)

# Get mean expression across WT early and WT late conditions
MPC$baseMean = rowMeans(MPC[,grepl("WT_Early_c|WT_Late_c",colnames(MPC))])
MPC = MPC[,c("SYMBOL","baseMean")]

# Import DEG file of WT Late vs Early
deg_path = file.path(work_dir, "Data/DESeq2/ERCCnorm_DE.WT_Late.vs.WT_Early.csv")
if (!file.exists(deg_path)) {stop(" No DESeq2 file exists.\n\tRun DESeq2.R script to generate.", call. = FALSE)} # check file exists
df = read.csv(deg_path,check.names = TRUE)
df = df[!grepl("ENSG",df$ENSEMBL) & !grepl("ERCC-",df$SYMBOL),c("SYMBOL","log2FoldChange","padj")] # get only EBV genes

# Merge mean MPC data
df = merge(x = df, y = MPC, by = "SYMBOL")
rm(MPC)

# Rename some EBV genes
df$SYMBOL = sub("BILF2-tPT2A-mScarlet-I", "BILF2-mSI", df$SYMBOL)
df$SYMBOL = sub("mGreenLantern", "mGL", df$SYMBOL)
df$SYMBOL = sub("BSLF2/BMLF1", "SM", df$SYMBOL)
df$SYMBOL = sub("BBLF2/BBLF3", "BBLF2/3", df$SYMBOL)
df$SYMBOL = sub("BGRF1/BDRF1", "BG/DRF1", df$SYMBOL)

# Remove genes
df = df[!grepl("EBER1|EBER2|BWRF1", df$SYMBOL),]

# get lists of EBV genes by their expression kinetics
Early = c("BZLF1","BRLF1","mGL", "BHLF1","BFLF2","BORF2","BaRF1","BMRF1","BMRF2","SM","BSLF1","BSRF1","BLLF3","BLLF2","BRRF1","BKRF3",
          "BBLF4","BBLF2/3","BGLF5","BGLF4","BGLF3.5","BGLF3","BDLF3.5","BVRF1","BVLF1","BALF2","BALF1","BARF1",
          "BCRF1","BFRF2","BG/DRF1","BTRF1","LF3","LF2","LF1","BILF1","BALF5","BNLF2a/b")
Leaky = c("BFLF1","BORF1","BLRF1","BLRF2","BRRF2","BKRF2","BBRF1","BBRF2","BBRF3","BBLF1","BDLF3","BcRF1","BXLF2","BdRF1",
          "BALF4","BALF3","BFRF0.5","BKRF4")
Late = c("BNRF1","BFRF1","BFRF3","BPLF1","BOLF1","BLLF1","BZLF2","BGLF2","BGLF1","BDLF2","BDLF1","BcLF1","BXRF1",
         "BVRF2","BILF2-mSI","BDLF4","LMP1_lytic")
Latent = c("LMP2", "LMP2A", "LMP2B","EBER1","EBER2","EBNALP","EBNA2","EBNA3A","EBNA3B","EBNA3C","EBNA1","BART","LMP1","BHRF1")
df$Kinetics =  ifelse(df$SYMBOL %in% Early, "Early",
                       ifelse(df$SYMBOL %in% Leaky, "Leaky",
                              ifelse(df$SYMBOL %in% Late, "Late",
                                     ifelse(df$SYMBOL %in% Latent, "Latent", NA))))
rm(list=setdiff(ls(), c("df", "out_dir")))

# Annotate df by DEG
df$Sig = ifelse(df$padj <= 0.05, "p < 0.05","NS") # significance
df$Kinetics = factor(df$Kinetics, levels = c("Latent","Early","Leaky","Late")) # factor kinetics
df$Shape = ifelse(abs(df$log2FoldChange) > 5, 17, 19) # if values are beyond the limits, change shape to triangle 

# Generate plot
ggplot(df, aes(y = log2FoldChange, x = reorder(SYMBOL, log2FoldChange), color = Sig, size = baseMean)) +
  geom_point(alpha = 0.6, shape = df$Shape) +
  #scale_color_gradient(low = "red", high = "blue") +
  scale_color_manual(values = c("#bebebe","blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), text = element_text(size = 10)) +
  ylab("") + 
  xlab("") + 
  ggtitle("WT Late versus Early") +
  geom_hline(yintercept = 1, color = "black", linetype = "dashed") +
  geom_hline(yintercept = -1, color = "black", linetype = "dashed") +
  facet_grid(~Kinetics,scales = "free", space='free') +
  scale_y_continuous(limits = c(-5,5), expand = c(0,0.25), oob = squish) +
  scale_size(breaks = c(10,1000,3000,6000,11000), range = c(3,12))

## Get point data
ggsave(file.path(out_dir,"Supplementary_Figure_3B.svg"), dpi = 300, width = 12, height = 3.5, units = "in")

## print session info ##
print("Session Info below: ")
sessionInfo()
