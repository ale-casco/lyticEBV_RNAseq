rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("reshape2")
#install.packages("ggpubr")
#install.packages("rstatix")
#install.packages("svglite")

## Setup the environment

# Import packages
library(rstudioapi)
library(reshape2)
library(ggpubr)
library(rstatix)
library(svglite)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Output directory
out_dir=file.path(work_dir, "Figures/Figure_5")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Load molecules per cell (MPC) data
MPC_path = file.path(work_dir, "Data/MPC/MPC.csv")
if (!file.exists(MPC_path)) {stop(" MPC.csv file does not exist.\n\tRun MPC.R script to generate.", call. = FALSE)} # check file exists
MPC=read.csv(MPC_path)
MPC = MPC[grepl("ENSG",MPC$ENSEMBL),grepl("ENSEMBL|SYMBOL|WT_Latent_c|WT_Early_c|WT_Late_c",colnames(MPC))]

# Subset genes of interest
MPC = MPC[grepl("OSM|CASTOR1|TBC1D10A",MPC$SYMBOL) & !grepl("MOSMO",MPC$SYMBOL),]

# Melt and prep data frame
MPC = melt(MPC)
MPC$Condition = sub("_c.*","",sub("WT_","",MPC$variable))
MPC = MPC[MPC$Condition != "Late",]

# Factor data
MPC$SYMBOL = factor(MPC$SYMBOL, levels = c("TBC1D10A","CASTOR1","OSM"))
MPC$Condition = factor(MPC$Condition, levels = c("Latent","Early","Late"))

# Perform Welch Two Sample t-test
stat.test <- MPC %>%
  group_by(SYMBOL) %>%
  t_test(value ~ Condition, var.equal = FALSE) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")
stat.test <- stat.test %>%
  add_xy_position(fun = "mean_se", x = "SYMBOL", dodge = 0.8)
stat.test

# Plot
plot_out = ggbarplot(MPC,"SYMBOL","value",fill = "Condition", palette = c("#BFBFBF","#00FF00"),
          label = FALSE, position = position_dodge(0.75),add = "mean_se",
          xlab = FALSE, ylab = "Molecules Per Cell") +
  stat_pvalue_manual(stat.test,  label = "p.adj.signif", tip.length = 0.01, bracket.nudge.y = -2)

# Export
plot_out
ggsave(file.path(out_dir,"Figure_5B.svg"), dpi = 300, width = 5.21, height = 5.21, units = "in")

## print session info ##
print("Session Info below: ")
sessionInfo()
