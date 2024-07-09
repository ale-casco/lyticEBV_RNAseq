rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("reshape2")
#install.packages("ggplot2")
#install.packages("svglite")

## Setup the environment

# Import packages
library(rstudioapi)
library(ggplot2)
library(reshape2)
library(svglite)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Output directory
out_dir=file.path(work_dir, "Figures/Figure_4")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# location of files
base=file.path(work_dir,"Data/metrics/host")

# Sample information file
sampleinfo_dir=file.path(work_dir,"Data/refs/sampleinfo.txt")

# read in files
files = list.files(file.path(base), pattern=".metrics",full.names = TRUE)
sampleinfo = read.delim(sampleinfo_dir)
if(exists("df_out")) { rm(df_out) }
for ( i in files ) {
  tmp_df = t(read.delim(i,comment.char="#")[1,17:20])
  
  # Sum CDS and UTR into exonic
  tmp_df[2,] = tmp_df[1,] + tmp_df[2,]
  rownames(tmp_df)[2] = "PCT_EXONIC_BASES"
  tmp_df = as.data.frame(tmp_df[-1,])
  
  ID = sub("\\..*","",basename(i))
  tmp_sampleinfo = sampleinfo[grepl(ID,sampleinfo$UNIQUE_ID),]
  colnames(tmp_df) = paste0(tmp_sampleinfo$Condition,"_",tmp_sampleinfo$Replicate)
  
  if (!exists("df_out")) {
    df_out = tmp_df
  } else {
    df_out = cbind(df_out,tmp_df)
  }
}

df = melt(as.matrix(df_out))
df$Condition = sub("_c.*","",df$Var2)
df$Replicate = sub(".*_c","c",df$Var2)
rm(tmp_df,df_out,base,files,sampleinfo,i,ID,tmp_sampleinfo)

df$Var1 = ifelse(df$Var1 == "PCT_CODING_BASES","CDS",
                 ifelse(df$Var1 == "PCT_EXONIC_BASES","Exonic",
                        ifelse(df$Var1 == "PCT_INTRONIC_BASES","Intronic",
                               ifelse(df$Var1 == "PCT_INTERGENIC_BASES","Intergenic","Other"))))
df$Fraction  = sub(".*_","",df$Condition)

df$Var1 = factor(df$Var1, levels = c("Exonic","Intronic","Intergenic"))
df$Fraction = factor(df$Fraction, levels = c("Latent","Early","Late"))

# Remove dBGLF5
df_plot = df[!grepl("dBGLF5",df$Var2),]

# Factor data
df_plot$Fraction =  factor(df_plot$Fraction, levels = c("Latent","Early","Late"))
df_plot$Var1 = factor(df_plot$Var1, levels = c("Exonic","Intronic","Intergenic"))

stderror <- function(x) sd(x)/sqrt(length(x))
if(exists("df_out")) { rm(df_out) }
for ( i in unique(df_plot$Condition)) {
  tmp = df_plot[df_plot$Condition == i,]
  tmp_mean = aggregate(value ~ Var1 + Condition, data=tmp, mean)
  tmp_SEM = aggregate(value ~ Var1 + Condition, data=tmp, sd)
  tmp_mean$SEM = tmp_SEM$value
  if(!exists("df_out")) {
    df_out = tmp_mean
  } else {
    df_out = rbind(df_out,tmp_mean)
  }
}

df_out$Condition = sub("WT_","",df_out$Condition)
df_out$Condition = factor(df_out$Condition, levels = c("Latent","Early","Late"))

# Pie graph
ggplot(df_out,aes("",value,fill=Var1)) +
  geom_bar(stat = "identity", color = "white", linewidth = 1) +
  geom_text(aes(
    #label = paste0(round(value * 100,1), "%")
    label = paste0(round(df_out$value * 100,1), " ± ", round(df_out$SEM * 100, 1))
  ), 
  position = position_stack(vjust = 0.5), 
  color = "white", size = 3) +
  coord_polar(theta = "y") +
  facet_wrap(~Condition, ncol = 3) +
  #scale_fill_manual(values = c("#0048cc", "#cc8400")) +
  scale_fill_brewer('', palette = 'Set1') +
  theme_void() +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = c("#008b00","#ff0000","#ff8c00"))

ggsave(file.path(out_dir,"Figure_4A.svg"), dpi = 300, width = 5.21, height = 5.21, units = "in")

## print session info ##
print("Session Info below: ")
sessionInfo()

