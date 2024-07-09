rm(list=setdiff(ls(), "")) # remove all variables
## Install R packages if not already installed
#install.packages("rstudioapi")
#install.packages("svglite")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("flowCore")
#BiocManager::install("ggcyto")
#BiocManager::install("flowAI")

## Setup the environment

# Import packages
library(rstudioapi)
library(svglite)
library(flowCore) # interpret .fcs files
library(ggcyto) # advanced visualization using the ggplot nomenclature
library(flowAI) # "clean" the data (essentially removes events that deviate from the statistical norm)

# Set working directory
work_dir=dirname(dirname(dirname(getSourceEditorContext()$path)))

# Output directory
out_dir=file.path(work_dir, "Figures/Supplementary_Figure_4")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# location of files
base=file.path(work_dir, "Data/Flow_Data/Supplementary_Figure_4B")

# read in files
fcs_files <- list.files(file.path(base),pattern = ".fcs", full.names = TRUE)
fset <- read.flowSet(fcs_files)

# Perform Comp 
for (i in 1:length(spillover(fset[[1]]))) {
  if (!is.null(spillover(fset[[1]])[[i]]) == "TRUE") {
    spillover_value <- i
  }
  if (i == length(spillover(fset[[1]]))) {
    if (exists("spillover_value") == "FALSE") {
      print("FCS file does not contain a pre-calculated spillover matrix.")
      print("Must computate spillover matrix from a set of compensation controls.")
    }
  }
}

# This can than be accessed using the spillover method.
comp <- fsApply(fset, function(x) spillover(x)[[spillover_value]], simplify=FALSE)
fset_comp <- compensate(fset, comp)

# Cleaning
fset_comp_clean = flow_auto_qc(fset_comp, second_fractionFR = 1)

# Change channel names for the entire flowSet
colnames(fset_comp_clean)[colnames(fset_comp_clean)=="FL1-A"] <- "GFP"
colnames(fset_comp_clean)[colnames(fset_comp_clean)=="FL3-A"] <- "RFP"
colnames(fset_comp_clean)[colnames(fset_comp_clean)=="FL6-A"] <- "DAPI"
colnames(fset_comp_clean)[colnames(fset_comp_clean)=="FL10-A"] <- "AF647"
colnames(fset_comp_clean) <- sub("\\-", ".", colnames(fset_comp_clean))

# to enable gating, transform flowSet object to GatingSete object
gs <- GatingSet(fset_comp_clean)

# Transform data using estimateLogicle command
trans <- estimateLogicle(gs[[3]], c("GFP","RFP","DAPI","AF647")) # estimate transformation
gs <- transform(gs, trans) # apply transformation
rm(list=setdiff(ls(), c("gs","out_dir")))

# set debris gate
g.debris <- polygonGate(filterId = "Debris","FSC.A"=c(2e5,6e5,7.5e5,7e5,3.5e5,2e5),"SSC.A"=c(1.25e5,2e5,5.5e5,8e5,5.5e5,2.5e5)) # define gate
#ggcyto(gs,aes(x=FSC.A,y=SSC.A),subset="root")+geom_hex(bins = 150)+geom_gate(g.debris) +
#  ggcyto_par_set(limits = list(x = c(0,1e6), y = c(0, 1e6)))+ geom_stats()
gs_pop_add(gs,g.debris,parent="root")
recompute(gs)

# set singlet gate
g.singlets <- polygonGate(filterId = "Singlets","FSC.A"=c(1.9e5,10e5,10e5,6e5,1.9e5),"FSC.H"=c(0.6e5,4.3e5,6e5,5e5,3e5)) # define gate
#ggcyto(gs,aes(x=FSC.A,y=FSC.H),subset="Debris")+geom_hex(bins = 200)+geom_gate(g.singlets) +
#  ggcyto_par_set(limits = list(x = c(2e5,9.5e5), y = c(0.5e5, 4.5e5))) + geom_stats()
gs_pop_add(gs,g.singlets,parent="Debris")
recompute(gs)

# set live gate
g.live <- polygonGate(filterId = "Live","FSC.A"=c(1.6e5,9.6e5,9.6e5,1.6e5),"DAPI"=c(-0.7,-0.7,1.6,1.6)) # define gate
#ggcyto(gs,aes(x=FSC.A,y=DAPI),subset="Singlets")+geom_hex(bins = 200)+geom_gate(g.live) +
#  ggcyto_par_set(limits = list(x = c(2e5,9.5e5), y = c(-1, 4))) + geom_stats()
gs_pop_add(gs,g.live,parent="Singlets")
recompute(gs)

# Set gates
RFP=1.75
GFP=2.1
gp350_neg=1.45

g.latent_RFP <- polygonGate(filterId = "Latent_RFP", "RFP" = c(0,RFP,RFP,0), "GFP" = c(0.5,0.5,1.25,1.25))
g.gfp <- polygonGate(filterId = "Early", "RFP" = c(0,RFP-0.2,RFP-0.2,0), "GFP" = c(GFP,GFP,4,4))
g.rfp <- polygonGate(filterId = "Late", "RFP" = c(RFP,4,4,RFP), "GFP" = c(GFP,GFP,4,4))
g.lytic <- polygonGate(filterId = "Lytic", "FSC.A" = c(1.9e5,7.5e5,7.5e5,1.9e5), "GFP" = c(GFP,GFP,4,4))
g.latent <- polygonGate(filterId = "Latent", "FSC.A" = c(1.9e5,7.5e5,7.5e5,1.9e5), "GFP" = c(0.5,0.5,1.5,1.5))

g.DN <- polygonGate(filterId = "DN", "RFP" = c(-3,-3,RFP,RFP), "AF647" = c(-3,gp350_neg,gp350_neg,-3))
g.RFP_pos <- polygonGate(filterId = "RFP_pos", "RFP" = c(RFP,RFP,6,6), "AF647" = c(-3,gp350_neg,gp350_neg,-3))
g.AF647_pos <- polygonGate(filterId = "AF647_pos", "RFP" = c(-3,-3,RFP,RFP), "AF647" = c(gp350_neg,6,6,gp350_neg))
g.DP <- polygonGate(filterId = "DP", "RFP" = c(RFP,RFP,6,6), "AF647" = c(gp350_neg,6,6,gp350_neg))

gs_pop_add(gs,g.latent,parent="Live")
gs_pop_add(gs,g.gfp,parent="Live")
gs_pop_add(gs,g.rfp,parent="Live")
gs_pop_add(gs,g.lytic,parent="Live")
recompute(gs)

ggcyto(gs, aes(x = "RFP", y = "GFP"),  subset = "Live") + geom_hex(bins=200)+
  geom_gate(g.latent_RFP) + geom_gate(g.gfp) + geom_gate(g.rfp) +
  ggcyto_par_set(limits = list(x = c(0,4), y = c(0.5,4))) + geom_stats()

ggcyto(gs, aes(x = "AF647", y = "RFP"),  subset = "Latent") + geom_hex(bins=200)+
  geom_gate(g.DN) + geom_gate(g.RFP_pos) + geom_gate(g.AF647_pos) + geom_gate(g.DP) +
  ggcyto_par_set(limits = list(x = c(0,4), y = c(0.5,4))) + geom_stats()

ggcyto(gs, aes(x = "AF647", y = "RFP"),  subset = "Lytic") + geom_hex(bins=200)+
  geom_gate(g.DN) + geom_gate(g.RFP_pos) + geom_gate(g.AF647_pos) + geom_gate(g.DP) +
  ggcyto_par_set(limits = list(x = c(0,4), y = c(0.5,4))) + geom_stats()

ggcyto(gs, aes(x = "FSC.A", y = "GFP"),  subset = "Live") + geom_hex(bins=200)+
  geom_gate(g.latent) + geom_gate(g.lytic) + 
  ggcyto_par_set(limits = list(x = c(2e5,8e5), y = c(0.5,4))) + geom_stats()

# Get theme functions
theme_white = theme(text = element_text(color = "white"), axis.text = element_text(color = "white"),
                    axis.line = element_line(color = "white"), axis.ticks = element_line(color = "white"),
                    strip.text.x = element_text(color = "white"),panel.border = element_rect(color = "white"))
theme_mine <- function(base_size = 18) {
  # Starts with theme_grey and then modify some parts
  theme_bw(base_size = base_size) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(size = 18),
      strip.text.y = element_text(size = 18),
      axis.text.x = element_text(size=14),
      axis.text.y = element_text(size=14,hjust=1),
      axis.ticks =  element_line(colour = "black"), 
      axis.title.x= element_text(size=16),
      axis.title.y= element_text(size=16,angle=90),
      panel.background = element_blank(), 
      #panel.border =element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.spacing = unit(1.0, "lines"), 
      plot.background = element_rect(fill = "white", color = "white"), 
      plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      legend.position = "none"
    )
}

# Plot mGL vs FSC-A
p1 = ggcyto(gs[[-3]], aes(x = "FSC.A", y = "GFP"),  subset = "Live") + 
  geom_hex(bins=100) + 
  #geom_gate(c("Latent","Lytic")) + 
  #geom_stats(adjust= c(0.2,0.1), gate = "Latent",digits = 2) + 
  #geom_stats(adjust= c(0.76,0.9), gate = "Lytic",digits = 2) + # 4 column
  theme_mine() + theme(legend.position = "none", plot.title = element_blank()) +
  #ggcyto_par_set(limits = list(x = c(0,1.2e6), y = c(-0.15,3.85))) +
  #axis_y_inverse_trans() +
  facet_wrap(~name, ncol = 4) +
  labs(x = "FSC-A",y = "mGreenLantern")

# Plot mSI vs gp350
p2 = ggcyto(gs[[-3]], aes(x = "AF647", y = "RFP"),  subset = "Lytic") + 
  geom_hex(bins=100) +
  #geom_gate(g.DN) + geom_gate(g.RFP_pos) + geom_gate(g.AF647_pos) + geom_gate(g.DP) + 
  #geom_stats(adjust= c(0.5,0.5)) +
  theme_mine() + theme(legend.position = "none", plot.title = element_blank()) +
  ggcyto_par_set(limits = list(x = c(0,3.5), y = c(-0.25,4))) +
  #axis_y_inverse_trans() + axis_x_inverse_trans() +
  facet_wrap(~name, ncol = 4) +
  labs(y = "mScarlet-I", x = "Alexa Fluor 647-anti-gp350")

p3 = ggcyto(gs[[-3]], aes(x = "AF647", y = "RFP"),  subset = "Latent") + 
  geom_hex(bins=100) +
  #geom_gate(g.DN) + geom_gate(g.RFP_pos) + geom_gate(g.AF647_pos) + geom_gate(g.DP)  +
  #geom_stats(adjust= c(0.5,0.5)) +
  theme_mine() + theme(legend.position = "none", plot.title = element_blank()) +
  ggcyto_par_set(limits = list(x = c(0,3.5), y = c(-0.25,4))) +
  #axis_y_inverse_trans() + axis_x_inverse_trans() +
  facet_wrap(~name, ncol = 4) +
  labs(y = "mScarlet-I", x = "Alexa Fluor 647-anti-gp350")

p2 + geom_gate(g.DN) + geom_gate(g.RFP_pos) + geom_gate(g.AF647_pos) + geom_gate(g.DP) + geom_stats(adjust= c(0.5,0.5))
p3 + geom_gate(g.DN) + geom_gate(g.RFP_pos) + geom_gate(g.AF647_pos) + geom_gate(g.DP) + geom_stats(adjust= c(0.5,0.5))

## plot data points without borders+labels and save as PNG
# FSC-A versus mGreenLantern
p1 + theme_white
ggsave(file.path(out_dir,"Supplementary_Figure_4B.mGL_vs_FSC-A.png"), dpi = 300, width = 10.24, height = 5.21, units = "in")

# gp350 versus mScarlet-I (Lytic)
p2 + theme_white
ggsave(file.path(out_dir,"Supplementary_Figure_4B.mSI_vs_gp350.Lytic.png"), dpi = 300, width = 10.24, height = 5.21, units = "in")

# gp350 versus mScarlet-I (Latent)
p3 + theme_white
ggsave(file.path(out_dir,"Supplementary_Figure_4B.mSI_vs_gp350.Latent.png"), dpi = 300, width = 10.24, height = 5.21, units = "in")

## plot borders+labels and save as SVG
# FSC-A versus mGreenLantern
bg_plot_1 = p1
bg_plot_1$layers <- NULL
bg_plot_1+ geom_gate(c("Latent","Lytic")) + geom_stats(adjust= c(0.5,0.5))
ggsave(file.path(out_dir,"Supplementary_Figure_4B.mGL_vs_FSC-A.Background.svg"), dpi = 300, width = 10.24, height = 5.21, units = "in")

# gp350 versus mScarlet-I (Lytic)
bg_plot_2 = p2
bg_plot_2$layers <- NULL
bg_plot_2 + geom_gate(g.DN) + geom_gate(g.RFP_pos) + geom_gate(g.AF647_pos) + geom_gate(g.DP) + geom_stats(adjust= c(0.5,0.5)) 
ggsave(file.path(out_dir,"Supplementary_Figure_4B.mSI_vs_gp350.Lytic.Background.svg"), dpi = 300, width = 10.24, height = 5.21, units = "in")

# gp350 versus mScarlet-I (Latent)
bg_plot_3 = p3
bg_plot_3$layers <- NULL
bg_plot_3 + geom_gate(g.DN) + geom_gate(g.RFP_pos) + geom_gate(g.AF647_pos) + geom_gate(g.DP) + geom_stats(adjust= c(0.5,0.5)) 
ggsave(file.path(out_dir,"Supplementary_Figure_4B.mSI_vs_gp350.Latent.Background.svg"), dpi = 300, width = 10.24, height = 5.21, units = "in")

# Get statistics
gs_pop_add(gs,g.DN,parent="Latent")
gs_pop_add(gs,g.RFP_pos,parent="Latent")
gs_pop_add(gs,g.AF647_pos,parent="Latent")
gs_pop_add(gs,g.DP,parent="Latent")
gs_pop_add(gs,g.DN,parent="Lytic")
gs_pop_add(gs,g.RFP_pos,parent="Lytic")
gs_pop_add(gs,g.AF647_pos,parent="Lytic")
gs_pop_add(gs,g.DP,parent="Lytic")
recompute(gs)

pData(gs)
stats_df = data.frame(gs_pop_get_stats(gs[[1]]))
stats_df = stats_df[grepl("Live",stats_df$pop),]
stats_df$pop = sub(".*\\/Live/","",stats_df$pop)
stats_df$pop = sub(".*\\/Live","Live",stats_df$pop)

lytic =  stats_df$count[stats_df$pop == "Lytic"]/stats_df$count[stats_df$pop == "Live"]
latent =  stats_df$count[stats_df$pop == "Latent"]/stats_df$count[stats_df$pop == "Live"]

DN = stats_df$count[stats_df$pop == "Lytic/DN"]
DP = stats_df$count[stats_df$pop == "Lytic/DP"]
RFP = stats_df$count[stats_df$pop == "Lytic/RFP_pos"]
AF647 = stats_df$count[stats_df$pop == "Lytic/AF647_pos"]

# Create statistics table
if (exists("stats_df")) {rm(stats_df)}
stats_df = data.frame(matrix(nrow = 1, ncol = 5))
colnames(stats_df) = c("Sample","Condition","Pop","Counts","Percent")
for (i in 1:(length(gs))) {
  counts = data.frame(gs_pop_get_stats(gs[[i]])) # Get counts
  counts = counts[grepl("Live",counts$pop),]
  counts$pop = sub(".*\\/Live/","",counts$pop)
  counts$pop = sub(".*\\/Live","Live",counts$pop)
  
  # Get lytic counts
  Live_lytic = counts$count[grepl("root",counts$pop)]
  Total_lytic = counts$count[counts$pop == "Lytic"]
  DN_lytic = counts$count[counts$pop == "Lytic/DN"]
  RFP_lytic = counts$count[counts$pop == "Lytic/RFP_pos"]
  DP_lytic = counts$count[counts$pop == "Lytic/DP"]
  AF647_lytic = counts$count[counts$pop == "Lytic/AF647_pos"]
  
  # Get latent counts
  Live_latent = counts$count[grepl("root",counts$pop)]
  Total_latent = counts$count[counts$pop == "Latent"]
  DN_latent = counts$count[counts$pop == "Latent/DN"]
  RFP_latent = counts$count[counts$pop == "Latent/RFP_pos"]
  DP_latent = counts$count[counts$pop == "Latent/DP"]
  AF647_latent = counts$count[counts$pop == "Latent/AF647_pos"]
  
  # Get percentages relative to total latent and lytic live cells
  stats_df[NROW(stats_df)+1,] = c(sub("\\.fcs","",unique(counts$sample)), "Latent","Total",Total_latent,round((Total_latent/(Total_latent + Total_lytic))*100,2))
  stats_df[NROW(stats_df)+1,] = c(sub("\\.fcs","",unique(counts$sample)), "Latent","DN",DN_latent,round((DN_latent/(Total_latent + Total_lytic))*100,2))
  stats_df[NROW(stats_df)+1,] = c(sub("\\.fcs","",unique(counts$sample)), "Latent","RFP+",RFP_latent,round((RFP_latent/(Total_latent + Total_lytic))*100,2))
  stats_df[NROW(stats_df)+1,] = c(sub("\\.fcs","",unique(counts$sample)), "Latent","DP",DP_latent,round((DP_latent/(Total_latent + Total_lytic))*100,2))
  stats_df[NROW(stats_df)+1,] = c(sub("\\.fcs","",unique(counts$sample)), "Latent","AF647+",AF647_latent,round((AF647_latent/(Total_latent + Total_lytic))*100,2))
  
  stats_df[NROW(stats_df)+1,] = c(sub("\\.fcs","",unique(counts$sample)), "Lytic","Total",Total_lytic,round((Total_lytic/(Total_latent + Total_lytic))*100,2))
  stats_df[NROW(stats_df)+1,] = c(sub("\\.fcs","",unique(counts$sample)), "Lytic","DN",DN_lytic,round((DN_lytic/(Total_latent + Total_lytic))*100,2))
  stats_df[NROW(stats_df)+1,] = c(sub("\\.fcs","",unique(counts$sample)), "Lytic","RFP+",RFP_lytic,round((RFP_lytic/(Total_latent + Total_lytic))*100,2))
  stats_df[NROW(stats_df)+1,] = c(sub("\\.fcs","",unique(counts$sample)), "Lytic","DP",DP_lytic,round((DP_lytic/(Total_latent + Total_lytic))*100,2))
  stats_df[NROW(stats_df)+1,] = c(sub("\\.fcs","",unique(counts$sample)), "Lytic","AF647+",AF647_lytic,round((AF647_lytic/(Total_latent + Total_lytic))*100,2))
}
stats_df = stats_df[-1,]

# Export
write.csv(stats_df, file.path(out_dir,"Supplementary_Figure_4B.Stats.csv"), row.names = FALSE)

## print session info ##
print("Session Info below: ")
sessionInfo()
