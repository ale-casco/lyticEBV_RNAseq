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
out_dir=file.path(work_dir, "Figures/Figure_1C")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# location of files
base = file.path(work_dir, "Data/Flow_Data/Figure_1C")

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
colnames(fset_comp_clean) <- sub("\\-", ".", colnames(fset_comp_clean))

# to enable gating, transform flowSet object to GatingSete object
gs <- GatingSet(fset_comp_clean)

# Transform data using estimateLogicle command
trans <- estimateLogicle(gs[[2]], c("GFP","RFP","DAPI")) # estimate transformation
gs <- transform(gs, trans) # apply transformation
rm(list=setdiff(ls(), c("gs","out_dir")))

# set debris gate
g.debris <- polygonGate(filterId = "Debris","FSC.A"=c(3.5e5,6.5e5,7.5e5,7e5,4e5,2e5),"SSC.A"=c(0.1e5,0.5e5,2.5e5,5.5e5,6e5,2.5e5))
#ggcyto(gs,aes(x=FSC.A,y=SSC.A),subset="root")+geom_hex(bins = 256)+geom_gate(g.debris) +
#  ggcyto_par_set(limits = list(x = c(0,1e6), y = c(0, 1e6))) + geom_stats()
gs_pop_add(gs,g.debris,parent="root")
recompute(gs)

# set singlet gate
g.singlets <- polygonGate(filterId = "Singlets","FSC.A"=c(1.9e5,7.5e5,7.5e5,6e5,1.9e5),"FSC.H"=c(0.5e5,3.3e5,4.6e5,4.4e5,3e5))
#ggcyto(gs,aes(x=FSC.A,y=FSC.H),subset="Debris")+geom_hex(bins = 200)+geom_gate(g.singlets) +
#  ggcyto_par_set(limits = list(x = c(2e5,7.5e5), y = c(0.5e5, 4.5e5))) + geom_stats()
gs_pop_add(gs,g.singlets,parent="Debris")
recompute(gs)

# set live gate
g.live <- polygonGate(filterId = "Live","FSC.A"=c(1.9e5,7.5e5,7.5e5,1.9e5),"DAPI"=c(-1,-1,1.35,1.35))
#ggcyto(gs,aes(x=FSC.A,y=DAPI),subset="Singlets")+geom_hex(bins = 200)+geom_gate(g.live) +
#  ggcyto_par_set(limits = list(x = c(2e5,7.5e5), y = c(-1, 4))) + geom_stats()
gs_pop_add(gs,g.live,parent="Singlets")
recompute(gs)

# Set GFP-/RFP- (latent), GFP+/RFP- (early), and GFP+/RFP+ (late) gates
g.latent <- polygonGate(filterId = "Latent", "RFP" = c(0.2,1.5,1.85,1.85,1.6,0.2), "GFP" = c(0.4,0.4,1.5,1.83,1.83,0.96))
g.early <- polygonGate(filterId = "Early", "RFP" = c(0.2,1.33,1.33,0.2), "GFP" = c(1.13,1.86221428571429,3.8,3.8))
g.late <- polygonGate(filterId = "Late", "RFP" = c(1.53,3.9,3.9,1.53), "GFP" = c(2.03,2.03,3.8,3.8))
#ggcyto(gs, aes(x = "RFP", y = "GFP"),  subset = "Live") + geom_hex(bins=200)+
#  geom_gate(g.latent) + geom_gate(g.early) + geom_gate(g.late) +
#  ggcyto_par_set(limits = list(x = c(0,4), y = c(0.5,4))) + geom_stats()
gs_pop_add(gs,g.latent,parent="Live")
gs_pop_add(gs,g.early,parent="Live")
gs_pop_add(gs,g.late,parent="Live")
recompute(gs)

rm(list=setdiff(ls(), c("gs","out_dir")))

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

# Change names
pData(gs)$name = c(" Control LCL", "Stimulated +GCV", "Stimulated", " Unstimulated")

# plot data points without borders+labels and save as PNG
ggcyto(gs, aes(x = "RFP", y = "GFP"),  subset = "Live") + geom_hex(bins=100) + 
  #geom_gate(c("Latent","Early","Late")) + 
  #geom_stats(adjust= c(0.2,0.1), gate = "Latent",digits = 2) + 
  #geom_stats(adjust= c(0.23,0.93), gate = "Early",digits = 2) + 
  #geom_stats(adjust= c(0.8,0.9), gate = "Late",digits = 1) + # 2 column
  #geom_stats(adjust= c(0.76,0.9), gate = "Late",digits = 1) + # 4 column
  theme_mine() + theme(legend.position = "none", plot.title = element_blank()) +
  ggcyto_par_set(limits = list(x = c(0,4), y = c(0.5,4))) +
  #axis_y_inverse_trans() + axis_x_inverse_trans() +
  facet_wrap(~name, ncol = 4) +
  labs(x = "mScarlet-I",y = "mGreenLantern") +
  theme_white

ggsave(file.path(out_dir,"Figure_1C.png"), dpi = 300, width = 10, height = 3.3, units = "in")

# plot borders+labels and save as SVG
ggcyto(gs, aes(x = "RFP", y = "GFP"),  subset = "Live") + geom_gate(c("Latent","Early","Late")) + 
  geom_stats(adjust= c(0.2,0.1), gate = "Latent") + 
  geom_stats(adjust= c(0.23,0.93), gate = "Early") + 
  #geom_stats(adjust= c(0.8,0.9), gate = "Late",digits = 1) + # 2 column
  geom_stats(adjust= c(0.76,0.9), gate = "Late") + # 4 column
  theme_mine() + theme(legend.position = "none", plot.title = element_blank()) +
  ggcyto_par_set(limits = list(x = c(0,4), y = c(0.5,4))) +
  #axis_y_inverse_trans() + axis_x_inverse_trans() +
  facet_wrap(~name, ncol = 4) +
  labs(x = "mScarlet-I",y = "mGreenLantern")

ggsave(file.path(out_dir,"Figure_1C.Background.svg"), dpi = 300, width = 10, height = 3.3, units = "in")

# to make final figures, overlay plots

# Count number of cells per sample
cs <- gs_pop_get_data(gs, "Live", inverse.transform = FALSE)
fsApply(cs,nrow)

## print session info ##
print("Session Info below: ")
sessionInfo()
