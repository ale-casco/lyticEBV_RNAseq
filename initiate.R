## Install rstudioapi package to get project directory
if (!require("rstudioapi", quietly = TRUE)) {
  install.packages("rstudioapi")
}

## Set project directory
project_dir = dirname(rstudioapi::getSourceEditorContext()$path)

## Open project
# Checks if project is open or not. If not, it will proceed to open.
if (!is.null(rstudioapi::getActiveProject())) {
  project_files <- list.files(rstudioapi::getActiveProject(), pattern = "\\.Rproj$")
  if(length(project_files)==1){
    project_file_name <- project_files[1]
    if (project_file_name != "lyticEBV_RNAseq.Rproj") {
      openProject(file.path(project_dir,"lyticEBV_RNAseq.Rproj"))
    }
  } else {
    stop("Multiple or no .Rproj files found in the project directory.")
  }
} else {
  openProject(file.path(project_dir,"lyticEBV_RNAseq.Rproj"))
}

## Install project manager package renv
if (!require("renv", quietly = TRUE)) {
  install.packages("renv")
}

# lockfile with the current state of dependencies in the project library
#renv::snapshot()
#renv::status()

## Open  all data analysis R scripts
data_analysis_scripts = paste(project_dir,"Rscripts/Data_analysis",c("gtf2gene_info.R","DESeq2.R","DESeq2_AltHypothesis.R","Raji_DESeq2.R","DESeq2_ARTDeco_readin.R",
                                                                     "DESeq2_ARTDeco_readthrough.R","DESeq2_ERCCnorm_SizeFactors.R","Expressed_genes.R","MPC.R",
                                                                     "SpliceWiz.R","IsoformSwitchAnalyzeR_external.R","IsoformSwitchAnalyzeR.R"),sep="/")
data_analysis_script_codes = lapply(data_analysis_scripts, rstudioapi::documentOpen)
#for (i in unlist(data_analysis_script_codes)) { .rs.api.documentClose(i) } # close data analysis scripts

## Open all main figure R scripts
figure_scripts = list.files(file.path(project_dir,"Rscripts/Figures"),pattern="\\.R$",full.names = TRUE)
figure_scripts = figure_scripts[!grepl("Supplementary_Figure_|Table_S",figure_scripts)]
figure_script_codes = lapply(figure_scripts, rstudioapi::documentOpen)
#for (i in unlist(figure_script_codes)) { .rs.api.documentClose(i) } # Main figure scripts

## Open all supplementary figure R scripts
supplementary_figure_scripts = list.files(file.path(project_dir,"Rscripts/Figures"),pattern="\\.R$",full.names = TRUE)
supplementary_figure_scripts = supplementary_figure_scripts[grepl("Supplementary_Figure_|Table_S",supplementary_figure_scripts)]
supplementary_figure_script_codes = lapply(supplementary_figure_scripts, rstudioapi::documentOpen)
#for (i in unlist(supplementary_figure_script_codes)) { .rs.api.documentClose(i) } # Supplementary figure scripts
