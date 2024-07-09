## Set scripts directory
scripts_dir = dirname(rstudioapi::getSourceEditorContext()$path)

## Open  all data analysis R scripts
data_analysis_scripts = paste(scripts_dir,"Data_analysis",c("gtf2gene_info.R","DESeq2.R","DESeq2_AltHypothesis.R","Raji_DESeq2.R","DESeq2_ARTDeco_readin.R",
                                                                     "DESeq2_ARTDeco_readthrough.R","DESeq2_ERCCnorm_SizeFactors.R","Expressed_genes.R","MPC.R",
                                                                     "SpliceWiz.R","IsoformSwitchAnalyzeR_external.R","IsoformSwitchAnalyzeR.R"),sep="/")
data_analysis_script_codes = lapply(data_analysis_scripts, rstudioapi::documentOpen)
#for (i in unlist(data_analysis_script_codes)) { .rs.api.documentClose(i) } # close data analysis scripts

## Open all main figure R scripts
figure_scripts = list.files(file.path(scripts_dir,"Figures"),pattern="\\.R$",full.names = TRUE)
figure_scripts = figure_scripts[!grepl("Supplementary_Figure_|Table_S",figure_scripts)]
figure_script_codes = lapply(figure_scripts, rstudioapi::documentOpen)
#for (i in unlist(figure_script_codes)) { .rs.api.documentClose(i) } # Main figure scripts

## Open all supplementary figure R scripts
supplementary_figure_scripts = list.files(file.path(scripts_dir,"Figures"),pattern="\\.R$",full.names = TRUE)
supplementary_figure_scripts = supplementary_figure_scripts[grepl("Supplementary_Figure_|Table_S",supplementary_figure_scripts)]
supplementary_figure_script_codes = lapply(supplementary_figure_scripts, rstudioapi::documentOpen)
#for (i in unlist(supplementary_figure_script_codes)) { .rs.api.documentClose(i) } # Supplementary figure scripts
