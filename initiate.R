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
      rstudioapi::openProject(file.path(project_dir,"lyticEBV_RNAseq.Rproj"))
    }
  } else {
    stop("Multiple or no .Rproj files found in the project directory.")
  }
} else {
  rstudioapi::openProject(file.path(project_dir,"lyticEBV_RNAseq.Rproj"))
}

## Install project manager package renv
if (!require("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Restore renv library
renv::load() # activate renv environment
renv::restore() # restore library from renv.lock
renv::status() # check if the project is in sync with lock file and report any inconsistencies

# Open data analysis and figure generation scripts
rstudioapi::navigateToFile(file.path(project_dir,"Rscripts/open_scripts.R"))
