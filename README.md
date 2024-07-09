# lyticEBV_RNAseq

## Overview
This repository contains the scripts and workflows used in our study on bulk RNA-seq data from FAC-sorted latent, early lytic, and late lytic fractions of lymphoblastoid cell lines (LCLs) transformed by a dual-fluorescent lytic reporter Epstein-Barr virus (EBV). Data processing was conducted using bash scripts within an Ubuntu Linux environment, with additional processing and visualization performed in R on a Windows platform.

## Data Availability

Raw (genome-aligned BAM files) and processed (RSEM output) RNA-seq data used in this study can be accessed through the NCBI GEO database under the accession number GSE271717. The full manuscript detailing our findings is available on Cell Reports <under revisions>.

## Repository Contents

- **Linux_workflow**: Contains the detailed Linux workflow and associated bash scripts employed for data processing.
- **Rscripts**: Contains the R scripts utilized for data analysis and generation of figures.
- **initiate.R**: This script facilitates the setup of the R environment. Executing it will launch the RStudio project **'lyticEBV_project.Rproj'**, which leverages the **'renv'** environment. This environment ensures reproducibility by encapsulating the specific R package versions used throughout the study. See **'renv.lock'** for a snapshot of the R environment.

## Prerequisites

### Ubuntu Linux

Ensure you have Ubuntu Linnux installed on your system. Processing of the raw data and bash scripts have been tested on Ubuntu 20.04 LTS. Bash is required for executing the data processing scripts.

### R and R Packages

This analysis relies on R for data processing and visualization. The R scripts manage dependencies using the **'renv'** package. Please note that the **'renv'** environment was created on a Windows system and may not be fully compatible with Ubuntu or MacOS.

## Initiate R environment

1. **Set up R environment**

   This script will:
   - Check if **'rstudioapi'** is installed.
   - Use **'rstudioapi'** to define the project directory.
   - Check if the project is open.
   - If the project is not open, it will open the project using **'openProject(file.path(project_dir, "lytic_EBV_project.Rproj"))'**.
   
By opening the project, the necessary **'renv'** project libraries will be loaded automatically, so there is no need to install or run **'renv::restore()'** manually.

## Contributing

We welcome contributions from the community. If you find any issues or have suggestions for improvements, please open an issue or submit a pull request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For any questions or further information, please contact casco@wisc.edu


Please cite accordingly if use this code in your work, including the original citations for any programs used.
