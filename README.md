# lyticEBV_RNAseq

## Overview
This repository contains the scripts and workflows used in our study on bulk RNA-seq data from FAC-sorted latent, early lytic, and late lytic fractions of lymphoblastoid cell lines (LCLs) transformed by a dual-fluorescent lytic reporter Epstein-Barr virus (EBV). Our analysis was performed using an Ubuntu Linux environment, utilizing bash for data processing and R for data visualization. 

## Data Availability

Raw and processed RNA-seq data used in this study can be accessed through the NCBI GEO data base under the accession number <to be provided>. The full manuscript detailing our findings is available on Cell Reports <insert website>.

## Repository Contents

- **'workflow/'**: Contains bash scripts for data processing.
- **'Rscripts/'**: Contains R scripts for data analysis and visualization.
- **'renv.lock/'**: Snapshop of the R environment to ensure reproducibility of the analysis.
- **'lyticEBV_project.Rproj/'**: RStudio project file for the analysis.
- **'initiate/'**: Script to set up the R environment.

## Prerequisites

### Ubuntu Linux

Ensure you have Ubuntu Linnux installed on your system. Processing of the raw data and bash scripts have been tested on Ubuntu 20.04 LTS. Bash is required for executing the data processing scripts.

### R and R Packages

This analysis relies on R for data processing and visualization. The R scripts manage dependencies using the **'renv'** package. Please note that the **'renv'** environment was created on a Windows system and may not be fully compatible with Ubuntu or MacOS.

## Installation

1. **Clone the repository:**

   ```bash
   git clone https://github.com/ale-casco/lyticEBV_RNAseq.git
   cd repository-name

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
