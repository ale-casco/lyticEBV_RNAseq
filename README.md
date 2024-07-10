# lyticEBV_RNAseq

## Overview
This repository contains the scripts and workflows used in our study on bulk RNA-seq data from FAC-sorted latent, early lytic, and late lytic fractions of lymphoblastoid cell lines (LCLs) transformed by a dual-fluorescent lytic reporter Epstein-Barr virus (EBV).

Data processing was primarily conducted in an Ubuntu Linux environment. Additional processing and visualization were performed in R.

## Data Availability

The raw (genome-aligned BAM files) and processed (RSEM output) RNA-seq data used in this study are available through the NCBI GEO database under the accession number GSE271717. The full manuscript detailing our findings is available on Cell Reports (currently under revision).

## Repository Contents

- 'Linux_workflow': Contains the detailed Linux workflow and associated bash scripts employed for data processing.
- **Rscripts**: Contains the R scripts utilized for data analysis and generation of figures.
- **initiate.R**: This script facilitates the setup of the R environment. Executing it will launch the RStudio project **'lyticEBV_project.Rproj'**, which leverages the **'renv'** environment. This environment ensures reproducibility by encapsulating the specific R package versions used throughout the study. See **'renv.lock'** for a snapshot of the R environment.

## Prerequisites

### Ubuntu Linux

Ensure you have Ubuntu Linnux installed on your system. Processing of the raw data and bash scripts have been tested on Ubuntu 20.04 LTS. Bash is required for executing the data processing scripts.

### R and R Packages

This analysis relies on R for data processing and visualization. The R scripts manage dependencies using the **'renv'** package.

## Start Guide

1. Run the **'initiate.R'** script.
   This script will:
   - Use **'rstudioapi'** to define the project directory and open the R project **'lytic_EBV_project.Rproj'**.
   - Use **'renv'** to activate and restore the renv library from the lock file.

2. Run data analysis or figure visualization scripts.
   - Note: due to size limitations, some files are missing. See below under 'Data Exclusions'.

## Data Exclusions

The following files and directories are not included in this repository, but are referenced in the workflow:
| Path                                                                            | Notes                                                                                  |
| ------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------- |
| ~/Data/ARTDeco/preprocess_files                                                 | Run ARTDeco separately to generate these files.                                        |
| ~/Data/SpliceWiz/NxtSE ,<br> ~/Data/SpliceWiz/output ,<br> ~/Data/SpliceWiz/ref | <be>Run SpliceWiz.R to generate these files.                                           |
| ~/Data/bams                                                                     | Download these files from the NCBI GEO database under accession number GSE271717.      |
| ~/Data/refs/GRCh38.p14.ERCC.M81_DFLR.M81.chrEBV.inverted.fa                     | See README file in the ~/Data/refs directory for instructions on generating this file. |

For convenience, a README file with a tree visualization of the file structure and disk usage for each excluded directory is included in its place.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions or further information about this project, please contact:
- **Alejandro Casco, PhD**: casco@wisc.edu
- **Eric Johannsen, MD**: ejohannsen@medicine.wisc.edu.
