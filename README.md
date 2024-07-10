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

This analysis relies on R for data processing and visualization. The R scripts manage dependencies using the **'renv'** package.

## Start Guide

1. Run the **'initiate.R'** script.
   This script will:
   - Use **'rstudioapi'** to define the project directory and open the R project **'lytic_EBV_project.Rproj'**.
   - Use **'renv'** to activate and restore the renv library from the lock file.

2. Run data analysis or figure visualization scripts.
   - Note: due to size limitations, some files are missing. See below under 'Data Exclusions'.

## Data Exclusions

Files are replaced with a README file with tree visualization of the file structure and disk usage in the directory for missing files
| Path                                            | Notes                               |
| ----------------------------------------------- |:-----------------------------------:|
| ~/Data/ARTDeco/preprocess_files                 | Run ARTDeco                         |
| ~/Data/SpliceWiz/NxtSE\                         |                                     |
  ~/Data/SpliceWiz/output\                          Run SpliceWiz.R                     |
  ~/Data/SpliceWiz/ref\                                                                 |
| ~/Data/bams                                     | NCBI GEO accession number GSE271717 |
| ~/Data/refs/\                                   | See README in ~/Data/refs for       |
| GRCh38.p14.ERCC.M81_DFLR.M81.chrEBV.inverted.fa | instructions to generate            |

|_. First Header |_. Second Header     |
| Content Cell   | Content Cell Line 1
                   Content Cell Line 2 |
| Content Cell   | Content Cell        |

   - ~/Data/ARTDeco/preprocess_files                               # Run ARTDeco separately
   - ~/Data/SpliceWiz/NxtSE                                        # Run SpliceWiz.R to generate
   - ~/Data/SpliceWiz/output                                       # Run SpliceWiz.R to generate
   - ~/Data/SpliceWiz/ref                                          # Run SpliceWiz.R to generate
   - ~/Data/bams                                                   # Download from the NCBI GEO database under accession number GSE271717
   - ~/Data/refs/GRCh38.p14.ERCC.M81_DFLR.M81.chrEBV.inverted.fa   # see README in ~/Data/refs for how to generate

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Authour

Alejandro Casco, PhD Candidate.
For any questions or further information, please contact Alejandro Casco casco@wisc.edu or Dr. Eric Johannsen ejohannsen@medicine.wisc.edu.


Please cite accordingly if use this code in your work, including the original citations for any programs used.
