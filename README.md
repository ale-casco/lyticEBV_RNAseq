# lyticEBV_RNAseq

## Overview
This repository contains the scripts and workflows used in our study on bulk RNA-seq data from FAC-sorted latent, early lytic, and late lytic fractions of lymphoblastoid cell lines (LCLs) transformed by a dual-fluorescent lytic reporter Epstein-Barr virus (EBV).

Data processing was primarily conducted in an Ubuntu Linux environment. Additional processing and visualization were performed in R.

## Data Availability

The raw (genome-aligned BAM files) and processed (RSEM output) RNA-seq data used in this study are available through the NCBI GEO database under the accession number GSE271717. The full manuscript detailing our findings is available on Cell Reports (currently under revision).

## Repository Contents

- **`Linux_workflow`**: Contains the detailed Linux workflow and associated bash scripts employed for data processing within the Ubuntu Linux environment.
- **`Rscripts`**: Contains the R scripts used for data analysis and figure generation.
- **`initiate.R`**: This script sets up the R environment, launching the RStudio project **`lyticEBV_project.Rproj`** and activating the **`renv`** environment to ensure reproducibility of the R package versions used in the study. See **`renv.lock`** for a snapshot of the R environment.

## Getting Started

### Prerequisites

- **Ubuntu Linux:** The data processing scripts are designed for Ubuntu Linux (tested on 20.04 LTS). Bash is required.
- **R:** The R environment is managed with the **`renv`** package to handle dependencies.

Ensure you have Ubuntu Linnux installed on your system. Processing of the raw data and bash scripts have been tested on Ubuntu 20.04 LTS. Bash is required for executing the data processing scripts.

## Instructions

1. Run **`initiate.R`:**
   - This script will open the RStudio project **`lyticEBV_project.Rproj`**.
   - It will also activate the **`renv`** environment, restoring the necessary R packages.

2. **Run Analysis Scripts:**
   - Execute the data analysis or figure generation scripts as needed from the **`Rscripts`** directory.
   -  **Note:** Some large files are not included due to size limitations. Refer to the "Data Exclusions" section below.

## Data Exclusions

The following files/directories are referenced in the workflow but not included in this repository:
| Path                                                                                  | Notes                                                    |
| ------------------------------------------------------------------------------------- | -------------------------------------------------------- |
| `~/Data/SpliceWiz/NxtSE` ,<br> `~/Data/SpliceWiz/output` ,<br> `~/Data/SpliceWiz/ref` | <be>Generate by running `SpliceWiz.R`                    |
| `~/Data/bams`                                                                         | Download from NCBI GEO (accession GSE271717).            |
| `~/Data/refs/GRCh38.p14.ERCC.M81_DFLR.M81.chrEBV.inverted.fa`                         | See README in `~/Data/refs` for generation instructions. |
| `~/Data/ARTDeco/preprocess_files`                                                     | Generate by running ARTDeco. Not required for analysis.  |

Placeholder README files with directory structure and size information are provided for these exclusions.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions or further information about this project, please contact:
- **Alejandro Casco, PhD**: casco@wisc.edu
- **Eric Johannsen, MD**:   ejohannsen@medicine.wisc.edu
