# Somatic Variant Calling Analysis for Acral Melanoma PDX Study - Figures

## Overview
This analysis pipeline processes somatic variant data from Patient-Derived Xenograft (PDX) samples and their matched patient tumors, focusing on mutation patterns and variant type comparisons

## Analysis Scripts

# `01_Figure_2a.R`
- **Purpose**: Generate oncoplots of mutational landscape
- **Key Features**:
  - Processes MAF files from tumor samples
  - Filters to one sample per patient
  - Visualizes top mutated COSMIC v97 cancer genes
  - Annotates with clinical metadata

# `02_Figure_4c.R`
- **Purpose**: Compare PDX and tumor samples mutation profiles
- **Key Features**:
  - Co-barplots of mutation frequencies
  - Variant type comparisons

## Project Structure
project/
├── all_samples/ # Analysis of combined PDX+tumor samples
├── analysis/ # Additional analysis outputs
├── COSMIC/ # Cancer Gene Census data (v97)
├── data/ # Input data files
│ └── METADATA_all_clinical_reviewed_dec_2024.tsv
├── environment_oncoplot.RData # Saved R environment
├── matched_samples/ # Matched PDX-tumor pairs analysis
├── patient_samples/ # Patient-only analyses
├── PDX/ # PDX sample data
│ └── all_samples/ # MAF and TMB files
├── pdx_cell_samples/ # PDX-cell line comparisons
├── PDX_Tumor_samples/ # Comparative analyses
├── scripts/ # Analysis scripts
│ ├── 01_Figure_2a.R # Oncoplot generation
│ └── 02_Figure_4c.R # PDX-tumor comparison
└── Tumour/ # Patient tumor data
└── all_samples/ # MAF and TMB files

## Setup & Reproducibility

### R Environment
We used `renv` for dependency management. To reproduce our environment:

```R
# Install renv if needed
if (!require("renv")) install.packages("renv")

# Restore the project environment
renv::restore()
```

**Environment Details:**
- R version: 4.4.1
- Dependencies: See `renv.lock` in each analysis folder

### File Path Management
We used the `here` package for robust path handling:
- Automatically detects project root
- Creates platform-independent paths
- Seamless RStudio integration

## Additional Notes
- All scripts are contained in the `scripts/` directory
- Each analysis folder contains its specific dependencies

## ✉️ Contact
For questions about this repository or analyses:  
**Annie Cristhine Moraes Sousa Squiavinato**  
📧 as81@sanger.ac.uk
