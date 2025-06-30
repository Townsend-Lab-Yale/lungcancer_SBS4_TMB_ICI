## lungcancer_SBS4_TMB_ICI
### Integrated LUAD Data
A comprehensive LUAD dataset of 9230 samples with whole-genome sequence data (WGS), whole-exome sequence data (WES), and targeted sequence data (TGS), including 1066 smokers and 656 never-smokers. 
### Code
All necessary codes for the manuscript are included here. 
#### System Requirements
All codes have been tested on MacOS with R version 4.3.0.
#### Installation Guide
Users should install the following packages to use the code, from an R terminal
```
install.packages(c('ggplot2', 'data.table', 'dplyr', 'rtracklayer', 'stringr', 'ggpubr', 'patchwork'))
remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR@v2.10.1", dependencies = T, force = T)
remotes::install_github("Townsend-Lab-Yale/ces.refset.hg19@*release", dependencies = T, force = T) 
```
#### Demo 
Expected output files can be found in R_data directory
#### Instructions for Use
First, run `01_CES_analysis.R` file to process all maf files, calculate the cancer effect size information for all samples, smokers, and never-smokers.
Then, run `02_model_plot.R` file to do the analysis and visualization of fig1 and fig2.
More detailed annotations can be found in the code files.
#### License
All codes are licensed under GNU General Public License v3.0





