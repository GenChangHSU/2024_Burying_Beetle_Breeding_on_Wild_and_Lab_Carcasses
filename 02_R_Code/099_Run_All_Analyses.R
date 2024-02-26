## -----------------------------------------------------------------------------
## Title: Run the entire analysis
##
## Author: Gen-Chang Hsu
##
## Date: 2024-02-17
##
## Description: Source and run all the R scripts in a logical order for complete 
##              analysis results
## 
## -----------------------------------------------------------------------------
source("./02_R_Code/01_Data_Cleaning.R")
source("./02_R_Code/02_Models.R")
source("./02_R_Code/03_Models_by_Taxon.R")
source("./02_R_Code/04_Figures.R")

