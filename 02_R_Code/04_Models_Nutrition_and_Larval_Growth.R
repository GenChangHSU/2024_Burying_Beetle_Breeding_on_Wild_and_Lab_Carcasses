## -----------------------------------------------------------------------------
## Title: Analysis of carcass nutritional composition and larval growth 
##
## Author: Syuan-Jyun Sun and Gen-Chang Hsu
##
## Date: 2024-05-25
##
## Description:
## 1. Model the relationship between protein content vs. carcass type
## 2. Model the relationship between protein content vs. wild carcass taxon
## 3. Model the relationship between fat content vs. carcass type
## 4. Model the relationship between fat content vs. wild carcass taxon
## 5. Model the relationship between larval growth vs. carcass type
## 6. Model the relationship between larval growth vs. wild carcass taxon
## 7.1. Model the relationship between larval growth vs. protein content (lab and wild carcasses)
## 7.2. Model the relationship between larval growth vs. protein content (wild carcasses only)
## 8.1. Model the relationship between larval growth vs. fat content (lab and wild carcasses)
## 8.2. Model the relationship between larval growth vs. fat content (wild carcasses only)
##
## -----------------------------------------------------------------------------
set.seed(123)


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(car)
library(glmmTMB)
library(DHARMa)
library(performance)
library(lmtest)
library(sjPlot)
library(broom)
library(broom.mixed)
library(emmeans)
library(multcomp)


# Import files -----------------------------------------------------------------
nutrition_data_clean <- read_csv("./03_Outputs/Data_Clean/Nutrition_Data_Clean.csv")
larval_growth_data_clean <- read_csv("./03_Outputs/Data_Clean/Larval_Growth_Data_Clean.csv")


############################### Code starts here ###############################

# 1. Protein content of lab vs. wild carcasses ---------------------------------






# 2. Fat content of lab vs. wild carcasses ----------------------------

# 3. Protein content of wild carcass taxa ----------------------------

# 4. Fat content of wild carcass taxa ----------------------------

# 5. Larval growth on lab vs. wild carcasses ------------------------------

# 6. Larval growth on wild carcass taxa ------------------------------

# 7.1. Larval growth vs. protein content for lab and wild carcasses -----------------------

# 7.2. Larval growth vs. protein content for wild carcasses -----------------------

# 8.1. Larval growth vs. fat content for lab and wild carcasses -----------------------

# 8.2. Larval growth vs. fat content for wild carcasses -----------------------



