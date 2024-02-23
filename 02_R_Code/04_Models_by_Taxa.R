## -----------------------------------------------------------------------------
## Title: Analysis of the relationships between carcass attributes and beetle breeding outcomes 
##        in different taxonic groups 
##
## Author: Gen-Chang Hsu
##
## Date: 2024-02-22
##
## Description:
##
##
##
## -----------------------------------------------------------------------------
set.seed(123)


# Libraries --------------------------------------------------------------------
library(tidyverse)


# Import files -----------------------------------------------------------------
carcass_data_clean <- read_csv("./03_Outputs/Data_Clean/Carcass_Data_Clean.csv")


############################### Code starts here ###############################

### Convert the variable "parent_generation" to a factor
carcass_data_clean <- carcass_data_clean %>% 
  mutate(parent_generation = as.factor(parent_generation))

### Exclude wild carcasses larger than 100 grams
carcass_data_clean <- carcass_data_clean %>% 
  filter(carcass_weight <= 100)

# 1. Data summary and visualization --------------------------------------------

















