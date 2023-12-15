## -----------------------------------------------------------------------------
## Title: Organize the burying beetle breeding experiment data 
##
## Author: Gen-Chang Hsu
##
## Date: 2023-12-14
##
## Description:
## 1. Organize the breeding experiment data
##
##
## -----------------------------------------------------------------------------
set.seed(123)


# Libraries --------------------------------------------------------------------
library(tidyverse)


# Import files -----------------------------------------------------------------
carcass_data_raw <- read_csv("./01_Data_Raw/Carcass_Data_20231214.csv")


############################### Code starts here ###############################

# 1. Organize the breeding experiment data -------------------------------------
carcass_data_clean <- carcass_data_raw %>% 
  select(date,
         carcass_sp_Chinese = sp,
         carcass_type = tr, 
         carcass_taxon = class,
         carcass_weight = carc_wt,
         parent_generation,
         pair_id = pair,
         male_size,
         female_size,
         clutch_size,
         n_larvae = larvae,
         total_larval_mass = tot_mass,
         carcass_weight_loss = carcass_used) %>% 
  mutate(average_larval_mass = if_else(n_larvae == 0, NA, total_larval_mass/n_larvae),
         larval_density = n_larvae/carcass_weight,
         efficiency = if_else(n_larvae == 0, NA, carcass_weight_loss/carcass_weight))

view(carcass_data_clean)



