## -----------------------------------------------------------------------------
## Title: Organize and explore the burying beetle breeding experiment data 
##
## Author: Gen-Chang Hsu
##
## Date: 2023-12-14
##
## Description:
## 1. Organize the breeding experiment data
## 2. Explore and summarize breeding experiment data
##
## -----------------------------------------------------------------------------
set.seed(123)


# Libraries --------------------------------------------------------------------
library(tidyverse)


# Import files -----------------------------------------------------------------
carcass_data_raw <- read_csv("./01_Data_Raw/Carcass_Data_20231214.csv")


############################### Code starts here ###############################

# 1. Data organization ---------------------------------------------------------
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
  mutate(breeding_success = if_else(n_larvae == 0, 0, 1),
         prop_eggs_developed = if_else(clutch_size == 0, NA, n_larvae/clutch_size),
         average_larval_mass = if_else(n_larvae == 0, NA, total_larval_mass/n_larvae),
         larval_density = if_else(n_larvae == 0, NA, n_larvae/carcass_weight),
         efficiency = if_else(n_larvae == 0, NA, carcass_weight_loss/carcass_weight)) %>% 
  mutate(date = if_else(parent_generation == 3, ymd(date), mdy(date)))

view(carcass_data_clean)

write_csv(carcass_data_clean, "./03_Outputs/Data_Clean/Carcass_Data_Clean.csv")


# 2. Data exploration and summary ----------------------------------------------
### A helper function for summarizing and visualizing the continuous variables
summary_visualize_fun <- function(var){
  wild <- carcass_data_clean %>% 
    filter(carcass_type == "wild") %>% 
    pull({{var}}) %>% 
    summary()
  
  lab <- carcass_data_clean %>% 
    filter(carcass_type == "lab") %>% 
    pull({{var}}) %>% 
    summary()
  
  p <- ggplot(carcass_data_clean) + 
    geom_histogram(aes(x = {{var}}, fill = carcass_type)) + 
    facet_wrap(~ carcass_type, nrow = 2) + 
    scale_fill_brewer(palette = "Set1")
  
  print(p)
  return(list(wild = wild, lab = lab))
}

### (1) Taxon compositions of the wild carcasses
carcass_data_clean %>% 
  filter(carcass_type == "wild") %>% 
  group_by(carcass_taxon) %>% 
  summarise(n = n()) %>% 
  mutate(prop = n/sum(n))

### (2) Weight of the wild and lab carcasses
summary_visualize_fun(var = carcass_weight)

### (3) Male size
summary_visualize_fun(var = male_size)

### (4) Female size
summary_visualize_fun(var = female_size)

### (5) Clutch size
summary_visualize_fun(var = clutch_size)  # quite a few zeros

### (6) Number of larvae
summary_visualize_fun(var = n_larvae)  # quite a few zeros

### (7) Total larval mass
summary_visualize_fun(var = total_larval_mass)

### (8) Carcass weight loss
summary_visualize_fun(var = carcass_weight_loss)  # a few extreme values

### (9) Breeding success










