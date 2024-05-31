## -----------------------------------------------------------------------------
## Title: Organize and explore the burying beetle experiment data
##
## Author: Gen-Chang Hsu
##
## Date: 2024-05-26
##
## Description:
## 1. Clean and organize the breeding experiment data
## 2. Explore and summarize the breeding experiment data
## 3. Clean and organize the carcass nutritional composition data
## 4. Explore and summarize the carcass nutritional composition data
## 5. Clean and organize the larval feeding experiments data
## 6. Explore and summarize the larval feeding experiments data
##
## -----------------------------------------------------------------------------
set.seed(123)


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(readxl)


# Import files -----------------------------------------------------------------
carcass_data_raw <- read_xls("./01_Data_Raw/Breeding_Data_All.xls", sheet = 1)
nutrition_data_raw <- read_xls("./01_Data_Raw/Nutrition_Data.xls", sheet = 1)
larval_growth_data_raw <- read_xls("./01_Data_Raw/Larval_Growth_Data.xls", sheet = 1)


############################### Code starts here ###############################

# 1. Breeding experiment data organization and cleaning ------------------------
carcass_data_clean <- carcass_data_raw %>% 
  dplyr::select(date,
         carcass_sp_Chinese = sp_chinese,
         carcass_type = tr, 
         carcass_taxon = class,
         carcass_weight = carc_wt,
         parent_generation,
         pair_id,
         male_size,
         female_size,
         clutch_size,
         n_larvae = larvae,
         total_larval_mass = tot_mass,
         carcass_weight_loss = carcass_used) %>% 
  mutate(generation_pair_id = str_c(parent_generation, pair_id, sep = "_"),
         breeding_success = if_else(n_larvae == 0, 0, 1),
         prop_eggs_developed = if_else(clutch_size == 0, NA, n_larvae/clutch_size),
         total_larval_mass = if_else(is.na(clutch_size), NA, total_larval_mass),
         average_larval_mass = if_else(n_larvae == 0, NA, total_larval_mass/n_larvae),
         larval_density = if_else(n_larvae == 0, NA, n_larvae/carcass_weight),
         prop_carcass_used = if_else(n_larvae == 0, NA, carcass_weight_loss/carcass_weight)) %>% 
  mutate(date = ymd(date)) %>%
  relocate(generation_pair_id, .after = pair_id)

write_csv(carcass_data_clean, "./03_Outputs/Data_Clean/Breeding_Data_Clean.csv")


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

### (1) Taxon composition of the wild carcasses
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
summary_visualize_fun(var = total_larval_mass)  # quite a few zeros

### (8) Carcass weight loss
summary_visualize_fun(var = carcass_weight_loss)  # a few extreme and impossible values

### (9) Breeding success
carcass_data_clean %>%
  group_by(carcass_type, breeding_success) %>% 
  summarise(n = n()) %>% 
  mutate(prop = n/sum(n))

### (10) Proportion of eggs developed
summary_visualize_fun(var = prop_eggs_developed)  # quite a few zeros and some impossible values

### (11) Average larval mass
summary_visualize_fun(var = average_larval_mass)

### (12) Larval density
summary_visualize_fun(var = larval_density)

### (13) Proportion of carcass used
summary_visualize_fun(var = prop_carcass_used)  # some impossible values


# 3. Nutritional composition data organization and cleaning --------------------
nutrition_data_clean <- nutrition_data_raw %>% 
  dplyr::select(block_id = bl,
                carcass_id = newid,
                carcass_source = carc.type,
                tissue_type = type,
                tissue_replication = rep,
                wet_mass_g = wetmass,
                water_mass_g = water,
                dry_mass_g = drymass,
                protein_mass_g = protein,
                fat_mass_g = oil,
                total_mass_g = sum) %>% 
  mutate(prop_protein = protein_mass_g/total_mass_g,
         prop_fat = fat_mass_g/total_mass_g) %>% 
  mutate(carcass_type = if_else(carcass_source == "lab", "lab", "wild"),
         carcass_taxon = if_else(carcass_source == "lab", "mouse", carcass_source),
         .after = carcass_source) %>% 
  mutate(carcass_taxon = case_when(carcass_taxon == "snake" ~ "reptile",
                                   carcass_taxon == "mouse" ~ "mammal",
                                   TRUE ~ carcass_taxon),
         tissue_type = case_when(tissue_type == "i" ~ "viscera",
                                 tissue_type == "m" ~ "muscle")) %>% 
  dplyr::select(-carcass_source)

write_csv(nutrition_data_clean, "./03_Outputs/Data_Clean/Nutrition_Data_Clean.csv")


# 4. Data exploration and summary ----------------------------------------------
### (1) Number of lab and wild carcasses analyzed
nutrition_data_clean %>% 
  group_by(carcass_type) %>% 
  summarise(n = length(unique(carcass_id)))
  
### (2) Number of carcasses for each wild taxon
nutrition_data_clean %>% 
  filter(carcass_type == "wild") %>% 
  group_by(carcass_taxon) %>% 
  summarise(n = length(unique(carcass_id)))

### (3) Mean protein and fat content of lab and wild carcasses
nutrition_data_clean %>% 
  group_by(carcass_type) %>%
  summarise(mean_prop_protein = mean(prop_protein),
            mean_prop_fat = mean(prop_fat))

### (4) Mean protein and fat content of wild carcass taxa
nutrition_data_clean %>%
  filter(carcass_type == "wild") %>% 
  group_by(carcass_taxon) %>%
  summarise(mean_prop_protein = mean(prop_protein),
            mean_prop_fat = mean(prop_fat))


# 5. Larval growth data organization and cleaning ------------------------------
larval_growth_data_clean <- larval_growth_data_raw %>% 
  dplyr::select(block_id = bl,
                carcass_id = id,
                carcass_type = tr,
                carcass_taxon = carc.type,
                tissue_type = type,
                larva_replication = rep,
                tissue_mass_g = mass,
                family_id = family,
                success = success,
                initial_larval_mass_g = mass_1st,
                end_larval_mass_g = mass_2nd,
                mean_prop_protein = avg_prop_protein,
                mean_prop_fat = avg_prop_oil) %>% 
  mutate(carcass_taxon = if_else(carcass_type == "lab", "mouse", carcass_taxon),
         larval_weight_gain_g = end_larval_mass_g - initial_larval_mass_g,
         .after = end_larval_mass_g) %>% 
  mutate(carcass_taxon = case_when(carcass_taxon == "snake" ~ "reptile",
                                   carcass_taxon == "mouse" ~ "mammal",
                                   TRUE ~ carcass_taxon),
         tissue_type = case_when(tissue_type == "i" ~ "viscera",
                                 tissue_type == "m" ~ "muscle"))
  
write_csv(larval_growth_data_clean, "./03_Outputs/Data_Clean/Larval_Growth_Data_Clean.csv")


# 6. Data exploration and summary ----------------------------------------------
### (1) Number of lab and wild carcasses analyzed
larval_growth_data_clean %>% 
  group_by(carcass_type, tissue_type) %>% 
  summarise(n = n())

### (2) Number of carcasses for each wild taxon
larval_growth_data_clean %>% 
  filter(carcass_type == "wild") %>% 
  group_by(carcass_taxon, tissue_type) %>% 
  summarise(n = n())

### (3) Larval survival rates and growth on lab and wild carcasses
larval_growth_data_clean %>% 
  group_by(carcass_type, tissue_type) %>% 
  summarise(n = n(),
            prop_survived = sum(success)/n,
            mean_larval_weight_gain_g = mean(larval_weight_gain_g, na.rm = T))

### (4) Larval survival rates and growth on wild carcass taxa
larval_growth_data_clean %>% 
  filter(carcass_type == "wild") %>% 
  group_by(carcass_taxon, tissue_type) %>% 
  summarise(n = n(),
            prop_survived = sum(success)/n,
            mean_larval_weight_gain_g = mean(larval_weight_gain_g, na.rm = T))


