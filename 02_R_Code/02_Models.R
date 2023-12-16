## -----------------------------------------------------------------------------
## Title: Analyze the relationships between carcass attributes and beetle breeding outcomes
##
## Author: Gen-Chang Hsu
##
## Date: 2023-12-15
##
## Description:
## 1. Explore and summarize the data
## 2. 
##
##
## -----------------------------------------------------------------------------
set.seed(123)


# Libraries --------------------------------------------------------------------
library(tidyverse)


# Import files -----------------------------------------------------------------
carcass_data_clean <- read_csv("./03_Outputs/Data_Clean/Carcass_Data_Clean.csv")


############################### Code starts here ###############################

# 1. Data exploration and summary ----------------------------------------------

### (1) Taxon compositions of the wild carcasses
carcass_data_clean %>% 
  filter(carcass_type == "wild") %>% 
  group_by(carcass_taxon) %>% 
  summarise(n = n()) %>% 
  mutate(prop = n/sum(n))

### (2) Weight of the wild and lab carcasses
carcass_data_clean %>% 
  filter(carcass_type == "wild") %>% 
  pull(carcass_weight) %>% 
  summary()

carcass_data_clean %>% 
  filter(carcass_type == "lab") %>% 
  pull(carcass_weight) %>% 
  summary()

ggplot(carcass_data_clean) + 
  geom_histogram(aes(x = carcass_weight, fill = carcass_type)) + 
  facet_wrap(~ carcass_type, nrow = 2) + 
  scale_fill_brewer(palette = "Set1")

### (3) Parent size





