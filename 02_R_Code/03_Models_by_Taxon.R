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
library(car)
library(emmeans)
library(multcomp)
library(glmmTMB)
library(DHARMa)
library(performance)
library(lmtest)
library(sjPlot)
library(broom)
library(broom.mixed)


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

### Number of observations for each taxon
carcass_data_clean %>% 
  group_by(carcass_taxon) %>% 
  summarise(n = n())

### Carcass weight distribution for each taxon
ggplot(carcass_data_clean, aes(x = carcass_taxon, y = carcass_weight)) + 
  geom_boxplot() + 
  geom_point() + 
  theme_classic()

### Test whether the carcass weight differed among the three taxa
carcass_weight_taxon_lm <- lm(carcass_weight ~ carcass_taxon, data = carcass_data_clean)
summary(carcass_weight_taxon_lm)
Anova(carcass_weight_taxon_lm, type = 2)
carcass_weight_taxon_emmeans <- emmeans(carcass_weight_taxon_lm, "carcass_taxon")
pairs(regrid(carcass_weight_taxon_emmeans))
cld(carcass_weight_taxon_emmeans, adjust = "Tukey", Letters = letters)

### A helper function for visualizing the bivariate relationship
plot_relationship_taxon <- function(yvar){
  ggplot(carcass_data_clean, aes(x = carcass_weight, y = {{yvar}}, color = carcass_taxon)) + 
    geom_point() + 
    geom_smooth(se = F) + 
    scale_color_brewer(palette = "Set1")
}

# 1. Clutch size vs. carcass weight and carcass taxon --------------------------
### Plot
plot_relationship_taxon(clutch_size)  # a quadratic relationship seems to exist

### Model
# (1) test quadratic term
clutch_size_poisson_linear_taxon <- glmmTMB(clutch_size ~ carcass_weight * carcass_taxon + male_size + female_size + parent_generation,
                                      data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                      family = "poisson",
                                      na.action = na.omit)

clutch_size_poisson_quadratic_taxon <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) * carcass_taxon + male_size + female_size + parent_generation,
                                         data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                         family = "poisson",
                                         na.action = na.omit)

lrtest(clutch_size_poisson_linear_taxon, clutch_size_poisson_quadratic_taxon)













