## -----------------------------------------------------------------------------
## Title: Analyze the relationships between carcass attributes and beetle breeding outcomes
##
## Author: Gen-Chang Hsu
##
## Date: 2023-12-15
##
## Description:
## 1. 
## 2. 
##
##
## -----------------------------------------------------------------------------
set.seed(123)


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(performance)
library(lmtest)


# Import files -----------------------------------------------------------------
carcass_data_clean <- read_csv("./03_Outputs/Data_Clean/Carcass_Data_Clean.csv")


############################### Code starts here ###############################

### A helper function for visualizing the bivariate relationship
plot_relationship <- function(yvar){
  ggplot(carcass_data_clean, aes(x = carcass_weight, y = {{yvar}}, color = carcass_type)) + 
    geom_point() + 
    geom_smooth(se = F) + 
    scale_color_brewer(palette = "Set1")
}

# 1. Carcass weight vs. clutch size by carcass type ----------------------------
### Plot
plot_relationship(clutch_size)  # a quadratic relationship seems to exist

### Model
# (1) Test the quadratic term
clutch_size_poisson_linear <- glmmTMB(clutch_size ~ carcass_weight * carcass_type + male_size + female_size + (1|generation_pair_id),
                                      data = carcass_data_clean,
                                      family = "poisson",
                                      na.action = na.omit)

clutch_size_poisson_quadratic <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + (1|generation_pair_id),
                                         data = carcass_data_clean,
                                         family = "poisson",
                                         na.action = na.omit)

lrtest(clutch_size_poisson_linear, clutch_size_poisson_quadratic)
AIC(clutch_size_poisson_linear, clutch_size_poisson_quadratic)

# (2) test overdispersion
clutch_size_poisson_quadratic <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + (1|generation_pair_id),
                                         data = carcass_data_clean,
                                         family = "poisson",
                                         na.action = na.omit)

clutch_size_nb_quadratic <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + (1|generation_pair_id),
                                    data = carcass_data_clean,
                                    family = "nbinom2",
                                    na.action = na.omit)

lrtest(clutch_size_poisson_quadratic, clutch_size_nb_quadratic)
AIC(clutch_size_poisson_quadratic, clutch_size_nb_quadratic)

# (3) test zero inflation
clutch_size_zi_nb_quadratic <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + (1|generation_pair_id),
                                       data = carcass_data_clean,
                                       ziformula = ~ 1,
                                       family = "nbinom2",
                                       na.action = na.omit)

testZeroInflation(clutch_size_nb_quadratic)
lrtest(clutch_size_nb_quadratic, clutch_size_zi_nb_quadratic)
AIC(clutch_size_nb_quadratic, clutch_size_zi_nb_quadratic)

# (4) test the interaction term
clutch_size_zi_nb_quadratic_wo_interaction <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + (1|generation_pair_id),
                                                      data = carcass_data_clean,
                                                      ziformula = ~ 1,
                                                      family = "nbinom2",
                                                      na.action = na.omit)

lrtest(clutch_size_zi_nb_quadratic, clutch_size_zi_nb_quadratic_wo_interaction)
AIC(clutch_size_zi_nb_quadratic, clutch_size_zi_nb_quadratic_wo_interaction)

# (5) model diagnostics
plot(simulateResiduals(clutch_size_zi_nb_quadratic))
check_model(clutch_size_zi_nb_quadratic)

# (6) model significance







