## -----------------------------------------------------------------------------
## Title: Analysis of the relationships between carcass attributes and beetle breeding outcomes 
##        in different taxonomic groups 
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

# 2. Clutch size vs. carcass weight and carcass taxon --------------------------
### Plot
plot_relationship_taxon(clutch_size)  # a quadratic relationship seems to exist

### Model (without reptiles)
# (1) test quadratic term
clutch_size_poisson_linear_taxon <- glmmTMB(clutch_size ~ carcass_weight * carcass_taxon + male_size + female_size + parent_generation,
                                            data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                            family = "poisson",
                                            na.action = na.omit)

clutch_size_poisson_quadratic_taxon <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) * carcass_taxon + male_size + female_size + parent_generation,
                                               data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                               family = "poisson",
                                               na.action = na.omit)

lrtest(clutch_size_poisson_linear_taxon, clutch_size_poisson_quadratic_taxon)  # quadratic term is significant
AIC(clutch_size_poisson_linear_taxon, clutch_size_poisson_quadratic_taxon)  # quadratic model is better

# (2) test overdispersion
clutch_size_nb_quadratic_taxon <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) * carcass_taxon + male_size + female_size + parent_generation,
                                          data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                          family = "nbinom2",
                                          na.action = na.omit)

lrtest(clutch_size_poisson_quadratic_taxon, clutch_size_nb_quadratic_taxon)  # overdispersion is significant
AIC(clutch_size_poisson_quadratic_taxon, clutch_size_nb_quadratic_taxon)  # negative binomial model is better

# (3) test zero inflation
clutch_size_zi_nb_quadratic_taxon <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) * carcass_taxon + male_size + female_size + parent_generation,
                                             data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                             ziformula = ~ 1,
                                             family = "nbinom2",
                                             na.action = na.omit)

testZeroInflation(clutch_size_nb_quadratic_taxon)
lrtest(clutch_size_nb_quadratic_taxon, clutch_size_zi_nb_quadratic_taxon)  # zero inflation is significant
AIC(clutch_size_nb_quadratic_taxon, clutch_size_zi_nb_quadratic_taxon)  # zero-inflated model is significant

# (4) test interaction term
clutch_size_zi_nb_quadratic_wo_interaction_taxon <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) + carcass_taxon + male_size + female_size + parent_generation,
                                                            data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                                            ziformula = ~ 1,
                                                            family = "nbinom2",
                                                            na.action = na.omit)

lrtest(clutch_size_zi_nb_quadratic_taxon, clutch_size_zi_nb_quadratic_wo_interaction_taxon)  # interaction is not significant
AIC(clutch_size_zi_nb_quadratic_taxon, clutch_size_zi_nb_quadratic_wo_interaction_taxon)  # model without interaction is better

# (5) model diagnostics
plot(simulateResiduals(clutch_size_zi_nb_quadratic_taxon))  # some patterns of heteroscedasticity
check_model(clutch_size_zi_nb_quadratic_taxon)  # some patterns of heteroscedasticity

# (6) model significance
clutch_size_zi_nb_quadratic_null_taxon <- glmmTMB(clutch_size ~ 1,
                                                  data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                                  ziformula = ~ 1,
                                                  family = "nbinom2",
                                                  na.action = na.omit)

lrtest(clutch_size_zi_nb_quadratic_null_taxon, clutch_size_zi_nb_quadratic_taxon)  # model is globally significant

# (7) model summary
summary(clutch_size_zi_nb_quadratic_taxon)
# tidy(clutch_size_zi_nb_quadratic_taxon) %>% view
model_summary(clutch_size_zi_nb_quadratic_taxon, model_name = "Clutch size", transform_estimate = "exp")
model_forest_plot(clutch_size_zi_nb_quadratic_taxon, model_name = "Clutch size", transform_estimate = "exp")
Anova(clutch_size_zi_nb_quadratic_taxon, type = 2)
# confint(profile(clutch_size_zi_nb_quadratic_taxon)) %>% view

# # (8) emmeans
# emmeans_carcass_type_clutch_size_taxon <- emmeans(clutch_size_zi_nb_quadratic_taxon, "carcass_taxon", type = "response")
# emmeans_parent_generation_clutch_size_taxon <- emmeans(clutch_size_zi_nb_quadratic_taxon, "parent_generation", type = "response")
# 
# pairs(regrid(emmeans_carcass_type_clutch_size_taxon))
# pairs(regrid(emmeans_parent_generation_clutch_size_taxon))
# 
# cld(emmeans_carcass_type_clutch_size_taxon, adjust = "Tukey", Letters = letters)
# cld(emmeans_parent_generation_clutch_size_taxon, adjust = "Tukey", Letters = letters)

# (9) model visualization
plot_model(clutch_size_zi_nb_quadratic_taxon, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_taxon"))

# (10) write the model results
write_rds(clutch_size_zi_nb_quadratic_taxon, "./03_Outputs/Data_Clean/clutch_size_zi_nb_quadratic_taxon.rds")


# 3. Breeding success vs. carcass weight and carcass type ----------------------
### Plot
plot_relationship_taxon(breeding_success)  # a quadratic relationship seems to exist

### Model (without reptiles)
# (1) Test quadratic term
breeding_success_logistic_linear_taxon <- glmmTMB(breeding_success ~ carcass_weight * carcass_taxon + male_size + female_size + parent_generation,
                                            data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                            family = "binomial",
                                            na.action = na.omit)

breeding_success_logistic_quadratic_taxon <- glmmTMB(breeding_success ~ poly(carcass_weight, 2) * carcass_taxon + male_size + female_size + parent_generation,
                                               data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                               family = "binomial",
                                               na.action = na.omit)

lrtest(breeding_success_logistic_linear_taxon, breeding_success_logistic_quadratic_taxon)  # quadratic term is significant
AIC(breeding_success_logistic_linear_taxon, breeding_success_logistic_quadratic_taxon)  # quadratic model is better

# (2) test interaction term
breeding_success_logistic_quadratic_wo_interaction_taxon <- glmmTMB(breeding_success ~ poly(carcass_weight, 2) + carcass_taxon + male_size + female_size + parent_generation,
                                                              data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                                              family = "binomial",
                                                              na.action = na.omit)

lrtest(breeding_success_logistic_quadratic_taxon, breeding_success_logistic_quadratic_wo_interaction_taxon)  # interaction is not significant
AIC(breeding_success_logistic_quadratic_taxon, breeding_success_logistic_quadratic_wo_interaction_taxon)  #  model with interaction is not better

# (3) model diagnostics
plot(simulateResiduals(breeding_success_logistic_quadratic_taxon))  # no obvious residual patterns
check_model(breeding_success_logistic_quadratic_taxon)  # residual patterns acceptable

# (4) model significance
breeding_success_logistic_quadratic_null_taxon <- glmmTMB(breeding_success ~ 1,
                                                    data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                                    family = "binomial",
                                                    na.action = na.omit)

# lrtest(breeding_success_logistic_quadratic_taxon, breeding_success_logistic_quadratic_null_taxon)  # non-comparable because of the difference in sample sizes

# (5) model summary
summary(breeding_success_logistic_quadratic_taxon)
# tidy(breeding_success_logistic_quadratic_taxon) %>% view
model_summary(breeding_success_logistic_quadratic_taxon, model_name = "Breeding success", transform_estimate = "exp")
model_forest_plot(breeding_success_logistic_quadratic_taxon, model_name = "Breeding success", transform_estimate = "exp")
Anova(breeding_success_logistic_quadratic_taxon, type = 2)
# confint(profile(breeding_success_logistic_quadratic_taxon)) %>% view

# (6) emmeans
# emmeans_carcass_type_breeding_success_taxon <- emmeans(breeding_success_logistic_quadratic_taxon, "carcass_taxon", type = "response")
# emmeans_parent_generation_breeding_success_taxon <- emmeans(breeding_success_logistic_quadratic_taxon, "parent_generation", type = "response")
# 
# pairs(regrid(emmeans_carcass_type_breeding_success_taxon))
# pairs(regrid(emmeans_parent_generation_breeding_success_taxon))
# 
# cld(emmeans_carcass_type_breeding_success_taxon, adjust = "Tukey", Letters = letters)
# cld(emmeans_parent_generation_breeding_success_taxon, adjust = "Tukey", Letters = letters)

# (7) model visualization
plot_model(breeding_success_logistic_quadratic_taxon, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_taxon"))

# (8) write the model results
write_rds(breeding_success_logistic_quadratic_taxon, "./03_Outputs/Data_Clean/breeding_success_logistic_quadratic_taxon.rds")














