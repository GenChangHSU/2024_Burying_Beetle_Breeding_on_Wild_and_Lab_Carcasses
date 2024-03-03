## -----------------------------------------------------------------------------
## Title: Analysis of the relationships between carcass attributes and beetle breeding outcomes 
##        in different taxonomic groups 
##
## Author: Gen-Chang Hsu
##
## Date: 2024-03-02
##
## Description:
## 1. Summarize and visualize the data by carcass taxon
## 2. Model the relationship between clutch size vs. carcass weight and carcass taxon
## 3. Model the relationship between breeding success vs. carcass weight and carcass taxon
## 4. Model the relationship between proportion of eggs developed vs. carcass weight and carcass taxon
## 5. Model the relationship between number of larvae vs. carcass weight and carcass taxon
## 6. Model the relationship between total larval mass vs. carcass weight and carcass taxon
## 7. Model the relationship between average larval mass vs. carcass weight and carcass taxon
## 8. Model the relationship between larval density vs. carcass weight and carcass taxon
## 9. Model the relationship between carcass weight loss vs. carcass weight and carcass taxon
## 10. Model the relationship between proportion of carcass used vs. carcass weight and carcass taxon
## 11. Model the relationship between average larval mass vs. larval density and carcass taxon
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


# 3. Breeding success vs. carcass weight and carcass taxon ---------------------
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


# 4. Proportion of eggs developed vs. carcass weight and carcass taxon ---------
### Plot
plot_relationship_taxon(prop_eggs_developed)  # there are some impossible values

### Convert the zeros to 0.001 and values larger than 1 to 0.999
carcass_data_clean_prop_eggs_developed <- carcass_data_clean %>% 
  mutate(prop_eggs_developed = case_when(prop_eggs_developed >= 1 ~ 0.999,
                                         prop_eggs_developed == 0 ~ 0.001,
                                         TRUE ~ prop_eggs_developed))

### Re-plot the modified data
ggplot(carcass_data_clean_prop_eggs_developed, aes(x = carcass_weight, y = prop_eggs_developed, color = carcass_taxon)) + 
  geom_point() + 
  geom_smooth(se = F) + 
  scale_color_brewer(palette = "Set1")  # a quadratic relationship seems to exist

### Model (without reptiles)
# (1) test quadratic term
prop_eggs_developed_beta_linear_taxon <- glmmTMB(prop_eggs_developed ~ carcass_weight * carcass_taxon + male_size + female_size + parent_generation,
                                           data = filter(carcass_data_clean_prop_eggs_developed, carcass_taxon != "reptiles"),
                                           family = beta_family("logit"),
                                           na.action = na.omit)

prop_eggs_developed_beta_quadratic_taxon <- glmmTMB(prop_eggs_developed ~ poly(carcass_weight, 2) * carcass_taxon + male_size + female_size + parent_generation,
                                              data = filter(carcass_data_clean_prop_eggs_developed, carcass_taxon != "reptiles"),
                                              family = beta_family("logit"),
                                              na.action = na.omit)

lrtest(prop_eggs_developed_beta_linear_taxon, prop_eggs_developed_beta_quadratic_taxon)  # quadratic model is better
AIC(prop_eggs_developed_beta_linear_taxon, prop_eggs_developed_beta_quadratic_taxon)  # quadratic model is better

# (2) test interaction term
prop_eggs_developed_beta_quadratic_wo_interaction_taxon <- glmmTMB(prop_eggs_developed ~ poly(carcass_weight, 2) + carcass_taxon + male_size + female_size + parent_generation,
                                                             data = filter(carcass_data_clean_prop_eggs_developed, carcass_taxon != "reptiles"),
                                                             family = beta_family("logit"),
                                                             na.action = na.omit)

lrtest(prop_eggs_developed_beta_quadratic_taxon, prop_eggs_developed_beta_quadratic_wo_interaction_taxon)  # interaction is not significant
AIC(prop_eggs_developed_beta_quadratic_taxon, prop_eggs_developed_beta_quadratic_wo_interaction_taxon)  # model without interaction is better

# (3) model diagnostics
plot(simulateResiduals(prop_eggs_developed_beta_quadratic_taxon))  # some residual patterns
check_model(prop_eggs_developed_beta_quadratic_taxon)

# (4) model significance
prop_eggs_developed_beta_quadratic_null_taxon <- glmmTMB(prop_eggs_developed ~ 1,
                                                   data = filter(carcass_data_clean_prop_eggs_developed, carcass_taxon != "reptiles"),
                                                   family = beta_family("logit"),
                                                   na.action = na.omit)

lrtest(prop_eggs_developed_beta_quadratic_taxon, prop_eggs_developed_beta_quadratic_null_taxon)  # model is globally significant

# (5) model summary
summary(prop_eggs_developed_beta_quadratic_taxon)
# tidy(prop_eggs_developed_beta_quadratic_taxon) %>% view
model_summary(prop_eggs_developed_beta_quadratic_taxon, model_name = "Proportion of eggs developed", transform_estimate = NULL)
model_forest_plot(prop_eggs_developed_beta_quadratic_taxon, model_name = "Proportion of eggs developed", transform_estimate = NULL)
Anova(prop_eggs_developed_beta_quadratic_taxon, type = 2)
# confint(profile(prop_eggs_developed_beta_quadratic_taxon)) %>% view

# (6) emmeans
# emmeans_carcass_type_prop_eggs_developed_taxon <- emmeans(prop_eggs_developed_beta_quadratic_taxon, "carcass_taxon", type = "response")
# emmeans_parent_generation_prop_eggs_developed_taxon <- emmeans(prop_eggs_developed_beta_quadratic_taxon, "parent_generation", type = "response")
# 
# pairs(regrid(emmeans_carcass_type_prop_eggs_developed_taxon))
# pairs(regrid(emmeans_parent_generation_prop_eggs_developed_taxon))
# 
# cld(emmeans_carcass_type_prop_eggs_developed_taxon, adjust = "Tukey", Letters = letters)
# cld(emmeans_parent_generation_prop_eggs_developed_taxon, adjust = "Tukey", Letters = letters)

# (7) model visualization
plot_model(prop_eggs_developed_beta_quadratic_taxon, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_taxon"))

# (8) write the model results
write_rds(prop_eggs_developed_beta_quadratic_taxon, "./03_Outputs/Data_Clean/prop_eggs_developed_beta_quadratic_taxon.rds")


# 5. Number of larvae vs. carcass weight and carcass taxon ---------------------
### Plot
plot_relationship_taxon(n_larvae)  # a quadratic relationship seems to exist

### Model (without reptiles)
# (1) test quadratic term
n_larvae_poisson_linear_taxon <- glmmTMB(n_larvae ~ carcass_weight * carcass_taxon + male_size + female_size + parent_generation,
                                   data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                   family = "poisson",
                                   na.action = na.omit)

n_larvae_poisson_quadratic_taxon <- glmmTMB(n_larvae ~ poly(carcass_weight, 2) * carcass_taxon + male_size + female_size + parent_generation,
                                      data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                      family = "poisson",
                                      na.action = na.omit)

lrtest(n_larvae_poisson_linear_taxon, n_larvae_poisson_quadratic_taxon)  # quadratic term is significant
AIC(n_larvae_poisson_linear_taxon, n_larvae_poisson_quadratic_taxon)  # quadratic model is better

# (2) test overdispersion
n_larvae_poisson_quadratic_taxon <- glmmTMB(n_larvae ~ poly(carcass_weight, 2) * carcass_taxon + male_size + female_size + parent_generation,
                                      data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                      family = "poisson",
                                      na.action = na.omit)

n_larvae_nb_quadratic_taxon <- glmmTMB(n_larvae ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                 data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                 family = "nbinom2",
                                 na.action = na.omit)

lrtest(n_larvae_poisson_quadratic_taxon, n_larvae_nb_quadratic_taxon)  # overdispersion is significant
AIC(n_larvae_poisson_quadratic_taxon, n_larvae_nb_quadratic_taxon)  # negative binomial model is better

# (3) test zero inflation
n_larvae_zi_nb_quadratic_taxon <- glmmTMB(n_larvae ~ poly(carcass_weight, 2) * carcass_taxon + male_size + female_size + parent_generation,
                                    data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                    ziformula = ~ 1,
                                    family = "nbinom2",
                                    na.action = na.omit)

testZeroInflation(n_larvae_nb_quadratic_taxon)
lrtest(n_larvae_nb_quadratic_taxon, n_larvae_zi_nb_quadratic_taxon)  # zero inflation is significant
AIC(n_larvae_nb_quadratic_taxon, n_larvae_zi_nb_quadratic_taxon)  # zero-inflated model is better

# (4) test interaction term
n_larvae_zi_nb_quadratic_wo_interaction_taxon <- glmmTMB(n_larvae ~ poly(carcass_weight, 2) + carcass_taxon + male_size + female_size + parent_generation,
                                                   data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                                   ziformula = ~ 1,
                                                   family = "nbinom2",
                                                   na.action = na.omit)

lrtest(n_larvae_zi_nb_quadratic_taxon, n_larvae_zi_nb_quadratic_wo_interaction_taxon)  # interaction is significant
AIC(n_larvae_zi_nb_quadratic_taxon, n_larvae_zi_nb_quadratic_wo_interaction_taxon)  # model with interaction is better

# (5) model diagnostics
plot(simulateResiduals(n_larvae_zi_nb_quadratic_taxon))  # some patterns of heteroscedasticity
check_model(n_larvae_zi_nb_quadratic_taxon)  # some patterns of heteroscedasticity

# (6) model significance
n_larvae_zi_nb_quadratic_null_taxon <- glmmTMB(n_larvae ~ 1,
                                         data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                         ziformula = ~ 1,
                                         family = "nbinom2",
                                         na.action = na.omit)

# lrtest(n_larvae_zi_nb_quadratic_null_taxon, n_larvae_zi_nb_quadratic_taxon)  # non-comparable because of the difference in sample sizes

# (7) model summary
summary(n_larvae_zi_nb_quadratic_taxon)
# tidy(n_larvae_zi_nb_quadratic_taxon) %>% view
model_summary(n_larvae_zi_nb_quadratic_taxon, model_name = "Number of larvae", transform_estimate = "exp")
model_forest_plot(n_larvae_zi_nb_quadratic_taxon, model_name = "Number of larvae", transform_estimate = "exp")
Anova(n_larvae_zi_nb_quadratic_taxon, type = 3)
# confint(profile(n_larvae_zi_nb_quadratic)) %>% view

# (8) emmeans
# emmeans_carcass_type_n_larvae_taxon <- emmeans(n_larvae_zi_nb_quadratic_taxon, "carcass_taxon", type = "response")
# emmeans_parent_generation_n_larvae_taxon <- emmeans(n_larvae_zi_nb_quadratic_taxon, "parent_generation", type = "response")
# 
# pairs(regrid(emmeans_carcass_type_n_larvae_taxon))
# pairs(regrid(emmeans_parent_generation_n_larvae_taxon))
# 
# cld(emmeans_carcass_type_n_larvae_taxon, adjust = "Tukey", Letters = letters)
# cld(emmeans_parent_generation_n_larvae_taxon, adjust = "Tukey", Letters = letters)

# (9) model visualization
plot_model(n_larvae_zi_nb_quadratic_taxon, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_taxon"))

# (10) write the model results
write_rds(n_larvae_zi_nb_quadratic_taxon, "./03_Outputs/Data_Clean/n_larvae_zi_nb_quadratic_taxon.rds")


# 6. Total larval mass vs. carcass weight and carcass taxon --------------------
### Plot
plot_relationship_taxon(total_larval_mass)  # a quadratic relationship seems to exist

### Model (without reptiles)
# (1) test quadratic term
total_larval_mass_gaussian_linear_taxon <- glmmTMB(total_larval_mass ~ carcass_weight * carcass_taxon + male_size + female_size + parent_generation,
                                             data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                             family = "gaussian",
                                             na.action = na.omit)

total_larval_mass_gaussian_quadratic_taxon <- glmmTMB(total_larval_mass ~ poly(carcass_weight, 2) * carcass_taxon + male_size + female_size + parent_generation,
                                                data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                                family = "gaussian",
                                                na.action = na.omit)

lrtest(total_larval_mass_gaussian_linear_taxon, total_larval_mass_gaussian_quadratic_taxon)  # quadratic term is significant
AIC(total_larval_mass_gaussian_linear_taxon, total_larval_mass_gaussian_quadratic_taxon)  # quadratic model is better

# (2) test interaction term
total_larval_mass_gaussian_quadratic_wo_interaction_taxon <- glmmTMB(total_larval_mass ~ poly(carcass_weight, 2) + carcass_taxon + male_size + female_size + parent_generation,
                                                               data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                                               family = "gaussian",
                                                               na.action = na.omit)

lrtest(total_larval_mass_gaussian_quadratic_taxon, total_larval_mass_gaussian_quadratic_wo_interaction_taxon)  # interaction is not significant
AIC(total_larval_mass_gaussian_quadratic_taxon, total_larval_mass_gaussian_quadratic_wo_interaction_taxon)  # model without interaction is better

# (3) model diagnostics
plot(simulateResiduals(total_larval_mass_gaussian_quadratic_taxon))  # some patterns of heteroscedasticity
check_model(total_larval_mass_gaussian_quadratic_taxon)  # some patterns of heteroscedasticity

# (4) model significance
total_larval_mass_gaussian_quadratic_null_taxon <- glmmTMB(total_larval_mass ~ 1,
                                                     data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                                     family = "gaussian",
                                                     na.action = na.omit)

lrtest(total_larval_mass_gaussian_quadratic_taxon, total_larval_mass_gaussian_quadratic_null_taxon)  # model is globally significant

# (5) model summary
summary(total_larval_mass_gaussian_quadratic_taxon)
# tidy(total_larval_mass_gaussian_quadratic_taxon) %>% view
model_summary(total_larval_mass_gaussian_quadratic_taxon, model_name = "Total larval mass", transform_estimate = NULL)
model_forest_plot(total_larval_mass_gaussian_quadratic_taxon, model_name = "Total larval mass", transform_estimate = NULL)
Anova(total_larval_mass_gaussian_quadratic_taxon, type = 2)
# confint(profile(total_larval_mass_gaussian_quadratic_taxon)) %>% view

# (6) emmeans
emmeans_carcass_type_total_larval_mass_taxon <- emmeans(total_larval_mass_gaussian_quadratic_taxon, "carcass_taxon")
emmeans_parent_generation_total_larval_mass_taxon <- emmeans(total_larval_mass_gaussian_quadratic_taxon, "parent_generation")

pairs(regrid(emmeans_carcass_type_total_larval_mass_taxon))
pairs(regrid(emmeans_parent_generation_total_larval_mass_taxon))

cld(emmeans_carcass_type_total_larval_mass_taxon, adjust = "Tukey", Letters = letters)
cld(emmeans_parent_generation_total_larval_mass_taxon, adjust = "Tukey", Letters = letters)

# (7) model visualization
plot_model(total_larval_mass_gaussian_quadratic_taxon, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_taxon"))

# (8) write the model results
write_rds(total_larval_mass_gaussian_quadratic_taxon, "./03_Outputs/Data_Clean/total_larval_mass_gaussian_quadratic_taxon.rds")


# 6. Average larval mass vs. carcass weight and carcass taxon ------------------
### Plot
plot_relationship_taxon(average_larval_mass)

### Model (without reptiles)
# (1) test quadratic term
average_larval_mass_gaussian_linear_taxon <- glmmTMB(average_larval_mass ~ carcass_weight * carcass_taxon + male_size + female_size + parent_generation,
                                               data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                               family = "gaussian",
                                               na.action = na.omit)

average_larval_mass_gaussian_quadratic_taxon <- glmmTMB(average_larval_mass ~ poly(carcass_weight, 2) * carcass_taxon + male_size + female_size + parent_generation,
                                                  data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                                  family = "gaussian",
                                                  na.action = na.omit)

lrtest(average_larval_mass_gaussian_linear_taxon, average_larval_mass_gaussian_quadratic_taxon)  # quadratic term is significantly better
AIC(average_larval_mass_gaussian_linear_taxon, average_larval_mass_gaussian_quadratic_taxon)  # quadratic model is better

# (2) test interaction term
average_larval_mass_gaussian_quadratic_wo_interaction_taxon <- glmmTMB(average_larval_mass ~ poly(carcass_weight, 2) + carcass_taxon + male_size + female_size + parent_generation,
                                                                 data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                                                 family = "gaussian",
                                                                 na.action = na.omit)

lrtest(average_larval_mass_gaussian_quadratic_taxon, average_larval_mass_gaussian_quadratic_wo_interaction_taxon)  # interaction is not significant
AIC(average_larval_mass_gaussian_quadratic_taxon, average_larval_mass_gaussian_quadratic_wo_interaction_taxon)  # model without interaction is slightly better

# (3) model diagnostics
plot(simulateResiduals(average_larval_mass_gaussian_quadratic_taxon))  # slight residual patterns
check_model(average_larval_mass_gaussian_quadratic_taxon)  # no obvious residual patterns

# (4) model significance
average_larval_mass_gaussian_quadratic_null_taxon <- glmmTMB(average_larval_mass ~ 1,
                                                       data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                                       family = "gaussian",
                                                       na.action = na.omit)

lrtest(average_larval_mass_gaussian_quadratic_taxon, average_larval_mass_gaussian_quadratic_null_taxon)  # model is globally significant

# (5) model summary
summary(average_larval_mass_gaussian_quadratic_taxon)
# tidy(average_larval_mass_gaussian_quadratic_taxon) %>% view
model_summary(average_larval_mass_gaussian_quadratic_taxon, model_name = "Average larval mass", transform_estimate = NULL)
model_forest_plot(average_larval_mass_gaussian_quadratic_taxon, model_name = "Average larval mass", transform_estimate = NULL)
Anova(average_larval_mass_gaussian_quadratic_taxon, type = 2)
# confint(profile(average_larval_mass_gaussian_quadratic_taxon)) %>% view

# (6) emmeans
emmeans_carcass_type_average_larval_mass_taxon <- emmeans(average_larval_mass_gaussian_quadratic_taxon, "carcass_taxon")
emmeans_parent_generation_average_larval_mass_taxon <- emmeans(average_larval_mass_gaussian_quadratic_taxon, "parent_generation")

pairs(regrid(emmeans_carcass_type_average_larval_mass_taxon))
pairs(regrid(emmeans_parent_generation_average_larval_mass_taxon))

cld(emmeans_carcass_type_average_larval_mass_taxon, adjust = "Tukey", Letters = letters)
cld(emmeans_parent_generation_average_larval_mass_taxon, adjust = "Tukey", Letters = letters)

# (7) model visualization
plot_model(average_larval_mass_gaussian_quadratic_taxon, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_taxon"))

# (8) write the model results
write_rds(average_larval_mass_gaussian_quadratic_taxon, "./03_Outputs/Data_Clean/average_larval_mass_gaussian_quadratic_taxon.rds")


# 8. Larval density vs. carcass weight and carcass taxon -----------------------
### Plot
plot_relationship_taxon(larval_density)

### Model (without reptiles)
# (1) test quadratic term
larval_density_gaussian_linear_taxon <- glmmTMB(larval_density ~ carcass_weight * carcass_taxon + male_size + female_size + parent_generation,
                                          data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                          family = "gaussian",
                                          na.action = na.omit)

larval_density_gaussian_quadratic_taxon <- glmmTMB(larval_density ~ poly(carcass_weight, 2) * carcass_taxon + male_size + female_size + parent_generation,
                                             data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                             family = "gaussian",
                                             na.action = na.omit)

lrtest(larval_density_gaussian_linear_taxon, larval_density_gaussian_quadratic_taxon)  # quadratic term is not significant
AIC(larval_density_gaussian_linear_taxon, larval_density_gaussian_quadratic_taxon)  # linear model is better

# (2) test interaction term
larval_density_gaussian_linear_wo_interaction_taxon <- glmmTMB(larval_density ~ carcass_weight + carcass_taxon + male_size + female_size + parent_generation,
                                                         data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                                         family = "gaussian",
                                                         na.action = na.omit)

lrtest(larval_density_gaussian_linear_taxon, larval_density_gaussian_linear_wo_interaction_taxon)  # interaction is not significant
AIC(larval_density_gaussian_linear_taxon, larval_density_gaussian_linear_wo_interaction_taxon)  # model without interaction is better

# (3) model diagnostics
plot(simulateResiduals(larval_density_gaussian_linear_taxon))  # residuals acceptable
check_model(larval_density_gaussian_linear_taxon)  # residuals acceptable

# (4) model significance
larval_density_gaussian_linear_null_taxon <- glmmTMB(larval_density ~ 1,
                                               data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                               family = "gaussian",
                                               na.action = na.omit)

# lrtest(larval_density_gaussian_linear_taxon, larval_density_gaussian_linear_null_taxon)  # non-comparable because of the difference in sample sizes

# (5) model summary
summary(larval_density_gaussian_linear_taxon)
# tidy(larval_density_gaussian_linear_taxon) %>% view
model_summary(larval_density_gaussian_linear_taxon, model_name = "Larval density", transform_estimate = NULL)
model_forest_plot(larval_density_gaussian_linear_taxon, model_name = "Larval density", transform_estimate = NULL)
Anova(larval_density_gaussian_linear_taxon, type = 2)
# confint(profile(larval_density_gaussian_linear_taxon)) %>% view

# (6) emmeans
emmeans_carcass_type_larval_density_taxon <- emmeans(larval_density_gaussian_linear_taxon, "carcass_taxon")
emmeans_parent_generation_larval_density_taxon <- emmeans(larval_density_gaussian_linear_taxon, "parent_generation")

pairs(regrid(emmeans_carcass_type_larval_density_taxon))
pairs(regrid(emmeans_parent_generation_larval_density_taxon))

cld(emmeans_carcass_type_larval_density_taxon, adjust = "Tukey", Letters = letters)
cld(emmeans_parent_generation_larval_density_taxon, adjust = "Tukey", Letters = letters)

# (7) model visualization
plot_model(larval_density_gaussian_linear_taxon, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_taxon"))

# (8) write the model results
write_rds(larval_density_gaussian_linear_taxon, "./03_Outputs/Data_Clean/larval_density_gaussian_linear_taxon.rds")


# 9. Carcass weight loss vs. carcass weight and carcass taxon ------------------
### Plot
plot_relationship_taxon(carcass_weight_loss)  # one impossible value and two outliers

### Exclude the impossible value, two outliers, and the observations without any larva
carcass_data_clean_carcass_weight_loss <- carcass_data_clean %>% 
  filter(carcass_weight_loss < 25 & carcass_weight_loss > 0) %>% 
  filter(breeding_success == 1)

### Replot the data
ggplot(carcass_data_clean_carcass_weight_loss, aes(x = carcass_weight, y = carcass_weight_loss, color = carcass_taxon)) + 
  geom_point() + 
  geom_smooth(se = F) + 
  scale_color_brewer(palette = "Set1")  # a quadratic relationship seems to exist

### Model (without reptiles)
# (1) test quadratic term
carcass_weight_loss_gaussian_linear_taxon <- glmmTMB(carcass_weight_loss ~ carcass_weight * carcass_taxon + male_size + female_size + parent_generation,
                                               data = filter(carcass_data_clean_carcass_weight_loss, carcass_taxon != "reptiles"),
                                               family = "gaussian",
                                               na.action = na.omit)

carcass_weight_loss_gaussian_quadratic_taxon <- glmmTMB(carcass_weight_loss ~ poly(carcass_weight, 2) * carcass_taxon + male_size + female_size + parent_generation,
                                                  data = filter(carcass_data_clean_carcass_weight_loss, carcass_taxon != "reptiles"),
                                                  family = "gaussian",
                                                  na.action = na.omit)

lrtest(carcass_weight_loss_gaussian_linear_taxon, carcass_weight_loss_gaussian_quadratic_taxon)  # quadratic term is significant
AIC(carcass_weight_loss_gaussian_linear_taxon, carcass_weight_loss_gaussian_quadratic_taxon)  # quadratic model is better

# (2) test interaction term
carcass_weight_loss_gaussian_quadratic_wo_interaction_taxon <- glmmTMB(carcass_weight_loss ~ poly(carcass_weight, 2) + carcass_taxon + male_size + female_size + parent_generation,
                                                                 data = filter(carcass_data_clean_carcass_weight_loss, carcass_taxon != "reptiles"),
                                                                 family = "gaussian",
                                                                 na.action = na.omit)

lrtest(carcass_weight_loss_gaussian_quadratic_taxon, carcass_weight_loss_gaussian_quadratic_wo_interaction_taxon)  # interaction term is not significant
AIC(carcass_weight_loss_gaussian_quadratic_taxon, carcass_weight_loss_gaussian_quadratic_wo_interaction_taxon)  # model without interaction is slightly better

# (3) model diagnostics
plot(simulateResiduals(carcass_weight_loss_gaussian_quadratic_taxon))  # residuals acceptable
check_model(carcass_weight_loss_gaussian_quadratic_taxon)  # residuals acceptable

# (4) model significance
carcass_weight_loss_gaussian_quadratic_null_taxon <- glmmTMB(carcass_weight_loss ~ 1,
                                                       data = filter(carcass_data_clean_carcass_weight_loss, carcass_taxon != "reptiles"),
                                                       family = "gaussian",
                                                       na.action = na.omit)

lrtest(carcass_weight_loss_gaussian_quadratic_taxon, carcass_weight_loss_gaussian_quadratic_null_taxon)  # model is globally significant

# (5) model summary
summary(carcass_weight_loss_gaussian_quadratic_taxon)
# tidy(carcass_weight_loss_gaussian_quadratic_taxon) %>% view
model_summary(carcass_weight_loss_gaussian_quadratic_taxon, model_name = "Carcass weight loss", transform_estimate = NULL)
model_forest_plot(carcass_weight_loss_gaussian_quadratic_taxon, model_name = "Carcass weight loss", transform_estimate = NULL)
Anova(carcass_weight_loss_gaussian_quadratic_taxon, type = 2)
# confint(profile(carcass_weight_loss_gaussian_quadratic_taxon)) %>% view

# (6) emmeans
emmeans_carcass_type_carcass_weight_loss_taxon <- emmeans(carcass_weight_loss_gaussian_linear_taxon, "carcass_taxon")
emmeans_parent_generation_carcass_weight_loss_taxon <- emmeans(carcass_weight_loss_gaussian_linear_taxon, "parent_generation")

pairs(regrid(emmeans_carcass_type_carcass_weight_loss_taxon))
pairs(regrid(emmeans_parent_generation_carcass_weight_loss_taxon))

cld(emmeans_carcass_type_carcass_weight_loss_taxon, adjust = "Tukey", Letters = letters)
cld(emmeans_parent_generation_carcass_weight_loss_taxon, adjust = "Tukey", Letters = letters)

# (7) model visualization
plot_model(carcass_weight_loss_gaussian_quadratic_taxon, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_taxon"))

# (8) write the model results
write_rds(carcass_weight_loss_gaussian_quadratic_taxon, "./03_Outputs/Data_Clean/carcass_weight_loss_gaussian_quadratic_taxon.rds")


# 9. Proportion of carcass used vs. carcass weight and carcass taxon -----------
### Plot
plot_relationship_taxon(prop_carcass_used)  # one impossible value and two outliers

### Remove the impossible values and two outliers
carcass_data_clean_prop_carcass_used <- carcass_data_clean %>% 
  filter(prop_carcass_used < 0.7 & prop_carcass_used > 0)

### Replot the data
ggplot(carcass_data_clean_prop_carcass_used, aes(x = carcass_weight, y = prop_carcass_used, color = carcass_taxon)) + 
  geom_point() + 
  geom_smooth(se = F) + 
  scale_color_brewer(palette = "Set1") 

### Model (without reptiles)
# (1) test quadratic term
prop_carcass_used_beta_linear_taxon <- glmmTMB(prop_carcass_used ~ carcass_weight * carcass_taxon + male_size + female_size + parent_generation,
                                         data = filter(carcass_data_clean_prop_carcass_used, carcass_taxon != "reptiles"),
                                         family = beta_family("logit"),
                                         na.action = na.omit)

prop_carcass_used_beta_quadratic_taxon <- glmmTMB(prop_carcass_used ~ poly(carcass_weight, 2) * carcass_taxon + male_size + female_size + parent_generation,
                                            data = filter(carcass_data_clean_prop_carcass_used, carcass_taxon != "reptiles"),
                                            family = beta_family("logit"),
                                            na.action = na.omit)

lrtest(prop_carcass_used_beta_linear_taxon, prop_carcass_used_beta_quadratic_taxon)  # quadratic term not significant
AIC(prop_carcass_used_beta_linear_taxon, prop_carcass_used_beta_quadratic_taxon)  # linear model is slightly better

# (2) test interaction term
prop_carcass_used_beta_linear_wo_interaction_taxon <- glmmTMB(prop_carcass_used ~ carcass_weight + carcass_taxon + male_size + female_size + parent_generation,
                                                        data = filter(carcass_data_clean_prop_carcass_used, carcass_taxon != "reptiles"),
                                                        family = beta_family("logit"),
                                                        na.action = na.omit)

lrtest(prop_carcass_used_beta_linear_taxon, prop_carcass_used_beta_linear_wo_interaction_taxon)  # interaction is not significant
AIC(prop_carcass_used_beta_linear_taxon, prop_carcass_used_beta_linear_wo_interaction_taxon)  # model without interaction is better

# (3) model diagnostics
plot(simulateResiduals(prop_carcass_used_beta_linear_taxon))  # no obvious residual patterns
check_model(prop_carcass_used_beta_linear_taxon)

# (4) model significance
prop_carcass_used_beta_linear_null_taxon <- glmmTMB(prop_carcass_used ~ 1,
                                              data = filter(carcass_data_clean_prop_carcass_used, carcass_taxon != "reptiles"),
                                              family = beta_family("logit"),
                                              na.action = na.omit)

lrtest(prop_carcass_used_beta_linear_taxon, prop_carcass_used_beta_linear_null_taxon)  # model is globally significant

# (5) model summary
summary(prop_carcass_used_beta_linear_taxon)
# tidy(prop_carcass_used_beta_linear_taxon) %>% view
model_summary(prop_carcass_used_beta_linear_taxon, model_name = "Proportion of carcass used", transform_estimate = "exp")
model_forest_plot(prop_carcass_used_beta_linear_taxon, model_name = "Proportion of carcass used", transform_estimate = "exp")
Anova(prop_carcass_used_beta_linear_taxon, type = 2)
# confint(profile(prop_carcass_used_beta_linear_taxon)) %>% view

# (6) emmeans
emmeans_carcass_type_prop_carcass_used_taxon <- emmeans(prop_carcass_used_beta_linear_taxon, "carcass_taxon")
emmeans_parent_generation_prop_carcass_used_taxon <- emmeans(prop_carcass_used_beta_linear_taxon, "parent_generation")

pairs(regrid(emmeans_carcass_type_prop_carcass_used_taxon))
pairs(regrid(emmeans_parent_generation_prop_carcass_used_taxon))

cld(emmeans_carcass_type_prop_carcass_used_taxon, adjust = "Tukey", Letters = letters)
cld(emmeans_parent_generation_prop_carcass_used_taxon, adjust = "Tukey", Letters = letters)

# (7) model visualization
plot_model(prop_carcass_used_beta_linear_taxon, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_taxon"))

# (8) write the model results
write_rds(prop_carcass_used_beta_linear_taxon, "./03_Outputs/Data_Clean/prop_carcass_used_beta_linear_taxon.rds")


# 10. Average larval mass vs. larval density -----------------------------------
### Plot
ggplot(carcass_data_clean, aes(x = larval_density, y = average_larval_mass, color = carcass_taxon)) + 
  geom_point() + 
  geom_smooth(se = F) + 
  scale_color_brewer(palette = "Set1")  # a linear relationship

### Model (without reptiles)
# (1) test interaction term
average_larval_mass_larval_density_gaussian_linear_taxon <- glmmTMB(average_larval_mass ~ larval_density * carcass_taxon + male_size + female_size + parent_generation,
                                                              data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                                              family = "gaussian",
                                                              na.action = na.omit)

average_larval_mass_larval_density_gaussian_linear_wo_interaction_taxon <- glmmTMB(average_larval_mass ~ larval_density + carcass_taxon + male_size + female_size + parent_generation,
                                                                             data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                                                             family = "gaussian",
                                                                             na.action = na.omit)

lrtest(average_larval_mass_larval_density_gaussian_linear_taxon, average_larval_mass_larval_density_gaussian_linear_wo_interaction_taxon)  # interaction is not significant
AIC(average_larval_mass_larval_density_gaussian_linear_taxon, average_larval_mass_larval_density_gaussian_linear_wo_interaction_taxon)  # model without interaction is slightly better

# (2) model diagnostics
plot(simulateResiduals(average_larval_mass_larval_density_gaussian_linear_taxon))  # no obvious residual patterns
check_model(average_larval_mass_larval_density_gaussian_linear_taxon)  # no obvious residual patterns

# (3) model significance
average_larval_mass_larval_density_gaussian_linear_null_taxon <- glmmTMB(average_larval_mass ~ 1,
                                                                   data = filter(carcass_data_clean, carcass_taxon != "reptiles"),
                                                                   family = "gaussian",
                                                                   na.action = na.omit)

lrtest(average_larval_mass_larval_density_gaussian_linear_taxon, average_larval_mass_larval_density_gaussian_linear_null_taxon)  # model is globally significant

# (4) model summary
summary(average_larval_mass_larval_density_gaussian_linear_taxon)
# tidy(average_larval_mass_larval_density_gaussian_linear_taxon) %>% view
model_summary(average_larval_mass_larval_density_gaussian_linear_taxon, model_name = "Average larval mass vs. Larval density", transform_estimate = NULL)
model_forest_plot(average_larval_mass_larval_density_gaussian_linear_taxon, model_name = "Average larval mass vs. Larval density", transform_estimate = NULL)
Anova(average_larval_mass_larval_density_gaussian_linear_taxon, type = 2)
# confint(profile(average_larval_mass_larval_density_gaussian_linear_taxon)) %>% view

# (5) emmeans
emmeans_carcass_type_average_larval_mass_larval_density_taxon <- emmeans(average_larval_mass_larval_density_gaussian_linear_taxon, "carcass_taxon")
emmeans_parent_generation_average_larval_mass_larval_density_taxon <- emmeans(average_larval_mass_larval_density_gaussian_linear_taxon, "parent_generation")

pairs(regrid(emmeans_carcass_type_average_larval_mass_larval_density_taxon))
pairs(regrid(emmeans_parent_generation_average_larval_mass_larval_density_taxon))

cld(emmeans_carcass_type_average_larval_mass_larval_density_taxon, adjust = "Tukey", Letters = letters)
cld(emmeans_parent_generation_average_larval_mass_larval_density_taxon, adjust = "Tukey", Letters = letters)

# (6) model visualization
plot_model(average_larval_mass_larval_density_gaussian_linear_taxon, 
           type = "pred", 
           terms = c("larval_density [0:2]", "carcass_taxon"))

# (7) write the model results
write_rds(average_larval_mass_larval_density_gaussian_linear_taxon, "./03_Outputs/Data_Clean/average_larval_mass_larval_density_gaussian_linear_taxon.rds")




















