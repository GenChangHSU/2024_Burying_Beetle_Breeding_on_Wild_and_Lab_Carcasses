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
library(car)
library(broom)
library(broom.mixed)
library(emmeans)


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

### Convert the variable "parent_generation" to a factor
carcass_data_clean <- carcass_data_clean %>% 
  mutate(parent_generation = as.factor(parent_generation))

# 1. Clutch size vs. carcass weight and carcass type ---------------------------
### Plot
plot_relationship(clutch_size)  # a quadratic relationship seems to exist

### Model
# (1) Test the quadratic term
clutch_size_poisson_linear <- glmmTMB(clutch_size ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                      data = carcass_data_clean,
                                      family = "poisson",
                                      na.action = na.omit)

clutch_size_poisson_quadratic <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                         data = carcass_data_clean,
                                         family = "poisson",
                                         na.action = na.omit)

lrtest(clutch_size_poisson_linear, clutch_size_poisson_quadratic)  # quadratic model is better
AIC(clutch_size_poisson_linear, clutch_size_poisson_quadratic)  # quadratic model is better

# (2) test overdispersion
clutch_size_poisson_quadratic <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                         data = carcass_data_clean,
                                         family = "poisson",
                                         na.action = na.omit)

clutch_size_nb_quadratic <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                    data = carcass_data_clean,
                                    family = "nbinom2",
                                    na.action = na.omit)

lrtest(clutch_size_poisson_quadratic, clutch_size_nb_quadratic)  # overdispersion is significant
AIC(clutch_size_poisson_quadratic, clutch_size_nb_quadratic)

# (3) test zero inflation
clutch_size_zi_nb_quadratic <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                       data = carcass_data_clean,
                                       ziformula = ~ 1,
                                       family = "nbinom2",
                                       na.action = na.omit)

testZeroInflation(clutch_size_nb_quadratic)
lrtest(clutch_size_nb_quadratic, clutch_size_zi_nb_quadratic)  # zero inflation is significant
AIC(clutch_size_nb_quadratic, clutch_size_zi_nb_quadratic)

# (4) test the interaction term
clutch_size_zi_nb_quadratic_wo_interaction <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                      data = carcass_data_clean,
                                                      ziformula = ~ 1,
                                                      family = "nbinom2",
                                                      na.action = na.omit)

lrtest(clutch_size_zi_nb_quadratic, clutch_size_zi_nb_quadratic_wo_interaction)  # the interaction is significant
AIC(clutch_size_zi_nb_quadratic, clutch_size_zi_nb_quadratic_wo_interaction)

# (5) model diagnostics
plot(simulateResiduals(clutch_size_zi_nb_quadratic))
check_model(clutch_size_zi_nb_quadratic)  # residual plot looks acceptable

# (6) model significance
clutch_size_zi_nb_quadratic_null <- glmmTMB(clutch_size ~ 1,
                                            data = carcass_data_clean,
                                            ziformula = ~ 1,
                                            family = "nbinom2",
                                            na.action = na.omit)

lrtest(clutch_size_zi_nb_quadratic_null, clutch_size_zi_nb_quadratic)  # model is globally significant

# (7) coefficient significance
summary(clutch_size_zi_nb_quadratic)
tidy(clutch_size_zi_nb_quadratic) %>% view
Anova(clutch_size_zi_nb_quadratic, type = 3)
confint(profile(clutch_size_zi_nb_quadratic)) %>% view

# (8) emmeans
emmeans_clutch_size <- emmeans(clutch_size_zi_nb_quadratic, "carcass_type", type = "response")
emmeans_parent_generation <- emmeans(clutch_size_zi_nb_quadratic, "parent_generation", type = "response")

pairs(regrid(emmeans_clutch_size))
pairs(regrid(emmeans_parent_generation))


# 2. Breeding success vs. carcass weight and carcass type ----------------------
### Plot
plot_relationship(breeding_success)  # a quadratic relationship seems to exist

### Model
# (1) Test the quadratic term
breeding_success_logistic_linear <- glmmTMB(breeding_success ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                            data = carcass_data_clean,
                                            family = "binomial",
                                            na.action = na.omit)

breeding_success_logistic_quadratic <- glmmTMB(breeding_success ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                               data = carcass_data_clean,
                                               family = "binomial",
                                               na.action = na.omit)  # there is a model convergence issue

lrtest(breeding_success_logistic_linear, breeding_success_logistic_quadratic)  # can't compare the two models
AIC(breeding_success_logistic_linear, breeding_success_logistic_quadratic)  # can't compare the two models

# (2) test the interaction term
breeding_success_logistic_quadratic_wo_interaction <- glmmTMB(breeding_success ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                              data = carcass_data_clean,
                                                              family = "binomial",
                                                              na.action = na.omit)

lrtest(breeding_success_logistic_quadratic, breeding_success_logistic_quadratic_wo_interaction)  # can't compare the two models
AIC(breeding_success_logistic_quadratic, breeding_success_logistic_quadratic_wo_interaction)  # can't compare the two models

# (3) model diagnostics
plot(simulateResiduals(breeding_success_logistic_quadratic))  # residual plot looks acceptable
check_model(breeding_success_logistic_quadratic)

# (4) model significance
breeding_success_logistic_quadratic_null <- glmmTMB(breeding_success ~ 1,
                                                    data = carcass_data_clean,
                                                    family = "binomial",
                                                    na.action = na.omit)

lrtest(breeding_success_logistic_quadratic, breeding_success_logistic_quadratic_null)  # can't compare the two models

# (5) coefficient significance
summary(breeding_success_logistic_quadratic)
tidy(breeding_success_logistic_quadratic) %>% view
Anova(breeding_success_logistic_quadratic, type = 3)
confint(profile(breeding_success_logistic_quadratic)) %>% view

# (6) emmeans
emmeans_breeding_success <- emmeans(breeding_success_logistic_quadratic, "carcass_type", type = "response")
emmeans_parent_generation <- emmeans(breeding_success_logistic_quadratic, "parent_generation", type = "response")

pairs(regrid(emmeans_breeding_success))
pairs(regrid(emmeans_parent_generation))


# 3. Proportion of eggs developed vs. carcass weight and carcass type ----------
### Plot
plot_relationship(prop_eggs_developed)  # there are some unreasonable values

### Convert the zeros to 0.001 and values larger than 1 to 0.999
carcass_data_clean_prop_eggs_developed <- carcass_data_clean %>% 
  mutate(prop_eggs_developed = case_when(prop_eggs_developed >= 1 ~ 0.999,
                                         prop_eggs_developed == 0 ~ 0.001,
                                         TRUE ~ prop_eggs_developed))

### Re-plot the modified data
ggplot(carcass_data_clean_prop_eggs_developed, aes(x = carcass_weight, y = prop_eggs_developed, color = carcass_type)) + 
  geom_point() + 
  geom_smooth(se = F) + 
  scale_color_brewer(palette = "Set1")  # a quadratic relationship seems to exist

### Model
# (1) Test the quadratic term
prop_eggs_developed_beta_linear <- glmmTMB(prop_eggs_developed ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                           data = carcass_data_clean_prop_eggs_developed,
                                           family = beta_family("logit"),
                                           na.action = na.omit)

prop_eggs_developed_beta_quadratic <- glmmTMB(prop_eggs_developed ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                              data = carcass_data_clean_prop_eggs_developed,
                                              family = beta_family("logit"),
                                              na.action = na.omit)

lrtest(prop_eggs_developed_beta_linear, prop_eggs_developed_beta_quadratic)  # quadratic model is not better
AIC(prop_eggs_developed_beta_linear, prop_eggs_developed_beta_quadratic)  # linear model is slightly better

# (2) test the interaction term
prop_eggs_developed_beta_linear_wo_interaction <- glmmTMB(prop_eggs_developed ~ carcass_weight + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                          data = carcass_data_clean_prop_eggs_developed,
                                                          family = beta_family("logit"),
                                                          na.action = na.omit)

lrtest(prop_eggs_developed_beta_linear, prop_eggs_developed_beta_linear_wo_interaction)  # interaction term near significant
AIC(prop_eggs_developed_beta_linear, prop_eggs_developed_beta_linear_wo_interaction) # model with interaction term is slightly better

# (3) model diagnostics
plot(simulateResiduals(prop_eggs_developed_beta_linear))  # residual plot looks fine
check_model(prop_eggs_developed_beta_linear)

# (4) model significance
prop_eggs_developed_beta_linear_null <- glmmTMB(prop_eggs_developed ~ 1,
                                                data = carcass_data_clean_prop_eggs_developed,
                                                family = beta_family("logit"),
                                                na.action = na.omit)

lrtest(prop_eggs_developed_beta_linear, prop_eggs_developed_beta_linear_null)  # model is marginally globally significant

# (5) coefficient significance
summary(prop_eggs_developed_beta_linear)
tidy(prop_eggs_developed_beta_linear) %>% view
Anova(prop_eggs_developed_beta_linear, type = 3)
confint(profile(prop_eggs_developed_beta_linear)) %>% view

# (6) emmeans
emmeans_prop_eggs_developed <- emmeans(prop_eggs_developed_beta_linear, "carcass_type", type = "response")
emmeans_parent_generation <- emmeans(prop_eggs_developed_beta_linear, "parent_generation", type = "response")

pairs(regrid(emmeans_prop_eggs_developed))
pairs(regrid(emmeans_parent_generation))



# 4. Number of larvae vs. carcass weight and carcass type ----------------------
### Plot
plot_relationship(n_larvae)  # a quadratic relationship seems to exist

### Model
# (1) Test the quadratic term
n_larvae_poisson_linear <- glmmTMB(n_larvae ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                   data = carcass_data_clean,
                                   family = "poisson",
                                   na.action = na.omit)

n_larvae_poisson_quadratic <- glmmTMB(n_larvae ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                      data = carcass_data_clean,
                                      family = "poisson",
                                      na.action = na.omit)

lrtest(n_larvae_poisson_linear, n_larvae_poisson_quadratic)  # the quadratic model is better
AIC(n_larvae_poisson_linear, n_larvae_poisson_quadratic)  # the quadratic model is better

# (2) test overdispersion
n_larvae_poisson_quadratic <- glmmTMB(n_larvae ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                      data = carcass_data_clean,
                                      family = "poisson",
                                      na.action = na.omit)

n_larvae_nb_quadratic <- glmmTMB(n_larvae ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                 data = carcass_data_clean,
                                 family = "nbinom2",
                                 na.action = na.omit)

lrtest(n_larvae_poisson_quadratic, n_larvae_nb_quadratic)  # overdispersion is significant
AIC(n_larvae_poisson_quadratic, n_larvae_nb_quadratic)  # negative binomial model is better

# (3) test zero inflation
n_larvae_zi_nb_quadratic <- glmmTMB(n_larvae ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                    data = carcass_data_clean,
                                    ziformula = ~ 1,
                                    family = "nbinom2",
                                    na.action = na.omit)

testZeroInflation(n_larvae_nb_quadratic)
lrtest(n_larvae_nb_quadratic, n_larvae_zi_nb_quadratic)  # zero inflation is significant
AIC(n_larvae_nb_quadratic, n_larvae_zi_nb_quadratic)  # zero-inflated model is better

# (4) test the interaction term
n_larvae_zi_nb_quadratic_wo_interaction <- glmmTMB(n_larvae ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                   data = carcass_data_clean,
                                                   ziformula = ~ 1,
                                                   family = "nbinom2",
                                                   na.action = na.omit)

lrtest(n_larvae_zi_nb_quadratic, n_larvae_zi_nb_quadratic_wo_interaction)  # the interaction term is marginally significant
AIC(n_larvae_zi_nb_quadratic, n_larvae_zi_nb_quadratic_wo_interaction)  # the model with the interaction term is slightly better

# (5) model diagnostics
plot(simulateResiduals(n_larvae_zi_nb_quadratic))
check_model(n_larvae_zi_nb_quadratic)  # there are some patterns in the residual plot

# (6) model significance
n_larvae_zi_nb_quadratic_null <- glmmTMB(n_larvae ~ 1,
                                         data = carcass_data_clean,
                                         ziformula = ~ 1,
                                         family = "nbinom2",
                                         na.action = na.omit)

lrtest(n_larvae_zi_nb_quadratic_null, n_larvae_zi_nb_quadratic)

# (7) coefficient significance
summary(n_larvae_zi_nb_quadratic)
tidy(n_larvae_zi_nb_quadratic) %>% view
Anova(n_larvae_zi_nb_quadratic, type = 3)
confint(profile(n_larvae_zi_nb_quadratic)) %>% view

# (8) emmeans
emmeans_n_larvae <- emmeans(n_larvae_zi_nb_quadratic, "carcass_type", type = "response")
emmeans_parent_generation <- emmeans(n_larvae_zi_nb_quadratic, "parent_generation", type = "response")

pairs(regrid(emmeans_n_larvae))
pairs(regrid(emmeans_parent_generation))


# 5. Total larval mass vs. carcass weight and carcass type ---------------------
### Plot
plot_relationship(total_larval_mass)  # a quadratic relationship seems to exist

### Model
# (1) Test the quadratic term
total_larval_mass_gaussian_linear <- glmmTMB(total_larval_mass ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                             data = carcass_data_clean,
                                             family = "gaussian",
                                             na.action = na.omit)

total_larval_mass_gaussian_quadratic <- glmmTMB(total_larval_mass ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                data = carcass_data_clean,
                                                family = "gaussian",
                                                na.action = na.omit)

lrtest(total_larval_mass_gaussian_linear, total_larval_mass_gaussian_quadratic)  # quadratic model is better
AIC(total_larval_mass_gaussian_linear, total_larval_mass_gaussian_quadratic)  # quadratic model is better

# (2) test the interaction term
total_larval_mass_gaussian_quadratic_wo_interaction <- glmmTMB(total_larval_mass ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                               data = carcass_data_clean,
                                                               family = "gaussian",
                                                               na.action = na.omit)

lrtest(total_larval_mass_gaussian_quadratic, total_larval_mass_gaussian_quadratic_wo_interaction)  # the interaction term is significant
AIC(total_larval_mass_gaussian_quadratic, total_larval_mass_gaussian_quadratic_wo_interaction)  # the model with the interaction term is better

# (3) model diagnostics
plot(simulateResiduals(total_larval_mass_gaussian_quadratic))
check_model(total_larval_mass_gaussian_quadratic)  # residual plot looks acceptable

# (4) model significance
total_larval_mass_gaussian_quadratic_null <- glmmTMB(total_larval_mass ~ 1,
                                                     data = carcass_data_clean,
                                                     family = "gaussian",
                                                     na.action = na.omit)

lrtest(total_larval_mass_gaussian_quadratic, total_larval_mass_gaussian_quadratic_null)

# (5) coefficient significance
summary(total_larval_mass_gaussian_quadratic)
tidy(total_larval_mass_gaussian_quadratic) %>% view
Anova(total_larval_mass_gaussian_quadratic, type = 3)
confint(profile(total_larval_mass_gaussian_quadratic)) %>% view

# (6) emmeans
emmeans_total_larval_mass <- emmeans(total_larval_mass_gaussian_quadratic, "carcass_type")
emmeans_parent_generation <- emmeans(total_larval_mass_gaussian_quadratic, "parent_generation")

pairs(regrid(emmeans_total_larval_mass))
pairs(regrid(emmeans_parent_generation))


# 6. Average larval mass vs. carcass weight and carcass type -------------------
### Plot
plot_relationship(average_larval_mass)

### Model
# (1) Test the quadratic term
average_larval_mass_gaussian_linear <- glmmTMB(average_larval_mass ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                               data = carcass_data_clean,
                                               family = "gaussian",
                                               na.action = na.omit)

average_larval_mass_gaussian_quadratic <- glmmTMB(average_larval_mass ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                  data = carcass_data_clean,
                                                  family = "gaussian",
                                                  na.action = na.omit)

lrtest(average_larval_mass_gaussian_linear, average_larval_mass_gaussian_quadratic)  # quadratic model is not significantly better
AIC(average_larval_mass_gaussian_linear, average_larval_mass_gaussian_quadratic)  # linear model is not better

# (2) test the interaction term
average_larval_mass_gaussian_linear_wo_interaction <- glmmTMB(average_larval_mass ~ carcass_weight + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                              data = carcass_data_clean,
                                                              family = "gaussian",
                                                              na.action = na.omit)

lrtest(average_larval_mass_gaussian_linear, average_larval_mass_gaussian_linear_wo_interaction)  # the interaction term is significant
AIC(average_larval_mass_gaussian_linear, average_larval_mass_gaussian_linear_wo_interaction)  # the model with the interaction term is better

# (3) model diagnostics
plot(simulateResiduals(average_larval_mass_gaussian_linear))
check_model(average_larval_mass_gaussian_linear)  # residual plot looks fine

# (4) model significance
average_larval_mass_gaussian_linear_null <- glmmTMB(average_larval_mass ~ 1,
                                                    data = carcass_data_clean,
                                                    family = "gaussian",
                                                    na.action = na.omit)

lrtest(average_larval_mass_gaussian_linear, average_larval_mass_gaussian_linear_null)

# (5) coefficient significance
summary(average_larval_mass_gaussian_linear)
tidy(average_larval_mass_gaussian_linear) %>% view
Anova(average_larval_mass_gaussian_linear, type = 3)
confint(profile(average_larval_mass_gaussian_linear)) %>% view

# (6) emmeans
emmeans_average_larval_mass <- emmeans(average_larval_mass_gaussian_linear, "carcass_type")
emmeans_parent_generation <- emmeans(average_larval_mass_gaussian_linear, "parent_generation")

pairs(regrid(emmeans_average_larval_mass))
pairs(regrid(emmeans_parent_generation))


# 7. Larval density vs. carcass weight and carcass type ------------------------
### Plot
plot_relationship(larval_density)

### Model
# (1) Test the quadratic term
larval_density_gaussian_linear <- glmmTMB(larval_density ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                          data = carcass_data_clean,
                                          family = "gaussian",
                                          na.action = na.omit)

larval_density_gaussian_quadratic <- glmmTMB(larval_density ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                             data = carcass_data_clean,
                                             family = "gaussian",
                                             na.action = na.omit)

lrtest(larval_density_gaussian_linear, larval_density_gaussian_quadratic)  # the quadratic model is not significantly better
AIC(larval_density_gaussian_linear, larval_density_gaussian_quadratic)  # the linear model is better

# (2) test the interaction term
larval_density_gaussian_linear_wo_interaction <- glmmTMB(larval_density ~ carcass_weight + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                         data = carcass_data_clean,
                                                         family = "gaussian",
                                                         na.action = na.omit)

lrtest(larval_density_gaussian_linear, larval_density_gaussian_linear_wo_interaction)  # the interaction term is marginally significant
AIC(larval_density_gaussian_linear, larval_density_gaussian_linear_wo_interaction)  # the model with interaction term is slightly better

# (3) model diagnostics
plot(simulateResiduals(larval_density_gaussian_linear))
check_model(larval_density_gaussian_linear)  # residual plot looks fine

# (4) model significance
larval_density_gaussian_linear_null <- glmmTMB(larval_density ~ 1,
                                               data = carcass_data_clean,
                                               family = "gaussian",
                                               na.action = na.omit)

lrtest(larval_density_gaussian_linear, larval_density_gaussian_linear_null)

# (5) coefficient significance
summary(larval_density_gaussian_linear)
tidy(larval_density_gaussian_linear) %>% view
Anova(larval_density_gaussian_linear, type = 3)
confint(profile(larval_density_gaussian_linear)) %>% view

# (7) emmeans
emmeans_larval_density <- emmeans(larval_density_gaussian_linear, "carcass_type")
emmeans_parent_generation <- emmeans(larval_density_gaussian_linear, "parent_generation")

pairs(regrid(emmeans_larval_density))
pairs(regrid(emmeans_parent_generation))


# 8. Carcass used vs. carcass weight and carcass type --------------------------
### Plot
plot_relationship(carcass_weight_loss)

### Model
# (1) Test the quadratic term (need to add parent generation to the model later)
carcass_weight_loss_gaussian_linear <- glmmTMB(carcass_weight_loss ~ carcass_weight * carcass_type + male_size + female_size + (1|generation_pair_id),
                                               data = carcass_data_clean,
                                               family = "gaussian",
                                               na.action = na.omit)

carcass_weight_loss_gaussian_quadratic <- glmmTMB(carcass_weight_loss ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + (1|generation_pair_id),
                                                  data = carcass_data_clean,
                                                  family = "gaussian",
                                                  na.action = na.omit)

lrtest(carcass_weight_loss_gaussian_linear, carcass_weight_loss_gaussian_quadratic)  # the quadratic model is not significantly better
AIC(carcass_weight_loss_gaussian_linear, carcass_weight_loss_gaussian_quadratic)  # the linear model is better

# (2) test the interaction term
carcass_weight_loss_gaussian_linear_wo_interaction <- glmmTMB(carcass_weight_loss ~ carcass_weight + carcass_type + male_size + female_size + (1|generation_pair_id),
                                                              data = carcass_data_clean,
                                                              family = "gaussian",
                                                              na.action = na.omit)

lrtest(carcass_weight_loss_gaussian_linear, carcass_weight_loss_gaussian_linear_wo_interaction)  # the interaction term is significant
AIC(carcass_weight_loss_gaussian_linear, carcass_weight_loss_gaussian_linear_wo_interaction)  # the model with interaction term is better

# (3) model diagnostics
plot(simulateResiduals(carcass_weight_loss_gaussian_linear))
check_model(carcass_weight_loss_gaussian_linear)

# (4) model significance
carcass_weight_loss_gaussian_linear_null <- glmmTMB(carcass_weight_loss ~ 1,
                                                    data = carcass_data_clean,
                                                    family = "gaussian",
                                                    na.action = na.omit)

lrtest(carcass_weight_loss_gaussian_linear, carcass_weight_loss_gaussian_linear_null)  # the model is globally significant

# (5) coefficient significance
summary(carcass_weight_loss_gaussian_linear)
tidy(carcass_weight_loss_gaussian_linear) %>% view
Anova(carcass_weight_loss_gaussian_linear, type = 3)
confint(profile(carcass_weight_loss_gaussian_linear)) %>% view

# (6) emmeans
emmeans_carcass_weight_loss <- emmeans(carcass_weight_loss_gaussian_linear, "carcass_type")
emmeans_parent_generation <- emmeans(carcass_weight_loss_gaussian_linear, "parent_generation")

pairs(regrid(emmeans_carcass_weight_loss))
pairs(regrid(emmeans_parent_generation))


# 9. Carcass use efficiency vs. carcass weight and carcass type ----------------
### Plot
plot_relationship(efficiency)  # there seem to be two outliers

### Remove the outliers
carcass_data_clean_efficiency <- carcass_data_clean %>% 
  filter(efficiency < 0.5)

### Model
# (1) Test the quadratic term (need to add parent generation as a fixed-effects variable later)
efficiency_beta_linear <- glmmTMB(efficiency ~ carcass_weight * carcass_type + male_size + female_size + (1|generation_pair_id),
                                          data = carcass_data_clean_efficiency,
                                          family = beta_family("logit"),
                                          na.action = na.omit)

efficiency_beta_quadratic <- glmmTMB(efficiency ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + (1|generation_pair_id),
                                             data = carcass_data_clean_efficiency,
                                             family = beta_family("logit"),
                                             na.action = na.omit)

lrtest(efficiency_beta_linear, efficiency_beta_quadratic)
AIC(efficiency_beta_linear, efficiency_beta_quadratic)

# (2) test the interaction term
efficiency_beta_linear_wo_interaction <- glmmTMB(efficiency ~ carcass_weight + carcass_type + male_size + female_size + (1|generation_pair_id),
                                                         data = carcass_data_clean_efficiency,
                                                         family = beta_family("logit"),
                                                         na.action = na.omit)

lrtest(efficiency_beta_linear, efficiency_beta_linear_wo_interaction)
AIC(efficiency_beta_linear, efficiency_beta_linear_wo_interaction)

# (3) model diagnostics
plot(simulateResiduals(efficiency_beta_linear))
check_model(efficiency_beta_linear)

# (4) model significance
efficiency_beta_linear_null <- glmmTMB(efficiency ~ 1,
                                               data = carcass_data_clean_efficiency,
                                               family = beta_family("logit"),
                                               na.action = na.omit)

lrtest(efficiency_beta_linear, efficiency_beta_linear_null)

# (5) refit the final model using REML
efficiency_beta_linear <- update(efficiency_beta_linear, REML = T)  # refit the

# (6) coefficient significance
summary(efficiency_beta_linear)
tidy(efficiency_beta_linear) %>% view
Anova(efficiency_beta_linear, type = 3)
confint(profile(efficiency_beta_linear)) %>% view

# (7) emmeans
emmeans_efficiency <- emmeans(efficiency_beta_linear, "carcass_type")
emmeans_parent_generation <- emmeans(efficiency_beta_linear, "parent_generation")

pairs(regrid(emmeans_efficiency))
pairs(regrid(emmeans_parent_generation))


# 10. Total larval mass vs. carcass weight loss --------------------------------
### Plot
ggplot(carcass_data_clean, aes(x = total_larval_mass, y = carcass_weight_loss, color = carcass_type)) + 
  geom_point() + 
  geom_smooth(se = F) + 
  scale_color_brewer(palette = "Set1")

### Model
# (1) test the interaction term (need to add parent generation to the model later)
total_larval_mass_carcass_weight_loss_gaussian_linear <- glmmTMB(total_larval_mass ~ carcass_weight_loss * carcass_type + (1|generation_pair_id),
                                                              data = carcass_data_clean,
                                                              family = "gaussian",
                                                              na.action = na.omit)

total_larval_mass_carcass_weight_loss_gaussian_linear_wo_interaction <- glmmTMB(total_larval_mass ~ carcass_weight_loss + carcass_type + male_size + female_size + (1|generation_pair_id),
                                                                             data = carcass_data_clean,
                                                                             family = "gaussian",
                                                                             na.action = na.omit)

lrtest(total_larval_mass_carcass_weight_loss_gaussian_linear, total_larval_mass_carcass_weight_loss_gaussian_linear_wo_interaction)
AIC(total_larval_mass_carcass_weight_loss_gaussian_linear, total_larval_mass_carcass_weight_loss_gaussian_linear_wo_interaction)

# (2) model diagnostics
plot(simulateResiduals(total_larval_mass_carcass_weight_loss_gaussian_linear))
check_model(total_larval_mass_carcass_weight_loss_gaussian_linear)

# (3) model significance
total_larval_mass_carcass_weight_loss_gaussian_linear_null <- glmmTMB(total_larval_mass ~ 1,
                                                                   data = carcass_data_clean,
                                                                   family = "gaussian",
                                                                   na.action = na.omit)

lrtest(total_larval_mass_carcass_weight_loss_gaussian_linear, total_larval_mass_carcass_weight_loss_gaussian_linear_null)

# (4) refit the final model using REML
total_larval_mass_carcass_weight_loss_gaussian_linear <- update(total_larval_mass_carcass_weight_loss_gaussian_linear, REML = T)  # refit the

# (5) coefficient significance
summary(total_larval_mass_carcass_weight_loss_gaussian_linear)
tidy(total_larval_mass_carcass_weight_loss_gaussian_linear) %>% view
Anova(total_larval_mass_carcass_weight_loss_gaussian_linear, type = 3)
confint(profile(total_larval_mass_carcass_weight_loss_gaussian_linear)) %>% view

# (6) emmeans
emmeans_carcass_type <- emmeans(total_larval_mass_carcass_weight_loss_gaussian_linear, "carcass_type")
emmeans_parent_generation <- emmeans(total_larval_mass_carcass_weight_loss_gaussian_linear, "parent_generation")

pairs(regrid(emmeans_carcass_type))
pairs(regrid(emmeans_parent_generation))


# 11. Average larval mass vs. larval density -----------------------------------
### Plot
ggplot(carcass_data_clean, aes(x = average_larval_mass, y = larval_density, color = carcass_type)) + 
  geom_point() + 
  geom_smooth(se = F) + 
  scale_color_brewer(palette = "Set1")

### Model
# (1) test the interaction term
average_larval_mass_larval_density_gaussian_linear <- glmmTMB(average_larval_mass ~ larval_density * carcass_type + male_size + female_size + (1|generation_pair_id),
                                               data = carcass_data_clean,
                                               family = "gaussian",
                                               na.action = na.omit)

average_larval_mass_larval_density_gaussian_linear_wo_interaction <- glmmTMB(average_larval_mass ~ larval_density + carcass_type + male_size + female_size + (1|generation_pair_id),
                                                              data = carcass_data_clean,
                                                              family = "gaussian",
                                                              na.action = na.omit)

lrtest(average_larval_mass_larval_density_gaussian_linear, average_larval_mass_larval_density_gaussian_linear_wo_interaction)
AIC(average_larval_mass_larval_density_gaussian_linear, average_larval_mass_larval_density_gaussian_linear_wo_interaction)

# (2) model diagnostics
plot(simulateResiduals(average_larval_mass_larval_density_gaussian_linear))
check_model(average_larval_mass_larval_density_gaussian_linear)

# (3) model significance
average_larval_mass_larval_density_gaussian_linear_null <- glmmTMB(average_larval_mass ~ 1,
                                                    data = carcass_data_clean,
                                                    family = "gaussian",
                                                    na.action = na.omit)

lrtest(average_larval_mass_larval_density_gaussian_linear, average_larval_mass_larval_density_gaussian_linear_null)

# (4) refit the final model using REML
average_larval_mass_larval_density_gaussian_linear <- update(average_larval_mass_larval_density_gaussian_linear, REML = T)  # refit the

# (5) coefficient significance
summary(average_larval_mass_larval_density_gaussian_linear)
tidy(average_larval_mass_larval_density_gaussian_linear) %>% view
Anova(average_larval_mass_larval_density_gaussian_linear, type = 3)
confint(profile(average_larval_mass_larval_density_gaussian_linear)) %>% view

# (6) emmeans
emmeans_carcass_type <- emmeans(average_larval_mass_larval_density_gaussian_linear, "carcass_type")
emmeans_parent_generation <- emmeans(average_larval_mass_larval_density_gaussian_linear, "parent_generation")

pairs(regrid(emmeans_carcass_type))
pairs(regrid(emmeans_parent_generation))



















