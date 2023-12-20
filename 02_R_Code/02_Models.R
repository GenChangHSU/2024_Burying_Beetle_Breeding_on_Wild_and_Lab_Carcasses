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

lrtest(clutch_size_poisson_linear, clutch_size_poisson_quadratic)
AIC(clutch_size_poisson_linear, clutch_size_poisson_quadratic)

# (2) test overdispersion
clutch_size_poisson_quadratic <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                         data = carcass_data_clean,
                                         family = "poisson",
                                         na.action = na.omit)

clutch_size_nb_quadratic <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                    data = carcass_data_clean,
                                    family = "nbinom2",
                                    na.action = na.omit)

lrtest(clutch_size_poisson_quadratic, clutch_size_nb_quadratic)
AIC(clutch_size_poisson_quadratic, clutch_size_nb_quadratic)

# (3) test zero inflation
clutch_size_zi_nb_quadratic <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                       data = carcass_data_clean,
                                       ziformula = ~ 1,
                                       family = "nbinom2",
                                       na.action = na.omit)

testZeroInflation(clutch_size_nb_quadratic)
lrtest(clutch_size_nb_quadratic, clutch_size_zi_nb_quadratic)
AIC(clutch_size_nb_quadratic, clutch_size_zi_nb_quadratic)

# (4) test the interaction term
clutch_size_zi_nb_quadratic_wo_interaction <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
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
clutch_size_zi_nb_quadratic_null <- glmmTMB(clutch_size ~ 1,
                                            data = carcass_data_clean,
                                            ziformula = ~ 1,
                                            family = "nbinom2",
                                            na.action = na.omit)

lrtest(clutch_size_zi_nb_quadratic_null, clutch_size_zi_nb_quadratic)

# (7) refit the final model using REML
clutch_size_zi_nb_quadratic <- update(clutch_size_zi_nb_quadratic, REML = T)  # refit the

# (8) coefficient significance
summary(clutch_size_zi_nb_quadratic)
tidy(clutch_size_zi_nb_quadratic) %>% view
Anova(clutch_size_zi_nb_quadratic, type = 3)
confint(profile(clutch_size_zi_nb_quadratic)) %>% view

# (8) emmeans
emmeans_clutch_size <- emmeans(clutch_size_zi_nb_quadratic, "carcass_type", type = "response")
emmeans_parent_generation <- emmeans(clutch_size_zi_nb_quadratic, "parent_generation", type = "response")

pairs(regrid(emmeans_clutch_size))
pairs(regrid(emmeans_parent_generation))


# 2. Number of larvae vs. carcass weight and carcass type ----------------------
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

lrtest(n_larvae_poisson_linear, n_larvae_poisson_quadratic)
AIC(n_larvae_poisson_linear, n_larvae_poisson_quadratic)

# (2) test overdispersion
n_larvae_poisson_quadratic <- glmmTMB(n_larvae ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                      data = carcass_data_clean,
                                      family = "poisson",
                                      na.action = na.omit)

n_larvae_nb_quadratic <- glmmTMB(n_larvae ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                 data = carcass_data_clean,
                                 family = "nbinom2",
                                 na.action = na.omit)

lrtest(n_larvae_poisson_quadratic, n_larvae_nb_quadratic)
AIC(n_larvae_poisson_quadratic, n_larvae_nb_quadratic)

# (3) test zero inflation
n_larvae_zi_nb_quadratic <- glmmTMB(n_larvae ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                    data = carcass_data_clean,
                                    ziformula = ~ 1,
                                    family = "nbinom2",
                                    na.action = na.omit)

testZeroInflation(n_larvae_nb_quadratic)
lrtest(n_larvae_nb_quadratic, n_larvae_zi_nb_quadratic)
AIC(n_larvae_nb_quadratic, n_larvae_zi_nb_quadratic)

# (4) test the interaction term
n_larvae_zi_nb_quadratic_wo_interaction <- glmmTMB(n_larvae ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                   data = carcass_data_clean,
                                                   ziformula = ~ 1,
                                                   family = "nbinom2",
                                                   na.action = na.omit)

lrtest(n_larvae_zi_nb_quadratic, n_larvae_zi_nb_quadratic_wo_interaction)
AIC(n_larvae_zi_nb_quadratic, n_larvae_zi_nb_quadratic_wo_interaction)

# (5) model diagnostics
plot(simulateResiduals(n_larvae_zi_nb_quadratic))
check_model(n_larvae_zi_nb_quadratic)

# (6) model significance
n_larvae_zi_nb_quadratic_null <- glmmTMB(n_larvae ~ 1,
                                         data = carcass_data_clean,
                                         ziformula = ~ 1,
                                         family = "nbinom2",
                                         na.action = na.omit)

lrtest(n_larvae_zi_nb_quadratic_null, n_larvae_zi_nb_quadratic)

# (7) refit the final model using REML
n_larvae_zi_nb_quadratic <- update(n_larvae_zi_nb_quadratic, REML = T)  # refit the

# (8) coefficient significance
summary(n_larvae_zi_nb_quadratic)
tidy(n_larvae_zi_nb_quadratic) %>% view
Anova(n_larvae_zi_nb_quadratic, type = 3)
confint(profile(n_larvae_zi_nb_quadratic)) %>% view

# (8) emmeans
emmeans_n_larvae <- emmeans(n_larvae_zi_nb_quadratic, "carcass_type", type = "response")
emmeans_parent_generation <- emmeans(n_larvae_zi_nb_quadratic, "parent_generation", type = "response")

pairs(regrid(emmeans_n_larvae))
pairs(regrid(emmeans_parent_generation))


# 3. Breeding success vs. carcass weight and carcass type ----------------------
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
                                      na.action = na.omit)

lrtest(breeding_success_logistic_linear, breeding_success_logistic_quadratic)
AIC(breeding_success_logistic_linear, breeding_success_logistic_quadratic)

# (2) test the interaction term
breeding_success_logistic_quadratic_wo_interaction <- glmmTMB(breeding_success ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                   data = carcass_data_clean,
                                                   family = "binomial",
                                                   na.action = na.omit)

lrtest(breeding_success_logistic_quadratic, breeding_success_logistic_quadratic_wo_interaction)
AIC(breeding_success_logistic_quadratic, breeding_success_logistic_quadratic_wo_interaction)

# (3) model diagnostics
plot(simulateResiduals(breeding_success_logistic_quadratic))
check_model(breeding_success_logistic_quadratic)

# (4) model significance
breeding_success_logistic_quadratic_null <- glmmTMB(breeding_success ~ 1,
                                         data = carcass_data_clean,
                                         family = "binomial",
                                         na.action = na.omit)

lrtest(breeding_success_logistic_quadratic, breeding_success_logistic_quadratic_null)

# (5) refit the final model using REML
breeding_success_logistic_quadratic <- update(breeding_success_logistic_quadratic, REML = T)  # refit the

# (6) coefficient significance
summary(breeding_success_logistic_quadratic)
tidy(breeding_success_logistic_quadratic) %>% view
Anova(breeding_success_logistic_quadratic, type = 3)
confint(profile(breeding_success_logistic_quadratic)) %>% view

# (7) emmeans
emmeans_breeding_success <- emmeans(breeding_success_logistic_quadratic, "carcass_type", type = "response")
emmeans_parent_generation <- emmeans(breeding_success_logistic_quadratic, "parent_generation", type = "response")

pairs(regrid(emmeans_breeding_success))
pairs(regrid(emmeans_parent_generation))


# 4. Proportion of eggs developed vs. carcass weight and carcass type ----------
### Plot
plot_relationship(prop_eggs_developed)  # a quadratic relationship seems to exist

### Convert the zeros to 0.001 and values larger than 1 to 0.999
carcass_data_clean_prop_eggs_developed <- carcass_data_clean %>% 
  mutate(prop_eggs_developed = case_when(prop_eggs_developed >= 1 ~ 0.999,
                                         prop_eggs_developed == 0 ~ 0.001,
                                         TRUE ~ prop_eggs_developed))

### Re-plot the modified data
ggplot(carcass_data_clean_prop_eggs_developed, aes(x = carcass_weight, y = prop_eggs_developed, color = carcass_type)) + 
  geom_point() + 
  geom_smooth(se = F, method = "glm") + 
  scale_color_brewer(palette = "Set1")

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

lrtest(prop_eggs_developed_beta_linear, prop_eggs_developed_beta_quadratic)
AIC(prop_eggs_developed_beta_linear, prop_eggs_developed_beta_quadratic)

# (2) test the interaction term
prop_eggs_developed_beta_linear_wo_interaction <- glmmTMB(prop_eggs_developed ~ carcass_weight + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                          data = carcass_data_clean_prop_eggs_developed,
                                                          family = beta_family("logit"),
                                                          na.action = na.omit)

lrtest(prop_eggs_developed_beta_linear, prop_eggs_developed_beta_linear_wo_interaction)
AIC(prop_eggs_developed_beta_linear, prop_eggs_developed_beta_linear_wo_interaction)

# (3) model diagnostics
plot(simulateResiduals(prop_eggs_developed_beta_linear))
check_model(prop_eggs_developed_beta_linear)

# (4) model significance
prop_eggs_developed_beta_linear_null <- glmmTMB(prop_eggs_developed ~ 1,
                                                data = carcass_data_clean_prop_eggs_developed,
                                                family = beta_family("logit"),
                                                na.action = na.omit)

lrtest(prop_eggs_developed_beta_linear, prop_eggs_developed_beta_linear_null)

# (5) refit the final model using REML
prop_eggs_developed_beta_linear <- update(prop_eggs_developed_beta_linear, REML = T)  # refit the

# (6) coefficient significance
summary(prop_eggs_developed_beta_linear)
tidy(prop_eggs_developed_beta_linear) %>% view
Anova(prop_eggs_developed_beta_linear, type = 3)
confint(profile(prop_eggs_developed_beta_linear)) %>% view

# (7) emmeans
emmeans_prop_eggs_developed <- emmeans(prop_eggs_developed_beta_linear, "carcass_type", type = "response")
emmeans_parent_generation <- emmeans(prop_eggs_developed_beta_linear, "parent_generation", type = "response")

pairs(regrid(emmeans_prop_eggs_developed))
pairs(regrid(emmeans_parent_generation))


# 5. Average larval mass vs. carcass weight and carcass type -------------------
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

lrtest(average_larval_mass_gaussian_linear, average_larval_mass_gaussian_quadratic)
AIC(average_larval_mass_gaussian_linear, average_larval_mass_gaussian_quadratic)

# (2) test the interaction term
average_larval_mass_gaussian_linear_wo_interaction <- glmmTMB(average_larval_mass ~ carcass_weight + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                              data = carcass_data_clean,
                                                              family = "gaussian",
                                                              na.action = na.omit)

lrtest(average_larval_mass_gaussian_linear, average_larval_mass_gaussian_linear_wo_interaction)
AIC(average_larval_mass_gaussian_linear, average_larval_mass_gaussian_linear_wo_interaction)

# (3) model diagnostics
plot(simulateResiduals(average_larval_mass_gaussian_linear))
check_model(average_larval_mass_gaussian_linear)

# (4) model significance
average_larval_mass_gaussian_linear_null <- glmmTMB(average_larval_mass ~ 1,
                                                    data = carcass_data_clean,
                                                    family = "gaussian",
                                                    na.action = na.omit)

lrtest(average_larval_mass_gaussian_linear, average_larval_mass_gaussian_linear_null)

# (5) refit the final model using REML
average_larval_mass_gaussian_linear <- update(average_larval_mass_gaussian_linear, REML = T)  # refit the

# (6) coefficient significance
summary(average_larval_mass_gaussian_linear)
tidy(average_larval_mass_gaussian_linear) %>% view
Anova(average_larval_mass_gaussian_linear, type = 3)
confint(profile(average_larval_mass_gaussian_linear)) %>% view

# (7) emmeans
emmeans_average_larval_mass <- emmeans(average_larval_mass_gaussian_linear, "carcass_type")
emmeans_parent_generation <- emmeans(average_larval_mass_gaussian_linear, "parent_generation")

pairs(regrid(emmeans_average_larval_mass))
pairs(regrid(emmeans_parent_generation))


# 6. Larval density vs. carcass weight and carcass type ------------------------
### Plot
plot_relationship(larval_density)

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

lrtest(average_larval_mass_gaussian_linear, average_larval_mass_gaussian_quadratic)
AIC(average_larval_mass_gaussian_linear, average_larval_mass_gaussian_quadratic)

# (2) test the interaction term
average_larval_mass_gaussian_linear_wo_interaction <- glmmTMB(average_larval_mass ~ carcass_weight + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                              data = carcass_data_clean,
                                                              family = "gaussian",
                                                              na.action = na.omit)

lrtest(average_larval_mass_gaussian_linear, average_larval_mass_gaussian_linear_wo_interaction)
AIC(average_larval_mass_gaussian_linear, average_larval_mass_gaussian_linear_wo_interaction)

# (3) model diagnostics
plot(simulateResiduals(average_larval_mass_gaussian_linear))
check_model(average_larval_mass_gaussian_linear)

# (4) model significance
average_larval_mass_gaussian_linear_null <- glmmTMB(average_larval_mass ~ 1,
                                                    data = carcass_data_clean,
                                                    family = "gaussian",
                                                    na.action = na.omit)

lrtest(average_larval_mass_gaussian_linear, average_larval_mass_gaussian_linear_null)

# (5) refit the final model using REML
average_larval_mass_gaussian_linear <- update(average_larval_mass_gaussian_linear, REML = T)  # refit the

# (6) coefficient significance
summary(average_larval_mass_gaussian_linear)
tidy(average_larval_mass_gaussian_linear) %>% view
Anova(average_larval_mass_gaussian_linear, type = 3)
confint(profile(average_larval_mass_gaussian_linear)) %>% view

# (7) emmeans
emmeans_average_larval_mass <- emmeans(average_larval_mass_gaussian_linear, "carcass_type")
emmeans_parent_generation <- emmeans(average_larval_mass_gaussian_linear, "parent_generation")

pairs(regrid(emmeans_average_larval_mass))
pairs(regrid(emmeans_parent_generation))



















