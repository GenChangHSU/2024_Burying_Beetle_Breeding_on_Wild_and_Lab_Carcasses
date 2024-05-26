## -----------------------------------------------------------------------------
## Title: Analysis of the relationships between carcass attributes and beetle breeding outcomes
##
## Author: Gen-Chang Hsu
##
## Date: 2024-05-03
##
## Description:
## 1. Model the relationship between clutch size vs. carcass weight and carcass type
## 2. Model the relationship between breeding success vs. carcass weight and carcass type
## 3. Model the relationship between proportion of eggs developed vs. carcass weight and carcass type
## 4. Model the relationship between number of larvae vs. carcass weight and carcass type
## 5. Model the relationship between total larval mass vs. carcass weight and carcass type
## 6. Model the relationship between average larval mass vs. carcass weight and carcass type
## 7. Model the relationship between larval density vs. carcass weight and carcass type
## 8. Model the relationship between total carcass use vs. carcass weight and carcass type
## 9. Model the relationship between proportion of carcass used vs. carcass weight and carcass type
## 10. Model the relationship between average larval mass vs. larval density and carcass type
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
library(sjPlot)
library(emmeans)
library(multcomp)


# Model summary and plot functions from the package "sjPlot" -------------------
model_summary <- function(model, model_name, transform_estimate) {
  tab_model(model,
            dv.labels = model_name,
            auto.label = T,
            show.est = T,
            show.se = T,
            show.ci = 0.95,
            show.stat = T,
            show.p = T,
            show.reflvl = F,
            col.order = c("est", "se", "ci", "stat", "p"),
            transform = transform_estimate,
            show.zeroinf = T,
            string.ci = "95% CI",
            string.se = "SE",
            string.stat = "Test statistic",
            string.p = "P")
}

model_forest_plot <- function(model, model_name, transform_estimate) {
  plot_model(model, 
             sort.est = T,
             transform = transform_estimate,
             show.values = T,
             show.p = T,
             value.offset = 0.3,
             vline.color = "red",
             title = model_name)
}

# Import files -----------------------------------------------------------------
carcass_data_clean <- read_csv("./03_Outputs/Data_Clean/Breeding_Data_Clean.csv")


############################### Code starts here ###############################

### A helper function for visualizing the bivariate relationships
plot_relationship <- function(yvar){
  ggplot(carcass_data_clean, aes(x = carcass_weight, y = {{yvar}}, color = carcass_type)) + 
    geom_point() + 
    geom_smooth(se = F) + 
    scale_color_brewer(palette = "Set1")
}

### Convert the variable "parent_generation" to a factor
carcass_data_clean <- carcass_data_clean %>% 
  mutate(parent_generation = as.factor(parent_generation))

### Exclude the carcass pairs with wild carcasses larger than 100 grams
large_wild_carcasses_id <- carcass_data_clean %>% 
  filter(carcass_weight > 100) %>% 
  pull(generation_pair_id)
  
carcass_data_clean <- carcass_data_clean %>% 
  filter(!generation_pair_id %in% large_wild_carcasses_id)
  
# 1. Clutch size vs. carcass weight and carcass type ---------------------------
### Plot
plot_relationship(clutch_size)  # a quadratic relationship seems to exist

### Model
# (1) test quadratic term
clutch_size_poisson_linear <- glmmTMB(clutch_size ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                      data = carcass_data_clean,
                                      family = "poisson",
                                      na.action = na.omit)

clutch_size_poisson_quadratic <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                         data = carcass_data_clean,
                                         family = "poisson",
                                         na.action = na.omit)

lrtest(clutch_size_poisson_linear, clutch_size_poisson_quadratic)  # quadratic term is significant
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
AIC(clutch_size_poisson_quadratic, clutch_size_nb_quadratic)  # negative binomial model is better

# (3) test zero inflation
clutch_size_zi_nb_quadratic <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                       data = carcass_data_clean,
                                       ziformula = ~ poly(carcass_weight, 2),
                                       family = "nbinom2",
                                       na.action = na.omit)

testZeroInflation(clutch_size_nb_quadratic)  # zero inflation is significant
lrtest(clutch_size_nb_quadratic, clutch_size_zi_nb_quadratic)  # zero inflation is significant
AIC(clutch_size_nb_quadratic, clutch_size_zi_nb_quadratic)  # zero-inflated model is significant

# (4) test interaction term
clutch_size_zi_nb_quadratic_wo_interaction <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                      data = carcass_data_clean,
                                                      ziformula = ~ poly(carcass_weight, 2),
                                                      family = "nbinom2",
                                                      na.action = na.omit)

lrtest(clutch_size_zi_nb_quadratic, clutch_size_zi_nb_quadratic_wo_interaction)  # interaction is not significant
AIC(clutch_size_zi_nb_quadratic, clutch_size_zi_nb_quadratic_wo_interaction)  # model without interaction is better

# (5) model diagnostics
plot(simulateResiduals(clutch_size_zi_nb_quadratic))  # acceptable
check_model(clutch_size_zi_nb_quadratic)  # acceptable

# (6) model significance
clutch_size_zi_nb_quadratic_null <- glmmTMB(clutch_size ~ 1,
                                            data = carcass_data_clean,
                                            ziformula = ~ 1,
                                            family = "nbinom2",
                                            na.action = na.omit)

lrtest(clutch_size_zi_nb_quadratic_null, clutch_size_zi_nb_quadratic)  # model is globally significant

# (7) model summary
summary(clutch_size_zi_nb_quadratic)
model_summary(clutch_size_zi_nb_quadratic, model_name = "Clutch size", transform_estimate = "exp")
model_forest_plot(clutch_size_zi_nb_quadratic, model_name = "Clutch size", transform_estimate = "exp")
Anova(clutch_size_zi_nb_quadratic, type = 2)
# confint(profile(clutch_size_zi_nb_quadratic)) %>% view

# (8) emmeans
emmeans_carcass_type_clutch_size <- emmeans(clutch_size_zi_nb_quadratic, "carcass_type", type = "response")
pairs(regrid(emmeans_carcass_type_clutch_size))
cld(emmeans_carcass_type_clutch_size, Letters = letters)

# (9) model visualization
plot_model(clutch_size_zi_nb_quadratic, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_type"))

# (10) write the model results
write_rds(clutch_size_zi_nb_quadratic, "./03_Outputs/Data_Clean/clutch_size_zi_nb_quadratic.rds")


# 2. Breeding success vs. carcass weight and carcass type ----------------------
### Plot
plot_relationship(breeding_success)  # a quadratic relationship seems to exist

### Model
# (1) test quadratic term
breeding_success_binomial_linear <- glmmTMB(breeding_success ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                            data = carcass_data_clean,
                                            family = "binomial",
                                            na.action = na.omit)

breeding_success_binomial_quadratic <- glmmTMB(breeding_success ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                               data = carcass_data_clean,
                                               family = "binomial",
                                               na.action = na.omit)

lrtest(breeding_success_binomial_linear, breeding_success_binomial_quadratic)  # quadratic term is significant
AIC(breeding_success_binomial_linear, breeding_success_binomial_quadratic)  # quadratic model is better

# (2) test overdispersion
breeding_success_betabinomial_quadratic <- glmmTMB(breeding_success ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                   data = carcass_data_clean,
                                                   family = "betabinomial",
                                                   na.action = na.omit)

lrtest(breeding_success_binomial_quadratic, breeding_success_betabinomial_quadratic)  # can't compare the two models  
AIC(breeding_success_binomial_quadratic, breeding_success_betabinomial_quadratic)  # AIC not available for the beta-binomial model

# (3) test interaction term
breeding_success_binomial_quadratic_wo_interaction <- glmmTMB(breeding_success ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                              data = carcass_data_clean,
                                                              family = "binomial",
                                                              na.action = na.omit)

lrtest(breeding_success_binomial_quadratic, breeding_success_binomial_quadratic_wo_interaction)  # interaction is not significant
AIC(breeding_success_binomial_quadratic, breeding_success_binomial_quadratic_wo_interaction)  #  model with interaction is not better

# (4) model diagnostics
plot(simulateResiduals(breeding_success_binomial_quadratic))  # slight pattern
check_model(breeding_success_binomial_quadratic)  # acceptable

# (5) model significance
breeding_success_binomial_quadratic_null <- glmmTMB(breeding_success ~ 1,
                                                    data = carcass_data_clean,
                                                    family = "binomial",
                                                    na.action = na.omit)

# lrtest(breeding_success_binomial_quadratic, breeding_success_binomial_quadratic_null)  # non-comparable because of the difference in sample sizes

# (6) model summary
summary(breeding_success_binomial_quadratic)
model_summary(breeding_success_binomial_quadratic, model_name = "Breeding success", transform_estimate = "exp")
model_forest_plot(breeding_success_binomial_quadratic, model_name = "Breeding success", transform_estimate = "exp")
Anova(breeding_success_binomial_quadratic, type = 2)
# confint(profile(breeding_success_binomial_quadratic)) %>% view

# (7) emmeans
emmeans_carcass_type_breeding_success <- emmeans(breeding_success_binomial_quadratic, "carcass_type", type = "response")
pairs(regrid(emmeans_carcass_type_breeding_success))
cld(emmeans_carcass_type_breeding_success, Letters = letters)

# (8) model visualization
plot_model(breeding_success_binomial_quadratic, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_type"))

# (9) write the model results
write_rds(breeding_success_binomial_quadratic, "./03_Outputs/Data_Clean/breeding_success_binomial_quadratic.rds")


# 3. Proportion of eggs developed vs. carcass weight and carcass type ----------
### Plot
plot_relationship(prop_eggs_developed)  # there are some impossible values

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
# (1) test quadratic term
prop_eggs_developed_beta_linear <- glmmTMB(prop_eggs_developed ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                           data = carcass_data_clean_prop_eggs_developed,
                                           family = beta_family("logit"),
                                           na.action = na.omit)

prop_eggs_developed_beta_quadratic <- glmmTMB(prop_eggs_developed ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                              data = carcass_data_clean_prop_eggs_developed,
                                              family = beta_family("logit"),
                                              na.action = na.omit)

lrtest(prop_eggs_developed_beta_linear, prop_eggs_developed_beta_quadratic)  # quadratic model is better
AIC(prop_eggs_developed_beta_linear, prop_eggs_developed_beta_quadratic)  # quadratic model is better

# (2) test interaction term
prop_eggs_developed_beta_quadratic_wo_interaction <- glmmTMB(prop_eggs_developed ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                             data = carcass_data_clean_prop_eggs_developed,
                                                             family = beta_family("logit"),
                                                             na.action = na.omit)

lrtest(prop_eggs_developed_beta_quadratic, prop_eggs_developed_beta_quadratic_wo_interaction)  # interaction is not significant
AIC(prop_eggs_developed_beta_quadratic, prop_eggs_developed_beta_quadratic_wo_interaction)  # model without interaction is better

# (3) model diagnostics
plot(simulateResiduals(prop_eggs_developed_beta_quadratic))  # some residual patterns
check_model(prop_eggs_developed_beta_quadratic)

# (4) model significance
prop_eggs_developed_beta_quadratic_null <- glmmTMB(prop_eggs_developed ~ 1,
                                                   data = carcass_data_clean_prop_eggs_developed,
                                                   family = beta_family("logit"),
                                                   na.action = na.omit)

lrtest(prop_eggs_developed_beta_quadratic, prop_eggs_developed_beta_quadratic_null)  # model is globally significant

# (5) model summary
summary(prop_eggs_developed_beta_quadratic)
model_summary(prop_eggs_developed_beta_quadratic, model_name = "Proportion of eggs developed", transform_estimate = "exp")
model_forest_plot(prop_eggs_developed_beta_quadratic, model_name = "Proportion of eggs developed", transform_estimate = "exp")
Anova(prop_eggs_developed_beta_quadratic, type = 2)
# confint(profile(prop_eggs_developed_beta_quadratic)) %>% view

# (6) emmeans
emmeans_carcass_type_prop_eggs_developed <- emmeans(prop_eggs_developed_beta_quadratic, "carcass_type", type = "response")
pairs(regrid(emmeans_carcass_type_prop_eggs_developed))
cld(emmeans_carcass_type_prop_eggs_developed, Letters = letters)

# (7) model visualization
plot_model(prop_eggs_developed_beta_quadratic, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_type"))

# (8) write the model results
write_rds(prop_eggs_developed_beta_quadratic, "./03_Outputs/Data_Clean/prop_eggs_developed_beta_quadratic.rds")


# 4. Number of larvae vs. carcass weight and carcass type ----------------------
### Plot
plot_relationship(n_larvae)  # a quadratic relationship seems to exist

### Model
# (1) test quadratic term
n_larvae_poisson_linear <- glmmTMB(n_larvae ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                   data = carcass_data_clean,
                                   family = "poisson",
                                   na.action = na.omit)

n_larvae_poisson_quadratic <- glmmTMB(n_larvae ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                      data = carcass_data_clean,
                                      family = "poisson",
                                      na.action = na.omit)

lrtest(n_larvae_poisson_linear, n_larvae_poisson_quadratic)  # quadratic term is significant
AIC(n_larvae_poisson_linear, n_larvae_poisson_quadratic)  # quadratic model is better

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
                                    ziformula = ~ poly(carcass_weight, 2),
                                    family = "nbinom2",
                                    na.action = na.omit)

testZeroInflation(n_larvae_nb_quadratic)  # zero inflation is significant
lrtest(n_larvae_nb_quadratic, n_larvae_zi_nb_quadratic)  # zero inflation is significant
AIC(n_larvae_nb_quadratic, n_larvae_zi_nb_quadratic)  # zero-inflated model is better

# (4) test interaction term
n_larvae_zi_nb_quadratic_wo_interaction <- glmmTMB(n_larvae ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                   data = carcass_data_clean,
                                                   ziformula = ~ poly(carcass_weight, 2),
                                                   family = "nbinom2",
                                                   na.action = na.omit)

lrtest(n_larvae_zi_nb_quadratic, n_larvae_zi_nb_quadratic_wo_interaction)  # interaction is significant
AIC(n_larvae_zi_nb_quadratic, n_larvae_zi_nb_quadratic_wo_interaction)  # model with interaction is better

# (5) model diagnostics
plot(simulateResiduals(n_larvae_zi_nb_quadratic))  # acceptable
check_model(n_larvae_zi_nb_quadratic)  # some patterns of heteroscedasticity

# (6) model significance
n_larvae_zi_nb_quadratic_null <- glmmTMB(n_larvae ~ 1,
                                         data = carcass_data_clean,
                                         ziformula = ~ 1,
                                         family = "nbinom2",
                                         na.action = na.omit)

# lrtest(n_larvae_zi_nb_quadratic_null, n_larvae_zi_nb_quadratic)  # non-comparable because of the difference in sample sizes

# (7) model summary
summary(n_larvae_zi_nb_quadratic)
model_summary(n_larvae_zi_nb_quadratic, model_name = "Number of larvae", transform_estimate = "exp")
model_forest_plot(n_larvae_zi_nb_quadratic, model_name = "Number of larvae", transform_estimate = "exp")
Anova(n_larvae_zi_nb_quadratic, type = 2)
# confint(profile(n_larvae_zi_nb_quadratic)) %>% view

# (8) emmeans
emmeans_carcass_type_n_larvae <- emmeans(n_larvae_zi_nb_quadratic, "carcass_type", type = "response")
pairs(regrid(emmeans_carcass_type_n_larvae))
cld(emmeans_carcass_type_n_larvae, Letters = letters)

# (9) model visualization
plot_model(n_larvae_zi_nb_quadratic, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_type"))

# (10) write the model results
write_rds(n_larvae_zi_nb_quadratic, "./03_Outputs/Data_Clean/n_larvae_zi_nb_quadratic.rds")


# 5.1 Total larval mass vs. carcass weight and carcass type (with zeros) -------
### Plot
plot_relationship(total_larval_mass)  # a quadratic relationship seems to exist

### Model
# (1) test quadratic term
total_larval_mass_tweedie_linear <- glmmTMB(total_larval_mass ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                            data = carcass_data_clean,
                                            family = "tweedie",
                                            na.action = na.omit)

total_larval_mass_tweedie_quadratic <- glmmTMB(total_larval_mass ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                               data = carcass_data_clean,
                                               family = "tweedie",
                                               na.action = na.omit)

lrtest(total_larval_mass_tweedie_linear, total_larval_mass_tweedie_quadratic)  # quadratic term is significant
AIC(total_larval_mass_tweedie_linear, total_larval_mass_tweedie_quadratic)  # quadratic model is better

# (2) test interaction term
total_larval_mass_tweedie_quadratic_wo_interaction <- glmmTMB(total_larval_mass ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                              data = carcass_data_clean,
                                                              family = "tweedie",
                                                              na.action = na.omit)

lrtest(total_larval_mass_tweedie_quadratic, total_larval_mass_tweedie_quadratic_wo_interaction)  # interaction is not significant
AIC(total_larval_mass_tweedie_quadratic, total_larval_mass_tweedie_quadratic_wo_interaction)  # model without interaction is better

# (3) model diagnostics
plot(simulateResiduals(total_larval_mass_tweedie_quadratic))  # some patterns of heteroscedasticity
check_model(total_larval_mass_tweedie_quadratic)  # some patterns of heteroscedasticity

# (4) model significance
total_larval_mass_tweedie_quadratic_null <- glmmTMB(total_larval_mass ~ 1,
                                                    data = carcass_data_clean,
                                                    family = "tweedie",
                                                    na.action = na.omit)

lrtest(total_larval_mass_tweedie_quadratic, total_larval_mass_tweedie_quadratic_null)  # model is globally significant

# (5) model summary
summary(total_larval_mass_tweedie_quadratic)
model_summary(total_larval_mass_tweedie_quadratic, model_name = "Total larval mass", transform_estimate = NULL)
model_forest_plot(total_larval_mass_tweedie_quadratic, model_name = "Total larval mass", transform_estimate = NULL)
Anova(total_larval_mass_tweedie_quadratic, type = 2)
# confint(profile(total_larval_mass_tweedie_quadratic)) %>% view

# (6) emmeans
emmeans_carcass_type_total_larval_mass <- emmeans(total_larval_mass_tweedie_quadratic, "carcass_type")
pairs(regrid(emmeans_carcass_type_total_larval_mass))
cld(emmeans_carcass_type_total_larval_mass, Letters = letters)

# (7) model visualization
plot_model(total_larval_mass_tweedie_quadratic, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_type"))

# (8) write the model results
write_rds(total_larval_mass_tweedie_quadratic, "./03_Outputs/Data_Clean/total_larval_mass_tweedie_quadratic.rds")


# 5.2 Total larval mass vs. carcass weight and carcass type (without zeros) ----
### Remove zeros
carcass_data_clean_total_larval_mass <- carcass_data_clean %>% 
  filter(total_larval_mass > 0)

### Plot
ggplot(carcass_data_clean_total_larval_mass, aes(x = carcass_weight, y = total_larval_mass, color = carcass_type)) + 
  geom_point() + 
  geom_smooth(se = F) + 
  scale_color_brewer(palette = "Set1")  # a quadratic relationship seems to exist

### Model
# (1) test quadratic term
total_larval_mass_gaussian_linear <- glmmTMB(total_larval_mass ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                             data = carcass_data_clean_total_larval_mass,
                                             family = "gaussian",
                                             na.action = na.omit)

total_larval_mass_gaussian_quadratic <- glmmTMB(total_larval_mass ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                data = carcass_data_clean_total_larval_mass,
                                                family = "gaussian",
                                                na.action = na.omit)

lrtest(total_larval_mass_gaussian_linear, total_larval_mass_gaussian_quadratic)  # quadratic term is significant
AIC(total_larval_mass_gaussian_linear, total_larval_mass_gaussian_quadratic)  # quadratic model is better

# (2) test interaction term
total_larval_mass_gaussian_quadratic_wo_interaction <- glmmTMB(total_larval_mass ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                               data = carcass_data_clean_total_larval_mass,
                                                               family = "gaussian",
                                                               na.action = na.omit)

lrtest(total_larval_mass_gaussian_quadratic, total_larval_mass_gaussian_quadratic_wo_interaction)  # interaction is not significant
AIC(total_larval_mass_gaussian_quadratic, total_larval_mass_gaussian_quadratic_wo_interaction)  # model without interaction is better

# (3) model diagnostics
plot(simulateResiduals(total_larval_mass_gaussian_quadratic))  # some patterns of heteroscedasticity
check_model(total_larval_mass_gaussian_quadratic)  # some patterns of heteroscedasticity

# (4) model significance
total_larval_mass_gaussian_quadratic_null <- glmmTMB(total_larval_mass ~ 1,
                                                     data = carcass_data_clean_total_larval_mass,
                                                     family = "gaussian",
                                                     na.action = na.omit)

lrtest(total_larval_mass_gaussian_quadratic, total_larval_mass_gaussian_quadratic_null)  # model is globally significant

# (5) model summary
summary(total_larval_mass_gaussian_quadratic)
model_summary(total_larval_mass_gaussian_quadratic, model_name = "Total larval mass", transform_estimate = NULL)
model_forest_plot(total_larval_mass_gaussian_quadratic, model_name = "Total larval mass", transform_estimate = NULL)
Anova(total_larval_mass_gaussian_quadratic, type = 2)
# confint(profile(total_larval_mass_gaussian_quadratic)) %>% view

# (6) emmeans
emmeans_carcass_type_total_larval_mass <- emmeans(total_larval_mass_gaussian_quadratic, "carcass_type")
pairs(regrid(emmeans_carcass_type_total_larval_mass))
cld(emmeans_carcass_type_total_larval_mass, Letters = letters)

# (7) model visualization
plot_model(total_larval_mass_gaussian_quadratic, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_type"))

# (8) write the model results
write_rds(total_larval_mass_gaussian_quadratic, "./03_Outputs/Data_Clean/total_larval_mass_gaussian_quadratic.rds")


# 5.3 Total larval mass vs. carcass weight and carcass type (large wild carcasses removed) ----
### Remove zeros and two large wild carcasses
carcass_data_clean_total_larval_mass_large_removed <- carcass_data_clean %>% 
  filter(total_larval_mass > 0) %>% 
  filter(carcass_weight < 90)

### Plot
ggplot(carcass_data_clean_total_larval_mass_large_removed, aes(x = carcass_weight, y = total_larval_mass, color = carcass_type)) + 
  geom_point() + 
  geom_smooth(se = F) + 
  scale_color_brewer(palette = "Set1")  # a quadratic relationship seems to exist

### Model
# (1) test quadratic term
total_larval_mass_gaussian_linear_large_removed <- glmmTMB(total_larval_mass ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                             data = carcass_data_clean_total_larval_mass_large_removed,
                                             family = "gaussian",
                                             na.action = na.omit)

total_larval_mass_gaussian_quadratic_large_removed <- glmmTMB(total_larval_mass ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                data = carcass_data_clean_total_larval_mass_large_removed,
                                                family = "gaussian",
                                                na.action = na.omit)

lrtest(total_larval_mass_gaussian_linear_large_removed, total_larval_mass_gaussian_quadratic_large_removed)  # quadratic term is significant
AIC(total_larval_mass_gaussian_linear_large_removed, total_larval_mass_gaussian_quadratic_large_removed)  # quadratic model is better

# (2) test interaction term
total_larval_mass_gaussian_quadratic_large_removed_wo_interaction <- glmmTMB(total_larval_mass ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                               data = carcass_data_clean_total_larval_mass_large_removed,
                                                               family = "gaussian",
                                                               na.action = na.omit)

lrtest(total_larval_mass_gaussian_quadratic_large_removed, total_larval_mass_gaussian_quadratic_large_removed_wo_interaction)  # interaction is not significant
AIC(total_larval_mass_gaussian_quadratic_large_removed, total_larval_mass_gaussian_quadratic_large_removed_wo_interaction)  # model without interaction is better

# (3) model diagnostics
plot(simulateResiduals(total_larval_mass_gaussian_quadratic_large_removed))  # some patterns of heteroscedasticity
check_model(total_larval_mass_gaussian_quadratic_large_removed)  # some patterns of heteroscedasticity

# (4) model significance
total_larval_mass_gaussian_quadratic_large_removed_null <- glmmTMB(total_larval_mass ~ 1,
                                                     data = carcass_data_clean_total_larval_mass_large_removed,
                                                     family = "gaussian",
                                                     na.action = na.omit)

lrtest(total_larval_mass_gaussian_quadratic_large_removed, total_larval_mass_gaussian_quadratic_large_removed_null)  # model is globally significant

# (5) model summary
summary(total_larval_mass_gaussian_quadratic_large_removed)
model_summary(total_larval_mass_gaussian_quadratic_large_removed, model_name = "Total larval mass", transform_estimate = NULL)
model_forest_plot(total_larval_mass_gaussian_quadratic_large_removed, model_name = "Total larval mass", transform_estimate = NULL)
Anova(total_larval_mass_gaussian_quadratic_large_removed, type = 2)
# confint(profile(total_larval_mass_gaussian_quadratic_large_removed)) %>% view

# (6) emmeans
emmeans_carcass_type_total_larval_mass_large_removed <- emmeans(total_larval_mass_gaussian_quadratic_large_removed, "carcass_type")
pairs(regrid(emmeans_carcass_type_total_larval_mass_large_removed))
cld(emmeans_carcass_type_total_larval_mass_large_removed, Letters = letters)

# (7) model visualization
plot_model(total_larval_mass_gaussian_quadratic_large_removed, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_type"))

# (8) write the model results
write_rds(total_larval_mass_gaussian_quadratic_large_removed, "./03_Outputs/Data_Clean/total_larval_mass_gaussian_quadratic_large_removed.rds")


# 6. Average larval mass vs. carcass weight and carcass type -------------------
### Plot
plot_relationship(average_larval_mass)

### Model
# (1) Test quadratic term
average_larval_mass_gaussian_linear <- glmmTMB(average_larval_mass ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                               data = carcass_data_clean,
                                               family = "gaussian",
                                               na.action = na.omit)

average_larval_mass_gaussian_quadratic <- glmmTMB(average_larval_mass ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                  data = carcass_data_clean,
                                                  family = "gaussian",
                                                  na.action = na.omit)

lrtest(average_larval_mass_gaussian_linear, average_larval_mass_gaussian_quadratic)  # quadratic term is significantly better
AIC(average_larval_mass_gaussian_linear, average_larval_mass_gaussian_quadratic)  # quadratic model is better

# (2) test interaction term
average_larval_mass_gaussian_quadratic_wo_interaction <- glmmTMB(average_larval_mass ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                                 data = carcass_data_clean,
                                                                 family = "gaussian",
                                                                 na.action = na.omit)

lrtest(average_larval_mass_gaussian_quadratic, average_larval_mass_gaussian_quadratic_wo_interaction)  # interaction is not significant
AIC(average_larval_mass_gaussian_quadratic, average_larval_mass_gaussian_quadratic_wo_interaction)  # model without interaction is slightly better

# (3) model diagnostics
plot(simulateResiduals(average_larval_mass_gaussian_quadratic))  # no pattern
check_model(average_larval_mass_gaussian_quadratic)  # no pattern

# (4) model significance
average_larval_mass_gaussian_quadratic_null <- glmmTMB(average_larval_mass ~ 1,
                                                       data = carcass_data_clean,
                                                       family = "gaussian",
                                                       na.action = na.omit)

lrtest(average_larval_mass_gaussian_quadratic, average_larval_mass_gaussian_quadratic_null)  # model is globally significant

# (5) model summary
summary(average_larval_mass_gaussian_quadratic)
model_summary(average_larval_mass_gaussian_quadratic, model_name = "Average larval mass", transform_estimate = NULL)
model_forest_plot(average_larval_mass_gaussian_quadratic, model_name = "Average larval mass", transform_estimate = NULL)
Anova(average_larval_mass_gaussian_quadratic, type = 2)
# confint(profile(average_larval_mass_gaussian_quadratic)) %>% view

# (6) emmeans
emmeans_carcass_type_average_larval_mass <- emmeans(average_larval_mass_gaussian_quadratic, "carcass_type")
pairs(regrid(emmeans_carcass_type_average_larval_mass))
cld(emmeans_carcass_type_average_larval_mass, Letters = letters)

# (7) model visualization
plot_model(average_larval_mass_gaussian_quadratic, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_type"))

# (8) write the model results
write_rds(average_larval_mass_gaussian_quadratic, "./03_Outputs/Data_Clean/average_larval_mass_gaussian_quadratic.rds")


# 7. Larval density vs. carcass weight and carcass type ------------------------
### Plot
plot_relationship(larval_density)

### Model
# (1) test quadratic term
larval_density_gaussian_linear <- glmmTMB(larval_density ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                          data = carcass_data_clean,
                                          family = "gaussian",
                                          na.action = na.omit)

larval_density_gaussian_quadratic <- glmmTMB(larval_density ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                             data = carcass_data_clean,
                                             family = "gaussian",
                                             na.action = na.omit)

lrtest(larval_density_gaussian_linear, larval_density_gaussian_quadratic)  # quadratic term is not significant
AIC(larval_density_gaussian_linear, larval_density_gaussian_quadratic)  # linear model is better

# (2) test interaction term
larval_density_gaussian_linear_wo_interaction <- glmmTMB(larval_density ~ carcass_weight + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                         data = carcass_data_clean,
                                                         family = "gaussian",
                                                         na.action = na.omit)

lrtest(larval_density_gaussian_linear, larval_density_gaussian_linear_wo_interaction)  # interaction is not significant
AIC(larval_density_gaussian_linear, larval_density_gaussian_linear_wo_interaction)  # model without interaction is better

# (3) model diagnostics
plot(simulateResiduals(larval_density_gaussian_linear))  # acceptable
check_model(larval_density_gaussian_linear)  # acceptable

# (4) model significance
larval_density_gaussian_linear_null <- glmmTMB(larval_density ~ 1,
                                               data = carcass_data_clean,
                                               family = "gaussian",
                                               na.action = na.omit)

# lrtest(larval_density_gaussian_linear, larval_density_gaussian_linear_null)  # non-comparable because of the difference in sample sizes

# (5) model summary
summary(larval_density_gaussian_linear)
model_summary(larval_density_gaussian_linear, model_name = "Larval density", transform_estimate = NULL)
model_forest_plot(larval_density_gaussian_linear, model_name = "Larval density", transform_estimate = NULL)
Anova(larval_density_gaussian_linear, type = 2)
# confint(profile(larval_density_gaussian_linear)) %>% view

# (6) emmeans
emmeans_carcass_type_larval_density <- emmeans(larval_density_gaussian_linear, "carcass_type")
pairs(regrid(emmeans_carcass_type_larval_density))
cld(emmeans_carcass_type_larval_density, Letters = letters)

# (7) model visualization
plot_model(larval_density_gaussian_linear, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_type"))

# (8) write the model results
write_rds(larval_density_gaussian_linear, "./03_Outputs/Data_Clean/larval_density_gaussian_linear.rds")


# 8. Carcass weight loss vs. carcass weight and carcass type -------------------
### Plot
plot_relationship(carcass_weight_loss)  # one impossible value and two outliers

### Exclude the impossible value, two outliers, and the observations without any larva
carcass_data_clean_carcass_weight_loss <- carcass_data_clean %>% 
  filter(carcass_weight_loss < 25 & carcass_weight_loss > 0) %>% 
  filter(breeding_success == 1)

### Replot the data
ggplot(carcass_data_clean_carcass_weight_loss, aes(x = carcass_weight, y = carcass_weight_loss, color = carcass_type)) + 
  geom_point() + 
  geom_smooth(se = F) + 
  scale_color_brewer(palette = "Set1")  # a quadratic relationship seems to exist

### Model
# (1) test quadratic term
carcass_weight_loss_gaussian_linear <- glmmTMB(carcass_weight_loss ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                               data = carcass_data_clean_carcass_weight_loss,
                                               family = "gaussian",
                                               na.action = na.omit)

carcass_weight_loss_gaussian_quadratic <- glmmTMB(carcass_weight_loss ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                  data = carcass_data_clean_carcass_weight_loss,
                                                  family = "gaussian",
                                                  na.action = na.omit)

lrtest(carcass_weight_loss_gaussian_linear, carcass_weight_loss_gaussian_quadratic)  # quadratic term is significant
AIC(carcass_weight_loss_gaussian_linear, carcass_weight_loss_gaussian_quadratic)  # quadratic model is slightly better

# (2) test interaction term
carcass_weight_loss_gaussian_quadratic_wo_interaction <- glmmTMB(carcass_weight_loss ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                                 data = carcass_data_clean_carcass_weight_loss,
                                                                 family = "gaussian",
                                                                 na.action = na.omit)

lrtest(carcass_weight_loss_gaussian_quadratic, carcass_weight_loss_gaussian_quadratic_wo_interaction)  # interaction term is significant
AIC(carcass_weight_loss_gaussian_quadratic, carcass_weight_loss_gaussian_quadratic_wo_interaction)  # model with interaction is slightly better

# (3) model diagnostics
plot(simulateResiduals(carcass_weight_loss_gaussian_quadratic))  # acceptable
check_model(carcass_weight_loss_gaussian_quadratic)  # acceptable

# (4) model significance
carcass_weight_loss_gaussian_quadratic_null <- glmmTMB(carcass_weight_loss ~ 1,
                                                    data = carcass_data_clean_carcass_weight_loss,
                                                    family = "gaussian",
                                                    na.action = na.omit)

lrtest(carcass_weight_loss_gaussian_quadratic, carcass_weight_loss_gaussian_quadratic_null)  # model is globally significant

# (5) model summary
summary(carcass_weight_loss_gaussian_quadratic)
model_summary(carcass_weight_loss_gaussian_quadratic, model_name = "Carcass weight loss", transform_estimate = NULL)
model_forest_plot(carcass_weight_loss_gaussian_quadratic, model_name = "Carcass weight loss", transform_estimate = NULL)
Anova(carcass_weight_loss_gaussian_quadratic, type = 2)
# confint(profile(carcass_weight_loss_gaussian_quadratic)) %>% view

# (6) emmeans
emmeans_carcass_type_carcass_weight_loss <- emmeans(carcass_weight_loss_gaussian_linear, "carcass_type")
pairs(regrid(emmeans_carcass_type_carcass_weight_loss))
cld(emmeans_carcass_type_carcass_weight_loss, Letters = letters)

# (7) model visualization
plot_model(carcass_weight_loss_gaussian_quadratic, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_type"))

# (8) write the model results
write_rds(carcass_weight_loss_gaussian_quadratic, "./03_Outputs/Data_Clean/carcass_weight_loss_gaussian_quadratic.rds")


# 9. Proportion of carcass used vs. carcass weight and carcass type ------------
### Plot
plot_relationship(prop_carcass_used)  # one impossible value and two outliers

### Remove the impossible values and two outliers
carcass_data_clean_prop_carcass_used <- carcass_data_clean %>% 
  filter(prop_carcass_used < 0.7 & prop_carcass_used > 0)

### Replot the data
ggplot(carcass_data_clean_prop_carcass_used, aes(x = carcass_weight, y = prop_carcass_used, color = carcass_type)) + 
  geom_point() + 
  geom_smooth(se = F) + 
  scale_color_brewer(palette = "Set1") 

### Model
# (1) test quadratic term
prop_carcass_used_beta_linear <- glmmTMB(prop_carcass_used ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                         data = carcass_data_clean_prop_carcass_used,
                                         family = beta_family("logit"),
                                         na.action = na.omit)

prop_carcass_used_beta_quadratic <- glmmTMB(prop_carcass_used ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                            data = carcass_data_clean_prop_carcass_used,
                                            family = beta_family("logit"),
                                            na.action = na.omit)

lrtest(prop_carcass_used_beta_linear, prop_carcass_used_beta_quadratic)  # quadratic term not significant
AIC(prop_carcass_used_beta_linear, prop_carcass_used_beta_quadratic)  # quadratic model is not better

# (2) test interaction term
prop_carcass_used_beta_linear_wo_interaction <- glmmTMB(prop_carcass_used ~ carcass_weight + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                        data = carcass_data_clean_prop_carcass_used,
                                                        family = beta_family("logit"),
                                                        na.action = na.omit)

lrtest(prop_carcass_used_beta_linear, prop_carcass_used_beta_linear_wo_interaction)  # interaction is not significant
AIC(prop_carcass_used_beta_linear, prop_carcass_used_beta_linear_wo_interaction)  # model without interaction is better

# (3) model diagnostics
plot(simulateResiduals(prop_carcass_used_beta_linear))  # no pattern
check_model(prop_carcass_used_beta_linear)

# (4) model significance
prop_carcass_used_beta_linear_null <- glmmTMB(prop_carcass_used ~ 1,
                                              data = carcass_data_clean_prop_carcass_used,
                                              family = beta_family("logit"),
                                              na.action = na.omit)

lrtest(prop_carcass_used_beta_linear, prop_carcass_used_beta_linear_null)  # model is globally significant

# (5) model summary
summary(prop_carcass_used_beta_linear)
model_summary(prop_carcass_used_beta_linear, model_name = "Proportion of carcass used", transform_estimate = "exp")
model_forest_plot(prop_carcass_used_beta_linear, model_name = "Proportion of carcass used", transform_estimate = "exp")
Anova(prop_carcass_used_beta_linear, type = 2)
# confint(profile(prop_carcass_used_beta_linear)) %>% view

# (6) emmeans
emmeans_carcass_type_prop_carcass_used <- emmeans(prop_carcass_used_beta_linear, "carcass_type", type = "response")
pairs(regrid(emmeans_carcass_type_prop_carcass_used))
cld(emmeans_carcass_type_prop_carcass_used, Letters = letters)

# (7) model visualization
plot_model(prop_carcass_used_beta_linear, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_type"))

# (8) write the model results
write_rds(prop_carcass_used_beta_linear, "./03_Outputs/Data_Clean/prop_carcass_used_beta_linear.rds")


# 10. Average larval mass vs. larval density -----------------------------------
### Plot
ggplot(carcass_data_clean, aes(x = larval_density, y = average_larval_mass, color = carcass_type)) + 
  geom_point() + 
  geom_smooth(se = F) + 
  scale_color_brewer(palette = "Set1")  # a linear relationship

### Model
# (1) test interaction term
average_larval_mass_larval_density_gaussian_linear <- glmmTMB(average_larval_mass ~ larval_density * carcass_type,
                                                              data = carcass_data_clean,
                                                              family = "gaussian",
                                                              na.action = na.omit)

average_larval_mass_larval_density_gaussian_linear_wo_interaction <- glmmTMB(average_larval_mass ~ larval_density + carcass_type,
                                                                             data = carcass_data_clean,
                                                                             family = "gaussian",
                                                                             na.action = na.omit)

lrtest(average_larval_mass_larval_density_gaussian_linear, average_larval_mass_larval_density_gaussian_linear_wo_interaction)  # interaction is not significant
AIC(average_larval_mass_larval_density_gaussian_linear, average_larval_mass_larval_density_gaussian_linear_wo_interaction)  # model without interaction is slightly better

# (2) model diagnostics
plot(simulateResiduals(average_larval_mass_larval_density_gaussian_linear))  # acceptable
check_model(average_larval_mass_larval_density_gaussian_linear)  # acceptable

# (3) model significance
average_larval_mass_larval_density_gaussian_linear_null <- glmmTMB(average_larval_mass ~ 1,
                                                                   data = carcass_data_clean,
                                                                   family = "gaussian",
                                                                   na.action = na.omit)

lrtest(average_larval_mass_larval_density_gaussian_linear, average_larval_mass_larval_density_gaussian_linear_null)  # model is globally significant

# (4) model summary
summary(average_larval_mass_larval_density_gaussian_linear)
model_summary(average_larval_mass_larval_density_gaussian_linear, model_name = "Average larval mass vs. Larval density", transform_estimate = "exp")
model_forest_plot(average_larval_mass_larval_density_gaussian_linear, model_name = "Average larval mass vs. Larval density", transform_estimate = "exp")
Anova(average_larval_mass_larval_density_gaussian_linear, type = 2)
# confint(profile(average_larval_mass_larval_density_gaussian_linear)) %>% view

# (5) emmeans
emmeans_carcass_type_average_larval_mass_larval_density <- emmeans(average_larval_mass_larval_density_gaussian_linear, "carcass_type")
pairs(regrid(emmeans_carcass_type_average_larval_mass_larval_density))
cld(emmeans_carcass_type_average_larval_mass_larval_density, Letters = letters)

# (6) model visualization
plot_model(average_larval_mass_larval_density_gaussian_linear, 
           type = "pred", 
           terms = c("larval_density [0:2]", "carcass_type"))

# (7) write the model results
write_rds(average_larval_mass_larval_density_gaussian_linear, "./03_Outputs/Data_Clean/average_larval_mass_larval_density_gaussian_linear.rds")


# 11. Package citations --------------------------------------------------------
capture.output(utils:::print.bibentry(citation(), style = "Bibtex"),
               utils:::print.bibentry(citation("glmmTMB"), style = "Bibtex"),
               utils:::print.bibentry(citation("DHARMa"), style = "Bibtex"),
               utils:::print.bibentry(citation("car"), style = "Bibtex"),
               utils:::print.bibentry(citation("emmeans"), style = "Bibtex"),
               file = "./05_References/R_Pacakge_Citations.bib")










