## -----------------------------------------------------------------------------
## Title: Analysis of the relationships between carcass attributes and beetle breeding outcomes
##
## Author: Gen-Chang Hsu
##
## Date: 2024-02-06
##
## Description:
## 1. Model the relationship between clutch size vs. carcass weight and carcass type
## 2. Model the relationship between breeding success vs. carcass weight and carcass type
## 3. 
##
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
library(sjPlot)
library(broom)
library(broom.mixed)
library(emmeans)
library(multcomp)


# Model summary and plot functions from the package "sjPlot" -------------------
model_summary <- function(model, model_name, transform_estimate) {
  tab_model(model,
            dv.labels = model_name,
            auto.label = T,
            show.est = T,
            show.se = T,
            show.ci = T,
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

### Exclude wild carcasses larger than 100 grams
carcass_data_clean <- carcass_data_clean %>% 
  filter(carcass_weight <= 100)
  

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
                                       ziformula = ~ 1,
                                       family = "nbinom2",
                                       na.action = na.omit)

testZeroInflation(clutch_size_nb_quadratic)
lrtest(clutch_size_nb_quadratic, clutch_size_zi_nb_quadratic)  # zero inflation is significant
AIC(clutch_size_nb_quadratic, clutch_size_zi_nb_quadratic)  # zero-inflated model is significant

# (4) test interaction term
clutch_size_zi_nb_quadratic_wo_interaction <- glmmTMB(clutch_size ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                      data = carcass_data_clean,
                                                      ziformula = ~ 1,
                                                      family = "nbinom2",
                                                      na.action = na.omit)

lrtest(clutch_size_zi_nb_quadratic, clutch_size_zi_nb_quadratic_wo_interaction)  # interaction is not significant
AIC(clutch_size_zi_nb_quadratic, clutch_size_zi_nb_quadratic_wo_interaction)  # model without interaction is better

# (5) model diagnostics
plot(simulateResiduals(clutch_size_zi_nb_quadratic))  # some patterns of heteroscedasticity
check_model(clutch_size_zi_nb_quadratic)  # some patterns of heteroscedasticity

# (6) model significance
clutch_size_zi_nb_quadratic_null <- glmmTMB(clutch_size ~ 1,
                                            data = carcass_data_clean,
                                            ziformula = ~ 1,
                                            family = "nbinom2",
                                            na.action = na.omit)

lrtest(clutch_size_zi_nb_quadratic_null, clutch_size_zi_nb_quadratic)  # model is globally significant

# (7) model summary
summary(clutch_size_zi_nb_quadratic)
# tidy(clutch_size_zi_nb_quadratic) %>% view
model_summary(clutch_size_zi_nb_quadratic, model_name = "Clutch size", transform_estimate = "exp")
model_forest_plot(clutch_size_zi_nb_quadratic, model_name = "Clutch size", transform_estimate = "exp")
Anova(clutch_size_zi_nb_quadratic, type = 2)
# confint(profile(clutch_size_zi_nb_quadratic)) %>% view

# (8) emmeans
emmeans_carcass_type_clutch_size <- emmeans(clutch_size_zi_nb_quadratic, "carcass_type", type = "response")
emmeans_parent_generation_clutch_size <- emmeans(clutch_size_zi_nb_quadratic, "parent_generation", type = "response")

pairs(regrid(emmeans_carcass_type_clutch_size))
pairs(regrid(emmeans_parent_generation_clutch_size))

cld(emmeans_carcass_type_clutch_size, adjust = "Tukey", Letters = letters)
cld(emmeans_parent_generation_clutch_size, adjust = "Tukey", Letters = letters)

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
# (1) Test quadratic term
breeding_success_logistic_linear <- glmmTMB(breeding_success ~ carcass_weight * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                            data = carcass_data_clean,
                                            family = "binomial",
                                            na.action = na.omit)

breeding_success_logistic_quadratic <- glmmTMB(breeding_success ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                               data = carcass_data_clean,
                                               family = "binomial",
                                               na.action = na.omit)

lrtest(breeding_success_logistic_linear, breeding_success_logistic_quadratic)  # quadratic term is significant
AIC(breeding_success_logistic_linear, breeding_success_logistic_quadratic)  # quadratic model is better

# (2) test interaction term
breeding_success_logistic_quadratic_wo_interaction <- glmmTMB(breeding_success ~ poly(carcass_weight, 2) + carcass_type + male_size + female_size + parent_generation + (1|generation_pair_id),
                                                              data = carcass_data_clean,
                                                              family = "binomial",
                                                              na.action = na.omit)

lrtest(breeding_success_logistic_quadratic, breeding_success_logistic_quadratic_wo_interaction)  # interaction is not significant
AIC(breeding_success_logistic_quadratic, breeding_success_logistic_quadratic_wo_interaction)  #  model with interaction is not better

# (3) model diagnostics
plot(simulateResiduals(breeding_success_logistic_quadratic))  # no obvious residual patterns
check_model(breeding_success_logistic_quadratic)  # residual patterns acceptable

# (4) model significance
breeding_success_logistic_quadratic_null <- glmmTMB(breeding_success ~ 1,
                                                    data = carcass_data_clean,
                                                    family = "binomial",
                                                    na.action = na.omit)

lrtest(breeding_success_logistic_quadratic, breeding_success_logistic_quadratic_null)  # non-comparable because of the difference in sample sizes

# (5) model summary
summary(breeding_success_logistic_quadratic)
# tidy(clutch_size_zi_nb_quadratic) %>% view
model_summary(breeding_success_logistic_quadratic, model_name = "Breeding success", transform_estimate = "exp")
model_forest_plot(breeding_success_logistic_quadratic, model_name = "Breeding success", transform_estimate = "exp")
Anova(breeding_success_logistic_quadratic, type = 2)
# confint(profile(clutch_size_zi_nb_quadratic)) %>% view

# (6) emmeans
emmeans_carcass_type_breeding_success <- emmeans(breeding_success_logistic_quadratic, "carcass_type", type = "response")
emmeans_parent_generation_breeding_success <- emmeans(breeding_success_logistic_quadratic, "parent_generation", type = "response")

pairs(regrid(emmeans_carcass_type_breeding_success))
pairs(regrid(emmeans_parent_generation_breeding_success))

cld(emmeans_carcass_type_breeding_success, adjust = "Tukey", Letters = letters)
cld(emmeans_parent_generation_breeding_success, adjust = "Tukey", Letters = letters)

# (7) model visualization
plot_model(breeding_success_logistic_quadratic, 
           type = "pred", 
           terms = c("carcass_weight [0:100]", "carcass_type"))

# (8) write the model results
write_rds(breeding_success_logistic_quadratic, "./03_Outputs/Data_Clean/breeding_success_logistic_quadratic.rds")















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
emmeans_carcass_type_prop_eggs_developed <- emmeans(prop_eggs_developed_beta_linear, "carcass_type", type = "response")
emmeans_parent_generation_prop_eggs_developed <- emmeans(prop_eggs_developed_beta_linear, "parent_generation", type = "response")

pairs(regrid(emmeans_carcass_type_prop_eggs_developed))
pairs(regrid(emmeans_parent_generation_prop_eggs_developed))



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
emmeans_carcass_type_n_larvae <- emmeans(n_larvae_zi_nb_quadratic, "carcass_type", type = "response")
emmeans_parent_generation_n_larvae <- emmeans(n_larvae_zi_nb_quadratic, "parent_generation", type = "response")

pairs(regrid(emmeans_carcass_type_n_larvae))
pairs(regrid(emmeans_parent_generation_n_larvae))


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
emmeans_carcass_type_total_larval_mass <- emmeans(total_larval_mass_gaussian_quadratic, "carcass_type")
emmeans_parent_generation_total_larval_mass <- emmeans(total_larval_mass_gaussian_quadratic, "parent_generation")

pairs(regrid(emmeans_carcass_type_total_larval_mass))
pairs(regrid(emmeans_parent_generation_total_larval_mass))


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
emmeans_carcass_type_average_larval_mass <- emmeans(average_larval_mass_gaussian_linear, "carcass_type")
emmeans_parent_generation_average_larval_mass <- emmeans(average_larval_mass_gaussian_linear, "parent_generation")

pairs(regrid(emmeans_carcass_type_average_larval_mass))
pairs(regrid(emmeans_parent_generation_average_larval_mass))


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
emmeans_carcass_type_larval_density <- emmeans(larval_density_gaussian_linear, "carcass_type")
emmeans_parent_generation_larval_density <- emmeans(larval_density_gaussian_linear, "parent_generation")

pairs(regrid(emmeans_carcass_type_larval_density))
pairs(regrid(emmeans_parent_generation_larval_density))


# 8. Carcass used vs. carcass weight and carcass type --------------------------
### Plot
plot_relationship(carcass_weight_loss)

### Only include observations with successful breeding events
carcass_data_clean_carcass_weight_loss <- carcass_data_clean %>% 
  filter(breeding_success == 1)

### Model
# (1) Test the quadratic term (need to add parent generation to the model later)
carcass_weight_loss_gaussian_linear <- glmmTMB(carcass_weight_loss ~ carcass_weight * carcass_type + male_size + female_size + (1|generation_pair_id),
                                               data = carcass_data_clean_carcass_weight_loss,
                                               family = "gaussian",
                                               na.action = na.omit)

carcass_weight_loss_gaussian_quadratic <- glmmTMB(carcass_weight_loss ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + (1|generation_pair_id),
                                                  data = carcass_data_clean_carcass_weight_loss,
                                                  family = "gaussian",
                                                  na.action = na.omit)

lrtest(carcass_weight_loss_gaussian_linear, carcass_weight_loss_gaussian_quadratic)  # the quadratic model is not significantly better
AIC(carcass_weight_loss_gaussian_linear, carcass_weight_loss_gaussian_quadratic)  # the linear model is better

# (2) test the interaction term
carcass_weight_loss_gaussian_linear_wo_interaction <- glmmTMB(carcass_weight_loss ~ carcass_weight + carcass_type + male_size + female_size + (1|generation_pair_id),
                                                              data = carcass_data_clean_carcass_weight_loss,
                                                              family = "gaussian",
                                                              na.action = na.omit)

lrtest(carcass_weight_loss_gaussian_linear, carcass_weight_loss_gaussian_linear_wo_interaction)  # the interaction term is not significant
AIC(carcass_weight_loss_gaussian_linear, carcass_weight_loss_gaussian_linear_wo_interaction)  # the model without interaction term is better

# (3) model diagnostics
plot(simulateResiduals(carcass_weight_loss_gaussian_linear))
check_model(carcass_weight_loss_gaussian_linear)  # residual plot looks acceptable (expect for one outlier)

# (4) model significance
carcass_weight_loss_gaussian_linear_null <- glmmTMB(carcass_weight_loss ~ 1,
                                                    data = carcass_data_clean_carcass_weight_loss,
                                                    family = "gaussian",
                                                    na.action = na.omit)

lrtest(carcass_weight_loss_gaussian_linear, carcass_weight_loss_gaussian_linear_null)  # the model is globally significant

# (5) coefficient significance
summary(carcass_weight_loss_gaussian_linear)
tidy(carcass_weight_loss_gaussian_linear) %>% view
Anova(carcass_weight_loss_gaussian_linear, type = 2)
confint(profile(carcass_weight_loss_gaussian_linear)) %>% view

# (6) emmeans
emmeans_carcass_type_carcass_weight_loss <- emmeans(carcass_weight_loss_gaussian_linear, "carcass_type")
emmeans_parent_generation_carcass_weight_loss <- emmeans(carcass_weight_loss_gaussian_linear, "parent_generation")

pairs(regrid(emmeans_carcass_type_carcass_weight_loss))
pairs(regrid(emmeans_parent_generation_carcass_weight_loss))


# 9. Carcass use efficiency vs. carcass weight and carcass type ----------------
### Plot
plot_relationship(efficiency)  # there seem to be two outliers

### Remove the outliers
carcass_data_clean_efficiency <- carcass_data_clean %>% 
  filter(efficiency < 0.5)

### Model
# (1) Test the quadratic term (need to add parent generation to the model later)
efficiency_beta_linear <- glmmTMB(efficiency ~ carcass_weight * carcass_type + male_size + female_size + (1|generation_pair_id),
                                  data = carcass_data_clean_efficiency,
                                  family = beta_family("logit"),
                                  na.action = na.omit)

efficiency_beta_quadratic <- glmmTMB(efficiency ~ poly(carcass_weight, 2) * carcass_type + male_size + female_size + (1|generation_pair_id),
                                     data = carcass_data_clean_efficiency,
                                     family = beta_family("logit"),
                                     na.action = na.omit)

lrtest(efficiency_beta_linear, efficiency_beta_quadratic)  # quadratic model is not significantly better
AIC(efficiency_beta_linear, efficiency_beta_quadratic)  # linear model is better

# (2) test the interaction term
efficiency_beta_linear_wo_interaction <- glmmTMB(efficiency ~ carcass_weight + carcass_type + male_size + female_size + (1|generation_pair_id),
                                                 data = carcass_data_clean_efficiency,
                                                 family = beta_family("logit"),
                                                 na.action = na.omit)

lrtest(efficiency_beta_linear, efficiency_beta_linear_wo_interaction)  # interaction term is not significant
AIC(efficiency_beta_linear, efficiency_beta_linear_wo_interaction)  # model without the interaction term is better

# (3) model diagnostics
plot(simulateResiduals(efficiency_beta_linear))
check_model(efficiency_beta_linear)  # residual plot looks acceptable

# (4) model significance
efficiency_beta_linear_null <- glmmTMB(efficiency ~ 1,
                                       data = carcass_data_clean_efficiency,
                                       family = beta_family("logit"),
                                       na.action = na.omit)

lrtest(efficiency_beta_linear, efficiency_beta_linear_null)  # the model is globally significant

# (5) coefficient significance
summary(efficiency_beta_linear)
tidy(efficiency_beta_linear) %>% view
Anova(efficiency_beta_linear, type = 2)
confint(profile(efficiency_beta_linear)) %>% view

# (6) emmeans
emmeans_carcass_type_efficiency <- emmeans(efficiency_beta_linear, "carcass_type")
emmeans_parent_generation_efficiency <- emmeans(efficiency_beta_linear, "parent_generation")

pairs(regrid(emmeans_carcass_type_efficiency))
pairs(regrid(emmeans_parent_generation_efficiency))


# 10. Average larval mass vs. larval density -----------------------------------
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

lrtest(average_larval_mass_larval_density_gaussian_linear, average_larval_mass_larval_density_gaussian_linear_wo_interaction)  # interaction term is not significant
AIC(average_larval_mass_larval_density_gaussian_linear, average_larval_mass_larval_density_gaussian_linear_wo_interaction)  # model without interaction term is slightly better

# (2) model diagnostics
plot(simulateResiduals(average_larval_mass_larval_density_gaussian_linear))
check_model(average_larval_mass_larval_density_gaussian_linear)  # residual plot looks fine

# (3) model significance
average_larval_mass_larval_density_gaussian_linear_null <- glmmTMB(average_larval_mass ~ 1,
                                                    data = carcass_data_clean,
                                                    family = "gaussian",
                                                    na.action = na.omit)

lrtest(average_larval_mass_larval_density_gaussian_linear, average_larval_mass_larval_density_gaussian_linear_null)

# (4) coefficient significance
summary(average_larval_mass_larval_density_gaussian_linear)
tidy(average_larval_mass_larval_density_gaussian_linear) %>% view
Anova(average_larval_mass_larval_density_gaussian_linear, type = 3)
confint(profile(average_larval_mass_larval_density_gaussian_linear)) %>% view

# (5) emmeans
emmeans_carcass_type_average_larval_mass_larval_density <- emmeans(average_larval_mass_larval_density_gaussian_linear, "carcass_type")
emmeans_parent_generation_average_larval_mass_larval_density <- emmeans(average_larval_mass_larval_density_gaussian_linear, "parent_generation")

pairs(regrid(emmeans_carcass_type_average_larval_mass_larval_density))
pairs(regrid(emmeans_parent_generation_average_larval_mass_larval_density))



















