## -----------------------------------------------------------------------------
## Title: Analysis of the breeding outcomes among different carcass taxa
##
## Author: Gen-Chang Hsu
##
## Date: 2024-05-18
##
## Description:
## 1. Summarize and visualize the data by carcass taxon
## 2. Model the relationship between number of larvae vs. carcass taxon
## 3. Model the relationship between total larval mass (without zeros) vs. carcass taxon
## 4. Model the relationship between average larval mass vs. carcass taxon
## 5. Model the relationship between proportion of carcass used vs. carcass taxon
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
carcass_data_clean <- read_csv("./03_Outputs/Data_Clean/Breeding_Data_Clean.csv")


############################### Code starts here ###############################

# 1. Data summary and visualization --------------------------------------------
### Exclude the carcass pairs with wild carcasses larger than 100 grams
large_wild_carcasses_id <- carcass_data_clean %>% 
  filter(carcass_weight > 100) %>% 
  pull(generation_pair_id)

carcass_data_clean <- carcass_data_clean %>% 
  filter(!generation_pair_id %in% large_wild_carcasses_id)

### Observations from the wild carcasses
carcass_data_clean_wild <- carcass_data_clean %>% 
  filter(carcass_type == "wild")

### Number of observations for each taxon
carcass_data_clean_wild %>%   
  group_by(carcass_taxon) %>% 
  summarise(n = n())

### Carcass weight range for each taxon
carcass_data_clean_wild %>% 
  group_by(carcass_taxon) %>% 
  summarise(Min_carcass_weight = min(carcass_weight, na.rm = T),
            Max_carcass_weight = max(carcass_weight, na.rm = T))

### Restrict the carcass weight range to that of the reptiles
carcass_data_clean_wild_restricted <- carcass_data_clean_wild %>% 
  filter(carcass_weight <= max(filter(carcass_data_clean_wild, carcass_taxon == "reptile")$carcass_weight, na.rm = T))

### Number of observations for each taxon
carcass_data_clean_wild_restricted %>%   
  group_by(carcass_taxon) %>% 
  summarise(n = n())

### Carcass weight distribution for each taxon
ggplot(carcass_data_clean_wild_restricted, 
       aes(x = carcass_taxon, y = carcass_weight, fill = carcass_taxon)) + 
  geom_boxplot() + 
  geom_point(position = position_jitter(width = 0.05)) + 
  scale_fill_brewer(palette = "Set1") + 
  theme_classic()


### Test whether the carcass weight differed among the three taxa
carcass_weight_taxon_wild_restricted_lm <- lm(carcass_weight ~ carcass_taxon, data = carcass_data_clean_wild_restricted)
summary(carcass_weight_taxon_wild_restricted_lm)
Anova(carcass_weight_taxon_wild_restricted_lm, type = 2)
carcass_weight_taxon_wild_restricted_emmeans <- emmeans(carcass_weight_taxon_wild_restricted_lm, "carcass_taxon")
pairs(carcass_weight_taxon_wild_restricted_emmeans)
cld(carcass_weight_taxon_wild_restricted_emmeans, Letters = letters)


# 2. Number of larvae vs. carcass taxon ----------------------------------------
# (1) test overdispersion
n_larvae_poisson_taxon <- glmmTMB(n_larvae ~ carcass_taxon + carcass_weight + male_size + female_size + parent_generation,
                                      data = carcass_data_clean_wild_restricted,
                                      family = "poisson",
                                      na.action = na.omit)

n_larvae_nb_taxon <- glmmTMB(n_larvae ~ carcass_taxon + carcass_weight + male_size + female_size + parent_generation,
                                            data = carcass_data_clean_wild_restricted,
                                            family = "nbinom2",
                                            na.action = na.omit)

lrtest(n_larvae_poisson_taxon, n_larvae_nb_taxon)  # overdispersion is significant
AIC(n_larvae_poisson_taxon, n_larvae_nb_taxon)  # negative binomial model is better

# (2) test zero inflation
n_larvae_zi_nb_taxon <- glmmTMB(n_larvae ~ carcass_taxon + carcass_weight + male_size + female_size + parent_generation,
                                       data = carcass_data_clean_wild_restricted,
                                       ziformula = ~ 1,
                                       family = "nbinom2",
                                       na.action = na.omit)

testZeroInflation(n_larvae_nb_taxon)
lrtest(n_larvae_nb_taxon, n_larvae_zi_nb_taxon)  # zero inflation is significant
AIC(n_larvae_nb_taxon, n_larvae_zi_nb_taxon)  # zero-inflated model is better

# (3) model diagnostics
plot(simulateResiduals(n_larvae_zi_nb_taxon))  # some patterns of heteroscedasticity
check_model(n_larvae_zi_nb_taxon)  # some patterns of heteroscedasticity

# (4) model significance
n_larvae_zi_nb_null_taxon <- glmmTMB(n_larvae ~ 1,
                                         data = carcass_data_clean_wild_restricted,
                                         ziformula = ~ 1,
                                         family = "nbinom2",
                                         na.action = na.omit)

# lrtest(n_larvae_zi_nb_null_taxon, n_larvae_zi_nb_taxon)  # non-comparable because of the difference in sample sizes

# (5) model summary
summary(n_larvae_zi_nb_taxon)
model_summary(n_larvae_zi_nb_taxon, model_name = "Number of larvae", transform_estimate = "exp")
model_forest_plot(n_larvae_zi_nb_taxon, model_name = "Number of larvae", transform_estimate = "exp")
Anova(n_larvae_zi_nb_taxon, type = 2)
# confint(profile(n_larvae_zi_nb_taxon)) %>% view

# (6) emmeans
emmeans_carcass_taxon_n_larvae <- emmeans(n_larvae_zi_nb_taxon, "carcass_taxon", type = "response")
pairs(regrid(emmeans_carcass_taxon_n_larvae))
cld(emmeans_carcass_taxon_n_larvae, Letters = letters)

# (7) model visualization
plot_model(n_larvae_zi_nb_taxon, 
           type = "pred", 
           terms = c("carcass_taxon"))

# (8) write the model results
write_rds(n_larvae_zi_nb_taxon, "./03_Outputs/Data_Clean/n_larvae_zi_nb_taxon.rds")


# 3. Total larval mass vs. carcass taxon (without zeros) -----------------------
### Remove zeros
carcass_data_clean_wild_restricted_total_larval_mass <- carcass_data_clean_wild_restricted %>% 
  filter(total_larval_mass > 0)

### Model
# (1) fit the model
total_larval_mass_gaussian_taxon <- glmmTMB(total_larval_mass ~ carcass_taxon + carcass_weight + male_size + female_size + parent_generation,
                                             data = carcass_data_clean_wild_restricted_total_larval_mass,
                                             family = "gaussian",
                                             na.action = na.omit)

# (2) model diagnostics
plot(simulateResiduals(total_larval_mass_gaussian_taxon))  # no obvious patterns
check_model(total_larval_mass_gaussian_taxon)  # acceptable

# (3) model significance
total_larval_mass_gaussian_taxon_null <- glmmTMB(total_larval_mass ~ 1,
                                                     data = carcass_data_clean_wild_restricted_total_larval_mass,
                                                     family = "gaussian",
                                                     na.action = na.omit)

lrtest(total_larval_mass_gaussian_taxon, total_larval_mass_gaussian_taxon_null)  # model is globally significant

# (4) model summary
summary(total_larval_mass_gaussian_taxon)
model_summary(total_larval_mass_gaussian_taxon, model_name = "Total larval mass", transform_estimate = NULL)
model_forest_plot(total_larval_mass_gaussian_taxon, model_name = "Total larval mass", transform_estimate = NULL)
Anova(total_larval_mass_gaussian_taxon, type = 2)
# confint(profile(total_larval_mass_gaussian_taxon)) %>% view

# (5) emmeans
emmeans_carcass_taxon_total_larval_mass <- emmeans(total_larval_mass_gaussian_taxon, "carcass_taxon")
pairs(regrid(emmeans_carcass_taxon_total_larval_mass))
cld(emmeans_carcass_taxon_total_larval_mass, Letters = letters)

# (6) model visualization
plot_model(total_larval_mass_gaussian_taxon, 
           type = "pred", 
           terms = c("carcass_taxon"))

# (7) write the model results
write_rds(total_larval_mass_gaussian_taxon, "./03_Outputs/Data_Clean/total_larval_mass_gaussian_taxon.rds")


# 4. Average larval mass vs. carcass taxon -------------------------------------
# (1) fit the model
average_larval_mass_gaussian_taxon <- glmmTMB(average_larval_mass ~ carcass_taxon + carcass_weight + male_size + female_size + parent_generation,
                                            data = carcass_data_clean_wild_restricted,
                                            family = "gaussian",
                                            na.action = na.omit)

# (2) model diagnostics
plot(simulateResiduals(average_larval_mass_gaussian_taxon))  # no obvious patterns
check_model(average_larval_mass_gaussian_taxon)  # no obvious patterns

# (3) model significance
average_larval_mass_gaussian_taxon_null <- glmmTMB(average_larval_mass ~ 1,
                                                 data = carcass_data_clean_wild_restricted,
                                                 family = "gaussian",
                                                 na.action = na.omit)

lrtest(average_larval_mass_gaussian_taxon, average_larval_mass_gaussian_taxon_null)  # model is globally significant

# (4) model summary
summary(average_larval_mass_gaussian_taxon)
model_summary(average_larval_mass_gaussian_taxon, model_name = "average larval mass", transform_estimate = NULL)
model_forest_plot(average_larval_mass_gaussian_taxon, model_name = "average larval mass", transform_estimate = NULL)
Anova(average_larval_mass_gaussian_taxon, type = 2)
# confint(profile(average_larval_mass_gaussian_taxon)) %>% view

# (5) emmeans
emmeans_carcass_taxon_average_larval_mass <- emmeans(average_larval_mass_gaussian_taxon, "carcass_taxon")
pairs(regrid(emmeans_carcass_taxon_average_larval_mass))
cld(emmeans_carcass_taxon_average_larval_mass, Letters = letters)

# (6) model visualization
plot_model(average_larval_mass_gaussian_taxon, 
           type = "pred", 
           terms = c("carcass_taxon"))

# (7) write the model results
write_rds(average_larval_mass_gaussian_taxon, "./03_Outputs/Data_Clean/average_larval_mass_gaussian_taxon.rds")


# 5. Proportion of carcass used vs. carcass taxon ------------------------------
### Remove the impossible values
carcass_data_clean_wild_restricted_prop_carcass_used <- carcass_data_clean_wild_restricted %>% 
  filter(prop_carcass_used < 1 & prop_carcass_used > 0)

### Model
# (1) fit the model
prop_carcass_used_beta_taxon <- glmmTMB(prop_carcass_used ~ carcass_taxon + carcass_weight + male_size + female_size + parent_generation,
                                         data = carcass_data_clean_wild_restricted_prop_carcass_used,
                                         family = beta_family("logit"),
                                         na.action = na.omit)

# (2) model diagnostics
plot(simulateResiduals(prop_carcass_used_beta_taxon))  # no pattern
check_model(prop_carcass_used_beta_taxon)  # no pattern

# (3) model significance
prop_carcass_used_beta_taxon_null <- glmmTMB(prop_carcass_used ~ 1,
                                              data = carcass_data_clean_wild_restricted_prop_carcass_used,
                                              family = beta_family("logit"),
                                              na.action = na.omit)

lrtest(prop_carcass_used_beta_taxon, prop_carcass_used_beta_taxon_null)  # model is globally significant

# (5) model summary
summary(prop_carcass_used_beta_taxon)
model_summary(prop_carcass_used_beta_taxon, model_name = "Proportion of carcass used", transform_estimate = "exp")
model_forest_plot(prop_carcass_used_beta_taxon, model_name = "Proportion of carcass used", transform_estimate = "exp")
Anova(prop_carcass_used_beta_taxon, type = 2)
# confint(profile(prop_carcass_used_beta_taxon)) %>% view

# (6) emmeans
emmeans_carcass_taxon_prop_carcass_used <- emmeans(prop_carcass_used_beta_taxon, "carcass_taxon", type = "response")
pairs(regrid(emmeans_carcass_taxon_prop_carcass_used))
cld(emmeans_carcass_taxon_prop_carcass_used, Letters = letters)

# (7) model visualization
plot_model(prop_carcass_used_beta_taxon, 
           type = "pred", 
           terms = c("carcass_taxon"))

# (8) write the model results
write_rds(prop_carcass_used_beta_taxon, "./03_Outputs/Data_Clean/prop_carcass_used_beta_taxon.rds")



