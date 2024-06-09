## -----------------------------------------------------------------------------
## Title: Analysis of carcass nutritional composition and larval growth 
##
## Author: Syuan-Jyun Sun and Gen-Chang Hsu
##
## Date: 2024-05-25
##
## Description:
## 1. Model the relationship between protein content vs. carcass type
## 2. Model the relationship between protein content vs. wild carcass taxon
## 3. Model the relationship between fat content vs. carcass type
## 4. Model the relationship between fat content vs. wild carcass taxon
## 5. Model the relationship between larval growth vs. carcass type
## 6. Model the relationship between larval growth vs. wild carcass taxon
## 7. Model the relationship between larval growth vs. nutritional content (lab and wild carcasses)
## 8. Model the relationship between larval growth vs. nutritional content (wild carcasses only)
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
nutrition_data_clean <- read_csv("./03_Outputs/Data_Clean/Nutrition_Data_Clean.csv")
larval_growth_data_clean <- read_csv("./03_Outputs/Data_Clean/Larval_Growth_Data_Clean.csv")


############################### Code starts here ###############################

# 1. Protein content of lab vs. wild carcasses ---------------------------------
### Model
# (1) fit the model
prop_protein_beta_carcass_type <- glmmTMB(prop_protein ~ carcass_type + tissue_type + (1|carcass_id),
                                          data = nutrition_data_clean,
                                          family = beta_family("logit"),
                                          na.action = na.omit)

# (2) model diagnostics
plot(simulateResiduals(prop_protein_beta_carcass_type))  # acceptable
check_model(prop_protein_beta_carcass_type)

# (3) model significance
prop_protein_beta_carcass_type_null <- glmmTMB(prop_protein ~ 1,
                                               data = nutrition_data_clean,
                                               family = beta_family("logit"),
                                               na.action = na.omit)

lrtest(prop_protein_beta_carcass_type, prop_protein_beta_carcass_type_null)  # model is globally significant

# (5) model summary
summary(prop_protein_beta_carcass_type)
model_summary(prop_protein_beta_carcass_type, transform_estimate = "exp")
model_forest_plot(prop_protein_beta_carcass_type, transform_estimate = "exp")
Anova(prop_protein_beta_carcass_type, type = 2)

# (6) emmeans
emmeans_prop_protein_carcass_type <- emmeans(prop_protein_beta_carcass_type, "carcass_type", type = "response")
pairs(regrid(emmeans_prop_protein_carcass_type))
cld(emmeans_prop_protein_carcass_type, Letters = letters)

# (7) model visualization
plot_model(prop_protein_beta_carcass_type, 
           type = "pred", 
           terms = c("carcass_type"))

# (8) write the model results
write_rds(prop_protein_beta_carcass_type, "./03_Outputs/Data_Clean/prop_protein_beta_carcass_type.rds")


# 2. Fat content of lab vs. wild carcasses ----------------------------
### Add a small value to avoid zero issue
nutrition_data_clean_fat <- nutrition_data_clean %>% 
  mutate(prop_fat = if_else(prop_fat == 0, 0.00000001, prop_fat))

### Model
# (1) fit the model
prop_fat_beta_carcass_type <- glmmTMB(prop_fat ~ carcass_type + tissue_type + (1|carcass_id),
                                      data = nutrition_data_clean_fat,
                                      family = beta_family("logit"),
                                      na.action = na.omit)

# (2) model diagnostics
plot(simulateResiduals(prop_fat_beta_carcass_type))  # acceptable
check_model(prop_fat_beta_carcass_type)

# (3) model significance
prop_fat_beta_carcass_type_null <- glmmTMB(prop_fat ~ 1,
                                           data = nutrition_data_clean_fat,
                                           family = beta_family("logit"),
                                           na.action = na.omit)

lrtest(prop_fat_beta_carcass_type, prop_fat_beta_carcass_type_null)  # model is globally significant

# (5) model summary
summary(prop_fat_beta_carcass_type)
model_summary(prop_fat_beta_carcass_type, transform_estimate = "exp")
model_forest_plot(prop_fat_beta_carcass_type, transform_estimate = "exp")
Anova(prop_fat_beta_carcass_type, type = 2)

# (6) emmeans
emmeans_prop_fat_carcass_type <- emmeans(prop_fat_beta_carcass_type, "carcass_type", type = "response")
pairs(regrid(emmeans_prop_fat_carcass_type))
cld(emmeans_prop_fat_carcass_type, Letters = letters)

# (7) model visualization
plot_model(prop_fat_beta_carcass_type, 
           type = "pred", 
           terms = c("carcass_type"))

# (8) write the model results
write_rds(prop_fat_beta_carcass_type, "./03_Outputs/Data_Clean/prop_fat_beta_carcass_type.rds")


# 3. Protein content of wild carcass taxa ----------------------------
### Model
# (1) fit the model
prop_protein_beta_carcass_taxon <- glmmTMB(prop_protein ~ carcass_taxon + tissue_type + (1|carcass_id),
                                           data = filter(nutrition_data_clean, carcass_type == "wild"),
                                           family = beta_family("logit"),
                                           na.action = na.omit)

# (2) model diagnostics
plot(simulateResiduals(prop_protein_beta_carcass_taxon))  # no pattern
check_model(prop_protein_beta_carcass_taxon)

# (3) model significance
prop_protein_beta_carcass_taxon_null <- glmmTMB(prop_protein ~ 1,
                                                data = filter(nutrition_data_clean, carcass_type == "wild"),
                                                family = beta_family("logit"),
                                                na.action = na.omit)

lrtest(prop_protein_beta_carcass_taxon, prop_protein_beta_carcass_taxon_null)  # model is globally significant

# (5) model summary
summary(prop_protein_beta_carcass_taxon)
model_summary(prop_protein_beta_carcass_taxon, transform_estimate = "exp")
model_forest_plot(prop_protein_beta_carcass_taxon, transform_estimate = "exp")
Anova(prop_protein_beta_carcass_taxon, type = 2)

# (6) emmeans
emmeans_prop_protein_carcass_taxon <- emmeans(prop_protein_beta_carcass_taxon, "carcass_taxon", type = "response")
pairs(regrid(emmeans_prop_protein_carcass_taxon))
cld(emmeans_prop_protein_carcass_taxon, Letters = letters)

# (7) model visualization
plot_model(prop_protein_beta_carcass_taxon, 
           type = "pred", 
           terms = c("carcass_taxon"))

# (8) write the model results
write_rds(prop_protein_beta_carcass_taxon, "./03_Outputs/Data_Clean/prop_protein_beta_carcass_taxon.rds")


# 4. Fat content of wild carcass taxa ------------------------------------------
### Add a small value to avoid zero issue
nutrition_data_clean_fat <- nutrition_data_clean %>% 
  mutate(prop_fat = if_else(prop_fat == 0, 0.00000001, prop_fat)) %>% 
  mutate(carcass_id = as.factor(carcass_id))

### Model
# (1) fit the model
prop_fat_beta_carcass_taxon <- glmmTMB(prop_fat ~ carcass_taxon + tissue_type + (1|carcass_id),
                                       data = filter(nutrition_data_clean_fat, carcass_type == "wild"),
                                       family = beta_family("logit"),
                                       na.action = na.omit)

# (2) model diagnostics
plot(simulateResiduals(prop_fat_beta_carcass_taxon))  # some patterns
check_model(prop_fat_beta_carcass_taxon)

# (3) model significance
prop_fat_beta_carcass_taxon_null <- glmmTMB(prop_fat ~ 1,
                                            data = filter(nutrition_data_clean_fat, carcass_type == "wild"),
                                            family = beta_family("logit"),
                                            na.action = na.omit)

lrtest(prop_fat_beta_carcass_taxon, prop_fat_beta_carcass_taxon_null)  # model is globally significant

# (5) model summary
summary(prop_fat_beta_carcass_taxon)
model_summary(prop_fat_beta_carcass_taxon, transform_estimate = "exp")
model_forest_plot(prop_fat_beta_carcass_taxon, transform_estimate = "exp")
Anova(prop_fat_beta_carcass_taxon, type = 2)

# (6) emmeans
emmeans_prop_fat_carcass_taxon <- emmeans(prop_fat_beta_carcass_taxon, "carcass_taxon", type = "response")
pairs(regrid(emmeans_prop_fat_carcass_taxon))
cld(emmeans_prop_fat_carcass_taxon, Letters = letters)

# (7) model visualization
plot_model(prop_fat_beta_carcass_taxon, 
           type = "pred", 
           terms = c("carcass_taxon"))

# (8) write the model results
write_rds(prop_fat_beta_carcass_taxon, "./03_Outputs/Data_Clean/prop_fat_beta_carcass_taxon.rds")


# 5. Larval growth on lab vs. wild carcasses -----------------------------------
### Model
# (1) fit the model
larval_growth_gaussian_carcass_type <- glmmTMB(larval_weight_gain_g ~ carcass_type + tissue_type + initial_larval_mass_g + (1|block_id/carcass_id) + (1|family_id),
                                               data = larval_growth_data_clean,
                                               family = "gaussian",
                                               na.action = na.omit)

# (2) model diagnostics
plot(simulateResiduals(larval_growth_gaussian_carcass_type))  # no patterns
check_model(larval_growth_gaussian_carcass_type)

# (3) model significance
larval_growth_gaussian_carcass_type_null <- glmmTMB(larval_weight_gain_g ~ 1,
                                                    data = larval_growth_data_clean,
                                                    family = "gaussian",
                                                    na.action = na.omit)

lrtest(larval_growth_gaussian_carcass_type, larval_growth_gaussian_carcass_type_null)  # model is globally significant

# (5) model summary
summary(larval_growth_gaussian_carcass_type)
model_summary(larval_growth_gaussian_carcass_type, transform_estimate = NULL)
model_forest_plot(larval_growth_gaussian_carcass_type, transform_estimate = NULL)
Anova(larval_growth_gaussian_carcass_type, type = 2)

# (6) emmeans
emmeans_prop_protein_carcass_type <- emmeans(larval_growth_gaussian_carcass_type, "carcass_type", type = "response")
pairs(regrid(emmeans_prop_protein_carcass_type))
cld(emmeans_prop_protein_carcass_type, Letters = letters)

# (7) model visualization
plot_model(larval_growth_gaussian_carcass_type, 
           type = "pred", 
           terms = c("carcass_type"))

# (8) write the model results
write_rds(larval_growth_gaussian_carcass_type, "./03_Outputs/Data_Clean/larval_growth_gaussian_carcass_type.rds")


# 6. Larval growth on wild carcass taxa ----------------------------------------
### Model
# (1) fit the model
larval_growth_gaussian_carcass_taxon <- glmmTMB(larval_weight_gain_g ~ carcass_taxon + tissue_type + initial_larval_mass_g + (1|block_id/carcass_id) + (1|family_id),
                                                data = filter(larval_growth_data_clean, carcass_type == "wild"),
                                                family = "gaussian",
                                                na.action = na.omit)

# (2) model diagnostics
plot(simulateResiduals(larval_growth_gaussian_carcass_taxon))  # some patterns
check_model(larval_growth_gaussian_carcass_taxon)

# (3) model significance
larval_growth_gaussian_carcass_taxon_null <- glmmTMB(larval_weight_gain_g ~ 1,
                                                     data = filter(larval_growth_data_clean, carcass_type == "wild"),
                                                     family = "gaussian",
                                                     na.action = na.omit)

lrtest(larval_growth_gaussian_carcass_taxon, larval_growth_gaussian_carcass_taxon_null)  # model is globally significant

# (5) model summary
summary(larval_growth_gaussian_carcass_taxon)
model_summary(larval_growth_gaussian_carcass_taxon, transform_estimate = NULL)
model_forest_plot(larval_growth_gaussian_carcass_taxon, transform_estimate = NULL)
Anova(larval_growth_gaussian_carcass_taxon, type = 2)

# (6) emmeans
emmeans_prop_protein_carcass_taxon <- emmeans(larval_growth_gaussian_carcass_taxon, "carcass_taxon", type = "response")
pairs(regrid(emmeans_prop_protein_carcass_taxon))
cld(emmeans_prop_protein_carcass_taxon, Letters = letters)

# (7) model visualization
plot_model(larval_growth_gaussian_carcass_taxon, 
           type = "pred", 
           terms = c("carcass_taxon"))

# (8) write the model results
write_rds(larval_growth_gaussian_carcass_taxon, "./03_Outputs/Data_Clean/larval_growth_gaussian_carcass_taxon.rds")


# 7. Larval growth vs. nutritional content for lab and wild carcasses ----------
### Model
# (1) fit the model
larval_growth_gaussian_nutrition_all <- glmmTMB(larval_weight_gain_g ~ mean_prop_protein + mean_prop_fat + initial_larval_mass_g + (1|family_id),
                                                data = larval_growth_data_clean,
                                                family = "gaussian",
                                                na.action = na.omit)

# larval_growth_gaussian_nutrition_all_with_interaction <- glmmTMB(larval_weight_gain_g ~ mean_prop_protein * mean_prop_fat + initial_larval_mass_g + (1|family_id),
#                                                 data = larval_growth_data_clean,
#                                                 family = "gaussian",
#                                                 na.action = na.omit)
# 
# AIC(larval_growth_gaussian_nutrition_all, larval_growth_gaussian_nutrition_all_with_interaction)
# lrtest(larval_growth_gaussian_nutrition_all, larval_growth_gaussian_nutrition_all_with_interaction)

# (2) model diagnostics
plot(simulateResiduals(larval_growth_gaussian_nutrition_all))  # no pattern
check_model(larval_growth_gaussian_nutrition_all)

# (3) model significance
larval_growth_gaussian_nutrition_null_all <- glmmTMB(larval_weight_gain_g ~ 1,
                                                     data = larval_growth_data_clean,
                                                     family = "gaussian",
                                                     na.action = na.omit)

lrtest(larval_growth_gaussian_nutrition_all, larval_growth_gaussian_nutrition_null_all)  # model is not globally significant

# (5) model summary
summary(larval_growth_gaussian_nutrition_all)
model_summary(larval_growth_gaussian_nutrition_all, transform_estimate = NULL)
model_forest_plot(larval_growth_gaussian_nutrition_all, transform_estimate = NULL)
Anova(larval_growth_gaussian_nutrition_all, type = 2)

# (6) model visualization
plot_model(larval_growth_gaussian_nutrition_all, 
           type = "pred", 
           terms = c("mean_prop_protein"))

plot_model(larval_growth_gaussian_nutrition_all, 
           type = "pred", 
           terms = c("mean_prop_fat"))

# (7) write the model results
write_rds(larval_growth_gaussian_nutrition_all, "./03_Outputs/Data_Clean/larval_growth_gaussian_nutrition_all.rds")


# 8. Larval growth vs. nutritional content for wild carcasses ------------------
### Model
# (1) fit the model
larval_growth_gaussian_nutrition_wild <- glmmTMB(larval_weight_gain_g ~ mean_prop_protein + mean_prop_fat + tissue_type + initial_larval_mass_g + (1|family_id),
                                                 data = filter(larval_growth_data_clean, carcass_type == "wild"),
                                                 family = "gaussian",
                                                 na.action = na.omit)

larval_growth_gaussian_nutrition_wild_with_interaction <- glmmTMB(larval_weight_gain_g ~ mean_prop_protein * mean_prop_fat + tissue_type + initial_larval_mass_g + (1|family_id),
                                                 data = filter(larval_growth_data_clean, carcass_type == "wild"),
                                                 family = "gaussian",
                                                 na.action = na.omit)

lrtest(larval_growth_gaussian_nutrition_wild, larval_growth_gaussian_nutrition_wild_with_interaction)  # interaction term not significant
AIC(larval_growth_gaussian_nutrition_wild, larval_growth_gaussian_nutrition_wild_with_interaction)  # model without interaction is better

# (2) model diagnostics
plot(simulateResiduals(larval_growth_gaussian_nutrition_wild))  # no pattern
check_model(larval_growth_gaussian_nutrition_wild)

# (3) model significance
larval_growth_gaussian_nutrition_wild_null <- glmmTMB(larval_weight_gain_g ~ 1,
                                                      data = filter(larval_growth_data_clean, carcass_type == "wild"),
                                                      family = "gaussian",
                                                      na.action = na.omit)

lrtest(larval_growth_gaussian_nutrition_wild, larval_growth_gaussian_nutrition_wild_null)  # model is globally significant

# (5) model summary
summary(larval_growth_gaussian_nutrition_wild)
model_summary(larval_growth_gaussian_nutrition_wild, transform_estimate = NULL)
model_forest_plot(larval_growth_gaussian_nutrition_wild, transform_estimate = NULL)
Anova(larval_growth_gaussian_nutrition_wild, type = 2)

# (6) model visualization
plot_model(larval_growth_gaussian_nutrition_wild, 
           type = "pred", 
           terms = c("mean_prop_protein"))

plot_model(larval_growth_gaussian_nutrition_wild, 
           type = "pred", 
           terms = c("mean_prop_fat"))

# (7) write the model results
write_rds(larval_growth_gaussian_nutrition_wild, "./03_Outputs/Data_Clean/larval_growth_gaussian_nutrition_wild.rds")

