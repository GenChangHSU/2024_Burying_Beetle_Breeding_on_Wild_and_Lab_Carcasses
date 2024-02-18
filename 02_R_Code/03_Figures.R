## -----------------------------------------------------------------------------
## Title: Visualization of the relationships between carcass attributes and beetle breeding outcomes
##
## Author: Gen-Chang Hsu
##
## Date: 2024-02-17
##
## Description:
## 1. Plot the relationship between clutch size vs. carcass weight and carcass type
## 2. Plot the relationship between breeding success vs. carcass weight and carcass type
## 3. Plot the relationship between proportion of eggs developed vs. carcass weight and carcass type
## 4. Plot the relationship between number of larvae vs. carcass weight and carcass type
## 5. Plot the relationship between total larval mass vs. carcass weight and carcass type
## 6. Plot the relationship between average larval mass vs. carcass weight and carcass type
## 7. Plot the relationship between larval density vs. carcass weight and carcass type
## 8. Plot the relationship between carcass weight loss vs. carcass weight and carcass type
## 9. Plot the relationship between proportion of carcass used vs. carcass weight and carcass type
## 10. Plot the relationship between average larval mass vs. larval density by carcass type
##
## -----------------------------------------------------------------------------
set.seed(123)


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(sjPlot)
library(grid)
library(ggplotify)


# Import files -----------------------------------------------------------------
carcass_data_clean <- read_csv("./03_Outputs/Data_Clean/Carcass_Data_Clean.csv")
clutch_size_zi_nb_quadratic <- read_rds("./03_Outputs/Data_Clean/clutch_size_zi_nb_quadratic.rds")
breeding_success_logistic_quadratic <- read_rds("./03_Outputs/Data_Clean/breeding_success_logistic_quadratic.rds")
prop_eggs_developed_beta_quadratic <- read_rds("./03_Outputs/Data_Clean/prop_eggs_developed_beta_quadratic.rds")
n_larvae_zi_nb_quadratic <- read_rds("./03_Outputs/Data_Clean/n_larvae_zi_nb_quadratic.rds")
total_larval_mass_gaussian_quadratic <- read_rds("./03_Outputs/Data_Clean/total_larval_mass_gaussian_quadratic.rds")
average_larval_mass_gaussian_quadratic <- read_rds("./03_Outputs/Data_Clean/average_larval_mass_gaussian_quadratic.rds")
larval_density_gaussian_linear <- read_rds("./03_Outputs/Data_Clean/larval_density_gaussian_linear.rds")
carcass_weight_loss_gaussian_quadratic <- read_rds("./03_Outputs/Data_Clean/carcass_weight_loss_gaussian_quadratic.rds")
prop_carcass_used_beta_linear <- read_rds("./03_Outputs/Data_Clean/prop_carcass_used_beta_linear.rds")
average_larval_mass_larval_density_gaussian_linear <- read_rds("./03_Outputs/Data_Clean/average_larval_mass_larval_density_gaussian_linear.rds")


# ggplot theme -----------------------------------------------------------------
my_ggtheme <- 
  theme(# axis
        axis.text.x = element_text(size = 14, color = "black", margin = margin(t = 3)),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 16, margin = margin(t = 10)),
        axis.title.y = element_text(size = 16, margin = margin(r = 8)),
        axis.ticks.length.x = unit(0.18, "cm"),
        axis.ticks.length.y = unit(0.15, "cm"),
                
        # plot
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.background = element_rect(colour = "transparent"),
        
        # panel
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        
        # legend
        legend.position = "right",
        legend.spacing.x = unit(0.2, "cm"),
        legend.spacing.y = unit(0.2, "cm"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.size = unit(0.5, "line"),
        legend.key = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.box.just = "center",
        legend.justification = c(0.5, 0.5),
        legend.title.align = 0.5,
        legend.background = element_rect(fill = "transparent", size = 0.25, linetype = "solid", colour = "black"),
        
        # facet strip
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13, hjust = 0.5)
        )

        
############################### Code starts here ###############################

### Convert the variable "parent_generation" to a factor
carcass_data_clean <- carcass_data_clean %>% 
  mutate(parent_generation = as.factor(parent_generation))

### Exclude wild carcasses larger than 100 grams
carcass_data_clean <- carcass_data_clean %>% 
  filter(carcass_weight <= 100)

# 1. Clutch size vs. carcass weight and carcass type ---------------------------
### A scatterplot with model fitted lines
p_clutch_size <- plot_model(clutch_size_zi_nb_quadratic, 
                            type = "pred", 
                            terms = c("carcass_weight [0:100]", "carcass_type")) +
  geom_point(data = carcass_data_clean, aes(x = carcass_weight, y = clutch_size, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-1, 75), expand = c(0, 0)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Clutch size", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_clutch_size <- ggplotGrob(p_clutch_size)
gtable_clutch_size$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_clutch_size$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_clutch_size)

ggsave("./03_Outputs/Figures/Clutch_Size_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 2. Breeding success vs. carcass weight and carcass type ----------------------
p_breeding_success <- plot_model(breeding_success_logistic_quadratic, 
                                 type = "pred", 
                                 terms = c("carcass_weight [0:100]", "carcass_type")) +
  geom_point(data = carcass_data_clean, aes(x = carcass_weight, y = breeding_success, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) +
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.05, 1.05), expand = c(0, 0), labels = str_glue("{seq(0, 100, 25)}%")) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Probability of \n breeding success", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_breeding_success <- ggplotGrob(p_breeding_success)
gtable_breeding_success$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_breeding_success$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_breeding_success)

ggsave("./03_Outputs/Figures/Breeding_Success_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 3. Proportion of eggs developed vs. carcass weight and carcass ---------------
### Convert the zeros to 0.001 and values larger than 1 to 0.999
carcass_data_clean_prop_eggs_developed <- carcass_data_clean %>% 
  mutate(prop_eggs_developed = case_when(prop_eggs_developed >= 1 ~ 0.999,
                                         prop_eggs_developed == 0 ~ 0.001,
                                         TRUE ~ prop_eggs_developed))

p_prop_eggs_developed <- plot_model(prop_eggs_developed_beta_quadratic, 
                                 type = "pred", 
                                 terms = c("carcass_weight [0:100]", "carcass_type")) +
  geom_point(data = carcass_data_clean_prop_eggs_developed, aes(x = carcass_weight, y = prop_eggs_developed, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.05, 1.05), expand = c(0, 0)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Proportion of eggs developed", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_prop_eggs_developed <- ggplotGrob(p_prop_eggs_developed)
gtable_prop_eggs_developed$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_prop_eggs_developed$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_prop_eggs_developed)

ggsave("./03_Outputs/Figures/Prop_Eggs_Developed_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 4. Number of larvae vs. carcass weight and carcass type ----------------------
p_n_larvae <- plot_model(n_larvae_zi_nb_quadratic, 
                         type = "pred", 
                         terms = c("carcass_weight [0:100]", "carcass_type")) +
  geom_point(data = carcass_data_clean, aes(x = carcass_weight, y = n_larvae, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-1, 55), expand = c(0, 0)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Number of larvae", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_n_larvae <- ggplotGrob(p_n_larvae)
gtable_n_larvae$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_n_larvae$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_n_larvae)

ggsave("./03_Outputs/Figures/N_Larvae_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 5. Total larval mass vs. carcass weight and carcass type ---------------------
p_total_larval_mass <- plot_model(total_larval_mass_gaussian_quadratic, 
                         type = "pred", 
                         terms = c("carcass_weight [0:100]", "carcass_type")) +
  geom_point(data = carcass_data_clean, aes(x = carcass_weight, y = total_larval_mass, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.2, 12), expand = c(0, 0)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Total larval mass (g)", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_total_larval_mass <- ggplotGrob(p_total_larval_mass)
gtable_total_larval_mass$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_total_larval_mass$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_total_larval_mass)

ggsave("./03_Outputs/Figures/Total_Larval_Mass_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 6. Average larval mass vs. carcass weight and carcass type -------------------
p_average_larval_mass <- plot_model(average_larval_mass_gaussian_quadratic, 
                                    type = "pred", 
                                    terms = c("carcass_weight [0:100]", "carcass_type")) +
  geom_point(data = carcass_data_clean, aes(x = carcass_weight, y = average_larval_mass, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 0.5), expand = c(0, 0)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Average larval mass (g)", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_average_larval_mass <- ggplotGrob(p_average_larval_mass)
gtable_average_larval_mass$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_average_larval_mass$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_average_larval_mass)

ggsave("./03_Outputs/Figures/Average_Larval_Mass_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 7. Larval density vs. carcass weight and carcass type ------------------------
p_larval_density <- plot_model(larval_density_gaussian_linear, 
                               type = "pred", 
                               terms = c("carcass_weight [0:100]", "carcass_type")) +
  geom_point(data = carcass_data_clean, aes(x = carcass_weight, y = larval_density, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.03, 3.2), expand = c(0, 0)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Larval density \n (number per gram carcass)", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_larval_density <- ggplotGrob(p_larval_density)
gtable_larval_density$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_larval_density$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_larval_density)

ggsave("./03_Outputs/Figures/Larval_Density_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 8. Carcass weight loss vs. carcass weight and carcass type --------------------------
### Exclude the impossible value, two outliers, and the observations without any larva
carcass_data_clean_carcass_weight_loss <- carcass_data_clean %>% 
  filter(carcass_weight_loss < 25 & carcass_weight_loss > 0) %>% 
  filter(breeding_success == 1)

p_carcass_weight_loss <- plot_model(carcass_weight_loss_gaussian_quadratic, 
                               type = "pred", 
                               terms = c("carcass_weight [0:100]", "carcass_type")) +
  geom_point(data = carcass_data_clean_carcass_weight_loss, aes(x = carcass_weight, y = carcass_weight_loss, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.1, 15.2), expand = c(0, 0)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Carcass weight loss (g)", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_carcass_weight_loss <- ggplotGrob(p_carcass_weight_loss)
gtable_carcass_weight_loss$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_carcass_weight_loss$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_carcass_weight_loss)

ggsave("./03_Outputs/Figures/Carcass_Weight_Loss_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 9. Proportion of carcass used vs. carcass weight and carcass type ------------
### Remove the impossible values and two outliers
carcass_data_clean_prop_carcass_used <- carcass_data_clean %>% 
  filter(prop_carcass_used < 0.7 & prop_carcass_used > 0)

p_prop_carcass_used <- plot_model(prop_carcass_used_beta_linear, 
                                  type = "pred", 
                                  terms = c("carcass_weight [0:100]", "carcass_type")) +
  geom_point(data = carcass_data_clean_prop_carcass_used, aes(x = carcass_weight, y = prop_carcass_used, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.05, 0.61), expand = c(0, 0)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Proportion of carcass used", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_prop_carcass_used <- ggplotGrob(p_prop_carcass_used)
gtable_prop_carcass_used$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_prop_carcass_used$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_prop_carcass_used)

ggsave("./03_Outputs/Figures/Proportion_of_Carcass_Used.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 10. Average larval mass vs. larval density -----------------------------------
p_average_larval_mass_larval_density <- plot_model(average_larval_mass_larval_density_gaussian_linear, 
                                  type = "pred", 
                                  terms = c("larval_density [0:2]", "carcass_type")) +
  geom_point(data = carcass_data_clean, aes(x = larval_density, y = average_larval_mass, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-0.05, 3.15), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 0.45), expand = c(0, 0)) + 
  labs(title = NULL, x = "Larval density", y = "Average larval mass", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_average_larval_mass_larval_density <- ggplotGrob(p_average_larval_mass_larval_density)
gtable_average_larval_mass_larval_density$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_average_larval_mass_larval_density$grobs[[15]]$grobs$`99_3f99c453c8478605472e33507cf97de4`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_average_larval_mass_larval_density)

ggsave("./03_Outputs/Figures/Average_Larval_Mass_Larval_Density.tiff", width = 5, height = 4, dpi = 600, device = "tiff")

