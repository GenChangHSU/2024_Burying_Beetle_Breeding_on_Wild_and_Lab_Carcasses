## -----------------------------------------------------------------------------
## Title: Visualization of the relationships between carcass attributes and beetle breeding outcomes
##
## Author: Gen-Chang Hsu
##
## Date: 2024-04-20
##
## Description:
## 1. Plot the relationship between clutch size vs. carcass weight and carcass type
## 2. Plot the relationship between breeding success vs. carcass weight and carcass type
## 3. Plot the relationship between proportion of eggs developed vs. carcass weight and carcass type
## 4. Plot the relationship between number of larvae vs. carcass weight and carcass type
## 5.1 Plot the relationship between total larval mass vs. carcass weight and carcass type (with zeros)
## 5.2 Plot the relationship between total larval mass vs. carcass weight and carcass type (without zeros)
## 6. Plot the relationship between average larval mass vs. carcass weight and carcass type
## 7. Plot the relationship between larval density vs. carcass weight and carcass type
## 8. Plot the relationship between carcass weight loss vs. carcass weight and carcass type
## 9. Plot the relationship between proportion of carcass used vs. carcass weight and carcass type
## 10. Plot the relationship between average larval mass vs. larval density by carcass type
## 11. Create a multipanel figure for the breeding outcomes
## 12. Create a multipanel figure for the average larval mass and larval density
##
## -----------------------------------------------------------------------------
set.seed(123)


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(sjPlot)
library(grid)
library(ggplotify)
library(patchwork)


# Import files -----------------------------------------------------------------
carcass_data_clean <- read_csv("./03_Outputs/Data_Clean/Carcass_Data_Clean.csv")
clutch_size_zi_nb_quadratic <- read_rds("./03_Outputs/Data_Clean/clutch_size_zi_nb_quadratic.rds")
breeding_success_binomial_quadratic <- read_rds("./03_Outputs/Data_Clean/breeding_success_binomial_quadratic.rds")
prop_eggs_developed_beta_quadratic <- read_rds("./03_Outputs/Data_Clean/prop_eggs_developed_beta_quadratic.rds")
n_larvae_zi_nb_quadratic <- read_rds("./03_Outputs/Data_Clean/n_larvae_zi_nb_quadratic.rds")
total_larval_mass_tweedie_quadratic <- read_rds("./03_Outputs/Data_Clean/total_larval_mass_tweedie_quadratic.rds")
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
                            terms = c("carcass_weight [1:100]", "carcass_type")) +
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
gtable_clutch_size$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_clutch_size$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_clutch_size)

ggsave("./03_Outputs/Figures/Clutch_Size_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 2. Breeding success vs. carcass weight and carcass type ----------------------
p_breeding_success <- plot_model(breeding_success_binomial_quadratic, 
                                 type = "pred", 
                                 terms = c("carcass_weight [1:100]", "carcass_type")) +
  geom_point(data = carcass_data_clean, aes(x = carcass_weight, y = breeding_success, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) +
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.05, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Breeding success (%)", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_breeding_success <- ggplotGrob(p_breeding_success)
gtable_breeding_success$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_breeding_success$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[9]]$gp$col <- "#FFFFFF"
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
                                    terms = c("carcass_weight [1:100]", "carcass_type")) +
  geom_point(data = carcass_data_clean_prop_eggs_developed, aes(x = carcass_weight, y = prop_eggs_developed, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.05, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Hatching success", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_prop_eggs_developed <- ggplotGrob(p_prop_eggs_developed)
gtable_prop_eggs_developed$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_prop_eggs_developed$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_prop_eggs_developed)

ggsave("./03_Outputs/Figures/Prop_Eggs_Developed_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 4. Number of larvae vs. carcass weight and carcass type ----------------------
p_n_larvae <- plot_model(n_larvae_zi_nb_quadratic, 
                         type = "pred", 
                         terms = c("carcass_weight [1:100]", "carcass_type")) +
  geom_point(data = carcass_data_clean, aes(x = carcass_weight, y = n_larvae, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-1, 55), expand = c(0, 0)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Brood size", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_n_larvae <- ggplotGrob(p_n_larvae)
gtable_n_larvae$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_n_larvae$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_n_larvae)

ggsave("./03_Outputs/Figures/N_Larvae_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 5.1 Total larval mass vs. carcass weight and carcass type (with zeros) -------
p_total_larval_mass_with_zeros <- plot_model(total_larval_mass_tweedie_quadratic, 
                                             type = "pred", 
                                             terms = c("carcass_weight [1:100]", "carcass_type")) +
  geom_point(data = carcass_data_clean, aes(x = carcass_weight, y = total_larval_mass, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  # scale_y_continuous(limits = c(-0.2, 12), expand = c(0, 0)) + 
  coord_cartesian(ylim = c(0.45, 12)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Brood mass (g)", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_total_larval_mass_with_zeros <- ggplotGrob(p_total_larval_mass_with_zeros)
gtable_total_larval_mass_with_zeros$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_total_larval_mass_with_zeros$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_total_larval_mass_with_zeros)

ggsave("./03_Outputs/Figures/Total_Larval_Mass_Carcass_Weight_with_Zeros.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 5.2 Total larval mass vs. carcass weight and carcass type (without zeros) ----
### Remove zeros
p_total_larval_mass_without_zeros <- plot_model(total_larval_mass_gaussian_quadratic, 
                                                type = "pred", 
                                                terms = c("carcass_weight [1:100]", "carcass_type")) +
  geom_point(data = filter(carcass_data_clean, total_larval_mass > 0), aes(x = carcass_weight, y = total_larval_mass, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  # scale_y_continuous(limits = c(-0.2, 12), expand = c(0, 0)) + 
  coord_cartesian(ylim = c(0.45, 12)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Brood mass (g)", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_total_larval_mass_without_zeros <- ggplotGrob(p_total_larval_mass_without_zeros)
gtable_total_larval_mass_without_zeros$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_total_larval_mass_without_zeros$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_total_larval_mass_without_zeros)

ggsave("./03_Outputs/Figures/Total_Larval_Mass_Carcass_Weight_without_Zeros.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 6. Average larval mass vs. carcass weight and carcass type -------------------
p_average_larval_mass <- plot_model(average_larval_mass_gaussian_quadratic, 
                                    type = "pred", 
                                    terms = c("carcass_weight [3:80]", "carcass_type")) +
  geom_point(data = carcass_data_clean, aes(x = carcass_weight, y = average_larval_mass, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  # scale_y_continuous(limits = c(0, 0.5), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.01, 0.5)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Average larval mass (g)", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_average_larval_mass <- ggplotGrob(p_average_larval_mass)
gtable_average_larval_mass$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_average_larval_mass$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_average_larval_mass)

ggsave("./03_Outputs/Figures/Average_Larval_Mass_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 7. Larval density vs. carcass weight and carcass type ------------------------
### Remove the point with a larval density of > 3
p_larval_density <- plot_model(larval_density_gaussian_linear, 
                               type = "pred", 
                               terms = c("carcass_weight [3:80]", "carcass_type")) +
  geom_point(data = filter(carcass_data_clean, larval_density < 3), aes(x = carcass_weight, y = larval_density, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  # scale_y_continuous(limits = c(-0.03, 3.2), expand = c(0, 0)) + 
  coord_cartesian(ylim = c(0.1, 2.2)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = expression(paste("Larval density (g"^-1, "carcass)")), color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_larval_density <- ggplotGrob(p_larval_density)
gtable_larval_density$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_larval_density$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_larval_density)

ggsave("./03_Outputs/Figures/Larval_Density_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 8. Carcass weight loss vs. carcass weight and carcass type --------------------------
### Exclude the impossible value, two outliers, and the observations without any larva
carcass_data_clean_carcass_weight_loss <- carcass_data_clean %>% 
  filter(carcass_weight_loss < 25 & carcass_weight_loss > 0) %>% 
  filter(breeding_success == 1)

p_carcass_weight_loss <- plot_model(carcass_weight_loss_gaussian_quadratic, 
                                    type = "pred", 
                                    terms = c("carcass_weight [3:80]", "carcass_type")) +
  geom_point(data = carcass_data_clean_carcass_weight_loss, aes(x = carcass_weight, y = carcass_weight_loss, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.1, 15.2), expand = c(0, 0)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Carcass consumed (g)", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_carcass_weight_loss <- ggplotGrob(p_carcass_weight_loss)
gtable_carcass_weight_loss$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_carcass_weight_loss$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_carcass_weight_loss)

ggsave("./03_Outputs/Figures/Carcass_Weight_Loss_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 9. Proportion of carcass used vs. carcass weight and carcass type ------------
### Remove the impossible values and two outliers
carcass_data_clean_prop_carcass_used <- carcass_data_clean %>% 
  filter(prop_carcass_used < 0.7 & prop_carcass_used > 0)

p_prop_carcass_used <- plot_model(prop_carcass_used_beta_linear, 
                                  type = "pred", 
                                  terms = c("carcass_weight [2:100]", "carcass_type")) +
  geom_point(data = carcass_data_clean_prop_carcass_used, aes(x = carcass_weight, y = prop_carcass_used, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.05, 0.61), expand = c(0, 0),  breaks = seq(0, 0.6, 0.2), labels = seq(0, 60, 20)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Carcass use efficiency (%)", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_prop_carcass_used <- ggplotGrob(p_prop_carcass_used)
gtable_prop_carcass_used$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_prop_carcass_used$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_prop_carcass_used)

ggsave("./03_Outputs/Figures/Proportion_of_Carcass_Used.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 10. Average larval mass vs. larval density -----------------------------------
### Remove the point with a larval density of > 3
p_average_larval_mass_larval_density <- plot_model(average_larval_mass_larval_density_gaussian_linear, 
                                                   type = "pred", 
                                                   terms = c("larval_density [0:2]", "carcass_type")) +
  geom_point(data = filter(carcass_data_clean, larval_density < 3), aes(x = larval_density, y = average_larval_mass, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-0.05, 2.15), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 0.45), expand = c(0, 0)) + 
  labs(title = NULL, x = expression(paste("Larval density (g"^-1, "carcass)")), y = "Average larval mass (g)", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

### Remove the legend key borders
gtable_average_larval_mass_larval_density <- ggplotGrob(p_average_larval_mass_larval_density)
gtable_average_larval_mass_larval_density$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_average_larval_mass_larval_density$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_average_larval_mass_larval_density)

ggsave("./03_Outputs/Figures/Average_Larval_Mass_Larval_Density.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 11. Multipanel figure for the breeding outcomes ------------------------------
### Panel (a)
p_clutch_size_multipanel <- plot_model(clutch_size_zi_nb_quadratic, 
                                       type = "pred", 
                                       terms = c("carcass_weight [1:100]", "carcass_type")) +
  geom_point(data = carcass_data_clean, aes(x = carcass_weight, y = clutch_size, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-1, 82), expand = c(0, 0), breaks = c(0, 20, 40, 60, 80)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Clutch size", color = NULL, subtitle = "(a)") +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"),
        plot.subtitle = element_text(size = 16),
        plot.margin = margin(r = 20, b = 5))

gtable_clutch_size_multipanel <- ggplotGrob(p_clutch_size_multipanel)
gtable_clutch_size_multipanel$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_clutch_size_multipanel$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_clutch_size_multipanel)

### Panel (b)
carcass_data_clean_prop_eggs_developed <- carcass_data_clean %>% 
  mutate(prop_eggs_developed = case_when(prop_eggs_developed >= 1 ~ 0.999,
                                         prop_eggs_developed == 0 ~ 0.001,
                                         TRUE ~ prop_eggs_developed))

p_prop_eggs_developed_multipanel <- plot_model(prop_eggs_developed_beta_quadratic, 
                                               type = "pred", 
                                               terms = c("carcass_weight [1:100]", "carcass_type")) +
  geom_point(data = carcass_data_clean_prop_eggs_developed, aes(x = carcass_weight, y = prop_eggs_developed, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.05, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Hatching success (%)", color = NULL, subtitle = "(b)") +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"),
        axis.title.y = element_text(size = 16, margin = margin(r = 5)),
        plot.subtitle = element_text(size = 16),
        plot.margin = margin(r = 20, b = 5))

### Remove the legend key borders
gtable_prop_eggs_developed_multipanel <- ggplotGrob(p_prop_eggs_developed_multipanel)
gtable_prop_eggs_developed_multipanel$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_prop_eggs_developed_multipanel$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_prop_eggs_developed_multipanel)

### Panel (c)
p_n_larvae_multipanel <- plot_model(n_larvae_zi_nb_quadratic, 
                                    type = "pred", 
                                    terms = c("carcass_weight [1:100]", "carcass_type")) +
  geom_point(data = carcass_data_clean, aes(x = carcass_weight, y = n_larvae, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-1, 52), expand = c(0, 0)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Brood size", color = NULL, subtitle = "(c)") +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"),
        plot.subtitle = element_text(size = 16),
        plot.margin = margin(r = 20, b = 5))

gtable_n_larvae_multipanel <- ggplotGrob(p_n_larvae_multipanel)
gtable_n_larvae_multipanel$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_n_larvae_multipanel$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_n_larvae_multipanel)

### Panel (d)
### Remove zeros
p_total_larval_mass_without_zeros_multipanel <- plot_model(total_larval_mass_gaussian_quadratic, 
                                                           type = "pred", 
                                                           terms = c("carcass_weight [1:100]", "carcass_type")) +
  geom_point(data = filter(carcass_data_clean, total_larval_mass > 0), aes(x = carcass_weight, y = total_larval_mass, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  scale_y_continuous(labels = c(" 0", " 3", " 6", " 9", " 12")) +
  coord_cartesian(ylim = c(0.45, 11.8)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Brood mass (g)", color = NULL, subtitle = "(d)") +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"),
        plot.subtitle = element_text(size = 16),
        axis.title.y = element_text(size = 16, margin = margin(r = 5, l = 3)),
        plot.margin = margin(r = 20, b = 5))

### Remove the legend key borders
gtable_total_larval_mass_without_zeros_multipanel <- ggplotGrob(p_total_larval_mass_without_zeros_multipanel)
gtable_total_larval_mass_without_zeros_multipanel$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_total_larval_mass_without_zeros_multipanel$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_total_larval_mass_without_zeros_multipanel)

(as.ggplot(gtable_clutch_size_multipanel) + as.ggplot(gtable_prop_eggs_developed_multipanel))/
  (as.ggplot(gtable_n_larvae_multipanel) + as.ggplot(gtable_total_larval_mass_without_zeros_multipanel))

ggsave("./03_Outputs/Figures/Breeding_Outcomes.tiff", width = 9.5, height = 8, dpi = 600, device = "tiff")


# 12. Multipanel figure for the average larval mass and larval density ---------
### Panel (a)
p_average_larval_mass_multipanel <- plot_model(average_larval_mass_gaussian_quadratic, 
                                    type = "pred", 
                                    terms = c("carcass_weight [3:80]", "carcass_type")) +
  geom_point(data = carcass_data_clean, aes(x = carcass_weight, y = average_larval_mass, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  # scale_y_continuous(limits = c(0, 0.5), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.01, 0.5)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = "Average larval mass (g)", color = NULL, subtitle = "(a)") +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"),
        plot.subtitle = element_text(size = 16))

### Remove the legend key borders
gtable_average_larval_mass_multipanel <- ggplotGrob(p_average_larval_mass_multipanel)
gtable_average_larval_mass_multipanel$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_average_larval_mass_multipanel$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_average_larval_mass_multipanel)

### Panel (b)
p_larval_density_multipanel <- plot_model(larval_density_gaussian_linear, 
                               type = "pred", 
                               terms = c("carcass_weight [3:80]", "carcass_type")) +
  geom_point(data = filter(carcass_data_clean, larval_density < 3), aes(x = carcass_weight, y = larval_density, color = carcass_type), inherit.aes = F) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 102), expand = c(0, 0)) + 
  # scale_y_continuous(limits = c(-0.03, 3.2), expand = c(0, 0)) + 
  coord_cartesian(ylim = c(0.1, 2.2)) + 
  labs(title = NULL, x = "Carcass weight (g)", y = expression(paste("Larval density (g"^-1, "carcass)")), color = NULL, subtitle = "(b)") +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 1.5, fill = "white"))) + 
  my_ggtheme + 
  theme(legend.position = c(0.85, 0.87),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"),
        plot.subtitle = element_text(size = 16))

### Remove the legend key borders
gtable_larval_density_multipanel <- ggplotGrob(p_larval_density_multipanel)
gtable_larval_density_multipanel$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[5]]$gp$col <- "#FFFFFF"
gtable_larval_density_multipanel$grobs[[15]]$grobs$`99_d874966217baa2d1a56f8468aa8e76ce`$grobs[[9]]$gp$col <- "#FFFFFF"
as.ggplot(gtable_larval_density_multipanel)

as.ggplot(gtable_average_larval_mass_multipanel) + as.ggplot(gtable_larval_density_multipanel)
ggsave("./03_Outputs/Figures/Average_Larval_Mass_Larval_Density.tiff", width = 9.5, height = 4, dpi = 600, device = "tiff")




