## -----------------------------------------------------------------------------
## Title: Visualize the relationships between carcass attributes and beetle breeding outcomes
##
## Author: Gen-Chang Hsu
##
## Date: 2023-12-20
##
## Description:
## 1. 
## 2. 
##
##
##
## -----------------------------------------------------------------------------
set.seed(123)


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(MASS)
library(betareg)
library(mgcv)


# Import files -----------------------------------------------------------------
carcass_data_clean <- read_csv("./03_Outputs/Data_Clean/Carcass_Data_Clean.csv")


# ggplot theme -----------------------------------------------------------------
my_theme <- 
  theme(# Axis
        axis.text.x = element_text(size = 14, color = "black", margin = margin(t = 3)),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 16, margin = margin(t = 10)),
        axis.title.y = element_text(size = 16, margin = margin(r = 8)),
        axis.ticks.length.x = unit(0.18, "cm"),
        axis.ticks.length.y = unit(0.15, "cm"),
                
        # Plot
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        plot.background = element_rect(colour = "transparent"),
        
        # Panel
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        
        # Legend
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
        
        # Facet strip
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 13, hjust = 0.5)
        )
        

############################### Code starts here ###############################

# 1. Clutch size vs. carcass weight and carcass type ---------------------------
ggplot(carcass_data_clean, aes(x = carcass_weight, y = clutch_size)) + 
  geom_point(aes(color = carcass_type)) + 
  geom_smooth(aes(group = carcass_type), color = NA, method = "glm.nb", formula = y ~ poly(x, 2), se = T, show.legend = F) +
  geom_smooth(aes(color = carcass_type), method = "glm.nb", formula = y ~ poly(x, 2), se = F) +
  scale_color_brewer(palette = "Set1", label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 128), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-1, 70), expand = c(0, 0)) + 
  labs(x = "Carcass weight (g)", y = "Clutch size", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 2, fill = "red"))) + 
  my_theme + 
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))
  
ggsave("./03_Outputs/Figures/Clutch_Size_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 2. Breeding success vs. carcass weight and carcass type ----------------------
ggplot(carcass_data_clean, aes(x = carcass_weight, y = breeding_success)) + 
  geom_point(aes(color = carcass_type)) + 
  geom_smooth(aes(group = carcass_type), color = NA, method = "glm", formula = y ~ poly(x, 2), method.args = list(family = "binomial"), se = T, show.legend = F) +
  geom_smooth(aes(color = carcass_type), method = "glm", formula = y ~ poly(x, 2), method.args = list(family = "binomial"), se = F) +
  scale_color_brewer(palette = "Set1", label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 128), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.05, 1.05), expand = c(0, 0)) + 
  labs(x = "Carcass weight (g)", y = "Breeding success", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 2, fill = "red"))) + 
  my_theme + 
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

ggsave("./03_Outputs/Figures/Breeding_Success_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 3. Proportion of eggs developed vs. carcass weight and carcass ---------------
### Convert the zeros to 0.001 and values larger than 1 to 0.999
carcass_data_clean_prop_eggs_developed <- carcass_data_clean %>% 
  mutate(prop_eggs_developed = case_when(prop_eggs_developed >= 1 ~ 0.999,
                                         prop_eggs_developed == 0 ~ 0.001,
                                         TRUE ~ prop_eggs_developed))

ggplot(carcass_data_clean_prop_eggs_developed, aes(x = carcass_weight, y = prop_eggs_developed)) + 
  geom_point(aes(color = carcass_type)) + 
  geom_smooth(aes(group = carcass_type), color = NA, method = "gam", formula = y ~ x, se = T, method.args = list(family = betar(link = "logit")), show.legend = F) +
  geom_smooth(aes(color = carcass_type), method = "gam", formula = y ~ x, method.args = list(family = betar(link = "logit")), se = F) +
  scale_color_brewer(palette = "Set1", label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 128), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.05, 1.05), expand = c(0, 0)) + 
  labs(x = "Carcass weight (g)", y = "Proportion of eggs developed", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 2, fill = "red"))) + 
  my_theme + 
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

ggsave("./03_Outputs/Figures/Prop_Eggs_Developed_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 4. Number of larvae vs. carcass weight and carcass type ----------------------
ggplot(carcass_data_clean, aes(x = carcass_weight, y = n_larvae)) + 
  geom_point(aes(color = carcass_type)) + 
  geom_smooth(aes(group = carcass_type), color = NA, method = "glm.nb", formula = y ~ poly(x, 2), se = T, show.legend = F) +
  geom_smooth(aes(color = carcass_type), method = "glm.nb", formula = y ~ poly(x, 2), se = F) +
  scale_color_brewer(palette = "Set1", label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 128), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-1, 50), expand = c(0, 0)) + 
  labs(x = "Carcass weight (g)", y = "Number of larvae", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 2, fill = "red"))) + 
  my_theme + 
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

ggsave("./03_Outputs/Figures/N_Larvae_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 5. Total larval mass vs. carcass weight and carcass type ---------------------
ggplot(carcass_data_clean, aes(x = carcass_weight, y = total_larval_mass)) + 
  geom_point(aes(color = carcass_type)) + 
  geom_smooth(aes(group = carcass_type), color = NA, method = "lm", formula = y ~ poly(x, 2), se = T, show.legend = F) +
  geom_smooth(aes(color = carcass_type), method = "lm", formula = y ~ poly(x, 2), se = F) +
  scale_color_brewer(palette = "Set1", label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 128), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.1, 12), expand = c(0, 0)) + 
  labs(x = "Carcass weight (g)", y = "Total larval mass (g)", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 2, fill = "red"))) + 
  my_theme + 
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

ggsave("./03_Outputs/Figures/Total_Larval_Mass_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 6. Average larval mass vs. carcass weight and carcass type -------------------
ggplot(carcass_data_clean, aes(x = carcass_weight, y = average_larval_mass)) + 
  geom_point(aes(color = carcass_type)) + 
  geom_smooth(aes(group = carcass_type), color = NA, method = "lm", formula = y ~ x, se = T, show.legend = F) +
  geom_smooth(aes(color = carcass_type), method = "lm", formula = y ~ x, se = F) +
  scale_color_brewer(palette = "Set1", label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 128), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.01, 0.5), expand = c(0, 0)) + 
  labs(x = "Carcass weight (g)", y = "Average larval mass (g)", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 2, fill = "red"))) + 
  my_theme + 
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

ggsave("./03_Outputs/Figures/Average_Larval_Mass_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 7. Larval density vs. carcass weight and carcass type ------------------------
ggplot(carcass_data_clean, aes(x = carcass_weight, y = larval_density)) + 
  geom_point(aes(color = carcass_type)) + 
  geom_smooth(aes(group = carcass_type), color = NA, method = "lm", formula = y ~ x, se = T, show.legend = F) +
  geom_smooth(aes(color = carcass_type), method = "lm", formula = y ~ x, se = F) +
  scale_color_brewer(palette = "Set1", label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 128), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.01, 2.5), expand = c(0, 0)) + 
  labs(x = "Carcass weight (g)", y = "Larval density \n (No. of larvae/gram carcass)", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 2, fill = "red"))) + 
  my_theme + 
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

ggsave("./03_Outputs/Figures/Larval_Density_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 8. Carcass used vs. carcass weight and carcass type --------------------------
ggplot(carcass_data_clean, aes(x = carcass_weight, y = carcass_weight_loss)) + 
  geom_point(aes(color = carcass_type)) + 
  geom_smooth(aes(group = carcass_type), color = NA, method = "lm", formula = y ~ x, se = T, show.legend = F) +
  geom_smooth(aes(color = carcass_type), method = "lm", formula = y ~ x, se = F) +
  scale_color_brewer(palette = "Set1", label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 128), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-1, 65), expand = c(0, 0)) + 
  labs(x = "Carcass weight (g)", y = "Carcass used (g)", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 2, fill = "red"))) + 
  my_theme + 
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

ggsave("./03_Outputs/Figures/Carcass_Weight_Loss_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 9. Carcass use efficiency vs. carcass weight and carcass type ----------------
### Remove the outliers
carcass_data_clean_efficiency <- carcass_data_clean %>% 
  filter(efficiency < 0.5)

ggplot(carcass_data_clean_efficiency, aes(x = carcass_weight, y = efficiency)) + 
  geom_point(aes(color = carcass_type)) + 
  # geom_smooth(aes(group = carcass_type), color = NA, method = "lm", formula = y ~ x, se = T, show.legend = F) +
  geom_smooth(aes(color = carcass_type), method = "gam", formula = y ~ x, method.args = list(family = betar(link = "logit")), se = F, linetype = "dashed") +
  scale_color_brewer(palette = "Set1", label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(-1, 128), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.01, 0.5), expand = c(0, 0)) + 
  labs(x = "Carcass weight (g)", y = "Carcass use efficiency", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 2, fill = "red", linetype = "solid"))) + 
  my_theme + 
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

ggsave("./03_Outputs/Figures/Efficiency_Carcass_Weight.tiff", width = 5, height = 4, dpi = 600, device = "tiff")


# 10. Average larval mass vs. larval density -----------------------------------
ggplot(carcass_data_clean, aes(x = larval_density, y = average_larval_mass)) + 
  geom_point(aes(color = carcass_type)) + 
  geom_smooth(aes(group = carcass_type), color = NA, method = "lm", formula = y ~ x, se = T, show.legend = F) +
  geom_smooth(aes(color = carcass_type), method = "lm", formula = y ~ x, se = F) +
  scale_color_brewer(palette = "Set1", label = c("Lab", "Wild")) + 
  scale_x_continuous(limits = c(0, 2.2), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(-0.01, 0.5), expand = c(0, 0)) + 
  labs(x = "Larval density \n (No. of larvae/gram carcass)", y = "Average larval mass (g)", color = NULL) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 2, fill = "red"))) + 
  my_theme + 
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_blank(),
        legend.key.width = unit(0.3, "in"))

ggsave("./03_Outputs/Figures/Average_Larval_Mass_Larval_Density.tiff", width = 5, height = 4, dpi = 600, device = "tiff")
















