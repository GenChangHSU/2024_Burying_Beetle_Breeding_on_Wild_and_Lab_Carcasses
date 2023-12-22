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
  
ggsave("./03_Outputs/Figures/Carcass_Weight_Clutch_Size.tiff", width = 5, height = 4, dpi = 600, device = "tiff")







