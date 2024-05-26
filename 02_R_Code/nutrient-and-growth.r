####nutrient composition of carcasses and larval growth####
setwd("/Users/sun/Library/CloudStorage/Dropbox/paper/new lab paper/Sun Lab dropbox/data/wild carcass/protein content")

library(glmmTMB)
library(lme4)
library(car)
library(ggplot2)
library(emmeans)
library(multcomp)
library(tidyverse)
library(ggpubr)

# ggplot theme -----------------------------------------------------------------
my_ggtheme <- 
  theme(
    # axis
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
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5), # corrected here
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


# 1. Nutritional composition ----------------------------------------------

#nutrient composition
nutrient<-read.csv("nutrient_for_avg.csv")

#group mouse, bird, snake as wild carcass
nutrient$carc.type <- factor(nutrient$carc.type, levels=c("lab", "mouse", "bird", "snake"))
levels(nutrient$carc.type)[levels(nutrient$carc.type) == "mouse"] <- "mammal"
levels(nutrient$carc.type)[levels(nutrient$carc.type) == "snake"] <- "reptile"
nutrient$tr <- ifelse(nutrient$carc.type == 'lab', 'lab', 'wild')
nutrient_wild <- nutrient[nutrient$carc.type!="lab",]

model_protein <- glmmTMB(prop_protein ~ tr + type + (1|bl/id/rep),
                     data = nutrient,
                     family = beta_family("logit"),
                     na.action = na.omit)

Anova(model_protein,type=2)

a_model_protein <- ggplot(nutrient) + 
  geom_point(aes(x = tr, y = prop_protein*100), alpha = 0.2, position = position_jitter(width = 0.05)) + 
  stat_summary(aes(x = tr, y = prop_protein*100, color = tr), 
               fun = mean, geom = "point", size = 3.5) + 
  stat_summary(aes(x = tr, y = prop_protein*100, color = tr), 
               fun.data = mean_se, geom = "errorbar", width = 0, size = 0.5, show.legend = FALSE) + 
  scale_x_discrete(labels = c("Lab", "Wild")) + 
  scale_y_continuous(limits = c(0, 60), expand = c(0, 0), 
                     breaks = seq(0, 60, by = 10), labels = scales::number_format(scale = 1, suffix = "")) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) +  
  labs(x = "Carcass type", y = "Proportion of protein (%)", subtitle = "(a)") + 
  my_ggtheme + 
  theme(plot.subtitle = element_text(size = 16),
        plot.margin = margin(t = 5, r = 5, l = 10),legend.position = "none")

print(a_model_protein)

model_protein_wild <- glmmTMB(prop_protein ~ carc.type + type + (1|bl/id/rep),
                         data = nutrient_wild,
                         family = beta_family("logit"),
                         na.action = na.omit)
Anova(model_protein_wild,type=2)

model_protein_wild_post=emmeans (model_protein_wild,  ~ carc.type , adjust="tukey")
pairs(model_protein_wild_post)

color_pal <- c("#E69F00", "#009E73", "#8856a7") %>% 
  set_names(nm = c("mammal", "bird", "reptile"))

b_model_protein <- ggplot(nutrient_wild) + 
  geom_point(aes(x = carc.type, y = prop_protein*100), alpha = 0.2, position = position_jitter(width = 0.05)) + 
  stat_summary(aes(x = carc.type, y = prop_protein*100, color = carc.type), 
               fun.data = "mean_se", size = 0.7, linewidth = 1, show.legend = F) + 
  scale_x_discrete(labels = c("Mammal", "Bird", "Reptile")) + 
  scale_y_continuous(limits = c(0, 60), expand = c(0, 0), 
                     breaks = seq(0, 60, by = 10), labels = scales::number_format(scale = 1, suffix = "")) + 
  scale_color_manual(values = color_pal) + 
  labs(x = "Carcass taxon", y = "Proportion of protein (%)", subtitle = "(b)") + 
  my_ggtheme + 
  theme(plot.subtitle = element_text(size = 16),
        plot.margin = margin(t = 5, r = 5, l = 10),legend.position = "none")

# Print the plot
print(b_model_protein)


model_oil <- glmmTMB(prop_oil ~ tr + type + (1|bl/id/rep),
                                      data = nutrient,
                                      family = beta_family("logit"),
                                      ziformula = ~1,
                                      na.action = na.omit)

Anova(model_oil,type=2)

c_model_oil <- ggplot(nutrient) + 
  geom_point(aes(x = tr, y = prop_oil * 100), alpha = 0.2, position = position_jitter(width = 0.05)) + 
  stat_summary(aes(x = tr, y = prop_oil * 100, color = tr), 
               fun = mean, geom = "point", size = 3.5) + 
  stat_summary(aes(x = tr, y = prop_oil * 100, color = tr), 
               fun.data = mean_se, geom = "errorbar", width = 0, size = 0.5, show.legend = FALSE) + 
  scale_x_discrete(labels = c("Lab", "Wild")) + 
  scale_y_continuous(limits = c(0, 60), expand = c(0, 0), 
                     breaks = seq(0, 60, by = 10), labels = scales::number_format(scale = 1, suffix = "")) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), labels = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), labels = c("Lab", "Wild")) +  
  labs(x = "Carcass type", y = "Proportion of fat (%)", subtitle = "(c)") + 
  my_ggtheme + 
  theme(plot.subtitle = element_text(size = 16),
        plot.margin = margin(t = 5, r = 5, l = 10), legend.position = "none")

# Print the plot
print(c_model_oil)


model_oil_wild <- glmmTMB(prop_oil ~ carc.type + type + (1|bl/id/rep),
                     data = nutrient_wild,
                     family = beta_family("logit"),
                     ziformula = ~1,
                     control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")),
                     na.action = na.omit)

Anova(model_oil_wild,type=2)

color_pal <- c("#E69F00", "#009E73", "#8856a7") %>% 
  set_names(nm = c("mammal", "bird", "reptile"))

d_model_oil <- ggplot(nutrient_wild) + 
  geom_point(aes(x = carc.type, y = prop_oil*100), alpha = 0.2, position = position_jitter(width = 0.05)) + 
  stat_summary(aes(x = carc.type, y = prop_oil*100, color = carc.type), 
               fun.data = "mean_se", size = 0.7, linewidth = 1, show.legend = F) + 
  scale_x_discrete(labels = c("Mammal", "Bird", "Reptile")) + 
  scale_y_continuous(limits = c(0, 60), expand = c(0, 0), 
                     breaks = seq(0, 60, by = 10), labels = scales::number_format(scale = 1, suffix = "")) + 
  scale_color_manual(values = color_pal) + 
  labs(x = "Carcass taxon", y = "Proportion of fat (%)", subtitle = "(d)") + 
  my_ggtheme + 
  theme(plot.subtitle = element_text(size = 16),
        plot.margin = margin(t = 5, r = 5, l = 10),legend.position = "none")

# Print the plot
print(d_model_oil)









# 2. Larval growth --------------------------------------------------------



feeding<-read.csv("feeding_exp.csv")
#group mouse, bird, snake as wild carcass
feeding$carc.type <- factor(feeding$carc.type, levels=c("lab", "mouse", "bird", "snake"))
levels(feeding$carc.type)[levels(feeding$carc.type) == "mouse"] <- "mammal"
levels(feeding$carc.type)[levels(feeding$carc.type) == "snake"] <- "reptile"
feeding_wild <- feeding[feeding$carc.type!="lab",]

#lab v. wild
model_feeding <- glmmTMB(mass_2nd ~ tr + type + mass_1st+ (1|bl/family/rep),
                          data = feeding,
                          family = "gaussian",
                          na.action = na.omit)

Anova(model_feeding,type=2)

e_model_feeding <- ggplot(feeding) + 
  geom_point(aes(x = tr, y = mass_2nd), alpha = 0.2, position = position_jitter(width = 0.05)) + 
  stat_summary(aes(x = tr, y = mass_2nd, color = tr), 
               fun = mean, geom = "point", size = 3.5) + 
  stat_summary(aes(x = tr, y = mass_2nd, color = tr), 
               fun.data = mean_se, geom = "errorbar", width = 0, size = 0.5, show.legend = FALSE) + 
  scale_x_discrete(labels = c("Lab", "Wild")) + 
  scale_y_continuous(limits = c(0, 0.15), expand = c(0, 0),breaks = seq(0, 0.15, by = 0.03)) + 
  scale_color_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) + 
  scale_fill_brewer(palette = "Set1", limits = c("lab", "wild"), label = c("Lab", "Wild")) +  
  labs(x = "Carcass type", y = "Larval mass (g)", subtitle = "(e)") + 
  my_ggtheme + 
  theme(plot.subtitle = element_text(size = 16),
        plot.margin = margin(t = 5, r = 5, l = 10),legend.position = "none")

print(e_model_feeding)
#Within wild comparisons
model_feeding_wild <- glmmTMB(mass_2nd ~ carc.type + type + mass_1st + (1|bl/family/rep),
                         data = feeding_wild,
                         family = "gaussian",
                         na.action = na.omit)

Anova(model_feeding_wild,type=2)
model_model_feeding_wild=emmeans (model_feeding_wild,  ~ carc.type , adjust="tukey")
pairwise_comparisons <- pairs(model_model_feeding_wild, adjust = "tukey")
summary(pairwise_comparisons)
cld_results <- cld(model_model_feeding_wild, adjust = "sidak", Letters = letters)
print(cld_results)

color_pal <- c("#E69F00", "#009E73", "#8856a7") %>% 
  set_names(nm = c("mammal", "bird", "reptile"))

f_model_feeding_wild <- ggplot(feeding_wild) + 
  geom_point(aes(x = carc.type, y = mass_2nd), alpha = 0.2, position = position_jitter(width = 0.05)) + 
  stat_summary(aes(x = carc.type, y = mass_2nd, color = carc.type), 
               fun.data = "mean_se", size = 0.7, linewidth = 1, show.legend = F) + 
  scale_x_discrete(labels = c("Mammal", "Bird", "Reptile")) + 
  scale_y_continuous(limits = c(0, 0.15), expand = c(0, 0),breaks = seq(0, 0.15, by = 0.03)) + 
  scale_color_manual(values = color_pal) + 
  labs(x = "Carcass taxon", y = "Larval mass (g)", subtitle = "(f)") + 
  my_ggtheme + 
  theme(plot.subtitle = element_text(size = 16),
        plot.margin = margin(t = 5, r = 5, l = 10),legend.position = "none")

print(f_model_feeding_wild)

ggarrange(a_model_protein,b_model_protein,c_model_oil,d_model_oil,e_model_feeding,f_model_feeding_wild,ncol = 2, nrow = 3,widths=c(1,1),common.legend = FALSE,legend="none")
ggsave("./Fig4.tiff", width = 8, height = 12, dpi = 600, device = "tiff")

#using nutrient content to predict larval mass
#lab v. wild
model_feeding_content <- glmmTMB(mass_2nd ~ avg_prop_protein + avg_prop_oil + type + mass_1st+ (1|bl/family/rep),
                                 data = feeding,
                                 family = "gaussian",
                                 na.action = na.omit)

Anova(model_feeding_content,type=2)

#within wild carcasses comparisons
model_feeding_wild_content <- glmmTMB(mass_2nd ~ avg_prop_protein  + avg_prop_oil + type + mass_1st+ (1|bl/family/rep),
                                      data = feeding_wild,
                                      family = "gaussian",
                                      na.action = na.omit)

Anova(model_feeding_wild_content,type=2)

a_model_prop_protein <- ggplot(feeding) + 
  geom_point(aes(x = avg_prop_protein * 100, y = mass_2nd), alpha = 0.2) + 
  geom_smooth(aes(x = avg_prop_protein * 100, y = mass_2nd), method = "lm", se = FALSE, color = "black", size = 1,linetype="dashed") + 
  scale_x_continuous(limits = c(15, 35), expand = c(0, 0), 
                     breaks = seq(15, 35, by = 5), 
                     labels = function(x) format(x, nsmall = 0, trim = TRUE)) + 
  scale_y_continuous(limits = c(0, 0.15), expand = c(0, 0), 
                     breaks = seq(0, 0.15, by = 0.03)) + 
  labs(x = "Averaged proportion of protein (%)", y = "Larval mass (g)", subtitle = "(a)") + 
  my_ggtheme + 
  theme(plot.subtitle = element_text(size = 16),
        plot.margin = margin(t = 5, r = 15, l = 10),
        legend.position = "none")

print(a_model_prop_protein)


b_model_prop_protein <- ggplot(feeding_wild) + 
  geom_point(aes(x = avg_prop_protein * 100, y = mass_2nd), alpha = 0.2) + 
  geom_smooth(aes(x = avg_prop_protein * 100, y = mass_2nd), method = "lm", se = FALSE, color = "black", size = 1,linetype="dashed") + 
  scale_x_continuous(limits = c(15, 35), expand = c(0, 0), 
                     breaks = seq(15, 35, by = 5), 
                     labels = function(x) format(x, nsmall = 0, trim = TRUE)) + 
  scale_y_continuous(limits = c(0, 0.15), expand = c(0, 0), 
                     breaks = seq(0, 0.15, by = 0.03)) + 
  labs(x = "Averaged proportion of protein (%)", y = "Larval mass (g)", subtitle = "(b)") + 
  my_ggtheme + 
  theme(plot.subtitle = element_text(size = 16),
        plot.margin = margin(t = 5, r = 15, l = 10),
        legend.position = "none")

print(b_model_prop_protein)

#lab v. wild
c_model_prop_fat <- ggplot(feeding) + 
  geom_point(aes(x = avg_prop_oil * 100, y = mass_2nd), alpha = 0.2) + 
  geom_smooth(aes(x = avg_prop_oil * 100, y = mass_2nd), method = "lm", se = FALSE, color = "black", size = 1,linetype="dashed") + 
  scale_x_continuous(limits = c(-0.25, 12), expand = c(0, 0), 
                     breaks = seq(0, 12, by = 3), 
                     labels = function(x) format(x, nsmall = 0, trim = TRUE)) + 
  scale_y_continuous(limits = c(0, 0.15), expand = c(0, 0), 
                     breaks = seq(0, 0.15, by = 0.03)) + 
  labs(x = "Averaged proportion of fat (%)", y = "Larval mass (g)", subtitle = "(c)") + 
  my_ggtheme + 
  theme(plot.subtitle = element_text(size = 16),
        plot.margin = margin(t = 5, r = 15, l = 10),
        legend.position = "none")

print(c_model_prop_fat)

#within wild carcasses comparisons
d_model_prop_fat <- ggplot(feeding_wild) + 
  geom_point(aes(x = avg_prop_oil * 100, y = mass_2nd), alpha = 0.2) + 
  geom_smooth(aes(x = avg_prop_oil * 100, y = mass_2nd), method = "lm", se = FALSE, color = "black", size = 1) + 
  scale_x_continuous(limits = c(-0.25, 12), expand = c(0, 0), 
                     breaks = seq(0, 12, by = 3), 
                     labels = function(x) format(x, nsmall = 0, trim = TRUE)) + 
  scale_y_continuous(limits = c(0, 0.15), expand = c(0, 0), 
                     breaks = seq(0, 0.15, by = 0.03)) + 
  labs(x = "Averaged proportion of fat (%)", y = "Larval mass (g)", subtitle = "(d)") + 
  my_ggtheme + 
  theme(plot.subtitle = element_text(size = 16),
        plot.margin = margin(t = 5, r = 15, l = 10),
        legend.position = "none")

print(d_model_prop_fat)

ggarrange(a_model_prop_protein,b_model_prop_protein,c_model_prop_fat,d_model_prop_fat,ncol = 2, nrow = 2,widths=c(1,1),common.legend = FALSE,legend="none")
ggsave("./FigS2.tiff", width = 10, height = 10, dpi = 600, device = "tiff")
