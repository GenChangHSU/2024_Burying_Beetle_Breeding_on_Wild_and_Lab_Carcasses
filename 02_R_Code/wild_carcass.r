####this is the code for wild carcass analysis####
setwd("/Users/sun/Library/CloudStorage/Dropbox/paper/new lab paper/Sun Lab dropbox/data/wild carcass")
setwd("/Users/syuan-jyunsun/Library/CloudStorage/Dropbox/paper/new lab paper/Sun Lab dropbox/data/wild carcass")
library(lme4)
library(car)
library(emmeans)
library(MASS)
library(ggplot2)
library(dplyr)
library(ggpubr)

install.packages("export")
library(export)

graph2office(x=success, 
             file="fig", 
             type = c("PPT"),  
             width = 6.75, 
             height =5)

data=read.csv("wild_carcass.csv")
newdata=data[data$note=="1",]
newdata=newdata[newdata$class_bird!="others",]

newdata=newdata[newdata$efficiency_new_method<1,]
newdata=newdata[newdata$clutch_size>0,]

newdatawild=newdata[newdata$tr=="wild",]
newdatalab=newdata[newdata$tr=="lab",]

#test if wild carcass has less meat compared to the lab, 
model=glmer(carcass_used~poly(carc_wt,degree=2)*tr+(1|pair),gaussian,data=newdata)
Anova(model,type=3)

carcass_use1 <- ggplot(newdata, aes(x=carc_wt, y=carcass_used)) +
  geom_point(aes(color=tr), size=2, alpha=0.4) +
  geom_smooth(data=subset(newdata, tr=="lab"),
              aes(fill="lab", color="lab", y=carcass_used, x=carc_wt),
              formula=y ~ poly(x, 2),  # Include polynomial term
              method="glm",
              method.args=list(family="gaussian"),
              se=F,
              alpha=0.3,
              size=1) +
  geom_smooth(data=subset(newdata, tr=="wild"),
              aes(fill="wild", color="wild", y=carcass_used, x=carc_wt),
              formula=y ~ poly(x, 2),  # Include polynomial term
              method="glm",
              method.args=list(family="gaussian"),
              se=F,
              alpha=0.3,
              size=1) +
  scale_color_manual(values=c("lab"="#d53e4f", "wild"="#3288bd")) +
  scale_fill_manual(values=c("lab"="#d53e4f", "wild"="#3288bd")) +
  guides(color=guide_legend(override.aes=list(fill=NA)), fill=FALSE) +
  theme_classic() +
  scale_y_continuous(limits=c(0, 35)) +
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=14),
    legend.title=element_blank(),
    legend.key.size=unit(2,"lines"),
    legend.text=element_text(size=12),
    legend.background=element_blank(),
    legend.key=element_rect(fill=NA, colour=NA)
  )
carcass_use1

ggarrange(carcass_use1,carcass_use, labels= c("A", "B"),ncol = 2, nrow = 1,widths=c(1,1),common.legend = TRUE,legend="bottom")


model=glmer(tot_mass~carcass_used*tr+(1|pair),gaussian,data=newdata)
Anova(model,type=3)

model=glmer(tot_mass~carcass_used+tr+(1|pair),gaussian,data=newdata)
summary(model)

carcass_use <- ggplot(newdata, aes(x=carcass_used, y=tot_mass)) +
  geom_point(aes(color=tr), size=2, alpha=0.4) +
  geom_smooth(data=subset(newdata, tr=="lab"),
              aes(fill="lab", color="lab", y=tot_mass, x=carcass_used),
              formula=y ~ poly(x, 1),  # Include polynomial term
              method="glm",
              method.args=list(family="gaussian"),
              se=F,
              alpha=0.3,
              size=1) +
  geom_smooth(data=subset(newdata, tr=="wild"),
              aes(fill="wild", color="wild", y=tot_mass, x=carcass_used),
              formula=y ~ poly(x, 1),  # Include polynomial term
              method="glm",
              method.args=list(family="gaussian"),
              se=F,
              alpha=0.3,
              size=1) +
  scale_color_manual(values=c("lab"="#d53e4f", "wild"="#3288bd")) +
  scale_fill_manual(values=c("lab"="#d53e4f", "wild"="#3288bd")) +
  guides(color=guide_legend(override.aes=list(fill=NA)), fill=FALSE) +
  theme_classic() +
  scale_y_continuous(limits=c(0, 10)) +
  scale_x_continuous(limits=c(0, 35)) +
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=14),
    legend.title=element_blank(),
    legend.key.size=unit(2,"lines"),
    legend.text=element_text(size=12),
    legend.background=element_blank(),
    legend.key=element_rect(fill=NA, colour=NA)
  )
carcass_use

model=glmer(carc_wt~tr+(1|parent_generation/pair),gaussian,data=newdata)
Anova(model,type=3)

model=glmer(success~carc_wt*tr+(1|parent_generation/pair),binomial,data=newdata)
model2=glmer(success~poly(carc_wt,degree=2)*tr+(1|parent_generation/pair),binomial,data=newdata)
model3=glmer(success~poly(carc_wt,degree=3)*tr+(1|parent_generation/pair),binomial,data=newdata)

anova(model,model2)
Anova(model2,type=3)

modellab=glmer(success~carc_wt+(1|parent_generation/pair),binomial,data=newdata[newdata$tr=="lab",])
modelwild=glmer(success~carc_wt+(1|parent_generation/pair),binomial,data=newdata[newdata$tr=="wild",])

Anova(modellab,type=3)
Anova(modelwild,type=3)

modelwild=glmer(success~carc_wt+class+(1|parent_generation),binomial,data=newdata[newdata$tr=="wild",])
Anova(modelwild,type=3)

success <- ggplot(newdata, aes(x=carc_wt, y=success)) +
  geom_point(aes(color=tr), size=2, alpha=0.4) +
  geom_smooth(data=subset(newdata, tr=="lab"),
              aes(fill="lab", color="lab"),  # Ensure color is also set here for the line
              method=glm,
              formula=y ~ poly(x, 2),
              method.args=list(family=binomial(link="logit")),
              se=F,
              alpha=0.3,
              size=1) +
  geom_smooth(data=subset(newdata, tr=="wild"),
              aes(fill="wild", color="wild"),  # Ensure color is also set here for the line
              method=glm,
              formula=y ~ poly(x, 2),
              method.args=list(family=binomial(link="logit")),
              se=F,
              alpha=0.3,
              size=1) +
  scale_color_manual(values=c("lab"="#d53e4f", "wild"="#3288bd")) +
  scale_fill_manual(values=c("lab"="#d53e4f", "wild"="#3288bd")) +
  labs(x = "Carcass mass (g)", y = "Probability of breeding success") +
  guides(color=guide_legend(override.aes=list(fill=NA)), fill=FALSE) +  # Use fill=FALSE to remove the fill legend
  theme_classic() +
  scale_y_continuous(limits=c(0, 1)) +
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=14),
    legend.title=element_blank(),
    legend.key.size=unit(2,"lines"),
    legend.text=element_text(size=12),
    legend.background=element_blank(), 
    legend.key=element_rect(fill=NA, colour=NA)
  )
success

clutchsize <- ggplot(newdata, aes(x=carc_wt, y=clutch_size)) +
  geom_point(aes(color=tr), size=2, alpha=0.4) +
  geom_smooth(data=subset(newdata, tr=="lab"),
              aes(fill="lab", color="lab", y=larvae, x=carc_wt),
              formula=y ~ poly(x, 2),  # Include polynomial term
              method="glm",
              method.args=list(family="poisson"),
              se=F,
              alpha=0.3,
              size=1) +
  geom_smooth(data=subset(newdata, tr=="wild"),
              aes(fill="wild", color="wild", y=larvae, x=carc_wt),
              formula=y ~ poly(x, 2),  # Include polynomial term
              method="glm",
              method.args=list(family="poisson"),
              se=F,
              alpha=0.3,
              size=1) +
  scale_color_manual(values=c("lab"="#d53e4f", "wild"="#3288bd")) +
  scale_fill_manual(values=c("lab"="#d53e4f", "wild"="#3288bd")) +
  guides(color=guide_legend(override.aes=list(fill=NA)), fill=FALSE) +
  theme_classic() +
  scale_y_continuous(limits=c(0, 80)) +
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=14),
    legend.title=element_blank(),
    legend.key.size=unit(2,"lines"),
    legend.text=element_text(size=12),
    legend.background=element_blank(),
    legend.key=element_rect(fill=NA, colour=NA)
  )
clutchsize

# Fit a series of models with increasing polynomial degree
model1 = glmer(larvae ~ carc_wt * tr + (1 | parent_generation), family=poisson, data=newdata)
model2 = glmer(larvae ~ poly(carc_wt, 2) * tr + (1 | parent_generation), family=poisson, data=newdata)
model3 = glmer(larvae ~ poly(carc_wt, 3) * tr + (1 | parent_generation), family=poisson, data=newdata)
# ...and so on for higher degrees if necessary

# Compare models using AIC
AIC(model1, model2, model3)

# Compare models using BIC
BIC(model1, model2, model3)

# Perform likelihood ratio tests for nested models
anova(model1, model2)
anova(model2, model3)

# ...and so on for each pair of nested models

# Check the summary for each model to examine the significance of the coefficients
summary(model1)
summary(model2)
summary(model3)
# ...and so on for each model

# Plot residuals for each model to visually inspect the fit
par(mfrow=c(2, 2)) # Adjust the grid size as needed for the number of models
plot(resid(model1) ~ fitted(model1))
plot(resid(model2) ~ fitted(model2))
plot(resid(model3) ~ fitted(model3))
# ...and so on for each model

# Reset the plot settings to default
par(mfrow=c(1, 1))

Anova(model4,type=3)

model=glmer(larvae~poly(carc_wt,degree=2)*tr+male_size+female_size+(1|parent_generation/pair),poisson,data=newdata)

Anova(model,type=3)


modellab=glmer(larvae~carc_wt+(1|parent_generation/pair),poisson,data=newdata[newdata$tr=="lab",])
modelwild=glmer(larvae~carc_wt+(1|parent_generation/pair),poisson,data=newdata[newdata$tr=="wild",])

Anova(modellab)
Anova(modelwild)

modelwild2=glmer(larvae~poly(carc_wt,degree=2)+class+(1|parent_generation),poisson,data=newdata[newdata$tr=="wild",])
modelwild1=glmer(larvae~poly(carc_wt,degree=1)+class+(1|parent_generation),poisson,data=newdata[newdata$tr=="wild",])

Anova(modelwild2,type=3)

summary_data <- newdata %>%
  group_by(class) %>%
  summarise(
    mean_larvae = mean(larvae),
    se_larvae = sd(larvae)/sqrt(n())
  ) 

print(summary_data)


treatment_colors <- c("#479BD5", "#E4191C","black")

broodmass <- ggplot(newdata, aes(x = class, y = larvae,fill=class)) +
  geom_point(aes(color = class), position = position_dodge(width = 0.2), size = 1, alpha = 0.4) + # Jittered data points
  geom_point(data = summary_data, aes(x = class, y = mean_larvae, color = class), shape = 16, size = 4,position = position_dodge(width = 0.5)) +  # Mean points
  geom_errorbar(data = summary_data, aes(y = NULL, ymin = mean_larvae - se_larvae, ymax = mean_larvae + se_larvae, color = class), width = 0.15,position = position_dodge(width = 0.5)) + 
  geom_line(data = summary_data, aes(x = class, y = mean_larvae, group = class, color = class), position = position_dodge(width = 0.5)) + # Connecting lines
  theme_classic() +
  scale_fill_manual(values = treatment_colors) +  # Specifying fill colors
  scale_color_manual(values = treatment_colors) + 
  labs(title = "beetle reproductive success")+
  labs(y = "Brood size") +
  labs(x = "Carcass types")+ 
  scale_y_continuous(limits=c(0,50),breaks=seq(0,50,by=10))+
  theme(
    axis.text = element_text(size = 14),       
    axis.title = element_text(size = 16),
    plot.title=element_text(face = "bold", size = 16),
    legend.position="none")


print(broodmass)



a=emmeans (model,  ~ class , adjust="tukey")
pairs(a)

broodsize <- ggplot(newdata, aes(x=carc_wt, y=larvae)) +
  geom_point(aes(color=tr), size=2, alpha=0.4) +
  geom_smooth(data=subset(newdata, tr=="lab"),
              aes(fill="lab", color="lab", y=larvae, x=carc_wt),
              formula=y ~ poly(x, 2),  # Include polynomial term
              method="glm",
              method.args=list(family="poisson"),
              se=F,
              alpha=0.3,
              size=1) +
  geom_smooth(data=subset(newdata, tr=="wild"),
              aes(fill="wild", color="wild", y=larvae, x=carc_wt),
              formula=y ~ poly(x, 2),  # Include polynomial term
              method="glm",
              method.args=list(family="poisson"),
              se=F,
              alpha=0.3,
              size=1) +
  scale_color_manual(values=c("lab"="#d53e4f", "wild"="#3288bd")) +
  scale_fill_manual(values=c("lab"="#d53e4f", "wild"="#3288bd")) +
  labs(x = "Carcass mass (g)", y = "Brood size") +
  guides(color=guide_legend(override.aes=list(fill=NA)), fill=FALSE) +
  theme_classic() +
  scale_y_continuous(limits=c(0, 50)) +
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=14),
    legend.title=element_blank(),
    legend.key.size=unit(2,"lines"),
    legend.text=element_text(size=12),
    legend.background=element_blank(),
    legend.key=element_rect(fill=NA, colour=NA)
  )
broodsize

model=glmer(tot_mass~poly(carc_wt,degree=2)*tr+male_size+female_size+(1|parent_generation/pair),gaussian,data=newdata)
Anova(model,type=3)

modelwild=glmer(tot_mass~poly(carc_wt,degree=2)+class+(1|parent_generation),gaussian,data=newdata[newdata$tr=="wild",])
Anova(modelwild,type=3)

broodmass <- ggplot(newdata, aes(x=carc_wt, y=tot_mass)) +
  geom_point(aes(color=tr), size=2, alpha=0.4) +
  geom_smooth(data=subset(newdata, tr=="lab"),
              aes(fill="lab", color="lab"),
              method="glm",
              formula=y ~ poly(x, 2),  # Include polynomial term
              method.args=list(family="gaussian"),
              se=F,
              alpha=0.3,
              size=1) +
  geom_smooth(data=subset(newdata, tr=="wild"),
              aes(fill="wild", color="wild"),
              method="glm",
              formula=y ~ poly(x, 2),  # Include polynomial term
              method.args=list(family="gaussian"),
              se=F,
              alpha=0.3,
              size=1) +
  scale_color_manual(values=c("lab"="#d53e4f", "wild"="#3288bd")) +
  scale_fill_manual(values=c("lab"="#d53e4f", "wild"="#3288bd")) +
  labs(x = "Carcass mass (g)", y = "Brood mass (g)") +
  guides(color=guide_legend(override.aes=list(fill=NA)), fill=FALSE) +
  theme_classic() +
  scale_y_continuous(limits=c(0, 10)) +
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=14),
    legend.title=element_blank(),
    legend.key.size=unit(2,"lines"),
    legend.text=element_text(size=12),
    legend.background=element_blank(),
    legend.key=element_rect(fill=NA, colour=NA)
  )
broodmass

ggarrange(success,broodsize,broodmass,yield, labels= c("A", "B", "C","D"),ncol = 2, nrow = 2,widths=c(1,1),common.legend = TRUE,legend="bottom")

ggarrange(broodsize,broodmass, labels= c("A", "B"),ncol = 2, nrow = 1,widths=c(1,1),common.legend = TRUE,legend="bottom")

model1=glmer(yield~poly(carc_wt,degree=1)*tr+larval_density+(1|parent_generation/pair),gaussian,data=newdata)
model2=glmer(yield~poly(carc_wt,degree=2)*tr+larval_density+(1|parent_generation/pair),gaussian,data=newdata)
model3=glmer(yield~poly(carc_wt,degree=3)*tr+larval_density+(1|parent_generation/pair),gaussian,data=newdata)


model1=glmer(yield~poly(larval_density,degree=1)*tr+(1|parent_generation/pair),gaussian,data=newdata[is.na(newdata$yield)==FALSE,])
model2=glmer(yield~poly(larval_density,degree=2)*tr+(1|parent_generation/pair),gaussian,data=newdata[is.na(newdata$yield)==FALSE,])
model3=glmer(yield~poly(larval_density,degree=3)*tr+(1|parent_generation/pair),gaussian,data=newdata[is.na(newdata$yield)==FALSE,])

Anova(model1,type=3)

model1=glmer(yield~poly(carc_wt,degree=1)*tr+larvae+(1|parent_generation/pair),gaussian,data=newdata[is.na(newdata$yield)==FALSE,])
model2=glmer(yield~poly(carc_wt,degree=2)*tr+larvae+(1|parent_generation/pair),gaussian,data=newdata[is.na(newdata$yield)==FALSE,])
model3=glmer(yield~poly(carc_wt,degree=3)*tr+larvae+(1|parent_generation/pair),gaussian,data=newdata[is.na(newdata$yield)==FALSE,])

model1lab=glmer(yield~poly(carc_wt,degree=1)+larvae+(1|parent_generation),gaussian,data=newdatalab[is.na(newdatalab$yield)==FALSE,])
model1wild=glmer(yield~poly(carc_wt,degree=1)+larvae+class+(1|parent_generation),gaussian,data=newdatawild[is.na(newdatawild$yield)==FALSE,])

Anova(model1lab)
Anova(model1wild)

model1wild=glmer(yield~poly(larval_density,degree=2)+class+(1|parent_generation),gaussian,data=newdatawild[is.na(newdatawild$yield)==FALSE,])


yield <- ggplot(newdata, aes(x=carc_wt, y=yield)) +
  geom_point(aes(color=tr), size=2, alpha=0.4) +
  geom_smooth(data=subset(newdata, tr=="lab"),
              aes(fill="lab", color="lab"),
              method="glm",
              formula=y ~ poly(x, 1),  # Include polynomial term
              method.args=list(family="gaussian"),
              se=F,
              alpha=0.3,
              size=1) +
  geom_smooth(data=subset(newdata, tr=="wild"),
              aes(fill="wild", color="wild"),
              method="glm",
              formula=y ~ poly(x, 1),  # Include polynomial term
              method.args=list(family="gaussian"),
              se=F,
              alpha=0.3,
              size=1) +
  scale_color_manual(values=c("lab"="#d53e4f", "wild"="#3288bd")) +
  scale_fill_manual(values=c("lab"="#d53e4f", "wild"="#3288bd")) +
  labs(x = "Carcass mass (g)", y = "Averaged larval mass (g)") +
  guides(color=guide_legend(override.aes=list(fill=NA)), fill=FALSE) +
  theme_classic() +
  scale_y_continuous(limits=c(0, 0.4)) +
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=14),
    legend.title=element_blank(),
    legend.key.size=unit(2,"lines"),
    legend.text=element_text(size=12),
    legend.background=element_blank(),
    legend.key=element_rect(fill=NA, colour=NA)
  )
yield

yield <- ggplot(newdata, aes(x=larval_density, y=yield)) +
  geom_point(aes(color=tr), size=2, alpha=0.4) +
  geom_smooth(data=subset(newdata, tr=="lab"),
              aes(fill="lab", color="lab"),
              method="glm",
              formula=y ~ poly(x, 1),  # Include polynomial term
              method.args=list(family="gaussian"),
              se=F,
              alpha=0.3,
              size=1) +
  geom_smooth(data=subset(newdata, tr=="wild"),
              aes(fill="wild", color="wild"),
              method="glm",
              formula=y ~ poly(x, 1),  # Include polynomial term
              method.args=list(family="gaussian"),
              se=F,
              alpha=0.3,
              size=1) +
  scale_color_manual(values=c("lab"="#d53e4f", "wild"="#3288bd")) +
  scale_fill_manual(values=c("lab"="#d53e4f", "wild"="#3288bd")) +
  labs(x = "Larval density", y = "Averaged larval mass (g)") +
  guides(color=guide_legend(override.aes=list(fill=NA)), fill=FALSE) +
  theme_classic() +
  scale_y_continuous(limits=c(0, 0.4)) +
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=14),
    legend.title=element_blank(),
    legend.key.size=unit(2,"lines"),
    legend.text=element_text(size=12),
    legend.background=element_blank(),
    legend.key=element_rect(fill=NA, colour=NA)
  )
yield

model1=glmer(efficiency~poly(carc_wt,degree=1)*tr+(1|parent_generation/pair),gaussian,data=newdata)
model1=glmer(efficiency~poly(carc_wt,degree=1)+tr+(1|parent_generation/pair),gaussian,data=newdata)

model1=glmer(efficiency_new_method~poly(carc_wt,degree=1)*tr+(1|pair),gaussian,data=newdata[is.na(newdata$efficiency_new_method)==FALSE,])

model2=glmer(efficiency~poly(carc_wt,degree=2)*tr+(1|parent_generation/pair),gaussian,data=newdata)

Anova(model1,type=3)


efficiency <- ggplot(newdata, aes(x=carc_wt, y=efficiency)) +
  geom_point(aes(color=tr), size=2, alpha=0.4) +
  geom_smooth(data=subset(newdata, tr=="lab"),
              aes(fill="lab", color="lab"),
              method="glm",
              formula=y ~ poly(x, 1),  # Include polynomial term
              method.args=list(family="gaussian"),
              se=F,
              alpha=0.3,
              size=1) +
  geom_smooth(data=subset(newdata, tr=="wild"),
              aes(fill="wild", color="wild"),
              method="glm",
              formula=y ~ poly(x, 1),  # Include polynomial term
              method.args=list(family="gaussian"),
              se=F,
              alpha=0.3,
              size=1) +
  scale_color_manual(values=c("lab"="#d53e4f", "wild"="#3288bd")) +
  scale_fill_manual(values=c("lab"="#d53e4f", "wild"="#3288bd")) +
  labs(x = "Larval density", y = "Averaged larval mass (g)") +
  guides(color=guide_legend(override.aes=list(fill=NA)), fill=FALSE) +
  theme_classic() +
  scale_y_continuous(limits=c(0, 0.4)) +
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=14),
    legend.title=element_blank(),
    legend.key.size=unit(2,"lines"),
    legend.text=element_text(size=12),
    legend.background=element_blank(),
    legend.key=element_rect(fill=NA, colour=NA)
  )
efficiency

yield<-ggplot(newdata,aes(x=larval_density,y=yield,color=tr, linetype=tr))+
  geom_point()+
  geom_smooth(method=glm, formula = y~ poly(x,1),se=FALSE,alpha=0.3,size=1)+
  scale_linetype_manual(values=c("solid","solid"))+
  ylab(expression('Averaged larval mass (g)'))+
  xlab(expression('Larval density'))+
  scale_color_manual(values=c("#d53e4f","#3288bd"))+
  scale_fill_manual(values=c("#d53e4f","#3288bd"))+
  theme_classic()+
  scale_y_continuous(limits=c(0,0.4))+
  theme(axis.text=element_text(size=12),axis.title = element_text(size = 14),legend.title=element_blank(),legend.key.size=unit(2,"lines"),legend.text=element_text(size=12),
  )
yield

yield<-ggplot(newdata,aes(x=tot_mass,y=yield,color=tr, linetype=tr))+
  geom_point()+
  geom_smooth(method=glm, formula = y~ poly(x,1),se=FALSE,alpha=0.3,size=1)+
  scale_linetype_manual(values=c("solid","solid"))+
  ylab(expression('Averaged larval mass (g)'))+
  xlab(expression('Larval density'))+
  scale_color_manual(values=c("#d53e4f","#3288bd"))+
  scale_fill_manual(values=c("#d53e4f","#3288bd"))+
  theme_classic()+
  scale_y_continuous(limits=c(0,0.4))+
  theme(axis.text=element_text(size=12),axis.title = element_text(size = 14),legend.title=element_blank(),legend.key.size=unit(2,"lines"),legend.text=element_text(size=12),
  )
yield

efficiency<-ggplot(newdata,aes(x=carc_wt,y=efficiency,color=tr, linetype=tr))+
  geom_point()+
  geom_smooth(method=glm, formula = y~ poly(x,1),se=FALSE,alpha=0.3,size=1)+
  scale_linetype_manual(values=c("solid","solid"))+
  xlab(expression('carc_wt'))+
  ylab(expression('efficiency'))+
  scale_color_manual(values=c("#d53e4f","#3288bd"))+
  scale_fill_manual(values=c("#d53e4f","#3288bd"))+
  theme_classic()+
  scale_y_continuous(limits=c(0,1))+
  theme(axis.text=element_text(size=12),axis.title = element_text(size = 14),legend.title=element_blank(),legend.key.size=unit(2,"lines"),legend.text=element_text(size=12),
  )
efficiency

ratio<-ggplot(newdata,aes(x=carc_wt,y=tot_mass/carc_wt,color=tr, linetype=tr))+
  geom_point(alpha=0.7,size=2)+
  geom_smooth(method=glm, formula = y~ poly(x,1),se=FALSE,alpha=0.3,size=2)+
  scale_linetype_manual(values=c("solid","solid"))+
  xlab(expression('Carcass mass (g)'))+
  ylab(expression('ratio'))+
  scale_color_manual(values=c("#d53e4f","#3288bd"))+
  scale_fill_manual(values=c("#d53e4f","#3288bd"))+
  theme_classic()+
  scale_y_continuous(limits=c(0,0.5))+
  theme(axis.text=element_text(size=12),axis.title = element_text(size = 14),legend.title=element_blank(),legend.key.size=unit(2,"lines"),legend.text=element_text(size=12),
  )
ratio


             