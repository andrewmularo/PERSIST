library(tidyverse)
library(performance)
library(ggfortify)
library(redres)
library(lme4)





setwd("C:/Users/andre/OneDrive/Documents/FitzLab/AdaptOrDie/Persist_MainMS/CTMax")

#Focus on only July samples, as others had too much noise
ctmax <- read.csv("Mularo_2025_CTMax_h2_raw_version1_onlyJuly_Mularo.csv", header = T)



## Calculate number of days in captivity 
ctmax$Date_acquired<- as.Date(ctmax$Date_acquired)
ctmax$Date_measured<- as.Date(ctmax$Date_measured)


ctmax<- ctmax %>% mutate(daysintank = difftime(Date_measured, Date_acquired, units = 'days'))

ctmax$daysintank<- as.numeric(ctmax$daysintank)

ctmax$CTMax<- as.numeric(ctmax$CTMax)
ctmax$Mass<- as.numeric(ctmax$Mass)


## Final Figure for this experiment
ggplot(ctmax, aes(y=CTMax, x=Evo_treatment, fill = Heat_treatment))+geom_boxplot()+theme_classic()+
  xlab("Evolutionary History Treatment")+
  ylab("Critical Thermal Maximum")+
  ggtitle("Critical Thermal Maximum")+
  scale_fill_manual(values = c("darkred", "grey"))


### Scale mass with sex
mass <- as.data.frame(ctmax %>% group_by(Sex) %>% mutate(zMass=scale(Mass)))
ctmax$zMass<- mass$zMass




a<- lmer(CTMax ~ Evo_treatment+Heat_treatment+(1|zMass)+(1|Sex)+(1|Continuous_Trial_Number)+(1|CTMax_filter)+(1|daysintank), data=ctmax)
summary(a)
AIC(a)



plot_redres(a, type = "std_cond")
plot_resqq(a)
plot_ranef(a)

#### Big violation in homogeneity of variance

### Am going to partition into heat and non-heated treatments

ctmax_heat<- ctmax %>% filter(Heat_treatment == "Heated")
ctmax_norm<- ctmax %>% filter(Heat_treatment == "Unheated")


### Compare gene treatment to null for each type
heat<- lmer(CTMax ~ Evo_treatment+(1|zMass)+(1|Sex)+(1|Continuous_Trial_Number)+(1|CTMax_filter)+(1|daysintank), data=ctmax_heat)
summary(heat)
AIC(heat)

heat_null<- lmer(CTMax ~ (1|zMass)+(1|Sex)+(1|Continuous_Trial_Number)+(1|CTMax_filter)+(1|daysintank), data=ctmax_heat)
summary(heat_null)
AIC(heat_null)

anova(heat, heat_null)


norm<- lmer(CTMax ~ Evo_treatment+(1|zMass)+(1|Sex)+(1|Continuous_Trial_Number)+(1|CTMax_filter)+(1|daysintank), data=ctmax_norm)
summary(norm)
AIC(norm)

norm_null<- lmer(CTMax ~ (1|zMass)+(1|Sex)+(1|Continuous_Trial_Number)+(1|CTMax_filter)+(1|daysintank), data=ctmax_norm)
summary(norm_null)
AIC(norm_null)

anova(norm, norm_null)


### Comparing null models from both heated and unheated treatments

anova(heat_null, norm_null)

