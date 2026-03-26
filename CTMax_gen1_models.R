library(tidyverse)
library(performance)
library(ggfortify)
library(redres)
library(lme4)
library(glmmTMB)
library(emmeans)


setwd("C:/Users/andre/OneDrive/Documents/FitzLab/AdaptOrDie/Persist_MainMS/CTMax")
ctmax <- read.csv("CTMax_G1_Final_AJM.csv", header = T)
#prep factors for formatting into analysis
ctmax$Sex <- factor(ctmax$Sex, levels = c("J","F","M"))
ctmax$Date <- as.Date(ctmax$Date, format = "%m/%d/%y")
ctmax$Family <- factor(ctmax$Family)
ctmax$GenLin <- factor(ctmax$GenLin)
ctmax$OriginalLin <- factor(ctmax$OriginalLin)
ctmax$GenTrt <- factor(ctmax$GenTrt)
ctmax$Trial<- factor(ctmax$Trial)
ctmax$Cup <- factor(ctmax$Cup)
ctmax$YSI <- factor(ctmax$YSI)


### Create zlength variable for mass to standardize the mass of each sex and avoid interactive effects of those variables
mass <- as.data.frame(ctmax %>% group_by(Sex) %>% mutate(zMass=scale(Mass)))
ctmax$zMass<- mass$zMass

ggplot(ctmax, aes(x=zMass, color=Sex))+
  geom_boxplot()
## Females are larger, males middle, juveniles smallest - makes sense


### Explore potential random effects
#### Correlation between length and mass
### Correlated, but will use mass in models to compare to gen24 
cor(ctmax[,c(13:14)], use="complete.obs")
ggplot(ctmax, aes(x=Length, y=Mass, color=Sex))+
  geom_point()


## Correlation with mass and ctmax
## Sorted by sex
ggplot(ctmax, aes(x=zMass,y=CTMax_Temp, color=Sex))+
  geom_point()+
  geom_smooth(method = "lm", formula=y~x)
#general
ggplot(ctmax, aes(x=zMass,y=CTMax_Temp))+
  geom_point()+
  geom_smooth(method = "lm", formula=y~x)+
  ggtitle("Mass v CTMax")

### It seems like there is an interaction with z transformed mass and ctmax, 



##Family (Individual sibling groups)
ggplot(ctmax, aes(x=Family, y=CTMax_Temp))+
  geom_boxplot()+
  ggtitle("Family Group")


##Genetic Lineage (Black Tank)
ggplot(ctmax, aes(x=GenLin, y=CTMax_Temp))+
  geom_boxplot()

##Original Lineage (Original Black tank)
ggplot(ctmax, aes(x=OriginalLin, y=CTMax_Temp))+
  geom_boxplot()


##Genetic Treatment
ggplot(ctmax, aes(x=GenTrt, y=CTMax_Temp))+
  geom_boxplot()

# Trial
ggplot(ctmax, aes(x=Trial, y=CTMax_Temp))+
  geom_boxplot()+
  ggtitle("Trial")
## This effect may be due to the switching of the YSI

##Cup within cooler
ggplot(ctmax, aes(x=Cup, y=CTMax_Temp))+
  geom_boxplot()+
  ggtitle("Cup")


##YSI
ggplot(ctmax, aes(x=YSI, y=CTMax_Temp))+
  geom_boxplot()+
  ggtitle("YSI")
#definitely looks to be an effect of YSI.

##Sex: Effects probably due to size and generation effects.
ggplot(ctmax, aes(x=Sex, y=CTMax_Temp))+
  geom_boxplot()+
  ggtitle("Sex")


#Distribution of our dependent variable, CTmax: looks slightly left-skewed but mostly normal. 
ggplot(ctmax, aes(x=CTMax_Temp))+
  geom_bar()


##Normality tests
ggdensity(ctmax$CTMax_Temp, 
          main = "Density plot of CTmax",
          xlab = "CTmax")
ggqqplot(ctmax$CTMax_Temp)
shapiro.test(ctmax$CTMax_Temp)
#qqplot confirms left skew, shapiro test indicates non-normality. Will need to check model residuals.



## Final Figure for Gen 1 

ggplot(ctmax, aes(y=CTMax_Temp, x=GenTrt, fill = GenTrt))+geom_boxplot()+theme_classic()+
  xlab("Evolutionary History Treatment")+
  ylab("Critical Thermal Maximum")+
  ggtitle("Generation 1")+
  scale_fill_manual(values = c("black","black","black"))





########## General Linear Models

# Desired full model
m <- lmer(CTMax_Temp ~ GenTrt+(1|YSI)+(1|zMass)+(1|Sex)+(1|Trial)+(1|Cup)+(1|OriginalLin/Family), data=ctmax)
summary(m)

plot_redres(m, type = "std_cond")
plot_resqq(m)
plot_ranef(m)
shapiro.test(resid(m))
hist(resid(m))


null <- lmer(CTMax_Temp ~ (1|YSI)+(1|zMass)+(1|Sex)+(1|Trial)+(1|Cup)+(1|OriginalLin/Family), data=ctmax)
summary(null)

plot_redres(null, type = "std_cond")
plot_resqq(null)
plot_ranef(null)
shapiro.test(resid(null))
hist(resid(null))

#### Visual inspection indicate some slight deviations from normality/homogeneity of
#### variance, but nothing too extreme

#### shapiro test indicates the residuals are not normally distributed for each model


anova(m,null)


### Generalized linear mixed models 

m <- glmmTMB(CTMax_Temp ~ GenTrt+(1|YSI)+(1|zMass)+(1|Sex)+(1|Trial)+(1|Cup)+(1|OriginalLin/Family), data=ctmax)
summary(m)

emmeans(m, pairwise~GenTrt)


