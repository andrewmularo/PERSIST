library(tidyverse)
library(performance)
library(ggfortify)
library(redres)
library(lme4)


setwd("C:/Users/andre/OneDrive/Documents/FitzLab/AdaptOrDie/Persist_MainMS/CTMax")
ctmax <- read.csv("CTMax_G24_Final_Mularo.csv", header = T)



### Main Figure for generation 24
ggplot(ctmax, aes(y=CTMax_Temp, x=GenTrt, fill = Heat))+geom_boxplot()+theme_classic()+
  xlab("Evolutionary History Treatment")+
  ylab("Critical Thermal Maximum")+
  ggtitle("Generation 24")+
  scale_fill_manual(values = c("darkred", "grey"))


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
ctmax$Mass<- as.numeric(ctmax$Mass)

### Create zlength variable for mass to standardize the mass of each sex and avoid interactive effects of those variables
mass <- as.data.frame(ctmax %>% group_by(Sex) %>% mutate(zMass=scale(Mass)))
ctmax$zMass<- mass$zMass

### Remove outlier - it is unlikely that this is a true heat stress response

ctmax$CTMax_Temp[ctmax$CTMax_Temp == 33.4]<- NA

mass <- as.data.frame(ctmax %>% group_by(Sex) %>% mutate(zMass=scale(Mass)))
ctmax$zMass<- mass$zMass

ggplot(ctmax, aes(x=zMass, color=Sex))+
  geom_boxplot()
#### This is super weird - juveniles are huge 
#### Scale was really difficult to use at the field station - lots of problems

### Will conduct z score distribution analysis of mass based on sex to try to correct for this
ctmax_weights<- ctmax %>% select(c("Sex", "Mass"))
ctmax_weights_m<- ctmax_weights %>% filter(Sex == "M")
ctmax_weights_f<- ctmax_weights %>% filter(Sex == "F")
ctmax_weights_j<- ctmax_weights %>% filter(Sex == "J")


#Z score for males
mass_m<- as.numeric(na.omit(ctmax_weights_m$Mass))
mean_mass_m<- mean(mass_m)
sd_mass_m<- sd(mass_m)
zscore_mass_m<- (mass_m-mean_mass_m)/sd_mass_m
Massm_outlierID<- cbind(mass_m,zscore_mass_m)


# Z score for females
mass_f<- as.numeric(na.omit(ctmax_weights_f$Mass))
mean_mass_f<- mean(mass_f)
sd_mass_f<- sd(mass_f)
zscore_mass_f<- (mass_f-mean_mass_f)/sd_mass_f
Massf_outlierID<- cbind(mass_f,zscore_mass_f)

#Z score for juveniles
mass_j<- as.numeric(na.omit(ctmax_weights_j$Mass))
mean_mass_j<- mean(mass_j)
sd_mass_j<- sd(mass_j)
zscore_mass_j<- (mass_j-mean_mass_j)/sd_mass_j
Massj_outlierID<- cbind(mass_j,zscore_mass_j)

### Outlier mass values include:
# Males = 0.366, 
### Females had no outliers
## Juveniles = 1.1590, 1.0830, 1.0230, 0.8800

ctmax[4, "Mass"]<- NA
ctmax[37, "Mass"]<- NA
ctmax[25, "Mass"]<- NA
ctmax[22, "Mass"]<- NA
ctmax[135, "Mass"]<- NA

ctmax[4, "zMass"]<- NA
ctmax[37, "zMass"]<- NA
ctmax[25, "zMass"]<- NA
ctmax[22, "zMass"]<- NA
ctmax[135, "zMass"]<- NA

ggplot(ctmax, aes(x=zMass, color=Sex))+
  geom_boxplot()


### CTMax vs z transformed weight
ggplot(ctmax, aes(x=zMass,y=CTMax_Temp))+
  geom_point()+
  geom_smooth(method = "lm", formula=y~x)+
  ggtitle("Mass v CTMax")



##Family (Individual sibling groups)
ggplot(ctmax, aes(x=Family, y=CTMax_Temp))+
  geom_boxplot()+
  ggtitle("Family Group")


# Trial
ggplot(ctmax, aes(x=Trial, y=CTMax_Temp))+
  geom_boxplot()+
  ggtitle("Trial")

##Cup within cooler
ggplot(ctmax, aes(x=Cup, y=CTMax_Temp))+
  geom_boxplot()+
  ggtitle("Cup")

### No YSI Variability for this generation

##Sex: Effects probably due to size and generation effects.
ggplot(ctmax, aes(x=Sex, y=CTMax_Temp))+
  geom_boxplot()+
  ggtitle("Sex")



ggplot(ctmax, aes(x=OriginalLin, y=CTMax_Temp))+
  geom_boxplot()+
  ggtitle("Sex")


##### Generalized linear mixed effects model

a<- glmmTMB(CTMax_Temp ~ GenTrt+Heat+(GenTrt*Heat)+(1|zMass)+(1|Sex)+(1|Trial)+(1|Cup)+(1|OriginalLin/Family), data=ctmax)
        

summary(a)


hist(resid(a))

emmeans(a, pairwise~GenTrt)

emmeans(a, pairwise~Heat)


