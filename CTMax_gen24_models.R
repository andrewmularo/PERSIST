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

#### Model AIC comparisons of Heat, Gene Treatment, Interaction of Heat and Gene Treatment and Null

### All variables
a<- lmer(CTMax_Temp ~ GenTrt+Heat+(GenTrt*Heat)+(1|zMass)+(1|Sex)+(1|Trial)+(1|Cup)+(1|OriginalLin/Family), data=ctmax)
summary(a)
AIC(a)

emmeans(a, pairwise~GenTrt)
emmeans(a, pairwise~Heat)




plot_redres(a, type = "std_cond")
plot_resqq(a)
plot_ranef(a)

## Remove Interactive Effect
b<- lmer(CTMax_Temp ~ GenTrt+Heat+(1|zMass)+(1|Sex)+(1|Trial)+(1|Cup)+(1|OriginalLin/Family), data=ctmax)
summary(b)
AIC(b)



plot_redres(b, type = "std_cond")
plot_resqq(b)
plot_ranef(b)


## Remove Heat
c<- lmer(CTMax_Temp ~ GenTrt+(1|zMass)+(1|Sex)+(1|Trial)+(1|Cup)+(1|OriginalLin/Family), data=ctmax)
summary(c)
AIC(c)


plot_redres(b, type = "std_cond")
plot_resqq(b)
plot_ranef(b)

## Remove Gene Treatment
d<- lmer(CTMax_Temp ~ Heat+(1|zMass)+(1|Sex)+(1|Trial)+(1|Cup)+(1|OriginalLin/Family), data=ctmax)
summary(d)
AIC(d)


## Remove Gene Treatment but leave in interactive effect of heat
e<- lmer(CTMax_Temp ~ Heat+(GenTrt*Heat)+(1|zMass)+(1|Sex)+(1|Trial)+(1|Cup)+(1|OriginalLin/Family), data=ctmax)
summary(e)
AIC(e)

## Remove Heat but leave in interactive effect of heat
f<- lmer(CTMax_Temp ~ GenTrt+(GenTrt*Heat)+(1|zMass)+(1|Sex)+(1|Trial)+(1|Cup)+(1|OriginalLin/Family), data=ctmax)
summary(f)
AIC(f)


## Null Model - no variables of interest
null <- lmer(CTMax_Temp ~ (1|zMass)+(1|Sex)+(1|Trial)+(1|Cup)+(1|OriginalLin/Family), data=ctmax)
summary(null)
AIC(null)

#### Lowest AIC is null model

### Let's see if original lineage has an significant effect on model structure

ggplot(ctmax, aes(x=OriginalLin, y= CTMax_Temp, fill=Heat))+
  geom_boxplot()+
  scale_fill_manual(values = c( "grey","darkred"))


lineage <- lmer(CTMax_Temp ~ OriginalLin+(1|zMass)+(1|Sex)+(1|Trial)+(1|Cup), data=ctmax)
summary(lineage)
AIC(lineage)

lineageheat <- lmer(CTMax_Temp ~ OriginalLin+Heat+(1|zMass)+(1|Sex)+(1|Trial)+(1|Cup), data=ctmax)
summary(lineageheat)
AIC(lineageheat)

lineagegene<- lmer(CTMax_Temp ~ OriginalLin+(Heat*GenTrt)+(1|zMass)+(1|Sex)+(1|Trial)+(1|Cup), data=ctmax)
summary(lineagegene)
AIC(lineagegene)

null<- lmer(CTMax_Temp ~ (1|zMass)+(1|Sex)+(1|Trial)+(1|Cup), data=ctmax)
summary(null)
AIC(null)

#### It doesn't, lowest AIC model is the null model





##### Final Fig for publication

### bring in gen 1 and extract relevant variables
ctmax1<-read.csv("CTMax_G1_Final_AJM.csv")
ctmax1<- ctmax1 %>% mutate(generation = "1")
ctmax1<- ctmax1 %>% mutate(Heat = "NA")
ctmax1<- ctmax1 %>% select("GenTrt", "CTMax_Temp", "generation", "Heat")


### Add new column for gen 24
ctmax<- ctmax %>% mutate(generation = "24")
ctmax24<- ctmax %>% select("GenTrt", "CTMax_Temp", "generation", "Heat")

## Combine both
ctmax_all<- rbind(ctmax1, ctmax24)

### Plot both
ggplot(ctmax_all, aes(x=GenTrt, y=CTMax_Temp,  fill = Heat))+
  geom_boxplot()+
  ggtitle("CTMax across Generations")+
  ylab("CTMax Temperature (C)")+
  xlab("Evolutionary History Treatment")+
  theme_classic()+
  scale_fill_manual(values=c("black", "#607B8B","#8B1A1A"))






