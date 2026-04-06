library(tidyverse)
library(ggplot2)
library(gganimate)
library(gapminder)
library(ggpubr)
library(redres)
library(lme4)
library(Rmisc)
library(glmmTMB)
library(emmeans)



setwd("C://Users/andre/OneDrive/Documents/FitzLab/AdaptOrDie/Persist_MainMS/DemographicData")


data<- read.csv("MegaMesocosmCensus_Population.csv")


data<- data %>% mutate(total = Males+Females+Juv+Babies)
data<- data %>% mutate(total_adults = Males+Females)

data_heat<- data %>% filter(TrtID == "H")
data_norm<- data %>% filter(TrtID == "C")



ggplot(data_heat, aes(x=Census, y=total_adults, group = GenTrt, color = GenTrt, fill = GenTrt))+
  geom_smooth()+geom_point()+
  theme_classic()+scale_fill_manual(values = c("#0070C0", "#00B050","#C04F15" ))+
  scale_color_manual(values = c("#0070C0", "#00B050","#C04F15" ))+
  theme(axis.text.x = element_blank())+
  ylab("Number of Adults")+
  xlab("Time")+
  ggtitle("Heated Tanks")+
  ylim(0,21)

ggplot(data_norm, aes(x=Census, y=total_adults, group = GenTrt, color = GenTrt, fill = GenTrt))+
  geom_smooth()+geom_point()+
  theme_classic()+scale_fill_manual(values = c("#0070C0", "#00B050","#C04F15" ))+
  scale_color_manual(values = c("#0070C0", "#00B050","#C04F15" ))+
  theme(axis.text.x = element_blank())+
  ylab("Number of Adults")+
  xlab("Time")+
  ggtitle("Unheated Tanks")+
  ylim(0,21)

data_heat_summary<- summarySE(data_heat, measurevar = "total_adults", groupvars = c("GenTrt", "Census"))
data_norm_summary<- summarySE(data_norm, measurevar = "total_adults", groupvars = c("GenTrt", "Census"))


a<- ggplot(data_heat_summary, aes(x=Census, y=total_adults, group = GenTrt, color = GenTrt, fill = GenTrt))+
  geom_errorbar(aes(ymin=total_adults-se, ymax=total_adults+se))+
  geom_line()+
  geom_point(size=5)+
  theme_classic()+scale_fill_manual("Evolutionary History Treatment", values = c("#0070C0", "#00B050","#C04F15" ))+
  scale_color_manual("Evolutionary History Treatment", values = c("#0070C0", "#00B050","#C04F15" ))+
  theme(axis.text.x = element_blank())+
  ylab("Number of Adults")+
  xlab("Time")+
  ggtitle("Heated Tanks")+
  ylim(0,25)

b<- ggplot(data_norm_summary, aes(x=Census, y=total_adults, group = GenTrt, color = GenTrt, fill = GenTrt))+
  geom_errorbar(aes(ymin=total_adults-se, ymax=total_adults+se))+
  geom_line()+
  geom_point(size=5)+
  theme_classic()+scale_fill_manual("Evolutionary History Treatment", values = c("#0070C0", "#00B050","#C04F15" ))+
  scale_color_manual("Evolutionary History Treatment", values = c("#0070C0", "#00B050","#C04F15" ))+
  theme(axis.text.x = element_blank())+
  ylab("Number of Adults")+
  xlab("Time")+
  ggtitle("Unheated Tanks")+
  ylim(0,25)

  
ggarrange(b,a)


#### Model


### Went with glmmTMB - model more robust to non-normality
model<- glmmTMB(total_adults~GenTrt*Census+TrtID+(1|GenLin), data)

summary(model)


hist(resid(model))

emmeans(model, pairwise~GenTrt)

emmeans(model, pairwise~TrtID)



#### Differences between gene treatments for each heat treatment

heatmodel<- glmmTMB(total_adults~GenTrt*Census+(1|GenLin), data_heat)
summary(heatmodel)
hist(resid(heatmodel))
emmeans(heatmodel, pairwise~GenTrt)



normmodel<- glmmTMB(total_adults~GenTrt*Census+(1|GenLin), data_norm)
summary(normmodel)
hist(resid(normmodel))
emmeans(normmodel, pairwise~GenTrt)

#### Really breaks the rules for normality of residuals




#### Man Whitney U test - end census pop size

data_21<- data %>% filter(Census == 21) 


ggplot(data_21, aes(x=GenTrt, y=total, fill = TrtID))+
  geom_boxplot()+
  scale_fill_manual("Heat Treatment", values = c("Black", "Red"))+
  theme_classic()+ylab("Number of Individuals")+
  xlab("Evolutionary History Treatment")
  

data_21_c<- data_21 %>% filter(GenTrt == "Control")
data_21_i<- data_21 %>% filter(GenTrt == "Inbred")
data_21_gf<- data_21 %>% filter(GenTrt == "GeneFlow")

wilcox.test(total~TrtID, data_21_c)
wilcox.test(total~TrtID, data_21_i)
wilcox.test(total~TrtID, data_21_gf)



data_21_noextinct<- data_21 %>% filter(total > 0)

ggplot(data_21_noextinct, aes(x=GenTrt, y=total, fill = TrtID))+
  geom_boxplot()+
  scale_fill_manual(values = c("Black", "Red"))+
  theme_classic()+ylab("Number of Individuals")+
  xlab("Evolutionary History Treatment")

data_21_noextinct_c<- data_21_noextinct %>% filter(GenTrt == "Control")
data_21_noextinct_i<- data_21_noextinct %>% filter(GenTrt == "Inbred")
data_21_noextinct_gf<- data_21_noextinct %>% filter(GenTrt == "GeneFlow")

wilcox.test(total~TrtID, data_21_noextinct_c)
wilcox.test(total~TrtID, data_21_noextinct_i)
wilcox.test(total~TrtID, data_21_noextinct_gf)





### All individuals
data_heat_summary_all<- summarySE(data_heat, measurevar = "total", groupvars = c("GenTrt", "Census"))
data_norm_summary_all<- summarySE(data_norm, measurevar = "total", groupvars = c("GenTrt", "Census"))

c<- ggplot(data_heat_summary_all, aes(x=Census, y=total, group = GenTrt, color = GenTrt, fill = GenTrt))+
  geom_errorbar(aes(ymin=total-se, ymax=total+se))+
  geom_line()+
  geom_point(size=5)+
  theme_classic()+scale_fill_manual(values = c("#0070C0", "#00B050","#C04F15" ))+
  scale_color_manual(values = c("#0070C0", "#00B050","#C04F15" ))+
  theme(axis.text.x = element_blank())+
  ylab("Number of Individuals")+
  xlab("Time")+
  ggtitle("Heated Tanks")+
  ylim(0,55)

d<- ggplot(data_norm_summary_all, aes(x=Census, y=total, group = GenTrt, color = GenTrt, fill = GenTrt))+
  geom_errorbar(aes(ymin=total-se, ymax=total+se))+
  geom_line()+
  geom_point(size=5)+
  theme_classic()+scale_fill_manual(values = c("#0070C0", "#00B050","#C04F15" ))+
  scale_color_manual(values = c("#0070C0", "#00B050","#C04F15" ))+
  theme(axis.text.x = element_blank())+
  ylab("Number of Individuals")+
  xlab("Time")+
  ggtitle("Unheated Tanks")+
  ylim(0,55)



ggarrange(d,c)


### Color by lineage

data_heat_gf<- data_heat %>% filter(GenTrt == "GeneFlow")
data_heat_in<- data_heat %>% filter(GenTrt == "Inbred")
data_heat_out<- data_heat %>% filter(GenTrt == "Control")



a<- ggplot(data_heat_gf, aes(x=Census, y=total_adults, group = TankID, color = GenLin))+
  geom_line()+geom_point()+
  theme_classic()+
  theme(axis.text.x = element_blank())+
  ylab("Number of Adults")+
  xlab("Time")+
  ggtitle("Gene Flow")+
  ylim(0,21)

b<- ggplot(data_heat_in, aes(x=Census, y=total_adults, group = TankID, color = GenLin))+
  geom_line()+geom_point()+
  theme_classic()+
  theme(axis.text.x = element_blank())+
  ylab("Number of Adults")+
  xlab("Time")+
  ggtitle("Inbred")+
  ylim(0,21)


c<- ggplot(data_heat_out, aes(x=Census, y=total_adults, group = TankID, color = GenLin))+
  geom_line()+geom_point()+
  theme_classic()+
  theme(axis.text.x = element_blank())+
  ylab("Number of Adults")+
  xlab("Time")+
  ggtitle("Outbred")+
  ylim(0,21)



ggarrange(a,b,c)



data_norm_gf<- data_norm %>% filter(GenTrt == "GeneFlow")
data_norm_in<- data_norm %>% filter(GenTrt == "Inbred")
data_norm_out<- data_norm %>% filter(GenTrt == "Control")



d<- ggplot(data_norm_gf, aes(x=Census, y=total_adults, group = TankID, color = GenLin))+
  geom_line()+geom_point()+
  theme_classic()+
  theme(axis.text.x = element_blank())+
  ylab("Number of Adults")+
  xlab("Time")+
  ggtitle("Gene Flow")+
  ylim(0,21)

e<- ggplot(data_norm_in, aes(x=Census, y=total_adults, group = TankID, color = GenLin))+
  geom_line()+geom_point()+
  theme_classic()+
  theme(axis.text.x = element_blank())+
  ylab("Number of Adults")+
  xlab("Time")+
  ggtitle("Inbred")+
  ylim(0,21)


f<- ggplot(data_norm_out, aes(x=Census, y=total_adults, group = TankID, color = GenLin))+
  geom_line()+geom_point()+
  theme_classic()+
  theme(axis.text.x = element_blank())+
  ylab("Number of Adults")+
  xlab("Time")+
  ggtitle("Outbred")+
  ylim(0,21)



ggarrange(d,e,f)





