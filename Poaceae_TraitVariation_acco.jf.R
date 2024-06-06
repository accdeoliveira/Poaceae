
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggpubr)


#getting all data ####

alldatabases <- read.table("../TRY_DATA/alldata_final.csv", header=T, sep = ";")
head(alldatabases)
tail(alldatabases)
view(alldatabases)
colnames(alldatabases)


allobservations <- alldatabases %>%
  filter(plant.duration == "annual" | plant.duration=="perennial") %>%
  group_by(Database, plant.duration) %>%
  dplyr::summarise(across(Germination_p:Root.tissue.density_g.cm3, ~ sum(!is.na(.x)), .names = "n_{.col}"))

allobservations
View(allobservations)

#summarizing with averages among all databases
allsumm <- alldatabases %>%
  select(-c(Database, species)) %>%
  group_by(species_tpl, plant.duration) %>%
  summarise_each(funs(mean(., na.rm = TRUE))) #dealing with NAs

View(allsumm)

write.csv(allsumm, file="allsumm_updt4.csv")


#REMOVING THE CROPS ####

allsummupd <- read.table("../allsumm_updt4.csv", header=T, sep = ";")
head(allsummupd)
tail(allsummupd)


allsummupd3 <- allsummupd %>%
  filter(plant.duration=="annual" | plant.duration=="perennial") %>%
  filter(species_tpl != "Zea mays",
         species_tpl != "Avena sativa",
         species_tpl !=  "Hordeum vulgare",
         species_tpl !=   "Oryza sativa",
         species_tpl !=   "Secale cereale",
         species_tpl !=   "Sorghum bicolor",
         species_tpl !=   "Triticum dicoccoides",
         species_tpl !=   "Triticum durum",
         species_tpl !=   "Zea mexicana")                  #removendo crops

head(allsummupd3)
tail(allsummupd3)

#number of observations of each trait within plant duration
colnames(allsummupd3)
observations <- allsummupd3 %>%
  group_by(plant.duration) %>%
  dplyr::summarise(across(Germination_p:Root.tissue.density_g.cm3, ~ sum(!is.na(.x)), .names = "n_{.col}"))

observations
View(observations)

allsummupd3%>%
  group_by(plant.duration)%>%
  count() #499 species - 8 crops removed 

colnames(allsummupd3)

# after remov crops: comparison by trait####
library(ggpubr)
#root traits

#SRL
srl <- glm(SRL_m.g~plant.duration, data = allsummupd3, family = Gamma) 

par(mfrow=c(2,2))
plot(srl)

library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = srl, n = 250)
plot(simulationOutput)

summary(srl)
summary.glm(srl)


library(emmeans)
compar<-emmeans(srl, list(pairwise ~ plant.duration), adjust = "t")
compar
contrast(compar, method = "pairwise") #p-value: 0.0477

contrast(regrid(compar, method = "pairwise")) #p-value: 0.0477

summary.glm

ra <- ggdensity(allsummupd3, x = "SRL_m.g",
                add = "mean", rug = TRUE,
                color = "plant.duration", fill = "plant.duration", xlab = "SRL (m/g)", palette = c("#999999", "#E69F00"))+ 
  annotate(geom="text", x=450, y=0.005, label="p-value = 0.048",
          color="black")+
  ylim(0,0.0065)

ra 

raop2 <- ggdensity(allsummupd3, x = "SRL_m.g",
                add = "mean", rug = TRUE,
                color = "plant.duration", fill = "plant.duration", xlab = "SRL (m/g)", palette = c("#E69F00","#009E73"))+ 
  scale_y_continuous(limits = c(0, 0.006), breaks = c(0, 0.003, 0.006)) +
  annotate(geom="text", x=450, y=0.005, label="p-value = 0.048",
           color="black")

raop2 #liked more

raop4 <- ggdensity(allsummupd3, x = "SRL_m.g",
                   add = "mean", rug = TRUE,
                   color = "plant.duration", fill = "plant.duration", xlab = "SRL (m/g)", palette = c("#0072B2", "#E69F00"))+ 
  annotate(geom="text", x=450, y=0.005, label="p-value = 0.048",
           color="black")+
  ylim(0,0.0065)

raop4 

#Root nitrogen
rootn <- lm(sqrt(RootN.mass_mg.g)~plant.duration, data = allsummupd3) #p-value: 0.9699


par(mfrow=c(2,2))
plot(rootn)
summary(rootn)

library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = rootn, n = 250)
plot(simulationOutput)

rb <- ggdensity(allsummupd3, x = "RootN.mass_mg.g",
                add = "mean", rug = TRUE,
                color = "plant.duration",  fill = "plant.duration",
                xlab = "Root N (mg/g)", palette = c("#E69F00","#009E73")) + 
  scale_y_continuous(limits = c(0, 0.11), breaks = c(0, 0.05, 0.1)) +
  annotate(geom="text", x=18, y=0.098, label="p-value = 0.969",
           color="black")+
  theme(axis.title.y = element_blank())
rb

#Root diameter
rootdiam <- glm(Root.diameter_mm~plant.duration, data = allsummupd3, family = Gamma) #p-value: 0.9768

par(mfrow=c(2,2))
plot(rootdiam)
summary(rootdiam)

library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = rootdiam, n = 250)
plot(simulationOutput)

library(emmeans)
compar<-emmeans(rootdiam, list(pairwise ~ plant.duration), adjust = "tukey")
compar
contrast(compar, method = "pairwise") #p-value: 0.9768

rc <- ggdensity(allsummupd3, x = "Root.diameter_mm",
                add = "mean", rug = TRUE,
                color = "plant.duration",  fill = "plant.duration", 
                xlab = "Root diameter (mm)", palette = c("#E69F00","#009E73")) + 
  scale_y_continuous(limits = c(0, 4.5), breaks = c(0, 2, 4)) +
  annotate(geom="text", x=0.5, y=4, label="p-value = 0.977",
           color="black")+
  theme(axis.title.y = element_blank()) 

rc

#RMF
rmf <- lm(RMF_g.g~plant.duration, data = allsummupd3) #p-value: 0.2032
par(mfrow=c(2,2))
plot(rmf)
summary(rmf)
summary.lm(rmf)

library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = rmf, n = 250)
plot(simulationOutput)

rd <- ggdensity(allsummupd3, x = "RMF_g.g",
                add = "mean", rug = TRUE,
                color = "plant.duration",  fill = "plant.duration",
                xlab = "RMF (g/g)", palette = c("#E69F00","#009E73")) + 
  scale_y_continuous(limits = c(0, 3.5), breaks = c(0, 1.5, 3)) +
  annotate(geom="text", x=0.6, y=3.5, label="p-value = 0.203",
           color="black")+
  theme(axis.title.y = element_blank()) 

rd

#Root depth
rootdepth <- lm(Root.depth_cm~plant.duration, data = allsummupd3) #p-value: 0.04672
par(mfrow=c(2,2))
plot(rootdepth)
summary(rootdepth)

library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = rootdepth, n = 250)
plot(simulationOutput)

re <- ggdensity(allsummupd3, x = "Root.depth_cm",
                add = "mean", rug = TRUE,
                color = "plant.duration",  fill = "plant.duration", 
                xlab = "Root depth (cm)",palette = c("#E69F00","#009E73")) + 
  scale_y_continuous(limits = c(0, 0.0082), breaks = c(0, 0.004, 0.008)) +
  annotate(geom="text", x=250, y=0.0075, label="p-value = 0.047",
           color="black")+
  theme(axis.title.y = element_blank())

re

#Root tissue density
roottisden <- lm(Root.tissue.density_g.cm3~plant.duration, data = allsummupd3) #p-value: 0.02677
par(mfrow=c(2,2))
plot(roottisden)
summary(roottisden)

library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = roottisden, n = 250)
plot(simulationOutput)

rf <- ggdensity(allsummupd3, x = "Root.tissue.density_g.cm3",
                add = "mean", rug = TRUE,
                color = "plant.duration",  fill = "plant.duration",
                xlab = "RTD (g/cm)", palette = c("#E69F00","#009E73")) + 
  scale_y_continuous(limits = c(0, 2.2), breaks = c(0, 1.1, 2.2)) +
  annotate(geom="text", x=0.55, y=2.2, label="p-value = 0.268",
           color="black") 

rf

#ggarrange(ra, rb, rc, rd, re, rf, ncol = 3, nrow = 2,  common.legend = T, legend="bottom", labels="AUTO")

#seed traits
#germination
germout <- allsummupd3[-c(215,297,358,190),]

#germ <- glm(log(Germination_p)~plant.duration, data = germout, family = Gamma) #p-value: 
germ <- lm((Germination_p)^2~plant.duration, data = germout) #p-value: 


par(mfrow=c(2,2))
plot(germ)

library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = germ, n = 250)
plot(simulationOutput) #melhor modelo, MAS normalidade não foi atendida

summary(germ)

library(emmeans)
compar<-emmeans(germ, list(pairwise ~ plant.duration), adjust = "tukey")
compar
contrast(regrid(compar, method = "pairwise"))


sa <- ggdensity(germout, x = "Germination_p",
                add = "mean", rug = TRUE,
                color = "plant.duration",  fill = "plant.duration",
                xlab = "Germination (%)", palette = c("#E69F00","#009E73")) + 
  scale_y_continuous(limits = c(0, 0.07), breaks = c(0, 0.035, 0.07)) +
  annotate(geom="text", x=60, y=0.06, label="p-value = 0.128",
           color="black")+
  theme(axis.title.y = element_blank())
sa

#n obs germination after removing outliers
germobs <- germout %>%
  group_by(plant.duration) %>%
  summarise(across(Germination_p:Root.tissue.density_g.cm3, ~ sum(!is.na(.x)), .names = "n_{.col}"))
germobs
#View(germobs)

#seed mass
seedm <- lm(sqrt(Seed.mass_mg)~plant.duration, data = allsummupd3) #p-value: 0.8316
par(mfrow=c(2,2))
plot(seedm)
summary(seedm)

library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = seedm, n = 250)
plot(simulationOutput)


sb <- ggdensity(allsummupd3, x = "Seed.mass_mg",
                add = "mean", rug = TRUE,
                color = "plant.duration",  fill = "plant.duration",
                xlab = "Seed mass (mg)", palette = c("#E69F00","#009E73")) + 
  scale_y_continuous(limits = c(0, 0.28), breaks = c(0, 0.1, 0.2)) +
  annotate(geom="text", x=15, y=0.25, label="p-value = 0.832",
           color="black")+
  theme(axis.title.y = element_blank()) 
sb

#ggarrange(sa, sb, ncol = 2, nrow = 1,  common.legend = T, legend="bottom", labels="AUTO")


#shoot traits
#height
height <- lm(log(Height_cm)~plant.duration, data = allsummupd3) #p-value: 5.753e-06

par(mfrow=c(2,2))
plot(height)
summary(height)

simulationOutput <- simulateResiduals(fittedModel = height, n = 250)
plot(simulationOutput)

sha <- ggdensity(allsummupd3, x = "Height_cm",
                 add = "mean", rug = TRUE,
                 color = "plant.duration",  fill = "plant.duration",
                 xlab = "Height (cm)", palette = c("#E69F00","#009E73")) +
  scale_y_continuous(limits = c(0, 0.024), breaks = c(0, 0.01, 0.02)) +
  theme(plot.margin = margin(t = 30, r = 10, b = 10, l = 10)) +
  annotate("text", x=150, y=0.020, label= "p-value < 0.001") 
sha


#LDMC
ldmc <- lm(LDMC_mg.g~plant.duration, data = allsummupd3) #p-value: 0.0003235
par(mfrow=c(2,2))
plot(ldmc)
summary(ldmc)

simulationOutput <- simulateResiduals(fittedModel = ldmc, n = 250)
plot(simulationOutput)

shb <- ggdensity(allsummupd3, x = "LDMC_mg.g",
                 add = "mean", rug = TRUE,
                 color = "plant.duration",  fill = "plant.duration",
                 xlab = "LDMC (mg/g)", palette = c("#E69F00","#009E73")) +
  scale_y_continuous(limits = c(0, 0.0072), breaks = c(0, 0.003, 0.006)) +
  annotate("text", x=400, y=0.007, label= "p-value < 0.001") +
  theme(plot.margin = margin(t = 30, r = 10, b = 10, l = 10), axis.title.y = element_blank()) 

shb

#Leaf nitrogen
leafn <- lm(LeafN.mass_mg.g~plant.duration, data = allsummupd3) #p-value: 0.001437
par(mfrow=c(2,2))
plot(leafn)
summary(leafn)

simulationOutput <- simulateResiduals(fittedModel = leafn, n = 250)
plot(simulationOutput)

shc <- ggdensity(allsummupd3, x = "LeafN.mass_mg.g",
                 add = "mean", rug = TRUE,
                 color = "plant.duration",  fill = "plant.duration",
                 xlab = "Leaf N (mg/g)", palette = c("#E69F00","#009E73")) +
  scale_y_continuous(limits = c(0, 0.062), breaks = c(0, 0.03, 0.06)) +
  theme(plot.margin = margin(t = 30, r = 10, b = 10, l = 10)) +
  annotate("text", x=28, y=0.06, label= "p-value = 0.0014") 
shc 

#SLA
sla <- lm(sqrt(SLA_mm2.mg)~plant.duration, data = allsummupd3) #p-value: 2.051e-05

par(mfrow=c(2,2))
plot(sla) #437 = Leymus condesatus - perennial (SLA=98.76)
summary(sla)

library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = sla, n = 250)
plot(simulationOutput)

library(ggpmisc)

shd <- ggdensity(allsummupd3, x = "SLA_mm2.mg",
                 add = "mean", rug = TRUE,
                 color = "plant.duration",  fill = "plant.duration",
                 xlab = "SLA (mm²/mg)", 
                 palette = c("#E69F00","#009E73"), 
                 legend.title = "Lifespan") +
  scale_y_continuous(limits = c(0, 0.05), breaks = c(0, 0.02, 0.04)) +
  annotate("text", x=45, y=0.045, label= "p-value < 0.001") +
  theme(axis.title.y = element_blank()) 

shd

#ggarrange(sha, shb, shc, shd, ncol = 2, nrow = 2,  common.legend = T, legend="bottom", labels="AUTO")

# plot grid ##
leg <- get_legend(shd)
grid <- ggarrange(shc, shd, shb, raop2, rc, rb, rf, rd, re, sha, sa, sb, ncol = 3, nrow = 4, common.legend = F, legend.grob = leg, legend="bottom", labels=c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"), align = "hv", vjust=1)

grid 



# after remov crops: SMA regression ####
library(smatr)
#?`smatr-package`
#?sma #sma(y~x*groups)

#Collaboration gradient
collaboration <- allsummupd3%>%
  select(plant.duration,SRL_m.g,Root.diameter_mm) %>%
  na.omit() #%>%dplyr::count(plant.duration) #annual = 26; perennial = 103

coll.test <- sma(SRL_m.g~Root.diameter_mm*plant.duration, log="xy", data=allsummupd3, na.action = na.omit)

plot(coll.test,which="res") #checking assumptions - no pattern
plot(coll.test,which="qq") #Produce a normal quantile plot

#plot collaboration
plot(coll.test, col = c("#E69F00","#009E73"), pch=c(16,1), lty = c(1,2), xlab = "Root diameter (mm)", ylab = "SRL (m/g)") ##annuals = Filled circles and solid lines ; perennials = open circles and dashed lines represent native species

#resultados
coll.test

summary(coll.test)
##test for shift in elevation and shift along the slope
#elevation:
sma(SRL_m.g~Root.diameter_mm+plant.duration, type="elevation", log="xy", data=collaboration)

#shift along a commom axis
sma(SRL_m.g~Root.diameter_mm+plant.duration, type="shift", log="xy", data=collaboration)


#Conservation gradient - Leaf
conservationleaf <- allsummupd3 %>%
  select(plant.duration,LeafN.mass_mg.g,SLA_mm2.mg) %>%
  na.omit() #%>% dplyr::count(plant.duration) #annual = 44; perennial = 103

lcons.test <-sma(LeafN.mass_mg.g~SLA_mm2.mg*plant.duration, log="xy", data=allsummupd3, na.action = na.omit)

par(mfrow=c(1,2))
plot(lcons.test,which="res") #checking assumptions - no pattern
plot(lcons.test,which="qq") #Produce a normal quantile plot

#plot conservation leaf
plot(lcons.test, col=c("#E69F00","#009E73"), pch=c(16,1), lty = c(1,2), xlab = "SLA (mm²/mg)", ylab = "Leaf N (mg/g)") #annuals = Filled circles and solid lines ; perennials = open circles and dashed lines represent native species

#resultados
lcons.test

summary(lcons.test)
##test for shift in elevation and shift along the slope
#elevation:
sma(LeafN.mass_mg.g~SLA_mm2.mg+plant.duration, type="elevation", log="xy", data=conservationleaf)
#shift along a commom axis
sma(LeafN.mass_mg.g~SLA_mm2.mg+plant.duration, type="shift", log="xy", data=conservationleaf)


#Conservation gradient - Root
conservationroot <- allsummupd3 %>%
  select(plant.duration,RootN.mass_mg.g,Root.tissue.density_g.cm3) %>%
  na.omit() #%>% dplyr::count(plant.duration) #annual = 14; perennial = 106

rcons.test <-sma(RootN.mass_mg.g~Root.tissue.density_g.cm3*plant.duration, log="xy", data=allsummupd3, na.action = na.omit)

plot(rcons.test,which="res") #checking assumptions - no pattern
plot(rcons.test,which="qq") #Produce a normal quantile plot

#plot conservation root
plot(rcons.test, col=c("#E69F00","#009E73"), pch=c(16,1), lty = c(1,2), xlab = "RTD (g/cm?)", ylab = "Root N (mg/g)") #annuals = Filled circles and solid lines ; perennials = open circles and dashed lines represent native species

#resultados
rcons.test

summary(rcons.test)
##test for shift in elevation and shift along the slope (IF THEY SHARE A COMMON SLOPE?)
#elevation:
sma(RootN.mass_mg.g~Root.tissue.density_g.cm3+plant.duration, type="elevation", log="xy", data=conservationroot)
#shift along a commom axis
sma(RootN.mass_mg.g~Root.tissue.density_g.cm3+plant.duration, type="shift", log="xy", data=conservationroot)


#Size gradient 
sizegradient <- allsummupd3 %>%
  select(plant.duration,Height_cm,Root.depth_cm) %>%
  na.omit() #%>% dplyr::count(plant.duration) #annual = 16; perennial = 42

size.test <-sma(Height_cm~Root.depth_cm*plant.duration, log="xy", data=allsummupd3, na.action = na.omit)

plot(size.test,which="res") #checking assumptions - no pattern
plot(size.test,which="qq") #Produce a normal quantile plot

#plot size
plot(size.test, col=c("#E69F00","#009E73"), pch=c(16,1), lty = F, ylab = "Height (cm)", xlab = "Root depth (cm)") #annuals = Filled circles and solid lines ; perennials = open circles and dashed lines represent native species


#resultados
size.test

summary(size.test)
##test for shift in elevation and shift along the slope (IF THEY SHARE A COMMON SLOPE?)
#elevation:
sma(Height_cm~Root.depth_cm+plant.duration, type="elevation", log="xy", data=sizegradient)
#shift along a commom axis
sma(Height_cm~Root.depth_cm+plant.duration, type="shift", log="xy", data=sizegradient)



#Seedmass x SRL 
seedsrl <- allsummupd3 %>%
  select(plant.duration,Seed.mass_mg,SRL_m.g) %>%
  na.omit() #%>% dplyr::count(plant.duration) #annual = 21; perennial = 76

seedsrl.test <-sma(Seed.mass_mg~SRL_m.g*plant.duration, log="xy", data=allsummupd3, na.action = na.omit)

plot(size.test,which="res") #checking assumptions - no pattern
plot(size.test,which="qq") #Produce a normal quantile plot

#plot seed root 1
plot(seedsrl.test, col=c("#E69F00","#009E73"), pch=c(16,1), lty = c(F,2), ylab = "Seed.mass (mg)", xlab = "SRL (m/g)") #annuals = Filled circles and solid lines ; perennials = open circles and dashed lines represent native species

#resultados
seedsrl.test

summary(seedsrl.test)
##test for shift in elevation and shift along the slope (IF THEY SHARE A COMMON SLOPE?)
#elevation:
sma(Seed.mass_mg~SRL_m.g+plant.duration, type="elevation", log="xy", data=seedsrl)
#shift along a commom axis
sma(Seed.mass_mg~SRL_m.g+plant.duration, type="shift", log="xy", data=seedsrl)



#Seedmass x root diameter 
seedrd <- allsummupd3 %>%
  select(plant.duration,Seed.mass_mg,Root.diameter_mm) %>%
  na.omit()  #%>% dplyr::count(plant.duration) #annual = 14; perennial = 58

seedrd.test <-sma(Seed.mass_mg~Root.diameter_mm*plant.duration, log="xy", data=allsummupd3, na.action = na.omit)

plot(seedrd.test,which="res") #checking assumptions - no pattern
plot(seedrd.test,which="qq") #Produce a normal quantile plot

#plot seed root 2
plot(seedrd.test, col=c("#E69F00","#009E73"), pch=c(16,1), lty = c(1,2), ylab = "Seed mass (mg)", xlab = "Root diameter (mm)") #annuals = Filled circles and solid lines ; perennials = open circles and dashed lines represent native species

#resultados
seedrd.test

summary(seedrd.test)
##test for shift in elevation and shift along the slope (IF THEY SHARE A COMMON SLOPE?)
#elevation:
sma(Seed.mass_mg~Root.diameter_mm+plant.duration, type="elevation", log="xy", data=seedrd)
#shift along a commom axis
sma(Seed.mass_mg~Root.diameter_mm+plant.duration, type="shift", log="xy", data=seedrd)

#RTD x LDMC 
rtdldmc <- allsummupd3 %>%
  select(plant.duration,Root.tissue.density_g.cm3,LDMC_mg.g) %>%
  na.omit()  #%>% dplyr::count(plant.duration) #annual = 10; perennial = 44

rtdldmc.test <-sma(LDMC_mg.g~Root.tissue.density_g.cm3*plant.duration, log="xy", data=allsummupd3, na.action = na.omit)

plot(rtdldmc.test,which="res") #checking assumptions - no pattern
plot(rtdldmc.test,which="qq") #Produce a normal quantile plot

#plot seed root 2
plot(rtdldmc.test, col=c("#E69F00","#009E73"), pch=c(16,1), lty = F, xlab = "RTD (g/cm³)", ylab = "LDMC (mg/g)") 
#annuals = Filled circles and solid lines ; perennials = open circles and dashed lines represent native species

#resultados
rtdldmc.test

summary(rtdldmc.test)
##test for shift in elevation and shift along the slope (IF THEY SHARE A COMMON SLOPE?)
#elevation:
sma(LDMC_mg.g~Root.tissue.density_g.cm3+plant.duration, type="elevation", log="xy", data=rtdldmc)
#shift along a commom axis
sma(LDMC_mg.g~Root.tissue.density_g.cm3+plant.duration, type="shift", log="xy", data=rtdldmc)


#RTD and SLA -  conservation above and belowgroud
rtdsla <- allsummupd3 %>%
  select(plant.duration,Root.tissue.density_g.cm3,SLA_mm2.mg) %>%
  na.omit() # %>% dplyr::count(plant.duration) #annual = 15; perennial = 57

rtdsla.test <-sma(SLA_mm2.mg~Root.tissue.density_g.cm3*plant.duration, log="xy", data=allsummupd3, na.action = na.omit)

plot(rtdsla.test,which="res") #checking assumptions - no pattern
plot(rtdsla.test,which="qq") #Produce a normal quantile plot

#plot sla x rtd
plot(rtdsla.test, col=c("#E69F00","#009E73"), pch=c(16,1), lty = c(F,2), xlab = "RTD (g/cm³)", ylab = "SLA (mm².mg)") 
#annuals = Filled circles and solid lines ; perennials = open circles and dashed lines represent native species

#resultados
rtdsla.test

summary(rtdsla.test)
##test for shift in elevation and shift along the slope (IF THEY SHARE A COMMON SLOPE?)
#elevation:
sma(SLA_mm2.mg~Root.tissue.density_g.cm3+plant.duration, type="elevation", log="xy", data=rtdsla)
#shift along a commom axis
sma(SLA_mm2.mg~Root.tissue.density_g.cm3+plant.duration, type="shift", log="xy", data=rtdsla)




#plot sma ####
par(mfrow= c(2,4), 
    oma = c(3,1,1,1), #tamanho da area externa
    mar = c(4,4.5,3,1), # Altera margens da figura (inferior, esquerda, superior, direita)
    mgp = c(2.5,1,0), xpd=F)   # Posicionamento dos valores dos eixos (t?tulo, valores, marcadores) 

#conserv
plot(lcons.test, col=c("#E69F00","#009E73"), pch=c(16,1), lty = c(1,2), xlab = "SLA (mm²/mg)", ylab = "Leaf N (mg/g)", cex.lab=1.5, cex=1.5, cex.axis=1.4) #annuals = Filled circles and solid lines ; perennials = open circles and dashed lines represent native species
mtext("(a)", side=3, line=1, cex=1, adj=-0.05, font=2)

#mtext("A)", side=3, line=-0.8, adj=0.05, cex=1.1, font=2, outer=TRUE)

plot(rcons.test, col=c("#E69F00","#009E73"), pch=c(16,1), lty = c(1,2), xlab = "RTD (g/cm³)", ylab = "Root N (mg/g)", cex.lab=1.5, cex=1.5, cex.axis=1.4) #annuals = Filled circles and solid lines ; perennials = open circles and dashed lines represent native species
mtext("(b)", side=3, line=1, cex=1, adj=-0.05, font=2)

#mtext("B)", side=3, line=-0.8, adj=0.4, cex=1.1, font=2, outer=TRUE)


plot(rtdldmc.test, col=c("#E69F00","#009E73"), pch=c(16,1), lty = F, xlab = "RTD (g/cm³)", ylab = "LDMC (mg/g)", cex.lab=1.5, cex=1.5, cex.axis=1.4)
mtext("(c)", side=3, line=1, cex=1, adj=-0.05, font=2)

plot(rtdsla.test, col=c("#E69F00","#009E73"), pch=c(16,1), lty = c(F,2), xlab = "RTD (g/cm³)", ylab = "SLA (mm².mg)", cex.lab=1.5, cex=1.5, cex.axis=1.4)
mtext("(d)", side=3, line=1, cex=1, adj=-0.05, font=2)

#collab
plot(coll.test, col = c("#E69F00","#009E73"), pch=c(16,1), lty = c(1,2), xlab = "Root diameter (mm)", ylab = "SRL (m/g)", cex.lab=1.5, cex=1.5,cex.axis=1.4) ##annuals = Filled circles and solid lines ; perennials = open circles and dashed lines represent native species
mtext("(e)", side=3, line=1, cex=1, adj=-0.05, font=2)

#mtext("C)", side=3, line=-0.8, adj=0.75, cex=1.1, font=2, outer=TRUE)

#COLLAB X SEED MASS
plot(seedsrl.test, col=c("#E69F00","#009E73"), pch=c(16,1), lty = c(F,2), ylab = "Seed mass (mg)", xlab = "SRL (m/g)", cex.lab=1.5, cex=1.5, cex.axis=1.4) #annuals = Filled circles and solid lines ; perennials = open circles and dashed lines represent native species
mtext("(f)", side=3, line=1, cex=1, adj=-0.05, font=2)
#mtext("E)", side=3, line=-20.5, adj=0.4, cex=1.1, font=2, outer=TRUE)

plot(seedrd.test, col=c("#E69F00","#009E73"), pch=c(16,1), lty = F, ylab = "Seed mass (mg)", xlab = "Root diameter (mm)", cex.lab=1.5, cex=1.5, cex.axis=1.4) #annuals = Filled circles and solid lines ; perennials = open circles and dashed lines represent native species
mtext("(g)", side=3, line=1, cex=1, adj=-0.05, font=2)
#mtext("F)", side=3, line=-20.5, adj=0.75, cex=1.1, font=2, outer=TRUE)

#size
plot(size.test, col=c("#E69F00","#009E73"), pch=c(16,1), lty = F, ylab = "Height (cm)", xlab = "Root depth (cm)", cex.lab=1.5, cex=1.5, cex.axis=1.4) #annuals = Filled circles and solid lines ; perennials = open circles and dashed lines represent native species
mtext("(h)", side=3, line=1, cex=1, adj=-0.05, font=2)
#mtext("D)", side=3, line=-20.5, adj=0.05, cex=1.1, font=2, outer=TRUE)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)

plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend(x = "bottom", legend= c("annual", "perennial"),
       cex=1.4, pch = unique(c(16,1)), lty=unique(c(1,2)),
       col= unique(c("#E69F00","#009E73")), 
       bty = "n", pt.cex = 1.5, 
       text.col = "black", text.font=1, horiz = T, inset = c(0,0))

#citation("emmeans")

