#THIS SCRIPT CONTAINS ALL MODELS, CHECKS, AND PLOT FOR THE COMPONENTS OF THIS ANALYSIS
library(ggplot2)
library(dplyr)
library(frair) 
library(lmerTest)
library(lme4)
library(tidyverse)
library(lamW)
library(LambertW)
library(MASS)
library(visreg)
library(ggeffects)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(brglm)


#WORKING DIRECTORY AND FILES----------------------
#setwd("C:/Data/SFU/R")
#File for functional response experiment
clams2 <- read.csv("crab_consumption_raw2.csv",na.strings=c(''), stringsAsFactors = TRUE) %>%
  mutate(period = ifelse(Rep < 7, "first", "second"),
         Sediment = factor(Sediment)) %>%
  rename(Substrate = Sediment)
#Filtering functional response file to exclude a secondary temperature experiment
#and then created a new column specifying whether each individual crab had eaten a single clam in a trial or not
firstperiod <- clams2 %>% filter(period == "first") %>% 
  mutate(eaten_logistic = case_when(eaten == 0 ~ 0,
                                    TRUE ~ 1))
#Only trials with varnish clams
varnishclams <- filter(firstperiod, Clam == "VC")
#Only trials with Japanese littleneck clams
JLNclams <- filter(firstperiod, Clam == "JLN")
#Only trials with varnish clams with no substrate
VC0 <- filter(varnishclams, Substrate == "0")
#Only trials with varnish clams within substrate
VC1 <- filter(varnishclams, Substrate == "1")
#Only trials with Japanese littleneck clams with no substrate
JLN0 <- filter(JLNclams, Substrate == "0")
#Only trials with Japanese littleneck clams within substrate
JLN1 <- filter(JLNclams, Substrate == "1")

#Data set that contains the sizes of every clam used and whether they were consumed or not in a trial
umbo <- read.csv("clam_size.csv",na.strings=c('NA'), stringsAsFactors = TRUE) %>%
  mutate(period = ifelse(Rep < 7, "first", "second")) %>% #specify clams that were used for functional response experiment (first)
  rename(eaten = `Consumed..1..consumed.`, #renaming variables
         length = `Length..siphon.to.foot.`,
         width = `Width..hinge.up.`,
         substrate = Sed..1.yes.,) %>% 
  #Create a new column with whether that clam was eaten or not
  mutate(eaten = factor(eaten,levels = c("0", "1")),
         eaten_yn = case_when(eaten == "0" ~ "Unconsumed", #if it was not eaten = unconsumed
                              TRUE ~ "Consumed")) %>% #if it was eaten = consumed 
  mutate(substrate = factor(substrate,levels = c("0", "1")), #specify whether that clam was buried in substrate or not
         substrate_yn = case_when(substrate == "0" ~ "No substrate",
                                 TRUE ~ "Substrate")) %>%
  filter(!is.na(eaten))

#Only include clams that were used for the functional response experiment 
umbofirst <- filter(umbo, period == "first")
#Only include clams that were eaten during their trial 
umbosecond <- filter(umbofirst, eaten == "1")
#Only varnish clams
VCumbo <- filter(umbo, Species == "VC")
#only Japanese littleneck clams
JLNumbo <- filter(umbo, Species == "JLN")

#Data set for the clam burial depth experiment 
burial <- read.csv("clam_burial.csv",na.strings=c('NA'), stringsAsFactors = TRUE) 
#Varnish clam burial depth only
VCburial <- filter(burial, species == "VC")
#Japanese littleneck clam burial depth only
JLNburial <- filter(burial, species == "JLN")

#FUNCTIONAL RESPONSE MODELS---------------------------------
#Test for Type II/Type III FR curve

#Japanese littleneck clams with no substrate
frair_test(eaten ~ Density, data = JLN0) #no evidence of a FR curve
#Japanese littleneck clams with substrate
frair_test(eaten ~ Density, data = JLN1) #no evidence of a FR curve
#Varnish clams with no substrate
frair_test(eaten ~ Density, data = VC0) #evidence for a Type II curve
#Varnish clams with substrate
frair_test(eaten ~ Density, data = VC1) #no evidence of a FR curve


#Frair fit test for varnish clams with no substrate
fit_VC0<- frair_fit(eaten ~ Density, data = VC0, response = 'rogersII',
                      start = list(a = 0.2, h = 0.2), fixed = list(T=1))

# A linear fit
fitI_VC0 <- frair_fit(eaten ~ Density, data = VC0, response = 'typeI',
                     start = list(a = 0.2), fixed=list(T=1))
#Confirms that it is not a type I

#Visualise fits
plot(fit_VC0, pch=20, col=rgb(0,0,0,0.2), xlim=c(0,16))
lines(fit_VC0)
lines(fitI_VC0, lty=3)

#Parameter estimates (attack rate and handling time) for varnish clams with no substrate
a <- fit_VC0$coefficients[1] #coefficient for attack rate
h <- fit_VC0$coefficients[2] #coefficient for handling time 


#this creates a data frame of our predictions to compare to what we actually got
fit_varnish <- data.frame(x = VC0$Density) #Calculate 'a' for each value of x, where x is density of varnish clams

#calculate expected number of varnish clams eaten, based on frair flexnr equation, using lambert function
fit_varnish$Ne <- fit_varnish$x - lambertW0(a * h * fit_varnish$x * exp(-a * (1 - h * fit_varnish$x)))/(a * h) 
#residual fits
fit_varnish$actual <- VC0$eaten
fit_varnish$resid <- fit_varnish$Ne - fit_varnish$actual

plot(x = fit_varnish$x, y = fit_varnish$Ne)
plot(x = fit_varnish$Ne, y = fit_varnish$resid) #checking residual fits
abline(h = 0, lty = 'dotted')

# Have a look at original fits returned by mle2 
summary(fit_VC0$fit) #checking to see if attack and handling time is signficant or not
summary(fitI_VC0$fit)

# Compare models using AIC
AIC(fitI_VC0$fit,fit_VC0$fit) #type II is for sure better
#this is just to reconfirm that this is a Type II curve, which is confirmed with data, yay

# Bootstrap data 
set.seed(309331)
fit_VC0_boot <- frair_boot(fit_VC0, start = NULL, strata=VC0[,6], nboot=2000,
                            para=TRUE, ncores=NaN, WARN.ONLY=FALSE)

fit_VC0_boot 
confint(fit_VC0_boot)

# Illustrate bootlines
plot(fit_VC0_boot, xlim=c(0,16), ylim = c(0, 16), type='n', main='All bootstrapped lines')
lines(fit_VC0_boot, all_lines=TRUE)
points(fit_VC0_boot, pch=20)

# Illustrate bootpolys
plot(fit_VC0_boot, xlim=c(0,16), ylim = c(0, 16), type='n', main='Empirical 95 percent CI')
drawpoly(fit_VC0_boot, col=hsv(2/6,0.2, 0.8))
points(fit_VC0_boot, pch=20)
lines(fit_VC0_boot, all_lines=FALSE)

#jitter to points
fit_VC0_boot$x <- jitter(fit_VC0_boot$x, amount = 0.1)

#final figure for varnish clams with no substrate
png("Figure1_COLOUR.png",width = 9, height = 8,units = "in", res= 300)
par(bg = 'white', fg = 'black')
plot(fit_VC0_boot, xlim=c(0, 16), ylim = c(0, 16), type='n',
     xlab = "Initial Clam Density",
     ylab="Clams Consumed", 
     cex.lab = 1.5,
     font.lab = 2,
     cex.axis = 1.2,
     cex.main = 1.5)
lines(fit_VC0_boot, lwd = 3, all_lines=FALSE, col= "#02AE23", lty = 2)
drawpoly(fit_VC0_boot, border = NA, col=adjustcolor("#02AE23", alpha.f = 0.4))
points(fit_VC0_boot, pch=17, col=adjustcolor("#02AE23", alpha.f = 0.4), cex = 1.4)
dev.off()

#GENERALISED LINEAR MODEL FOR THE PROBABILITY OF A CRAB CONSUMING A CLAM OR NOT---------------

#As every crab in the varnish clams without substrate ate at least one clam during a trial, the probability of them consuming a clam is 100% across all trials
#Therefore we need to use bias reduction in binomial-response generalized liner model
#This allows us to use the maximum penalized likelihood, where penalization is done using Jeffreys invariant prior

#generalized logistic binomial model to examine the effects of Cheliped (cheliped height), Density (initial density of clams),
#the interaction between clam (clam species) and substrate (presence/absence of substrate for a trial), and Cheliped:Clam on the probability of a clam being consumed
hurdle_mod3.5 <- brglm(eaten_logistic ~ Cheliped + Density + Clam*Substrate + Cheliped:Clam, 
                       family = binomial(link = "logit"),
                       data = firstperiod)
summary(hurdle_mod3.5) #cheliped height, density, clam species, and the interaction between clam and substrate are signficant
#Examine residuals using DHARMa
simulateResiduals(hurdle_mod3.5, plot = TRUE) #looks good


#Plot of the probability of a clam being eaten model
hurdle3.5_predict <- ggpredict(hurdle_mod3.5, 
                             terms = c("Density[n = 100]", "Substrate","Clam")) %>% 
  rename(Density = x,
         Substrate = group,
         Clam = facet) 

#Visualise the model
colclam <- c("#25B3C6","#5A32FC","#02AE23","#EF3929")
clamlabel <- c("JLN, control","JLN, substrate", "VC, control", "VC, substrate")

# New facet label names for treatment combinations
trt.labs <- list(
  'JLN:0' = "Japanese littleneck clam, no substrate", 
  'JLN:1' = "Japanese littleneck clam, substrate", 
  'VC:0' = "Varnish clam, no substrate",
  'VC:1' = "Varnish clam, substrate")

trt_labeller <- function(variable,value){
  return(trt.labs[value])
}

#final plot
ggplot(hurdle3.5_predict, aes(Density, predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Clam:Substrate), show.legend = FALSE, alpha = 0.3) +
  geom_line(aes(color = Clam:Substrate), show.legend = FALSE,) +
  geom_point(data = firstperiod, 
             aes(Density, eaten_logistic, color = Clam:Substrate),show.legend = FALSE,
             size = 1, position = position_jitter(height = 0.01)) + 
  scale_colour_manual(values = colclam, labels = clamlabel) +
  scale_fill_manual(values = colclam, labels = clamlabel, drop = FALSE) +
  facet_wrap(~Clam:Substrate, labeller = trt_labeller) +
  labs(x = "Density",
       y = "Probability of consuming a clam") +
  theme_classic(base_size = 15) 
#ggsave("Figure2_COLOUR.png",width = 9, height = 8)

#GENERALISED LINEAR MODEL FOR THE PROPORTION OF CLAMS CONSUMED-------------------
#binomial logistic regression model and logit link 

#As previous studies have shown that cheliped height has shown to have a signficant effect on the prey that a crab can consume
#we wanted to see if cheliped height would effect the proportion of clam consumed during a trial
#Logistic regression model to examine the effects of Cheliped (cheliped height), Density (initial density of clams),
#Clam*Substrate (the interaction between clam species and substrate presence/absence of substrate for a trial), 
#and Cheliped:Clam on the probability of a clam being consumed
model4.5 <- glmmTMB(proportion_eaten ~ Cheliped + Density + Clam*Substrate+ Cheliped:Clam,
                    family = binomial(link = "logit"),
                    weights = Density, 
                    data = firstperiod)
summary(model4.5)
#Examine residuals using DHARMa
plot(simulateResiduals(model4.5))
testDispersion(model4.5 )
testOutliers(model4.5)
testZeroInflation(model4.5)

#Plot of the proportion of clams consumed model 
model_predict <- ggpredict(model4.5, 
                           terms = c("Cheliped[n = 100]", "Substrate", "Clam")) %>% 
  rename(Cheliped = x,
         Substrate = group,
         Clam = facet) 

colclam <- c("#25B3C6","#5A32FC","#02AE23","#EF3929")
clamlabel <- c("JLN, control","JLN, substrate", "VC, control", "VC, substrate")

# New facet label names for treatment combinations
trt.labs <- list(
  'JLN:0' = "Japanese littleneck clam, no substrate", 
  'JLN:1' = "Japanese littleneck clam, substrate", 
  'VC:0' = "Varnish clam, no substrate",
  'VC:1' = "Varnish clam, substrate")

trt_labeller <- function(variable,value){
  return(trt.labs[value])
}

#Final plot
ggplot(model_predict, aes(Cheliped, predicted)) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Clam:Substrate),show.legend = FALSE, alpha = 0.3) +
  geom_line(aes(colour = Clam:Substrate), show.legend = FALSE,) +
  geom_point(data = firstperiod, 
             aes(Cheliped, proportion_eaten, color = Clam:Substrate),
             show.legend = FALSE,
             size = 1, position = position_jitter(height = 0.01)) +
  scale_colour_manual(values = colclam, labels = clamlabel) +
  scale_fill_manual(values = colclam, labels = clamlabel) + 
  facet_wrap(~Clam:Substrate,
             labeller = trt_labeller) +
  labs(x = "Cheliped height (mm)", y = "Proportion consumed") +
  theme_classic(base_size = 15) 
#ggsave("Figure3.png",width = 9, height = 8)

#SUPPLEMENTAL FIGURES----------------

#BURIAL DEPTH FOR BOTH CLAM SPECIES-------------------
#Burial depth model
burialmod <- lm(distance_buried ~ species, data = burial)
summary(burialmod)
plot(simulateResiduals(burialmod))

#Function to get standard error
standard_error <- function(x) sd(x) / sqrt(length(x))

#Mean and SE for burial depth of varnish clams
summary(VCburial)
standard_error(VCburial$distance_buried)
#Mean and SE for burial depth of Japanese littleneck clams
summary(JLNburial)
standard_error(JLNburial$distance_buried)

#Plot
ggplot(data = burial, 
       aes(x = species, y = distance_buried),
       size = 3) +  
  stat_summary(data = burial, fun = mean, size = 1) +
  stat_summary(data = burial, fun.data = mean_cl_normal,
               geom = "errorbar", width = 0.1) +
  labs(x = "", y = "Burial depth of clams (cm)") +
  geom_jitter(data = burial, aes (x = species, y = distance_buried),
              alpha = 0.2, height=0) +
  scale_x_discrete(labels = c("Japanese littleneck clam", "Varnish clam")) +
  theme_classic(base_size = 17) +
  ylim(11.5,0) 
#ggsave("FigureS1.png",width = 9, height = 9)

#SIZE RANGE OF CLAMS THAT WERE/WERE NOT CONSUMED ACROSS BOTH SUBSTRATE TREATMENTS--------------

#Comparing VC with eaten vs not eaten in substrate vs no substrate
#clammod2 <- glmmTMB(length ~ substrate*eaten, dispformula = ~eaten,data = VCumbo)
clammod2 <- glmmTMB(length ~ substrate*eaten, data = VCumbo)
summary(clammod2)
#test Residuals
simulation <- simulateResiduals(clammod2)
plot(simulation) 

emmeans(clammod2, pairwise ~ substrate| eaten)$contrasts
emmeans(clammod2, pairwise ~ eaten| substrate)$contrasts

#Comparing JLN with eaten vs not eaten in substrate vs no substrate
clammod3 <- glmmTMB(length ~ substrate*eaten,  data = JLNumbo)
summary(clammod3)
#test Residuals
simulation2 <- simulateResiduals(clammod3)
plot(simulation2) 

emmeans(clammod3, pairwise ~ substrate| eaten)$contrasts
emmeans(clammod3, pairwise ~ eaten| substrate)$contrasts

# New facet label names for clam names
clam.labs <- c("Japanese littleneck clam", "Varnish clam")
names(clam.labs) <- c("JLN", "VC")

#Plot
ggplot(data = umbo, 
       aes(x = eaten_yn, y = length),show.legend = FALSE,
       size = 3) +  
  facet_grid(Species~substrate_yn, 
             labeller = labeller(Species = clam.labs)) +
  stat_summary(data = umbo, fun = mean, size = 1, show.legend = FALSE) +
  stat_summary(data = umbo, fun.data = mean_cl_normal,
               geom = "errorbar", width = 0.1, show.legend = FALSE) +
  labs(x = "", y = "Length of clam (mm)") +
  geom_jitter(data = umbo, aes (x = eaten_yn, y = length),
              show.legend = FALSE,
              alpha = 0.2, height=0) +
  theme_bw(base_size = 17) +
  theme(panel.grid = element_blank())
#ggsave("FigureS2.png",width = 9, height = 7)
#ggsave("FigureS2_bw.png",width = 9, height = 7)

#PLOT OF RAW DATA FOR ALL FUNCTIONAL RESPONSE EXPERIMENTS 
#This is to provide a visualisation for the reason behind a lack of evidence for a FR type
#for the three treatment combinations (varnish clams with substrate and Japanese littleneck clams with/without substrate)

clamlabel <- c("JLN, control","JLN, substrate", "VC, control", "VC, substrate")

# New facet label names for treatment combinations
trt.labs <- list(
  'JLN:0' = "Japanese littleneck clam, no substrate", 
  'JLN:1' = "Japanese littleneck clam, substrate", 
  'VC:0' = "Varnish clam, no substrate",
  'VC:1' = "Varnish clam, substrate")

trt_labeller <- function(variable,value){
  return(trt.labs[value])
}

#Plot
ggplot(firstperiod, aes(Density, eaten)) + 
  geom_point(data = firstperiod, 
             aes(Density, eaten),show.legend = FALSE,
             size = 2, position = position_jitter(height = 0.05), alpha = 0.2) + 
  #scale_colour_manual(values = colclam, labels = clamlabel) +
  #scale_fill_manual(labels = clamlabel) +
  facet_wrap(~Clam:Substrate, 
             labeller = trt_labeller) +
  labs(x = "Initial clam density",
       y = "Clams consumed") +
  theme_classic(base_size = 17)
#ggsave("FigureS3_bw.png",width = 9, height = 9)

###ADDITIONAL ANALYSIS 
#Number of crabs that did not consume a single clam
#Varnish clams, no substrate
sum(VC0$eaten == "0") #0
#Varnish clams, substrate
sum(VC1$eaten == "0") #24
#Japanese littleneck clams, no substrate
sum(JLN0$eaten == "0") #26
#Japanese littleneck clams, substrate
sum(JLN1$eaten == "0") #26

#Distance buried as a function of clam size and species
burialmod1 <- glmmTMB(distance_buried ~ clam_length*species, data = burial)
summary(burialmod1)

