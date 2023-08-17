#Graphs for chapter 1-----------------------------
library(ggplot2)
library(lmerTest)
library(lme4)
library(dplyr)
library(tidyverse)
library(lamW)
library(LambertW)
library(frair)
library(visreg)
library(ggeffects)
library(car)
library(DHARMa)
library(hrbrthemes)
library(viridis)
library(glmmTMB)
library(MASS)
library(tweedie)
library(statmod)
library(stringr)

#setwd("C:/Data/SFU/R")
clams2 <- read.csv("crab_consumption_raw2.csv",na.strings=c(''), stringsAsFactors = TRUE) %>%
  mutate(period = ifelse(Rep < 7, "first", "second"),
         Sediment = factor(Sediment)) 
clam_temp <- read.csv("crab_temperature.csv",na.strings=c(''), stringsAsFactors = TRUE)
clams <- merge(clams2, clam_temp, by="Date") %>% 
  mutate(mean_temp = (Temp_Start+Temp_End)/2,
         mean_salinity = (Salinity_Start..ppt.+Salinity_End )/2,
         period = ifelse(Rep < 7, "first", "second"),
         Sediment = factor(Sediment))
clamscale <- read.csv("crab_consumption_raw2.csv",na.strings=c(''), stringsAsFactors = TRUE) %>%
  mutate(period = ifelse(Rep < 7, "first", "second"),
         Sediment = factor(Sediment),
         scaled_cheliped = scale(Cheliped)[,1],
         scaled_density = scale(Density)[,1],
         scaled_CW = scale(CW)[,1])

period1 <- filter(clams, period == "first")
JLNprop <- filter(period1, Clam == "JLN")  %>% 
  mutate(scaled_cheliped = scale(Cheliped)[,1])
VCprop <- filter(period1, Clam == "VC") %>% 
  mutate(scaled_cheliped = scale(Cheliped)[,1],
         scaled_density = scale(Density)[,1],
         scaled_temp = scale(mean_temp)[,1])
umbo <- read.csv("clam_size.csv",na.strings=c('NA'), stringsAsFactors = TRUE) %>%
  mutate(period = ifelse(Rep < 7, "first", "second"))
varnish <- filter(umbo, Species == "VC", period == "first") 
varnish2 <- filter(umbo, Species == "VC", period == "second")
littleneck <- filter(umbo, Species == "JLN")

#per clam df
#always run the NULL line before the for loop!!!!!
clamview <- NULL
for(i in 1:nrow(VCprop)){
  dens <- VCprop$Density[i]
  dead <- VCprop$eaten[i]
  alive <- VCprop$alive[i]
  df <- tibble(eaten = c(rep(1, times = dead), rep(0, times = alive)),
               Cheliped = rep(VCprop$Cheliped[i], times = dens),
               mean_temp = rep(VCprop$mean_temp[i], times = dens),
               Density = rep(VCprop$Density[i], times = dens),
               Sediment = rep(VCprop$Sediment[i], times = dens),
               crab_num = rep(as.factor(i), times = dens))
  clamview <- bind_rows(clamview, df)
}
clamview <- clamview %>% 
  mutate(scaled_cheliped = c(scale(Cheliped)),
         scaled_density = c(scale(Density)),
         scaled_temp = c(scale(mean_temp)))

#Raw data as ggplot---------------------------
# %>% stands for then
ggplot(clams2) + geom_jitter(aes(x = Density, y = eaten, colour = Sediment))
#year3 <- lm(proportion_eaten ~ Cheliped + Sediment, data = VCprop)
#year3_glm <- glm(cbind(eaten,alive) ~ Cheliped + Sediment, data = VCprop, family = binomial())
clams_glm <- glm(proportion_eaten ~ scaled_cheliped * Sediment, 
                         weights = Density, 
                         data = VCprop, 
                         family = binomial(link = "logit"))

colallclam <- c("JLN_0" = "#C2EAFF",
             "JLN_1" = "#004266", 
             "VC_0" = "#ADEBD2", 
             "VC_1" = "#186345")
clamalllabel <- c("JLN_0" = "JLN, no sed", 
               "JLN_1" = "JLN, with sed", 
               "VC_0" = "VC, no sed", 
               "VC_1" = "VC, with sed")

ggplot(data = clams2, aes(x = Density, y = eaten, 
                                             colour = Sediment, fill = Sediment)) +
  geom_ribbon(aes(ymin = conf.low, ymax=conf.high), alpha = 0.3, colour = NA) +
  geom_line() +
  geom_jitter(data = clams2, aes(x = Density, y = eaten)) +
  #scale_colour_manual deals with points and lines
  #scale_fill_manial deals with ribbons and bars and other things to fill in
  scale_colour_manual(values = colallclam, labels = clamalllabel) +
  scale_fill_manual(values = colallclam, labels = clamalllabel) +
  theme_classic() +
  labs(x = "Density of prey", y = "Number of clams eaten")



#clams <- read_csv("crab_consumption_raw.csv") %>% 
  #unite("combo", c(Clam,Sediment), remove = FALSE) %>% 
  #mutate(Clam=factor(Clam),
        #Sediment = factor(Sediment), 
         #we used rev so we didnt have to re write the levels in the right order cause we lazy af
         #combo = factor(combo, levels = rev(c("JLN_0", "JLN_1", "VC_0", "VC_1")))) 
#used the coolers website for the colours cause we fancy
#colclam <- c("JLN_0" = "#C2EAFF", JLN_1" = "#004266","VC_0" = "#ADEBD2", "VC_1" = "#186345")
#clamlabel <- c("JLN_0" = "JLN, no sed", "JLN_1" = "JLN, with sed", "VC_0" = "VC, no sed", "VC_1" = "VC, with sed")

ggplot(clams, aes(Density, eaten, colour = combo, fill = combo)) +
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE) +
  #scale_colour_manual deals with points and lines
  #scale_fill_manial deals with ribbons and bars and other things to fill in
  scale_colour_manual(values = colclam, labels = clamlabel) +
  scale_fill_manual(values = colclam, labels = clamlabel) +
  theme_classic() +
  labs(x = 'Density (clams/bin)', 
       y = "Clams Consumed", 
       fill = "Treatment", colour = "Treatment")


#Clam burial depth-------------
clam_burial <- read.csv("clam_burial.csv",na.strings=c(''), stringsAsFactors = TRUE)

ggplot(clam_burial) + geom_boxplot(aes(x = species, y = distance_buried)) +
  ylab ("Distance clam buried (cm)")

burial1 <- t.test(distance_buried~species, data = clam_burial)
burial1

summary(clam_burial)

#Clam size range-----------------------


#VC length and width (first period)
mean(varnish$ Length..siphon.to.foot)
standard_error <- function(x) sd(x) / sqrt(length(x)) # Create own function
standard_error(varnish$ Length..siphon.to.foot)

mean(varnish$Width..hinge.up.)
standard_error(varnish$Width..hinge.up.)

range(varnish$Width..hinge.up.)
range(varnish$Width..hinge.up.)

#VC length and width (second period)
mean(varnish2$ Length..siphon.to.foot)
standard_error <- function(x) sd(x) / sqrt(length(x)) # Create own function
standard_error(varnish2$ Length..siphon.to.foot)

mean(varnish2$Width..hinge.up.)
standard_error(varnish2$Width..hinge.up.)

range(varnish2$Width..hinge.up.)
range(varnish2$Width..hinge.up.)

#JLN (length and width)
mean(littleneck$Length)
standard_error <- function(x) sd(x) / sqrt(length(x)) # Create own function
standard_error(littleneck$Length)

mean(littleneck$Width)
standard_error(littleneck$Width)

range(little$Length)
range(little$Width)

#Clam size range graph-------------------

#Varnish clams comparing consumption/not consumed sizes 
ggplot(varnish %>% filter(!is.na(Consumed..1..consumed.)), aes(x=as.factor(Consumed..1..consumed.), y=Length..siphon.to.foot.)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  xlab("Consumption of clam") + ylab("Length of clam (siphon to foot)")

ggplot(varnish %>% filter(!is.na(Consumed..1..consumed.)), aes(x=as.factor(Consumed..1..consumed.), y=Length..siphon.to.foot.)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)) +
  xlab("Consumption of clam") + ylab("Length of clam (siphon to foot)")

VCsize_summary <- t.test(Length..siphon.to.foot.~ Consumed..1..consumed., data = varnish)
VCsize_summary
#no signficant difference between sizes

#VC consumption AND sediment graph


#Littleneck clams comparing consumption sizes 
ggplot(littleneck %>% filter(!is.na(Consumed..1..consumed.)), aes(x=as.factor(Consumed..1..consumed.), y=Length..siphon.to.foot.)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  xlab("Consumption of clam") + ylab("Length of clam (siphon to foot)")

ggplot(littleneck %>% filter(!is.na(Consumed..1..consumed.)), aes(x=as.factor(Consumed..1..consumed.), y=Length..siphon.to.foot.)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)) +
  xlab("Consumption of clam") + ylab("Length of clam (siphon to foot)")

JLNsize_summary <- t.test(Length..siphon.to.foot.~ Consumed..1..consumed., data = littleneck)
JLNsize_summary
#significant difference between sizes

#plot(x = VC$Consumed..1..consumed., y = umbo$Length..siphon.to.foot.,xlab = "Clam species", ylab = "Length (siphon to foot)", names = c("VC", "JLN"), las = 1)
#consumed <- filter(umbo, Consumed..1..consumed. == "1")
#plot(x = consumed$Species, y = consumed$Length..siphon.to.foot., xlab = "Clam species", ylab = "Length (siphon to foot)", names = c("VC", "JLN"), las = 1)
#plot(x = consumed$Species, y = consumed$Width..hinge.up., xlab = "Clam species", ylab = "Width", names = c("VC", "JLN"), las = 1)

#notconsumed <- filter(umbo, Consumed..1..consumed. == "0")
#plot(x = notconsumed$Species, y = notconsumed$Length..siphon.to.foot., xlab = "Clam species", ylab = "Length (siphon to foot)", names = c("VC", "JLN"), las = 1)
#plot(x = notconsumed$Species, y = notconsumed$Width..hinge.up., xlab = "Clam species", ylab = "Width", names = c("VC", "JLN"), las = 1)

#plot(x = umbo$Species, y = umbo$Length..siphon.to.foot., xlab = "Clam species", ylab = "Length (siphon to foot)", names = c("VC", "JLN"), las = 1)
#plot(x = umbo$Species, y = umbo$Width..hinge.up., xlab = "Clam species", ylab = "Width", names = c("VC", "JLN"), las = 1)

#claw size in relation to clam predation ------------------

#Varnish glams GLM
ggplot(VCprop) + geom_jitter(aes(x = Cheliped, y = proportion_eaten, colour = Sediment))
#year3 <- lm(proportion_eaten ~ Cheliped + Sediment, data = VCprop)
#year3_glm <- glm(cbind(eaten,alive) ~ Cheliped + Sediment, data = VCprop, family = binomial())
year3_glm_weights <- glm(proportion_eaten ~ scaled_cheliped * Sediment, 
                         weights = Density, 
                         data = VCprop, 
                         family = binomial(link = "logit"))
#summary(year3)
#summary(year3_glm)
summary(year3_glm_weights) #sediment and interaction has a signficant effect
CHmod_resid <- simulateResiduals(year3_glm_weights)
plot(CHmod_resid)
visreg(year3_glm_weights, scale = "response")

#plot output
colclam <- c("0" = "#ADEBD2", 
             "1" = "#186345")
clamlabel <- c("0" = "VC, no sed", 
               "1" = "VC, with sed")
year3_glm_weights_predict <- ggpredict(year3_glm_weights, 
                                       terms = c("scaled_cheliped", "Sediment")) %>% 
  rename(scaled_cheliped = x,
         Sediment = group) %>% 
  mutate(Cheliped = (scaled_cheliped*sd(VCprop$Cheliped))+mean(VCprop$Cheliped))

ggplot(data = year3_glm_weights_predict, aes(x = Cheliped, y = predicted, 
                                             colour = Sediment, fill = Sediment)) +
  geom_ribbon(aes(ymin = conf.low, ymax=conf.high), alpha = 0.3, colour = NA) +
  geom_line() +
  geom_jitter(data = VCprop, aes(x = Cheliped, y = proportion_eaten)) +
  #scale_colour_manual deals with points and lines
  #scale_fill_manial deals with ribbons and bars and other things to fill in
  scale_colour_manual(values = colclam, labels = clamlabel) +
  scale_fill_manual(values = colclam, labels = clamlabel) +
  theme_classic() +
  labs(x = "Cheliped (cm)", y = "Proportion of clams eaten")



#JLN GLM
#ggplot(JLNprop) + geom_jitter(aes(x = Cheliped, y = proportion_eaten, colour = Sediment))
#year2 <- lm(proportion_eaten ~ Cheliped + Sediment, data = JLNprop)
#summary(year2)

ggplot(JLNprop) + geom_jitter(aes(x = Cheliped, y = proportion_eaten, colour = Sediment))
#year3 <- lm(proportion_eaten ~ Cheliped + Sediment, data = VCprop)
#year3_glm <- glm(cbind(eaten,alive) ~ Cheliped + Sediment, data = VCprop, family = binomial())
year2_glm_weights <- glm(proportion_eaten ~ scaled_cheliped * Sediment, 
                         weights = Density, 
                         data = JLNprop, 
                         family = binomial(link = "logit"))
#summary(year3)
#summary(year3_glm)
summary(year2_glm_weights) #cheliped has a signficant effect
visreg(year2_glm_weights, scale = "response")


#plot output
colJclam <- c("0" = "#C2EAFF",
              "1" = "#004266")
clamJlabel <- c("0" = "JLN, no sand", 
                "1" = "JLN, with sand")
year2_glm_weights_predict <- ggpredict(year2_glm_weights, 
                                       terms = c("scaled_cheliped", "Sediment")) %>% 
  rename(scaled_cheliped = x,
         Sediment = group) %>% 
  mutate(Cheliped = (scaled_cheliped*sd(JLNprop$Cheliped))+mean(JLNprop$Cheliped))

ggplot(data = year2_glm_weights_predict, aes(x = Cheliped, y = predicted, 
                                             colour = Sediment, fill = Sediment)) +
  geom_ribbon(aes(ymin = conf.low, ymax=conf.high), alpha = 0.3, colour = NA) +
  geom_line() +
  geom_jitter(data = JLNprop, aes(x = Cheliped, y = proportion_eaten)) +
  #scale_colour_manual deals with points and lines
  #scale_fill_manial deals with ribbons and bars and other things to fill in
  scale_colour_manual(values = colJclam, labels = clamJlabel) +
  scale_fill_manual(values = colJclam, labels = clamJlabel) +
  theme_classic() +
  labs(x = "Cheliped (cm)", y = "Proportion of clams eaten")

#claw size and proportion eaten
modelclaw <- lm(proportion_eaten ~ Cheliped, 
             data = VCprop)
plot(modelclaw)
summary(modelclaw)

modclaw_resid <- simulateResiduals(modelclaw)
plot(modclaw_resid)


#claw size and VC proportion eaten
#Varnish glams GLM
ggplot(VCprop) + geom_jitter(aes(x = Cheliped, y = proportion_eaten, colour = Sediment))
#year3 <- lm(proportion_eaten ~ Cheliped + Sediment, data = VCprop)
#year3_glm <- glm(cbind(eaten,alive) ~ Cheliped + Sediment, data = VCprop, family = binomial())
lm_claw <- glm(proportion_eaten ~ Cheliped*Sediment, 
                         weights = Density, 
                         data = VCprop, 
                         family = binomial(link = "logit"))
#summary(year3)
#summary(year3_glm)
summary(lm_claw) 
clawmod_resid <- simulateResiduals(lm_claw)
plot(clawmod_resid)
visreg(year3_glm_weights, scale = "response")

year3_glm_weights_predict <- ggpredict(lm_claw, 
                                       terms = c("Cheliped", "Sediment")) %>% 
  mutate(Cheliped = x,
         Sediment = group) 

#plot output
colclam <- c("0" = "#ADEBD2", 
             "1" = "#186345")
clamlabel <- c("0" = "VC, no sed", 
               "1" = "VC, with sed")

ggplot(data = year3_glm_weights_predict, aes(x = Cheliped, y = predict, 
                                             colour = Sediment, fill = Sediment)) +
  geom_ribbon(aes(ymin = conf.low, ymax=conf.high), alpha = 0.3, colour = NA) +
  geom_line() +
  geom_jitter(data = VCprop, aes(x = Cheliped, y = proportion_eaten)) +
  #scale_colour_manual deals with points and lines
  #scale_fill_manial deals with ribbons and bars and other things to fill in
  scale_colour_manual(values = colclam, labels = clamlabel) +
  scale_fill_manual(values = colclam, labels = clamlabel) +
  theme_classic() +
  labs(x = "Cheliped (cm)", y = "Proportion of clams eaten")




#proportion eaten --------------------

#VC proportion
ggplot(VCprop) + geom_jitter(aes(x = Density, y = proportion_eaten, colour = Sediment))
VCprop_glm_weights <- glm(proportion_eaten ~ Density * Sediment, 
                         weights = Density, 
                         data = VCprop, 
                         family = binomial(link = "logit"))
summary(VCprop_glm_weights) #Density,sediment, and Density:Sediment have a significant effect
visreg(year3_glm_weights, scale = "response")

#plot output
colclam <- c("0" = "#ADEBD2", 
             "1" = "#186345")
clamlabel <- c("0" = "VC, no sed", 
               "1" = "VC, with sed")
VCprop_glm_weights_predict <- ggpredict(VCprop_glm_weights, 
                                       terms = c("Density", "Sediment")) %>% 
  rename(Density = x,
         Sediment = group) 

ggplot(data = VCprop_glm_weights_predict, aes(x = Density, y = predicted, 
                                             colour = Sediment, fill = Sediment)) +
  geom_ribbon(aes(ymin = conf.low, ymax=conf.high), alpha = 0.3, colour = NA) +
  geom_line() +
  geom_jitter(data = VCprop, aes(x = Density, y = proportion_eaten)) +
  #scale_colour_manual deals with points and lines
  #scale_fill_manial deals with ribbons and bars and other things to fill in
  scale_colour_manual(values = colclam, labels = clamlabel) +
  scale_fill_manual(values = colclam, labels = clamlabel) +
  theme_classic() +
  labs(x = "Density of clams (clams/bin)", y = "Proportion of clams eaten")

#JLN proportion
ggplot(JLNprop) + geom_jitter(aes(x = Density, y = proportion_eaten, colour = Sediment))
JLNprop_glm_weights <- glm(proportion_eaten ~ Density * Sediment, 
                         weights = Density, 
                         data = JLNprop, 
                         family = binomial(link = "logit"))
summary(JLNprop_glm_weights) #no significant effects
visreg(JLNprop_glm_weights, scale = "response")
JLNpropmod <- simulateResiduals(JLNprop_glm_weights)
plot(JLNpropmod)


#plot output
colJclam <- c("0" = "#C2EAFF",
              "1" = "#004266")
clamJlabel <- c("0" = "JLN, no sand", 
                "1" = "JLN, with sand")
JLNprop_glm_weights_predict <- ggpredict(JLNprop_glm_weights, 
                                       terms = c("Density", "Sediment")) %>% 
  rename(Density = x,
         Sediment = group) 

ggplot(data = JLNprop_glm_weights_predict, aes(x = Density, y = predicted, 
                                             colour = Sediment, fill = Sediment)) +
  geom_ribbon(aes(ymin = conf.low, ymax=conf.high), alpha = 0.3, colour = NA) +
  geom_line() +
  geom_jitter(data = JLNprop, aes(x = Density, y = proportion_eaten)) +
  #scale_colour_manual deals with points and lines
  #scale_fill_manial deals with ribbons and bars and other things to fill in
  scale_colour_manual(values = colJclam, labels = clamJlabel) +
  scale_fill_manual(values = colJclam, labels = clamJlabel) +
  theme_classic() +
  labs(x = "Cheliped (cm)", y = "Proportion of clams eaten")

#proportion all
firstperiod <- clams %>% filter(period == "first")

firstperiod2 <- clams %>% filter(period == "first") %>% 
  unite("combo", c(Clam,Sediment), remove = FALSE) %>% 
  mutate(Clam=factor(Clam),
  Sediment = factor(Sediment), 
  combo = factor(combo, levels = rev(c("JLN_0", "JLN_1", "VC_0", "VC_1")))) 

colclam <- c("JLN_0" = "#C2EAFF", "JLN_1" = "#004266","VC_0" = "#ADEBD2", "VC_1" = "#186345")
clamlabel <- c("JLN_0" = "JLN, no sed", "JLN_1" = "JLN, with sed", "VC_0" = "VC, no sed", "VC_1" = "VC, with sed")
  

ggplot(firstperiod, aes(Density, proportion_eaten, colour = c, fill = combo)) +
  geom_point() + 
  geom_smooth(method = "lm", se = TRUE) +
  #scale_colour_manual deals with points and lines
  #scale_fill_manial deals with ribbons and bars and other things to fill in
  scale_colour_manual(values = colclam, labels = clamlabel) +
  scale_fill_manual(values = colclam, labels = clamlabel) +
  theme_classic() +
  labs(x = 'Density (clams/bin)', 
       y = "Clams Consumed", 
       fill = "Treatment", colour = "Treatment")

#proportion with hurdle model
firstperiod <- clams %>% filter(period == "first")
crab_consump_logistic <- firstperiod %>% 
  mutate(eaten_logistic = case_when(eaten == 0 ~ 0,
                                    TRUE ~ 1))
crab_consump_gamma <- firstperiod %>% 
  filter(eaten != 0)

hurdle_mod1 <- glm(eaten_logistic ~ Sediment + Cheliped + 
                     Clam + Density + mean_temp, 
                   family = binomial,
                   data = crab_consump_logistic)

crab_consump_logistic$combo <- paste(crab_consump_logistic$Clam, crab_consump_logistic$Sediment)

colclam <- c("JLN 0" = "#C2EAFF", "JLN 1" = "#004266","VC 0" = "#ADEBD2", "VC 1" = "#186345")
clamlabel <- c("JLN 0" = "JLN, no sed", "JLN 1" = "JLN, with sed", "VC 0" = "VC, no sed", "VC 1" = "VC, with sed")

hurdle_glm_predict <- ggpredict(hurdle_mod1, 
                                        terms = c("Density", "Sediment")) %>% 
  rename(Density = x,
         Sediment = group) 

ggplot(data = hurdle_mod1, aes(x = Density, y = eaten, 
                                              colour = combo, fill = combo)) +
  geom_ribbon(aes(ymin = conf.low, ymax=conf.high), alpha = 0.3, colour = NA) +
  geom_line() +
  geom_jitter(data = hurdle_mod1, aes(x = Density, y = eaten)) +
  scale_colour_manual(values = colclam, labels = clamlabel) +
  scale_fill_manual(values = colclam, labels = clamlabel) +
  theme_classic() +
  labs(x = 'Density (clams/bin)', 
       y = "Clams Consumed")

ggplot(crab_consump_logistic, aes(x = Density, y = eaten, 
                                      colour = combo, fill = combo)) + 
  geom_point() +
  
  
#Models--------------------------------
#VC
VC_full <- glm(cbind(eaten,alive) ~ Density + Cheliped + Sediment, 
               data = VCprop, family = binomial(link = "logit"))
den_che <- glm(cbind(eaten,alive) ~ Density + Cheliped,  data = VCprop, 
               family = binomial(link = logit))
den_sca_che <- den_che <- glm(cbind(eaten,alive) ~ Density + scaled_cheliped,  data = VCprop, 
                              family = binomial(link = logit))
den_sed <- glm(cbind(eaten,alive) ~ Density + Sediment,  data = VCprop, 
              family = binomial(link = logit))
den <- glm(cbind(eaten,alive) ~ Density,  data = VCprop, 
           family = binomial(link = logit))
int <- glm(cbind(eaten,alive) ~ 1,  data = VCprop, 
           family = binomial(link = logit))
cara <-  glm(cbind(eaten,alive) ~ Density + CW + Sediment,  data = VCprop, 
             family = binomial(link = logit))
AIC(VC_full, den_che, den_sca_che, den_sed, den, int, cara)

###The best model is the one that includes Density and Sediment, closely followed by Density, Sediment and CW
par(mfrow = c(2,2))
plot(den_sed)

#JLN
JLN_full <- glm(cbind(eaten,alive) ~ Density + Cheliped + Sediment, 
               data = JLNprop, family = binomial(link = "logit"))
den2_che <- glm(cbind(eaten,alive) ~ Density + Cheliped,  data = JLNprop, 
               family = binomial(link = logit))
den2_sca_che <- den_che <- glm(cbind(eaten,alive) ~ Density + scaled_cheliped,  
                               data = JLNprop, 
                              family = binomial(link = logit))
den2_sed <- glm(cbind(eaten,alive) ~ Density + Sediment,  data = JLNprop, 
               family = binomial(link = logit))
den2 <- glm(cbind(eaten,alive) ~ Density,  data = JLNprop, 
           family = binomial(link = logit))
int2 <- glm(cbind(eaten,alive) ~ 1,  data = JLNprop, 
           family = binomial(link = logit))
cara2 <-  glm(cbind(eaten,alive) ~ Density + CW + Sediment,  data = JLNprop, 
             family = binomial(link = logit))
AIC(JLN_full, den2_che, den2_sca_che, den2_sed, den2, int2, cara2)

###Best model includes Density and Cheliped, closely followed by Density, Carapace width, and Sediment
par(mfrow = c(2,2))
plot(den2_che)



#Functional response (all)---------------------
clams2 <- read.csv("crab_consumption_raw.csv",na.strings=c(''), stringsAsFactors = TRUE) %>%
  mutate(period = ifelse(Rep < 7, "first", "second"))

# VC functional response ---------------------------------
firstperiod <- clams2 %>% filter(period == "first")
varnishclams <- filter(firstperiod, Clam == "VC")
green1 <- filter(varnishclams, Sediment == "1")
green2 <- filter(varnishclams, Sediment == "0")

# Test for type II including zeros
frair_test(eaten ~ Density, data = green1)
frair_test(eaten ~ Density, data = green2) #Type II, no sand


### Frair fit
outII_g3 <- frair_fit(eaten ~ Density, data = green2, response = 'rogersII',
                      start = list(a = 0.2, h = 0.2), fixed = list(T=1))
#this forces the curve into reading as a Type II curve. As the first terms would be negative if the curve is a TYpe II and positive if a Type III
#frair test tells us what inflection of the curve is
#the frair fit fits a curve to your data using the Ne equation
#now that has this curve can pull data from this information including the coefficients (ie handling time)

# A linear fit
outI_g3 <- frair_fit(eaten ~ Density, data = green2, response = 'typeI',
                     start = list(a = 0.2), fixed=list(T=1))
#this is to confirm that data is not a Type I

# Visualise fits
plot(outII_g3, pch=20, col=rgb(0,0,0,0.2), xlim=c(0,16))
lines(outII_g3)
lines(outI_g3, lty=3)

#So now we know that Type II is better than Type III using frair and visually confirmed
#That it is not a linear using the best fit graph

#####Greens resid
a <- outII_g3$coefficients[1] # Get coeffs
h <- outII_g3$coefficients[2]
#this gives us our attack rate and handling time (ie 2.89 and 0.12)

#this creates a data frame of our predictions to compare to what we actually got
fitsgrn3 <- data.frame(x = green2$Density) # Calculate 'a' for each value of x, where x is density of varnish clams

fitsgrn3$Ne <- fitsgrn3$x - lambertW0(a * h * fitsgrn3$x * exp(-a * (1 - h * fitsgrn3$x)))/(a * h) # calculate expected number of varnish clams eaten, based on frair flexnr equation, using lambert function

fitsgrn3$actual <- green2$eaten
fitsgrn3$resid <- fitsgrn3$Ne - fitsgrn3$actual

plot(x = fitsgrn3$x, y = fitsgrn3$Ne)
plot(x = fitsgrn3$Ne, y = fitsgrn3$resid) #checking residual fits
abline(h = 0, lty = 'dotted')

# Have a look at original fits returned by mle2 (*highly* recommended)
summary(outII_g3$fit) #checking to see if attack and handling time is signficant or not
summary(outI_g3$fit)

# Compare models using AIC
AIC(outI_g3$fit,outII_g3$fit) #type II is for sure better
#this is just to reconfirm that this is a Type II curve, whcih is confirmed with data, yay

##### nlm to look at asymptote (analysis not in manuscript)
green_asymp <- nls(eaten ~ SSasymp(density, Asym, R0, lrc), 
                   data = green) # asymptotic curve with initial guesses
summary(green_asymp)
#this asymptote tells us the max 
#attack rate and handling time attributes to where our asymptote lays
#attack rate and handling time for the two clam species may differ which could be influenced by biology and physical characteristiics
#asymptote at max density what is the max they can eat, can compare asympottoes between the eating of the two species
#so far this looks at an individual treatment combination
#next section of code is comparing them 

# Bootstrap
set.seed(309331)
outII_g3_boot <- frair_boot(outII_g3, start = NULL, strata=varnish2[,6], nboot=2000,
                           para=TRUE, ncores=NaN, WARN.ONLY=FALSE)

outII_g3_boot #a = 2.89, h = 0.12, T = 1
confint(outII_g3_boot)

# Illustrate bootlines
plot(outII_g3_boot, xlim=c(0,16), ylim = c(0, 16), type='n', main='All bootstrapped lines')
lines(outII_g3_boot, all_lines=TRUE)
points(outII_g3_boot, pch=20)

# Illustrate bootpolys
plot(outII_g3_boot, xlim=c(0,16), ylim = c(0, 16), type='n', main='Empirical 95 percent CI')
drawpoly(outII_g3_boot, col=hsv(2/6,0.2, 0.8))
points(outII_g3_boot, pch=20)
lines(outII_g3_boot, all_lines=FALSE)

par(bg = 'white', fg = 'black')
plot(outII_g3_boot, xlim=c(0, 16), ylim = c(0, 16), type='n',
     xlab = "Initial Clam Density",
     ylab="Clams Consumed", 
     cex.lab = 1.5,
     font.lab = 2,
     cex.axis = 1.2,
     cex.main = 1.5)
lines(outII_g3_boot, lwd = 3, all_lines=FALSE, col= "#625a94", lty = 2)
drawpoly(outII_g3_boot, border = NA, col=adjustcolor("#625a94", alpha.f = 0.4))
points(outII_g3_boot, pch=17, col=adjustcolor("#625a94", alpha.f = 0.4), cex = 1.4)

# BOOTSTRAP NeS ----------------------------------------------------------------
#compare Ne at max experimental density
# Set clam density to the max prey density offered
clam_den <- c(16)

# get the bootstrapped coefficients
a_m <- outII_g3_boot$bootcoefs[, 1]
h_m <- outII_g3_boot$bootcoefs[, 2]
T_m <- outII_g3_boot$bootcoefs[, 3]

# use the coefficients to calculate expected number of clams eaten, based on 
#frair flexnr equation, using the rogersII function
male_boot <- rogersII(16, a_m, h_m, T_m)


#VC linear models comparison between July and August
green <- clams %>%
  filter(Clam == "VC" & Sediment == "1")

#Green crabs
#cheliped are signif diff
che <- lm(Cheliped ~ period, data = green)
summary(che)
che_summary <- t.test(Cheliped ~ period, data = green)
che_summary

ggplot(green) + geom_boxplot(aes(x = period, y = Cheliped))

# CW aren't signif different
carapacewidth <- lm(CW ~ period, data = green)
summary(carapacewidth)
che_summary <- t.test(CW ~ period, data = green)
che_summary
ggplot(green) + geom_boxplot(aes(x = period, y = CW))

# does location affect eaten? NO
ggplot(green) + geom_jitter(aes(x = Density, y = proportion_eaten, colour = period))
year3 <- lm(proportion_eaten ~ Density + period, data = green)
summary(year3)

simulationOutput_model <- simulateResiduals(fittedModel = year3, plot = F)
plot(simulationOutput_model)
testDispersion(simulationOutput_model)
testOutliers(simulationOutput_model)
testZeroInflation(simulationOutput_model)

#Does temp differ? YES
green2 <- clams %>%
  filter(Clam == "VC" & Sediment == "1")
greenVC <- green2 %>% filter(period == "first")
greenVC2 <- green2 %>% filter(period == "second")
summary(greenVC)
summary(greenVC2)

tempdif <- lm(mean_temp ~ period, data = green2)
temp_summary <- t.test(mean_temp ~ period, data = green)
temp_summary
summary(tempdif)
summary(green2)
ggplot(green2) + geom_boxplot(aes(x = period, y = mean_temp))
tempresult1 <- t.test(mean_temp~period, data = green2)
tempresult1

#prop eaten
prop_summary <- t.test(proportion_eaten ~ period, data = green)
prop_summary

#temp and prop
temp_prop <- glm(proportion_eaten ~ mean_temp*period, data = green)
summary(temp_prop)
simulationOutput_tp <- simulateResiduals(fittedModel = temp_prop, plot = F)
plot(simulationOutput_tp)
testDispersion(simulationOutput_model)
testOutliers(simulationOutput_model)
testZeroInflation(simulationOutput_model)

#temp and proportion eaten? NO
#Starting with yes/no if the crabs ate
green2 <- clams %>%
  filter(Clam == "VC" & Sediment == "1")
crab_temp_logistic <- green2 %>% 
  #create a new variables of zeros and ones for ate, where 0 means they didn't eat and 1 means they did. 
  mutate(eaten_logistic2 = case_when(eaten == 0 ~ 0,
                                    TRUE ~ 1))
crab_temp_gamma <- crab_temp_logistic %>% filter(eaten != 0)

hurdle_tempprop <- glmmTMB(eaten_logistic2 ~ Cheliped + Density + mean_temp, 
                family = binomial(link = "logit"),
                data = crab_temp_logistic)
summary(hurdle_tempprop) 
tempprop_resid <- simulateResiduals(hurdle_tempprop)
plot(tempprop_resid) #does not works
testDispersion(tempprop_resid)
testOutliers(tempprop_resid)
testZeroInflation(tempprop_resid)

visreg(hurdle_tempprop)

hurdle_tempgamma <- glmmTMB(eaten_logistic2 ~ mean_temp, 
                           family = Gamma(link = "log"),
                           data = crab_temp_gamma)
summary(hurdle_tempgamma) #mean_temp not signficant
tempprop_resid <- simulateResiduals(hurdle_tempgamma)
plot(tempprop_resid) #does not works
visreg(tempprop)

ggplot(data = crab_temp_logistic, aes(x = mean_temp, y = eaten_logistic)) +
  geom_point() +
  stat_smooth(method = 'lm') + 
  theme_classic()



#ggpredict version with standarization------------------
green2 <- clams %>%
  filter(Clam == "VC" & Sediment == "1")
crab_temp_logistic <- green2 %>% 
  #create a new variables of zeros and ones for ate, where 0 means they didn't eat and 1 means they did. 
  mutate(eaten_logistic = case_when(eaten == 0 ~ 0,
                                    TRUE ~ 1))
crab_temp_log <- crab_temp_logistic %>%
  mutate(scaled_temp = scale(mean_temp)[,1])
tempprop <- glm(eaten_logistic ~ mean_temp, 
                family = binomial,
                data = crab_temp_log)
summary(tempprop)


predict_nb <- ggpredict(tempprop, terms = "mean_temp[n=100]") %>% 
  #note that specifying n= in square brackets just
  #tells ggpredict that you want it to create a 
  #larger df, which can help smooth out the plot
  #you can also specify a range here instead
  #(e.g., [10,30])
  #now we'll rename some of these columns so they line up with our original df
  rename(mean_temp = x)

ggplot(crab_temp_log, aes(x = mean_temp, y = eaten_logistic)) +
  #plot the raw data
  geom_point() +
  #now plot the  model output
  geom_ribbon(data = predict_nb,
              aes(y = predicted, ymin = conf.low, ymax = conf.high), 
              alpha = 0.3) +
  geom_line(data = predict_nb,
            aes(y = predicted)) +
  theme_classic() 


#The only other trick to know here is that if you've standardized your response variable
#to put it into the model (which you should!), 
#you'll also have to backtransform the predictor column before you plot anything.

mod_nb_stand <- glmmTMB(eaten_logistic ~ scaled_temp,
                        family = nbinom2(link = 'log'),
                        data = crab_temp_log)

predict_nb_stand <- ggpredict(mod_nb_stand, 
                              terms = "scaled_temp[n=100]") %>% 
  #note that specifying n= in square brackets just
  #tells ggpredict that you want it to create a 
  #larger df, which can help smooth out the plot
  #you can also specify a range here instead
  #(e.g., [10,30])
  #now we'll rename some of these columns so they line up with our original df
  mutate(mean_temp = 
           x*sd(crab_temp_log$mean_temp)+mean(crab_temp_log$mean_temp))

ggplot(crab_temp_log, aes(x = mean_temp, y = eaten_logistic)) +
  #plot the raw data
  geom_point() +
  #now plot the  model output
  geom_ribbon(data = predict_nb_stand,
              aes(y = predicted, ymin = conf.low, ymax = conf.high), 
              alpha = 0.3) +
  geom_line(data = predict_nb_stand,
            aes(y = predicted)) +
  theme_classic() 

#proportion eaten and period?
periodprop <- lm(proportion_eaten ~ period, data = green)
summary(periodprop)
ggplot(green) + geom_boxplot(aes(x = period, y = proportion_eaten))

#cheli x cara graph
ggplot(green) + geom_point(aes(x = CW, y = Cheliped, colour = period)) + 
  geom_smooth(aes(x = CW, y = Cheliped, colour = period), method = lm)

#Wilcox test
#difference in the proportion of clams eaten by period? Yes
green.period <- green[green$Density == 16,]
green.period$period <- as.factor(green.period$period)
green.period$proportion_eatenlog <- log(green.period$proportion_eaten)
plot(x = green.period$period, y = green.period$proportion_eaten)

green.period.lm <- lm(proportion_eaten ~ period, data = green.period)
summary(green.period.lm)
plot(green.period.lm, which = c(1,2,4))

green.period.nls <- glm(proportion_eaten ~ period, data = green.period[green.period$proportion_eaten > 0,], family = Gamma)
summary(green.period.nls)
plot(green.period.nls)

green.period$species <- 'C. maenas'

wilcox.test(formula = green.period$proportion_eaten ~ green.period$period)


#Functional response (nonzero)------------------------------
#no zero data
firstperiod <- clams2 %>% filter(period == "first")
nozero <- firstperiod %>% filter(eaten > 0)
varnishclams <- filter(nozero, Clam == "VC")

varnish1 <- filter(varnishclams, Sediment == "1")
varnish2 <- filter(varnishclams, Sediment == "0")
ggplot(varnishclams, aes(x=Density, y=proportion_eaten)) + geom_point()
frair_test(eaten ~ Density, data = varnish1) #no FR; sand
frair_test(eaten ~ Density, data = varnish2) #Type II; no sand
### Frair fit
#outII_g <- frair_fit(eaten ~ Density, data = varnish1, response = 'rogersII', start = list(a = 0.2, h = 0.2), fixed = list(T=1))

outII_r <- frair_fit(eaten ~ Density, data = varnish2, response = 'rogersII',
                     start = list(a = 0.2, h = 0.2), fixed = list(T=1))

# A linear fit
#outI_g <- frair_fit(eaten ~ Density, data = varnish1, response = 'typeI', start = list(a = 0.2), fixed=list(T=1))

outI_r <- frair_fit(eaten ~ Density, data = varnish2, response = 'typeI',
                    start = list(a = 0.2), fixed=list(T=1))
# Visualise fits
plot(outII_r, pch=20, col=rgb(0,0,0,0.2), xlim=c(0,16), ylim=c(0,10), 
     xlab = "Prey Density",ylab = "Number of VC eaten", las = 1)
lines(outII_r) #no sand, plain black line
#lines(outII_g, lty=3) #sediment, dotted line
#lines(outI_g, lty=4)
#lines(outI_r, lty=5)

#####Greens resid for sediment-----------------
#a <- outII_g$coefficients[1] # 5.06753
#h <- outII_g$coefficients[2] # 0.4714

#fitsgrn <- data.frame(x = green2$Density) # Calculate 'a' for each value of x, where x is density of VC clams

#fitsgrn$Ne <- fitsgrn$x - lambertW0(a * h * fitsgrn$x * exp(-a * (1 - h * fitsgrn$x)))/(a * h) # calculate expected number of VC eaten, based on frair flexnr equation, using lambert function

#fitsgrn$actual <- green2$eaten
#fitsgrn$resid <- fitsgrn$Ne - fitsgrn$actual

#plot(x = fitsgrn$x, y = fitsgrn$Ne)
plot(x = fitsgrn$Ne, y = fitsgrn$resid)
abline(h = 0, lty = 'dotted')

#####Greens resid for no sediment--------------------
a2 <- outII_r$coefficients[1] # 2.890235
h2 <- outII_r$coefficients[2] # 0.1197529

fitsgrn <- data.frame(x = varnish2$Density) # Calculate 'a' for each value of x, where x is density of VC clams

fitsgrn$Ne <- fitsgrn$x - lambertW0(a2 * h2 * fitsgrn$x * exp(-a2 * (1 - h2 * fitsgrn$x)))/(a2 * h2) # calculate expected number of VC eaten, based on frair flexnr equation, using lambert function

fitsgrn$actual <- varnish2$eaten
fitsgrn$resid <- fitsgrn$Ne - fitsgrn$actual

plot(x = fitsgrn$x, y = fitsgrn$Ne)
plot(x = fitsgrn$Ne, y = fitsgrn$resid)
abline(h = 0, lty = 'dotted')

# Have a look at original fits returned by mle2 (*highly* recommended)
summary(outII_r$fit) #checking to see if attack and handling time is signficant or not
summary(outI_r$fit)
# Bootstrap
set.seed(309331)
outII_r_boot <- frair_boot(outII_r, start = NULL, strata=varnish2[,6], nboot=2000,
                              para=TRUE, ncores=NaN, WARN.ONLY=FALSE)

outII_r_boot #a = 2.89, h = 0.12, T = 1
confint(outII_r_boot)

# Illustrate bootlines
plot(outII_r_boot, xlim=c(0,16), ylim = c(0, 16), type='n', main='All bootstrapped lines')
lines(outII_r_boot, all_lines=TRUE)
points(outII_r_boot, pch=20)

# Illustrate bootpolys
plot(outII_r_boot, xlim=c(0,16), ylim = c(0, 16), type='n', main='Empirical 95 percent CI')
drawpoly(outII_r_boot, col=hsv(2/6,0.2, 0.8))
points(outII_r_boot, pch=20)
lines(outII_r_boot, all_lines=FALSE)

par(bg = 'white', fg = 'black')
plot(outII_r_boot, xlim=c(0, 16), ylim = c(0, 16), type='n',
     xlab = "Initial Clam Density",
     ylab="Clams Consumed", 
     cex.lab = 1.5,
     font.lab = 2,
     cex.axis = 1.2,
     cex.main = 1.5)
lines(outII_r_boot, lwd = 3, all_lines=FALSE, col= "#625a94", lty = 2)
drawpoly(outII_r_boot, border = NA, col=adjustcolor("#625a94", alpha.f = 0.4))
points(outII_r_boot, pch=17, col=adjustcolor("#625a94", alpha.f = 0.4), cex = 1.4)

#####VC comparing first experiment to second experiment--------------------
secondperiod <- clams2 %>% filter(period == "second")
nozerosecond <- secondperiod %>% filter(eaten > 0)

#no zeros
ggplot(nozerosecond, aes(x=Density, y=proportion_eaten)) + geom_point()
frair_test(eaten ~ Density, data = nozerosecond) #Type II; sand in second period 

#taking zeros into account
frair_test(eaten ~ Density, data = secondperiod) #No FR

#fits for no zero data (all cant invert hessian?)
outII_c <- frair_fit(eaten ~ Density, data = nozerosecond, response = 'rogersII',
                     start = list(a = 0.2, h = 0.2), fixed = list(T=1))
outII_c_3 <- frair_fit(eaten ~ Density, data = nozerosecond, response = 'flexpnr',
                      start = list(b=0.0029242, h = 0.0312941), fixed = list(T=1, q = 1.9391532))
outII_c_2 <- frair_fit(eaten ~ Density, data = nozerosecond, response = 'hassIIInr',
                        start = list(b=0.0029242, c = 1.9391532, h = 0.0312941), fixed = list(T=1))


#FR for JLN-------------------------------
#no zeros
firstperiod <- clams2 %>% filter(period == "first")
nozero <- firstperiod %>% filter(eaten > 0)
littleclams <- filter(nozero, Clam == "JLN")

little1 <- filter(littleclams, Sediment == "1")
little2 <- filter(littleclams, Sediment == "0")
ggplot(littleclams, aes(x=Density, y=proportion_eaten)) + geom_point()
frair_test(eaten ~ Density, data = little1) #Type II; sand
frair_test(eaten ~ Density, data = little2) #No FR; no sand

#including zeros
little <- filter(firstperiod, Clam == "JLN")
little3 <- filter(little, Sediment == "1")
little4 <- filter(little, Sediment == "0")

frair_test(eaten ~ Density, data = little3) #No FR
frair_test(eaten ~ Density, data = little4) #No FR

### Frair fit
outII_v <- frair_fit(eaten ~ Density, data = little4, response = 'rogersII',
                     start = list(a = 0.2, h = 0.2), fixed = list(T=1))

# A linear fit
outI_v <- frair_fit(eaten ~ Density, data = little4, response = 'typeI',
                    start = list(a = 0.2), fixed=list(T=1))


# Visualise fits
plot(outII_v, pch=20, col=rgb(0,0,0,0.2), xlim=c(0,16), ylim=c(0,10), 
     xlab = "Prey Density",ylab = "Number of JLN eaten", las = 1)
lines(outII_v) #sediment, plain black line


#coeff
a3 <- outII_v$coefficients[1] # 3.444869
h3 <- outII_v$coefficients[2] # 0.5347771

fitsgrn2 <- data.frame(x = little1$Density) # Calculate 'a' for each value of x, where x is density of VC clams

fitsgrn2$Ne <- fitsgrn2$x - lambertW0(a3 * h3 * fitsgrn2$x * exp(-a3 * (1 - h3 * fitsgrn2$x)))/(a3 * h3) # calculate expected number of VC eaten, based on frair flexnr equation, using lambert function

fitsgrn2$actual <- little1$eaten
fitsgrn2$resid <- fitsgrn2$Ne - fitsgrn2$actual

plot(x = fitsgrn2$x, y = fitsgrn2$Ne)
plot(x = fitsgrn2$Ne, y = fitsgrn2$resid)
abline(h = 0, lty = 'dotted')

# Have a look at original fits returned by mle2 (*highly* recommended)
summary(outII_v$fit) #attack rate is not signficant
summary(outI_v$fit) #attack rate is signficant
# Bootstrap
set.seed(309331)
outII_v_boot <- frair_boot(outII_v, start = NULL, strata=little1[,6], nboot=2000,
                           para=TRUE, ncores=NaN, WARN.ONLY=FALSE)

outII_v_boot #a = 3.445, h = 0.535, T = 1
confint(outII_v_boot)

# Illustrate bootlines
plot(outII_v_boot, xlim=c(0,16), ylim = c(0, 5), type='n', main='All bootstrapped lines')
lines(outII_v_boot, all_lines=TRUE)
points(outII_v_boot, pch=20)

# Illustrate bootpolys
plot(outII_v_boot, xlim=c(0,16), ylim = c(0, 5), type='n', main='Empirical 95 percent CI')
drawpoly(outII_v_boot, col=hsv(2/6,0.2, 0.8))
points(outII_v_boot, pch=20)
lines(outII_v_boot, all_lines=FALSE)

par(bg = 'white', fg = 'black')
plot(outII_v_boot, xlim=c(0, 16), ylim = c(0, 16), type='n',
     xlab = "Initial Clam Density",
     ylab="Clams Consumed", 
     cex.lab = 1.5,
     font.lab = 2,
     cex.axis = 1.2,
     cex.main = 1.5)
lines(outII_v_boot, lwd = 3, all_lines=FALSE, col= "#625a94", lty = 2)
drawpoly(outII_v_boot, border = NA, col=adjustcolor("#625a94", alpha.f = 0.4))
points(outII_v_boot, pch=17, col=adjustcolor("#625a94", alpha.f = 0.4), cex = 1.4)


#Temperature effect?------------------------
#temp <- lm(eaten ~ mean_temp, data = clams)
#summary(temp) #no effect
#temp2 <- lm(proportion_eaten ~ mean_temp, data = clams)
#summary(temp) #0.0492 so mabye?
#sal <- lm(eaten ~ mean_salinity, data = clams)
#summary(sal) #no effect
#sal <- lm(proportion_eaten ~ mean_salinity, data = clams)
#summary(sal) #no effect 

#VC
#VC_full <- glm(cbind(eaten,alive) ~ Sediment + Density + scaled_cheliped + mean_temp, data = VCprop, family = binomial(link = logit))

VC_glm <- glm(proportion_eaten ~ scaled_density + scaled_temp + scaled_cheliped + Sediment, 
                         weights = Density, 
                         data = VCprop, 
                         family = binomial(link = "logit"))

#vif(VC_full)
vif(VC_glm)

 #moderate correlation between predictor variables
#summary(VC_full)

summary(VC_glm) #density, scaled_cheliped, and sediment all have significant effect


#DHARMa
#Sediment
simulationOutput2 <- simulateResiduals(fittedModel = VC_glm, plot = F)
residuals(simulationOutput2)
plotResiduals(simulationOutput2, form = VCprop$Sediment)

plot(simulationOutput2)
plotResiduals(simulationOutput2, VCprop$Sediment)
plotResiduals(simulationOutput2, VCprop$Density)
plotResiduals(simulationOutput2, VCprop$Cheliped)
plotResiduals(simulationOutput2, VCprop$mean_temp)



testDispersion(simulationOutput2)
testOutliers(simulationOutput2)
testZeroInflation(simulationOutput2)


#VC_glm model outputs
VC_glm_outputs <- simulateResiduals(fittedModel = VC_glm, plot = F)
plot(VC_glm_outputs)

plotResiduals(simulationOutput$scaledResiduals, VCprop$Sediment)
testDispersion(VC_glm_outputs)
testOutliers(VC_glm_outputs)
testZeroInflation(VC_glm_outputs)

#trying instead with the long form data
VC_glmer <- glmer(eaten ~ scaled_density + scaled_temp + scaled_cheliped + 
                  Sediment + (1|crab_num), 
              data = clamview, 
              family = binomial(link = "logit"))
vif(VC_glmer)
summary(VC_glmer)
simulationOutput_glmer <- simulateResiduals(fittedModel = VC_glmer, plot = F)
plot(simulationOutput_glmer)
plotResiduals(simulationOutput_glmer, clamview$Sediment)
plotResiduals(simulationOutput_glmer, clamview$scaled_density)
plotResiduals(simulationOutput_glmer, clamview$scaled_cheliped)
plotResiduals(simulationOutput_glmer, clamview$scaled_temp)
testZeroInflation(simulationOutput_glmer)


VC_glmer2 <- glmer(proportion_eaten ~ scaled_density + scaled_temp + scaled_cheliped + 
                    Sediment + (1|crab_num), 
                  data = clamview, 
                  family = binomial(link = "logit"))

###



###

#Predation and crab claw size to see if it is influencing it--------------------
firstperiod <- clams2 %>% filter(period == "first")
green <- filter(clams2, Clam == "VC")
green2 <- filter(green, Sediment == "1")
green3 <- filter(green, Sediment == "0")

plot(green2$proportion_eaten ~ green2$Cheliped, xlab = "Cheliped height (mm)", ylab = "Proportion VC eaten in sand")
abline(lm(proportion_eaten ~ Cheliped, data = green2))
regression1 <- lm(proportion_eaten ~ Cheliped, data = green2)
summary(regression1) #significant

plot(green3$proportion_eaten ~ green3$Cheliped, xlab = "Cheliped height (mm)", ylab = "Proportion VC eaten with no sand")
regression2 <- lm(proportion_eaten ~ Cheliped, data = green3)
summary(regression2) #not significant

jclam <- filter(clams2, Clam == "JLN")
jclam2 <- filter(green, Sediment == "1")
jclam3 <- filter(green, Sediment == "0")

plot(jclam2$proportion_eaten ~ jclam2$Cheliped, xlab = "Cheliped height (mm)", ylab = "Proportion JLN eaten in sand")
abline(lm(proportion_eaten ~ Cheliped, data = jclam2))
regression3 <- lm(proportion_eaten ~ Cheliped, data = jclam2)
summary(regression3) #significant

plot(jclam3$proportion_eaten ~ jclam3$Cheliped, xlab = "Cheliped height (mm)", ylab = "Proportion JLN eaten with no sand")
regression4 <- lm(proportion_eaten ~ Cheliped, data = jclam3)
summary(regression4) #not significant

#number eaten VC
plot(green2$eaten ~ green2$Cheliped, xlab = "Cheliped height (mm)", ylab = "Number of VC eaten in sand")
abline(lm(eaten ~ Cheliped, data = green2))
regression5 <- lm(eaten ~ Cheliped, data = green2)
summary(regression5)

plot(green3$eaten ~ green3$Cheliped, xlab = "Cheliped height (mm)", ylab = "Number of VC eaten with no sand")
regression2 <- lm(eaten ~ Cheliped, data = green3)
summary(regression2) #not significant

#number eaten JLN
plot(jclam2$eaten ~ jclam2$Cheliped, xlab = "Cheliped height (mm)", ylab = "Number of JLN eaten in sand")
abline(lm(eaten ~ Cheliped, data = jclam2))
regression3 <- lm(eaten ~ Cheliped, data = jclam2)
summary(regression3)

plot(jclam3$eaten ~ jclam3$Cheliped, xlab = "Cheliped height (mm)", ylab = "Number of JLN eaten with no sand")
regression4 <- lm(eaten ~ Cheliped, data = jclam3)
summary(regression4) #not significant



# Number eaten per cm cheliped --------------------------------------------

cheliped <- read.csv('functional_response_cheliped.csv', header = T)
greenchel <- cheliped[cheliped$species == 'green',]
redchel <- cheliped[cheliped$species == 'red',]

par(mfrow = c(1, 1))
plot(cheliped$density_per_cm_round[cheliped$species == 'green'], cheliped$eaten_per_cm[cheliped$species == 'green'], pch = 19, col = 'forestgreen',
     xlab = 'Initial prey density', ylab = 'Number of prey eaten per cm cheliped height')
points(cheliped$density_per_cm_round[cheliped$species == 'red'], cheliped$eaten_per_cm[cheliped$species == 'red'], pch = 17, col = 'firebrick')
legend('topleft', col = c('forestgreen', 'firebrick'), pch = c(19, 17), legend = c('Carcinus maenas', 'Cancer productus'))



#Hurdle model------------------
#Starting with yes/no if the crabs ate
firstperiod <- clams %>% filter(period == "first")
crab_consump_logistic <- firstperiod %>% 
  #create a new variables of zeros and ones for ate, where 0 means they didn't eat 
  #and 1 means they did. 
  mutate(eaten_logistic = case_when(eaten == 0 ~ 0,
                                       TRUE ~ 1))
#this gets rid of the crabs that didnt eat so there will be no observations 
#that have a value of 0 for consumption if you look at the data set after this line
crab_consump_gamma <- firstperiod %>% 
  #get rid of all the zeros for the analysis of the non-zero values
  filter(eaten != 0)

#this hurdle model is a generalized linear model that is first looking at the 
#crab_consump_logistic that has crabs that ate and ones that didnt, 
#and it is seeing if variables influenced if the crabs ate or not. 
#Family binomial means that its categorical (ate? yes or no)
#crab_consump_logistic is eat or not eat
hurdle_mod1 <- glm(eaten_logistic ~ Sediment + Cheliped + 
                       Clam + Density + mean_temp, 
                     family = binomial,
                     data = crab_consump_logistic)
hurdle1_resid <- simulateResiduals(hurdle_mod1)
plot(hurdle1_resid) #YES
car::vif(hurdle_mod1)

# this is looking at if those that did eat were influenced by variables.
#Family gamma means its a continuous non zero scale for consumption 
#that we are looking at 
#of the ones that ate, how much did they eat
#a hurdle gamma model allows for zeroes, 
#but will use a different distribution to model the 0-1 portion of your data 
#(i.e., it will categorize your data in zero and not-zero) 
#and then a Gamma distribution on the continuous data (i.e., the not-zero side). 
#It will then give you coefficient estimates for both parts of the model 
#(i.e., you are essentially running two separate models).
#gamma family might not be the right family, double check
hurdle_mod2 <- glm(eaten ~ Sediment + Cheliped + 
                       Clam + Density + mean_temp,
                     family = Gamma(link = "log"),
                     data = crab_consump_gamma)
car::vif(hurdle_mod2)

hurdle2_resid <- simulateResiduals(hurdle_mod2)
plot(hurdle2_resid) #nope


summary(hurdle_mod1) #this will tell you which predictors are related to 
#whether a crab eats = ALL

plot_model(hurdle_mod1)

summary(hurdle_mod2) #this will tell you which predictors are related to 
#how much a crab eats = all but cheliped (temp just barely)

plot_model(hurdle_mod2)

#hurdle models with glmmTMB-------------
#Starting with yes/no if the crabs ate
firstperiod <- clams2 %>% filter(period == "first")
crab_consump_logistic <- firstperiod %>% 
  mutate(eaten_logistic = case_when(eaten == 0 ~ 0,
                                    TRUE ~ 1))

#this gets rid of the crabs that didnt eat so there will be no observations 
#that have a value of 0 for consumption if you look at the data set after this line
#Family gamma means its a continuous non zero scale for consumption 
#that we are looking at of the ones that ate, how much did they eat
crab_consump_gamma <- crab_consump_logistic %>% 
  #get rid of all the zeros for the analysis of the non-zero values
  filter(eaten != 0)

#this hurdle model is a generalized linear model that is first looking at the 
#crab_consump_logistic that has crabs that ate and ones that didnt, 
#and it is seeing if variables influenced if the crabs ate or not. 

hurdle_mod3 <- glmmTMB(eaten_logistic ~ Cheliped + Density + Clam + Sediment, 
                   family = binomial(link = "logit"),
                   data = crab_consump_logistic)
summary(hurdle_mod3) #all signficant
hurdle3_resid <- simulateResiduals(hurdle_mod3)
plot(hurdle3_resid)

hurdle_mod4 <- glmmTMB(eaten ~ Sediment + Cheliped + 
                     Clam + Density + mean_temp,
                   family = Gamma(link = "log"),
                   data = crab_consump_gamma)
summary(hurdle_mod4) #Cheliped not signficant
hurdle4_resid <- simulateResiduals(hurdle_mod4)
plot(hurdle4_resid)
testDispersion(hurdle4_resid)
testOutliers(hurdle4_resid)
testZeroInflation(hurdle4_resid)

#predict for logistic model
hurdle3_predict <- ggpredict(hurdle_mod3, 
                                terms = c("Density[n = 100]", "Sediment", "Clam")) %>% 
  rename(Density = x,
         Sediment = group,
         Clam = facet) 
#plot predictions on raw data
colclam <- c("#C2EAFF","#004266","#ADEBD2","#186345")
clamlabel <- c("JLN, no sed","JLN, with sed", "VC, no sed", "VC, with sed")

ggplot(hurdle3_predict, aes(Density, predicted)) + 
  #geom_point(size = 2, position = position_dodge(width = 0.8), aes(color = Clam:Sediment)) +
  #geom_errorbar(aes(ymax = conf.high, ymin = conf.low, color = Clam:Sediment), 
                #width = 0.1, position = position_dodge(width = 0.8)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Clam:Sediment), alpha = 0.3) +
  geom_line(aes(color = Clam:Sediment)) +
  geom_point(data = crab_consump_logistic, 
             aes(Density, eaten_logistic, color = Clam:Sediment),
             size = 0.5, position = position_jitter(height = 0.01)) + 
  scale_colour_manual(values = colclam, labels = clamlabel) +
  scale_fill_manual(values = colclam, labels = clamlabel) +
  labs(x = "Density",
       y = "Probability of consuming a clam") +
  theme(plot.title = element_text(hjust = -0.1, vjust = 2.5, size = 12))

#predict for gamma model
hurdle4_predict <- ggpredict(hurdle_mod4, 
                             terms = c("Density[n = 100]", "Sediment", "Clam")) %>% 
  rename(Density = x,
         Sediment = group,
         Clam = facet) 
#plot predictions on raw data

ggplot(hurdle4_predict, aes(Density, predicted)) + 
  #geom_point(size = 2, position = position_dodge(width = 0.8), aes(color = Clam:Sediment)) +
  #geom_errorbar(aes(ymax = conf.high, ymin = conf.low, color = Clam:Sediment), 
  #width = 0.1, position = position_dodge(width = 0.8)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Clam:Sediment), alpha = 0.3) +
  geom_line(aes(color = Clam:Sediment)) +
  geom_point(data = crab_consump_gamma, 
             aes(Density, eaten, color = Clam:Sediment),
             size = 0.5, position = position_jitter(height = 0.01)) + 
  scale_colour_manual(values = colclam, labels = clamlabel) +
  scale_fill_manual(values = colclam, labels = clamlabel) +
  labs(x = "Density",
       y = "Number of clams consumed") +
  theme(plot.title = element_text(hjust = -0.1, vjust = 2.5, size = 12))

#temp hurdle model----------------
#Starting with yes/no if the crabs ate
green2 <- clams %>%
  filter(Clam == "VC" & Sediment == "1")
crab_temp_logistic <- green2 %>% 
  #create a new variables of zeros and ones for ate, where 0 means they didn't eat and 1 means they did. 
  mutate(eaten_logistic2 = case_when(eaten == 0 ~ 0,
                                     TRUE ~ 1))
crab_temp_gamma <- crab_temp_logistic %>% filter(eaten != 0)

hurdle_tempprop <- glmmTMB(eaten_logistic2 ~ Cheliped + Density + mean_temp + period, 
                           family = binomial(link = "logit"),
                           data = crab_temp_logistic)
summary(hurdle_tempprop) 

tempprop_resid <- simulateResiduals(hurdle_tempprop)
plot(tempprop_resid) #does not works
testDispersion(tempprop_resid)
testOutliers(tempprop_resid)
testZeroInflation(tempprop_resid)

hurdle_gammatemp <- glmmTMB(eaten ~ Cheliped + Density + mean_temp + period,
                       family = Gamma(link = "log"),
                       data = crab_temp_gamma )

summary(hurdle_gammatemp) 

tempprop_resid2 <- simulateResiduals(hurdle_gammatemp)
plot(tempprop_resid2) #does not works
testDispersion(tempprop_resid2)
testOutliers(tempprop_resid2)
testZeroInflation(tempprop_resid2)

#predict for logistic model
hurdle_temp_predict <- ggpredict(hurdle_tempprop, 
                             terms = c("Density[n = 100]", "period")) %>% 
  rename(Density = x,
         period = group) 

ggplot(hurdle_temp_predict, aes(Density, predicted)) + 
  #geom_point(size = 2, position = position_dodge(width = 0.8), aes(color = Clam:Sediment)) +
  #geom_errorbar(aes(ymax = conf.high, ymin = conf.low, color = Clam:Sediment), 
  #width = 0.1, position = position_dodge(width = 0.8)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = period), alpha = 0.3) +
  geom_line(aes(color = period)) +
  geom_point(data = crab_temp_logistic, 
             aes(Density, eaten_logistic2, color = period),
             size = 0.5, position = position_jitter(height = 0.01)) + 
  #scale_colour_manual(values = colclam, labels = clamlabel) +
  #scale_fill_manual(values = colclam, labels = clamlabel) +
  labs(x = "Density",
       y = "Probability of consuming a clam") +
  theme(plot.title = element_text(hjust = -0.1, vjust = 2.5, size = 12))


visreg(tempprop)

#predict for gamma model
hurdle_temp_predict2 <- ggpredict(hurdle_gammatemp, 
                             terms = c("Density[n = 100]", "period")) %>% 
  rename(Density = x,
         period = group) 

#plot predictions on raw data
ggplot(hurdle_temp_predict2, aes(Density, predicted)) + 
  #geom_point(size = 2, position = position_dodge(width = 0.8), aes(color = Clam:Sediment)) +
  #geom_errorbar(aes(ymax = conf.high, ymin = conf.low, color = Clam:Sediment), 
  #width = 0.1, position = position_dodge(width = 0.8)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = period), alpha = 0.3) +
  geom_line(aes(color = period)) +
  geom_point(data = crab_temp_gamma , 
             aes(Density, eaten, color = period),
             size = 0.5, position = position_jitter(height = 0.01)) + 
  #scale_colour_manual(values = colclam, labels = clamlabel) +
  #scale_fill_manual(values = colclam, labels = clamlabel) +
  labs(x = "Density",
       y = "Number of clams consumed") +
  theme(plot.title = element_text(hjust = -0.1, vjust = 2.5, size = 12))

#Model with CH and temp-----------
firstperiod <- clams2 %>% filter(period == "first")
model6 <- lm(proportion_eaten ~  Cheliped + mean_temp, 
              data = firstperiod, 
              family = binomial(link = "logit"))

summary(model6)
modl6_resid <- simulateResiduals(model6)
plot(modl6_resid) #not good
testDispersion(modl6_resid)
testOutliers(modl6_resid)
testZeroInflation(modl6_resid)

model7 <- glmmTMB(proportion_eaten ~  Cheliped + mean_temp, 
                  data = firstperiod,
                  family = binomial(link = "logit"))
  modl7_resid <- simulateResiduals(model7)
plot(modl7_resid) #not good

model8 <- glm(proportion_eaten ~  Cheliped + mean_temp, 
             data = firstperiod)
summary(model8)
model8_resid <- simulateResiduals(model8)
plot(model8_resid) #not good
testDispersion(model8_resid)
testOutliers(model8_resid)
testZeroInflation(model8_resid)

model15 <- glmmTMB(eaten ~ mean_temp + Cheliped,
                   family = nbinom2(link = 'log'),
                   data = firstperiod)
summary(model15)
norm_resid <- simulateResiduals(model15)
plot(norm_resid) #best so far but still not good
testDispersion(norm_resid)
testOutliers(norm_resid)
testZeroInflation(norm_resid)

model16 <- glmmTMB(eaten ~ mean_temp + Cheliped,
                   family = poisson(link = 'log'),
                   data = firstperiod)
summary(model16)
mod16_resid <- simulateResiduals(model16)
plot(mod16_resid) #not good
testDispersion(mod16_resid)
testOutliers(mod16_resid)
testZeroInflation(mod16_resid)

modeltemp <- lm(eaten ~ mean_temp,
              data = firstperiod)
summary(modeltemp)
modeltemp_resid <- simulateResiduals(modeltemp)
plot(modeltemp_resid)

modelCH <- lm(eaten ~ Cheliped,
                data = firstperiod)
summary(modelCH)
modelCH_resid <- simulateResiduals(modelCH)
plot(modelCH_resid)

scaled_glm <- glm(proportion_eaten ~ scaled_cheliped, data = clamscale,
                  weights = Density, 
                  family = binomial(link = logit))
summary(scaled_glm)
simulateResiduals(scaled_glm)

#CH and temp per clam species
JLNfirst <- filter(clams, Clam == "JLN")
VCfirst <- filter(clams, Clam == "VC")

modelCH_JLN <- glm(proportion_eaten ~ Cheliped, data = JLNfirst,
                   family = binomial(link = "logit"))
summary(modelCH_JLN)
simulateResiduals(modelCH_JLN)

modelCH_lm_JLN <- lm(eaten ~ Cheliped, data = JLNfirst)
summary(modelCH_lm_JLN)
simulateResiduals(modelCH_lm_JLN)

modelCH_lm_JLN2 <- lm(proportion_eaten ~ Cheliped, data = JLNfirst)
summary(modelCH_lm_JLN2)
simulateResiduals(modelCH_lm_JLN2)

JLNscale <- filter(clamscale, Clam == "JLN")
scaled_JLN_glm <- glm(proportion_eaten ~ scaled_cheliped, data = JLNscale,
                  weights = Density, 
                  family = binomial(link = logit))
summary(scaled_JLN_glm)
simulateResiduals(scaled_JLN_glm)


hurdle_mod3 <- glmmTMB(eaten_logistic ~ Cheliped + mean_temp, 
                       family = binomial(link = "logit"),
                       data = crab_consump_logistic)
summary(hurdle_mod3) #all signficant
hurdle3_resid <- simulateResiduals(hurdle_mod3)
plot(hurdle3_resid)

#cbind version
modelcbind <- glm(cbind(eaten, alive) ~ Cheliped + Clam, family="binomial", data = firstperiod)
summary(modelcbind)
plot(simulateResiduals(modelcbind))

#% consumed figure 3------------------
firstperiod <- clams2 %>% filter(period == "first")
firstscale <- clamscale %>% filter(period == "first")
#JLN <- firstperiod %>% filter(Clam == "JLN")
#as.integer(firstperiod$proportion_eaten)

#% consumed is linked to predator so individual crab CH
#firstclam <- firstperiod %>% mutate(percentconsumed = (eaten/Density)*100)

#hist(firstperiod$proportion_eaten, breaks=16, xlab = "clamproportion", freq = FALSE)

model2 <- glmmTMB(proportion_eaten ~ Cheliped + Sediment + Clam,
                  family = binomial(link = "logit"),
                  weights = Density, 
                  data = firstperiod)
summary(model2)
plot(simulateResiduals(model2))
testDispersion(model2)
testOutliers(model2)
testZeroInflation(model2)

success <- matrix(firstperiod$eaten, ncol = 1)
failures <- matrix(firstperiod$alive, ncol = 1)
clammat <- cbind(success, failures)

model2.4 <- glmmTMB(clammat ~ Cheliped + Sediment + Clam,
              family = binomial(link = "logit"),
              data = firstperiod)
summary(model2.4)
plot(simulateResiduals(model2.4))
plot(model2.4)
testDispersion(model2.4)
testOutliers(model2.4)
testZeroInflation(model2.4)

modelscale1 <- glm(proportion_eaten ~ scaled_cheliped + Sediment + Clam,
                    family = binomial(link = "logit"),
                    weights = Density,
                    data = firstscale)
plot(simulateResiduals(modelscale1))

model4.5 <- glmmTMB(proportion_eaten ~ Cheliped + Clam*Sediment,
                  family = binomial(link = "logit"),
                  weights = Density, 
                  data = firstperiod)
summary(model4.5)
plot(simulateResiduals(model4.5))
testDispersion(model4)
testOutliers(model4.5)
testZeroInflation(model4)


model3 <- glmmTMB(proportion_eaten ~ Cheliped + Sediment + Clam,
                 family='beta_family', 
                 data = firstperiod)
summary(model3)
plot(simulateResiduals(model3))

model4 <- glm(proportion_eaten ~ Cheliped + Sediment + Clam, family = quasibinomial(link = "logit"),
              data = firstperiod)
summary(model4)
plot(simulateResiduals(model4))
plot(model4)

model4.1 <- glm(proportion_eaten ~ Cheliped + Sediment, family = binomial(link = "logit"), data = JLN)
summary(model4)
plot(simulateResiduals(model4.1))
plot(model4)

model5 <- glm(proportion_eaten ~ Cheliped + Sediment + Clam,  
              family = poisson(link = 'log'),
              weights = Density,
              data = firstperiod)
summary(model5)
plot(simulateResiduals(model5))
testDispersion(model5)
testOutliers(model5)
testZeroInflation(model5)

modelpercent <- lm(percentconsumed ~ Cheliped + Sediment + Clam, data = firstclam)
summary(modelpercent) 
plot(simulateResiduals(modelpercent))

modelpercent2 <- glm(percentconsumed ~ Cheliped + Sediment + Clam,
                     data = firstclam,
                     family = binomial(link = "logit"))
summary(modelpercent2) 
plot(simulateResiduals(modelpercent2))

firstclam$percentconsumed <- as.integer(firstclam$percentconsumed)
modelpercent3 <- glm(percentconsumed ~ Cheliped + Sediment + Clam,
                     family = poisson(link = 'log'),
                     data = firstclam)
summary(modelpercent3) 
plot(simulateResiduals(modelpercent3))

firstclamnozero <- firstclam %>% filter(percentconsumed > 0)
modelpercent4 <- lm(percentconsumed ~ Cheliped + Sediment + Clam,
                    data = firstclamnozero)
plot(simulateResiduals(modelpercent4))

#glm with interaction------------------------
model4.5 <- glmmTMB(proportion_eaten ~ Cheliped + Clam*Sediment,
                    family = binomial(link = "logit"),
                    weights = Density, 
                    data = firstperiod)
summary(model4.5)
plot(simulateResiduals(model4.5))
testDispersion(model4.5 )
testOutliers(model4.5)
testZeroInflation(model4.5)

model_predict <- ggpredict(model4.5, 
                             terms = c("Cheliped[n = 100]", "Sediment", "Clam")) %>% 
  rename(Cheliped = x,
         Sediment = group,
         Clam = facet) 

colclam <- c("#C2EAFF","#004266","#ADEBD2","#186345")
clamlabel <- c("JLN, no sed","JLN, with sed", "VC, no sed", "VC, with sed")

ggplot(firstperiod, aes(x = Cheliped, y = proportion_eaten)) + 
  geom_point(data = firstperiod, aes(x = Cheliped, y = proportion_eaten,
                                     color = Clam:Sediment), size = 0.8) + 
  geom_ribbon(data = model_predict,
              aes(y = predicted, ymin = conf.low, ymax = conf.high,
                  fill = Clam:Sediment), alpha = 0.3) + 
  geom_line(data = model_predict,
            aes(y = predicted, colour = Clam:Sediment)) +
  labs(x = "Cheliped height (mm)", y = "Proportion consumed") + 
  scale_colour_manual(values = colclam, labels = clamlabel) +
  scale_fill_manual(values = colclam, labels = clamlabel) + 
  theme(plot.title = element_text(hjust = -0.1, vjust = 2.5, size = 12))
  


#FR with no fitted curve--------------
#JLN>
#clamnoVC <- firstperiod %>% 
#JLN <- firstperiod %>% filter(Clam == "JLN")
#ggplot(JLN, aes(x = Density, y = eaten, col = Sediment)) + geom_point()
#VC
#VC <- firstperiod %>% filter(Clam == "VC")
#VC1 <- VC %>% filter(Sediment == "1")
#ggplot(VC1, aes(x = Density, y = eaten, col = Sediment)) + geom_point()

colclam <- c("#C2EAFF","#004266","#186345")
clamlabel <- c("JLN, no sed","JLN, with sed", "VC, with sed")

newdf <- firstperiod %>% 
  filter(Clam == "JLN" | Sediment == "1")

ggplot(newdf, aes(Density, eaten)) + 
  #geom_point(size = 2, position = position_dodge(width = 0.8), aes(color = Clam:Sediment)) +
  #geom_errorbar(aes(ymax = conf.high, ymin = conf.low, color = Clam:Sediment), 
  #width = 0.1, position = position_dodge(width = 0.8)) +
  geom_point(data = newdf, 
             aes(Density, eaten, color = Clam:Sediment),
             size = 2, position = position_jitter(height = 0.01)) + 
  scale_colour_manual(values = colclam, labels = clamlabel) +
  scale_fill_manual(values = colclam, labels = clamlabel) +
  #geom_point(aes(shape = Clam:Sediment, color = Clam:Sediment), size = 2) +
  labs(x = "Initial clam density",
       y = "Clams consumed") +
  theme(plot.title = element_text(hjust = -0.1, vjust = 2.5, size = 12))

#how many non-eaters ------------
JLN <- firstperiod %>% filter(Clam == "JLN")
JLN1 <- JLN %>% filter(Sediment == "1")       
length(which(JLN1$eaten == "0")) #26
JLN0 <- JLN %>% filter(Sediment == "0")
length(which(JLN0$eaten == "0")) #26
VC <- firstperiod %>% filter(Clam == "VC")
VC1 <- VC %>% filter(Sediment == "1")
length(which(VC1$eaten == "0")) #24

# COMPARISON OF CRAB SIZE IN FUNCTIONAL RESPONSE EXPERIMENT ------------
cara_compare <- lm(CW ~ Sediment + Clam, data = firstperiod)
summary(cara_compare) #carapace width is barely significantly different between clam species
simulateResiduals(cara_compare, plot = TRUE)
#residuals look good

claw_compare <- lm(Cheliped ~ Sediment + Clam, data = firstperiod)
summary(claw_compare) # claw size is significantly different between sediment and clam treatments and just barely not significant with combination
simulateResiduals(claw_compare, plot = TRUE)

#EXTRA CODE---------------------------------
#Probability of consumption summary stats
#differences between probability of consumption between species and substrate presence
prob3 <- firstperiod %>% group_by(Clam, Substrate) %>%
  summarize(mean_eaten = mean(eaten_logistic), groups = 'drop') 
prob3

#Proportion of clams consumed summary stats
#differences between proportion of clams consumed between species and substrate presence
prob4 <- firstperiod %>% group_by(Clam, Substrate) %>%
  summarize(mean_proportion = mean(proportion_eaten), groups = 'drop') 
prob4


#Linear models to assess how an increase in cheliped height results in a change in the proportion of clams consumed for each treatment combination
varnish0claw <- lm(proportion_eaten ~ Cheliped, data = VC0)
summary(varnish0claw)
varnish1claw <- lm(proportion_eaten ~ Cheliped, data = VC1)
summary(varnish1claw)
JLN0claw <- lm(proportion_eaten ~ Cheliped, data = JLN0)
summary(JLN0claw)
JLN1claw <- lm(proportion_eaten ~ Cheliped, data = JLN1)
summary(JLN1claw)