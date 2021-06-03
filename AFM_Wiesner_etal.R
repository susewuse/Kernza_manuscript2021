library(plyr)
library(dplyr)
library(ggeffects)
library(DataCombine)
library(ggplot2)
library(tidyr)
library(REddyProc)
library(data.table)
library(nlme)
library(emmeans)
library(lsmeans)
library(effects)
library(lubridate)
library(parameters)
library(performance)
library(psych)
library(car)
library(gganimate)
library(sjPlot)
library(sjmisc)
library(lattice)  
library(bigleaf)
library(cowplot)
library(wesanderson)
library(evobiR)
library(lognorm)
library(merTools)
library(AICcmodavg)
library(genefilter)

K_long <- read.csv("/Users/susannewiesner/Downloads/AFM/IWGm_IWGb_final_halfhourly_data.csv")

#### harvest factor 3 days prior and post harvest
K_long$phen_s <-   ifelse(K_long$DateTime>=as.POSIXct("2019-08-26 00:00:00")&K_long$DateTime<as.POSIXct("2019-08-29 00:00:00"), "harvest", 
                   ifelse(K_long$DateTime>=as.POSIXct("2019-08-29 00:00:00")&K_long$DateTime<as.POSIXct("2019-09-01 00:00:00"), "post-harvest1",
                   ifelse(K_long$DateTime>=as.POSIXct("2019-09-01 00:00:00")&K_long$DateTime<as.POSIXct("2019-09-04 00:00:00"), "post-harvest2",
                   ifelse(K_long$DateTime>=as.POSIXct("2019-09-04 00:00:00")&K_long$DateTime<as.POSIXct("2019-09-07 00:00:00"), "post-harvest3",
                   ifelse(K_long$DateTime>=as.POSIXct("2019-09-07 00:00:00")&K_long$DateTime<as.POSIXct("2019-09-10 00:00:00"), "post-harvest4",
                   ifelse(K_long$DateTime>=as.POSIXct("2019-09-10 00:00:00")&K_long$DateTime<as.POSIXct("2019-09-13 00:00:00"), "post-harvest5",
                   ifelse(K_long$DateTime>=as.POSIXct("2019-09-13 00:00:00")&K_long$DateTime<as.POSIXct("2019-09-16 00:00:00"), "post-harvest6",
                   ifelse(K_long$DateTime>=as.POSIXct("2019-09-16 00:00:00")&K_long$DateTime<as.POSIXct("2019-09-19 00:00:00"), "post-harvest7",
                   ifelse(K_long$DateTime>=as.POSIXct("2019-09-19 00:00:00")&K_long$DateTime<as.POSIXct("2019-09-22 00:00:00"), "post-harvest8",
                   ifelse(K_long$DateTime>=as.POSIXct("2019-09-22 00:00:00")&K_long$DateTime<as.POSIXct("2019-09-25 00:00:00"), "post-harvest9",
                   ifelse(K_long$DateTime>=as.POSIXct("2019-09-25 00:00:00")&K_long$DateTime<as.POSIXct("2019-09-28 00:00:00"), "post-harvest10",

                   ifelse(K_long$DateTime>=as.POSIXct("2020-08-04 00:00:00")&K_long$DateTime<as.POSIXct("2020-08-07 00:00:00"), "harvest", 
                   ifelse(K_long$DateTime>=as.POSIXct("2020-08-07 00:00:00")&K_long$DateTime<as.POSIXct("2020-08-10 00:00:00"), "post-harvest1", 
                   ifelse(K_long$DateTime>=as.POSIXct("2020-08-10 00:00:00")&K_long$DateTime<as.POSIXct("2020-08-13 00:00:00"), "post-harvest2", 
                   ifelse(K_long$DateTime>=as.POSIXct("2020-08-13 00:00:00")&K_long$DateTime<as.POSIXct("2020-08-16 00:00:00"), "post-harvest3",
                   ifelse(K_long$DateTime>=as.POSIXct("2020-08-16 00:00:00")&K_long$DateTime<as.POSIXct("2020-08-19 00:00:00"), "post-harvest4",
                   ifelse(K_long$DateTime>=as.POSIXct("2020-08-19 00:00:00")&K_long$DateTime<as.POSIXct("2020-08-22 00:00:00"), "post-harvest5",
                   ifelse(K_long$DateTime>=as.POSIXct("2020-08-22 00:00:00")&K_long$DateTime<as.POSIXct("2020-08-25 00:00:00"), "post-harvest6",
                   ifelse(K_long$DateTime>=as.POSIXct("2020-08-25 00:00:00")&K_long$DateTime<as.POSIXct("2020-08-28 00:00:00"), "post-harvest7",
                   ifelse(K_long$DateTime>=as.POSIXct("2020-08-28 00:00:00")&K_long$DateTime<as.POSIXct("2020-08-31 00:00:00"), "post-harvest8",
                   ifelse(K_long$DateTime>=as.POSIXct("2020-08-31 00:00:00")&K_long$DateTime<as.POSIXct("2020-09-03 00:00:00"), "post-harvest9",
                   ifelse(K_long$DateTime>=as.POSIXct("2020-09-03 00:00:00")&K_long$DateTime<as.POSIXct("2020-09-06 00:00:00"), "post-harvest10",
                   ifelse(K_long$Month%in%c("11", "12", "4", "3"), "late vegetative",
                   ifelse(K_long$Month%in%c("1", "2"), "dormant",
                   ifelse(K_long$Month%in%c("8", "9", "10"), "early vegetative", "reproductive")))))))))))))))))))))))))

K_long$phen_s <- as.factor(K_long$phen_s)

########################### Correct energy balance using Fluxnet metho ################################################
 # https://fluxnet.org/data/fluxnet2015-dataset/data-processing/
### Used G gapfilled using soil heatflux from tall tower for 2019 and beginning of 2020
#EBC_CF <- (NETRAD - G) / (H + LE)

############### calculate energy correction factor ##################
K_long2 <- subset(K_long, Year %in% c("2019", "2020"))
KC_long <- subset(K_long2, Type == "IWGb")
KO_long <- subset(K_long2, Type == "IWGm")

nday = length(KC_long$DoY)/48
hrs = c(22,22.5,23,23.5,0,0.5,1,1.5,2,2.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5)

### IWGb
i = 1
KC_long$EBC_CF_ad <- NA
KC_long$EBC_CF_med <- NA
KC_long$EBC_CF_se <- NA

  for (i in 1:nday) {
    i = i+1
    loc_start = i*48L
    loc_end = loc_start+47L
    loc_doy = KC_long$DoY[loc_start]
    loc_yr = KC_long$Year[loc_start]

    dataF <- subset(KC_long, DoY >= loc_doy-15)
    dataF <- subset(dataF, DoY < loc_doy+15)
    dataF <- subset(dataF, Hour%in%hrs)
    dataF <- subset(dataF, Year == loc_yr)
  
     if (length(dataF) == 0) next
     ebc = median((dataF$Rn_Fill-dataF$G_f)/(dataF$H_uStar_f+dataF$LE_Rn), na.rm = T)
     ebc_lm =  lm((dataF$H_uStar_f+dataF$LE_Rn)~0+(dataF$Rn_Fill-dataF$G_f), na.action = na.omit)
     EBC_CF <- summary(ebc_lm)$coefficients[1,1]
     EBC_se <- coef(summary(ebc_lm))[, "Std. Error"]

     KC_long$EBC_CF_ad[loc_start:loc_end] <- EBC_CF
     KC_long$EBC_CF_med[loc_start:loc_end] <- ebc
     KC_long$EBC_CF_se[loc_start:loc_end] <- EBC_se
   }

#### IWGm
i = 1
KO_long$EBC_CF_ad <- NA
KO_long$EBC_CF_med <- NA
KO_long$EBC_CF_se <- NA

  for (i in 1:nday) {
    i = i+1
    loc_start = i*48L
    loc_end = loc_start+47L
    loc_doy = KO_long$DoY[loc_start]
    loc_yr = KO_long$Year[loc_start]

    dataF <- subset(KO_long, DoY >= loc_doy-15)
    dataF <- subset(dataF, DoY < loc_doy+15)
    dataF <- subset(dataF, Hour%in%hrs)
    dataF <- subset(dataF, Year == loc_yr)
  
     if (length(dataF) == 0) next
     ebc = median((dataF$Rn_Fill-dataF$G_f)/(dataF$H_uStar_f+dataF$LE_Rn), na.rm = T)
     ebc_lm =  lm((dataF$H_uStar_f+dataF$LE_Rn)~0+(dataF$Rn_Fill-dataF$G_f), na.action = na.omit)
     EBC_CF <- summary(ebc_lm)$coefficients[1,1]
     EBC_se <- coef(summary(ebc_lm))[, "Std. Error"]
     
     KO_long$EBC_CF_ad[loc_start:loc_end] <- EBC_CF
     KO_long$EBC_CF_med[loc_start:loc_end] <- ebc
     KO_long$EBC_CF_se[loc_start:loc_end] <- EBC_se
  }

hist(KC_long$EBC_CF_ad)
hist(KO_long$EBC_CF_ad)

KC_long$EBC_CF_ad <- zoo::na.approx(KC_long$EBC_CF_ad, rule=2)
KO_long$EBC_CF_ad <- zoo::na.approx(KO_long$EBC_CF_ad, rule=2)

K_long <- rbind(KC_long, KO_long)

################################## fix Energy fluxes with Fluxnet factor ##########################################
K_long$LE_F <- K_long$LE_f_N/K_long$EBC_CF_ad
K_long$H_F <- K_long$H_uStar_f/K_long$EBC_CF_ad

################################# calculate ET ###############################
K_long$ET <- fCalcETfromLE(K_long$LE_F, K_long$Tair_f)

################################# calculate Time of Day factor ###########################3
K_long$ToD <- ifelse(K_long$Hour>=6&K_long$Hour<=10, "morning", 
                   ifelse(K_long$Hour>10&K_long$Hour<=14, "midday", 
                   ifelse(K_long$Hour>14&K_long$Hour<=18, "afternoon", 
                   ifelse(K_long$Hour>18&K_long$Hour<=22, "evening", "night"))))

####################### ++++++++++++++++++ growing season models ++++++++++++++++++ ###########################
########### NEE model ############
lme_NEE <- lme(NEE_uStar_f ~ phen_s*Type+
                  phen_s*VPD_f+
                  phen_s*PPFD_m+
                  PPFD_m*VPD_f+
                  PPFD_m*Type+
                  VPD_f*Type+
                  ToD ,  random = (~1 | Year/DoY), 
                  method="ML", control=lmeControl(msMaxIter=100, returnObject=TRUE),
                  data=K_long, na.action = na.omit)

AIC(lme_NEE)
car::Anova(lme_NEE)
model_performance(lme_NEE)

############################# Plots NEE ##############################################
NEE_df_phen_stage_Type <- data.frame(Effect(c("phen_s", "Type"), lme_NEE))
NEE_df_phen_stage_Type$phen_stage <- factor(NEE_df_phen_stage_Type$phen_s, levels = c("late vegetative", "reproductive", "harvest", "post-harvest1", "post-harvest2",
                                                                              "post-harvest3", "post-harvest4", "post-harvest5", "post-harvest6", "post-harvest7",
                                                                              "post-harvest8","post-harvest9",  "post-harvest10", 
                                                                              "early vegetative",  "dormant"),
                                          labels = c("lv", "r", "h", "ph1", "ph2", "ph3", "ph4", "ph5", "ph6", "ph7",
                                                                              "ph8","ph9",  "ph10", "ev",  "d"))

lsm <- emmeans(lme_NEE, ~ Type:phen_s, data = K_long)
 pairs(lsm, by = "phen_s", adjust = "mvt")
 pairs(lsm, by = "Type", adjust = "mvt")

p_phen_stage_Type_NEE <- ggplot(NEE_df_phen_stage_Type, aes(x=phen_stage, y=fit, group=Type)) +
  geom_errorbar(data = NEE_df_phen_stage_Type, aes(ymin = fit-se, ymax = fit+se), size=0.3, width=0.3, position = position_dodge(width = 0.01)) +
  geom_point(data = NEE_df_phen_stage_Type, mapping=aes(x=phen_stage, y=fit, colour=Type), size=3, position = position_dodge(width = 0.01))+
  geom_line(data = NEE_df_phen_stage_Type, mapping=aes(x=phen_stage, y=fit, colour=Type), size=1, position = position_dodge(width = 0.01))+
  ylab(expression(NEE~(mu*mol~CO[2]~m^-2~s^-1))) + xlab("") +
  theme(axis.ticks.length=unit(-0.1, "cm"),
        legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(r=0.2,unit = "cm")),
        legend.title=element_blank())+
    scale_color_manual(values = c("darkgoldenrod", "grey23"))


##################### PPFD_f + Type ###########################
NEE_df_PPFD_VPD_f <- data.frame(Effect(c("PPFD_m", "VPD_f"), lme_NEE, xlevels = list(VPD_f = c(0, 1, 2, 3))))
NEE_df_PPFD_VPD_f$VPD_f <- factor(NEE_df_PPFD_VPD_f$VPD_f) # , levels = c("0", "0.5", "1", "1.5", "2", "2.5", "3")

lsm <- emmeans(lme_NEE, ~ VPD_f:PPFD_m, at = list(VPD_f = c(0, 1, 2, 3), PPFD_m = c(0, 300, 500, 800, 1000, 2000)))
 pairs(lsm, by = "PPFD_m", adjust = "mvt")
 pairs(lsm, by = "VPD_f", adjust = "mvt")

p_PPFD_f_VPD_f_NEE <- ggplot(NEE_df_PPFD_VPD_f, aes(x=PPFD_m, y=fit, group=VPD_f)) +
  geom_errorbar(data = NEE_df_PPFD_VPD_f, aes(ymin = fit-se, ymax = fit+se), size=0.3, width=0.4, position = position_dodge(width = 0.01)) +
  geom_point(data = NEE_df_PPFD_VPD_f, mapping=aes(x=PPFD_m, y=fit, colour=VPD_f),size=2, position = position_dodge(width = 0.01))+
  geom_line(data = NEE_df_PPFD_VPD_f, mapping=aes(x=PPFD_m, y=fit, colour=VPD_f), size=1, position = position_dodge(width = 0.01))+
  ylab(expression(NEE~(mu*mol~CO[2]~m^-2~s^-1))) + xlab(expression(PPFD~(mu*mol~m^-2~s^-1)))+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(r=0.2,unit = "cm")),
        legend.title=element_blank())+
    scale_color_manual(values = c("grey23", "#66a182", "#edae49", "lightblue"))


##################### PPFD_f + phen_s ###########################
NEE_df_PPFD_f_phen_s <- data.frame(Effect(c("phen_s", "PPFD_m"), lme_NEE, xlevels = list(PPFD_m = c(0, 500, 1000, 1500, 2000))))
NEE_df_PPFD_f_phen_s$phen_stage <- factor(NEE_df_PPFD_f_phen_s$phen_s, levels = c("late vegetative", "reproductive", "harvest", "post-harvest1", "post-harvest2",
                                                                              "post-harvest3", "post-harvest4", "post-harvest5", "post-harvest6", "post-harvest7", "post-harvest8", 
                                                                              "early vegetative",  "dormant"))

lsm <- lsmeans(lme_NEE, ~ phen_s:PPFD_m, at = list(PPFD_m = c(0, 500, 1000, 1500, 2000)))
 pairs(lsm, by = "PPFD_m", adjust = "mvt")
 pairs(lsm, by = "phen_s", adjust = "mvt")

p_PPFD_f_phen_s_NEE <- ggplot(NEE_df_PPFD_f_phen_s, aes(x=PPFD_m, y=fit, group=phen_s)) +
  geom_errorbar(data = NEE_df_PPFD_f_phen_s, aes(ymin = fit-se, ymax = fit+se), size=0.3, width=300, position = position_dodge(width = 10)) +
  geom_line(data = NEE_df_PPFD_f_phen_s, mapping=aes(x=PPFD_m, y=fit, colour=phen_s), size=1, position = position_dodge(width = 10))+
 geom_point(data = NEE_df_PPFD_f_phen_s, mapping=aes(x=PPFD_m, y=fit, colour=phen_s), size=2, position = position_dodge(width = 10))+
 ylab(expression(NEE~(mu*mol~CO[2]~m^-2~s^-1))) +xlab(expression(PPFD~(mu*mol~m^-2~s^-1)))+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        legend.position="nobe",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(r=0.2,unit = "cm")),
        legend.title=element_blank())+
    scale_color_manual(values = c("#208eb7", "#9ecbf4", "#214a65", "#78d8bc", "#1c875c", "#74ee65", "#59a20c", "#bde267",
                                   "#7c8869", "#fe8f06", "#994440", "#c98d9f", "#ea1349", "#f1c039", "#4b481f"))


##################### PPFD_f + Type ###########################
NEE_df_PPFD_f_Type <- data.frame(Effect(c("Type", "PPFD_m"), lme_NEE, xlevels = list(PPFD_m = c(0, 500, 1000, 1500, 2000))))

lsm <- lsmeans(lme_NEE, ~ Type:PPFD_m, at = list(PPFD_m = c(0, 500, 1000, 1500, 2000)))
 pairs(lsm, by = "PPFD_m", adjust = "mvt")
 pairs(lsm, by = "Type", adjust = "mvt")

p_PPFD_f_Type_NEE <- ggplot(NEE_df_PPFD_f_Type, aes(x=PPFD_m, y=fit, group=Type)) +
  geom_errorbar(data = NEE_df_PPFD_f_Type, aes(ymin = fit-se, ymax = fit+se), size=0.3, width=300, position = position_dodge(width = 10)) +
  geom_line(data = NEE_df_PPFD_f_Type, mapping=aes(x=PPFD_m, y=fit, colour=Type), size=1, position = position_dodge(width = 10))+
  geom_point(data = NEE_df_PPFD_f_Type, mapping=aes(x=PPFD_m, y=fit, colour=Type), size=2, position = position_dodge(width = 10))+
  ylab(expression(NEE~(mu*mol~CO[2]~m^-2~s^-1))) + xlab(expression(PPFD~(mu*mol~m^-2~s^-1)))+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(r=0.2,unit = "cm")),
        legend.title=element_blank())+
    scale_color_manual(values = c("darkgoldenrod", "grey23"))

##################### VPD_f + Type ###########################
NEE_df_VPD_f_Type <- data.frame(Effect(c("Type", "VPD_f"), lme_NEE, xlevels = list(VPD_f = c(0, 1, 2, 3))))

lsm <- lsmeans(lme_NEE, ~ Type:VPD_f, at = list(VPD_f = c(0, 500, 1000, 1500, 2000)))
 pairs(lsm, by = "VPD_f", adjust = "mvt")
 pairs(lsm, by = "Type", adjust = "mvt")

p_VPD_f_Type_NEE <- ggplot(NEE_df_VPD_f_Type, aes(x=VPD_f, y=fit, group=Type)) +
  geom_errorbar(data = NEE_df_VPD_f_Type, aes(ymin = fit-se, ymax = fit+se), size=0.3, width=0.3, position = position_dodge(width = 0.01)) +
  geom_line(data = NEE_df_VPD_f_Type, mapping=aes(x=VPD_f, y=fit, colour=Type), size=1, position = position_dodge(width = 0.01))+
  geom_point(data = NEE_df_VPD_f_Type, mapping=aes(x=VPD_f, y=fit, colour=Type), size=2, position = position_dodge(width = 0.01))+
  ylab(expression(NEE~(mu*mol~CO[2]~m^-2~s^-1))) + xlab(expression(VPD~(kPa)))+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(r=0.2,unit = "cm")),
        legend.title=element_blank())+
    scale_color_manual(values = c("darkgoldenrod", "grey23"))

##################### VPD_f + phen_s ###########################
NEE_df_VPD_f_phen_s <- data.frame(Effect(c("phen_s", "VPD_f"), lme_NEE, xlevels = list(VPD_f = c(0, 1, 2, 3))))
NEE_df_VPD_f_phen_s$phen_stage <- factor(NEE_df_VPD_f_phen_s$phen_s, levels = c("late vegetative", "reproductive", "harvest", "post-harvest1", "post-harvest2",
                                                                              "post-harvest3", "post-harvest4", "post-harvest5", "post-harvest6", "post-harvest7", "post-harvest8", 
                                                                              "early vegetative",  "dormant"))
lsm <- lsmeans(lme_NEE, ~ phen_s:VPD_f, at = list(VPD_f = c(0, 1, 2, 3)))
 pairs(lsm, by = "VPD_f", adjust = "mvt")
 pairs(lsm, by = "phen_s", adjust = "mvt")

p_VPD_phen_s_NEE <- ggplot(NEE_df_VPD_f_phen_s, aes(x=VPD_f, y=fit, group=phen_s)) +
  geom_errorbar(data = NEE_df_VPD_f_phen_s, aes(ymin = fit-se, ymax = fit+se), size=0.3, width=0.3, position = position_dodge(width = 0.01)) +
  geom_line(data = NEE_df_VPD_f_phen_s, mapping=aes(x=VPD_f, y=fit, colour=phen_s), size=1,position = position_dodge(width = 0.01))+
  geom_point(data = NEE_df_VPD_f_phen_s, mapping=aes(x=VPD_f, y=fit, colour=phen_s), size=2, position = position_dodge(width = 0.01))+
  ylab(expression(NEE~(mu*mol~CO[2]~m^-2~s^-1))) + xlab(expression("VPD (kPa)")) +
  theme(axis.ticks.length=unit(-0.1, "cm"),
        legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(r=0.2,unit = "cm")),
        legend.title=element_blank())+
    scale_color_manual(values = c("#208eb7", "#9ecbf4", "#214a65", "#78d8bc", "#1c875c", "#74ee65", "#59a20c", "#bde267",
                                   "#7c8869", "#fe8f06", "#994440", "#c98d9f", "#ea1349", "#f1c039", "#4b481f"))

#########################################################################################
###################### ++++++++++++++ GEE model ++++++++++++++++ ########################
#########################################################################################
lme_GEE <- lme(GPP_uStar_f ~ phen_s*Type+
                  phen_s*VPD_f+
                  phen_s*PPFD_m+
                  PPFD_m*VPD_f+
                  VPD_f+
                  PPFD_m*Type+
                  VPD_f*Type+
                  ToD , random = (~1 | Year/DoY),
                  method="ML", control=lmeControl(msMaxIter=100, returnObject=TRUE),
                  data=K_long, na.action = na.omit)

plot(lme_GEE) 
car::Anova(lme_GEE)
AIC(lme_GEE)
model_performance(lme_GEE)

############################# Plots GEE ##############################################
GEE_df_phen_stage_Type <- data.frame(Effect(c("phen_s", "Type"), lme_GEE))
GEE_df_phen_stage_Type$phen_stage <- factor(GEE_df_phen_stage_Type$phen_s, levels = c("late vegetative", "reproductive", "harvest", "post-harvest1", "post-harvest2",
                                                                              "post-harvest3", "post-harvest4", "post-harvest5", "post-harvest6", "post-harvest7",
                                                                              "post-harvest8","post-harvest9",  "post-harvest10", 
                                                                              "early vegetative",  "dormant"),
                                          labels = c("lv", "r", "h", "ph1", "ph2", "ph3", "ph4", "ph5", "ph6", "ph7",
                                                                              "ph8","ph9",  "ph10", "ev",  "d"))

lsm <- emmeans(lme_GEE, ~ Type:phen_s)
 pairs(lsm, by = "phen_s", adjust = "mvt")
 pairs(lsm, by = "Type", adjust = "mvt")

p_phen_stage_Type_GEE <- ggplot(GEE_df_phen_stage_Type, aes(x=phen_stage, y=fit, group=Type)) +
  geom_errorbar(data = GEE_df_phen_stage_Type, aes(ymin = fit-se, ymax = fit+se), size=0.3, width=0.3, position = position_dodge(width = 0.01)) +
  geom_point(data = GEE_df_phen_stage_Type, mapping=aes(x=phen_stage, y=fit, colour=Type), size=3, position = position_dodge(width = 0.01))+
  geom_line(data = GEE_df_phen_stage_Type, mapping=aes(x=phen_stage, y=fit, colour=Type), size=1, position = position_dodge(width = 0.01))+
  ylab(expression(GEE~(mu*mol~CO[2]~m^-2~s^-1))) + xlab("") +
  theme(axis.ticks.length=unit(-0.1, "cm"),
        legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(r=0.2,unit = "cm")),
        legend.title=element_blank())+
    scale_color_manual(values = c("darkgoldenrod", "grey23"))


##################### PPFD_f + Type ###########################
GEE_df_PPFD_VPD_f <- data.frame(Effect(c("PPFD_m", "VPD_f"), lme_GEE, xlevels = list(VPD_f = c(0, 1, 2, 3))))
GEE_df_PPFD_VPD_f$VPD_f <- factor(GEE_df_PPFD_VPD_f$VPD_f) # , levels = c("0", "0.5", "1", "1.5", "2", "2.5", "3")

lsm <- emmeans(lme_GEE, ~ VPD_f:PPFD_m, at = list(VPD_f = c(0, 1, 2, 3), PPFD_m = c(0, 500, 1000, 1500, 2000)))
 pairs(lsm, by = "PPFD_m", adjust = "mvt")
 pairs(lsm, by = "VPD_f", adjust = "mvt")

p_PPFD_f_VPD_f_GEE <- ggplot(GEE_df_PPFD_VPD_f, aes(x=PPFD_m, y=fit, group=VPD_f)) +
  geom_errorbar(data = GEE_df_PPFD_VPD_f, aes(ymin = fit-se, ymax = fit+se), size=0.3, width=0.4, position = position_dodge(width = 0.01)) +
  geom_point(data = GEE_df_PPFD_VPD_f, mapping=aes(x=PPFD_m, y=fit, colour=VPD_f),size=2, position = position_dodge(width = 0.01))+
  geom_line(data = GEE_df_PPFD_VPD_f, mapping=aes(x=PPFD_m, y=fit, colour=VPD_f), size=1, position = position_dodge(width = 0.01))+
  ylab(expression(GEE~(mu*mol~CO[2]~m^-2~s^-1))) + xlab(expression(PPFD~(mu*mol~m^-2~s^-1)))+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(r=0.2,unit = "cm")),
        legend.title=element_blank())+
    scale_color_manual(values = c("grey23", "#66a182", "#edae49", "lightblue"))


##################### PPFD_f + phen_s ###########################
GEE_df_PPFD_f_phen_s <- data.frame(Effect(c("phen_s", "PPFD_m"), lme_GEE))
GEE_df_PPFD_f_phen_s$phen_stage <- factor(GEE_df_PPFD_f_phen_s$phen_s, levels = c("late vegetative", "reproductive", "harvest", "post-harvest1", "post-harvest2",
                                                                              "post-harvest3", "post-harvest4", "post-harvest5", "post-harvest6", "post-harvest7", "post-harvest8", 
                                                                              "early vegetative",  "dormant"))
lsm <- lsmeans(lme_GEE, ~ phen_s:PPFD_m, at = list(PPFD_m = c(0, 500, 1000, 1500, 2000)))
 pairs(lsm, by = "phen_s", adjust = "mvt")

p_PPFD_f_phen_s_GEE <- ggplot(GEE_df_PPFD_f_phen_s, aes(x=PPFD_m, y=fit, group=phen_s)) +
  geom_errorbar(data = GEE_df_PPFD_f_phen_s, aes(ymin = fit-se, ymax = fit+se), size=0.3, width=300, position = position_dodge(width = 10)) +
  geom_line(data = GEE_df_PPFD_f_phen_s, mapping=aes(x=PPFD_m, y=fit, colour=phen_s), size=1, position = position_dodge(width = 10))+
  geom_point(data = GEE_df_PPFD_f_phen_s, mapping=aes(x=PPFD_m, y=fit, colour=phen_s), size=2, position = position_dodge(width = 10))+
  ylab(expression(GEE~(mu*mol~CO[2]~m^-2~s^-1))) + xlab(expression(PPFD~(mu*mol~m^-2~s^-1)))+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(r=0.2,unit = "cm")),
        legend.title=element_blank())+
    scale_color_manual(values = c("#208eb7", "#9ecbf4", "#214a65", "#78d8bc", "#1c875c", "#74ee65", "#59a20c", "#bde267",
                                   "#7c8869", "#fe8f06", "#994440", "#c98d9f", "#ea1349", "#f1c039", "#4b481f"))


##################### PPFD_f + Type ###########################
GEE_df_PPFD_f_Type <- data.frame(Effect(c("Type", "PPFD_m"), lme_GEE))

lsm <- lsmeans(lme_GEE, ~ Type:PPFD_m, at = list(PPFD_m = c(0, 500, 1000, 1500, 2000)))
 pairs(lsm, by = "PPFD_m", adjust = "mvt")
 pairs(lsm, by = "Type", adjust = "mvt")

p_PPFD_f_Type_GEE <- ggplot(GEE_df_PPFD_f_Type, aes(x=PPFD_m, y=fit, group=Type)) +
  geom_errorbar(data = GEE_df_PPFD_f_Type, aes(ymin = fit-se, ymax = fit+se), size=0.3, width=300, position = position_dodge(width = 10)) +
  geom_line(data = GEE_df_PPFD_f_Type, mapping=aes(x=PPFD_m, y=fit, colour=Type), size=1, position = position_dodge(width = 10))+
  geom_point(data = GEE_df_PPFD_f_Type, mapping=aes(x=PPFD_m, y=fit, colour=Type), size=2, position = position_dodge(width = 10))+
  ylab(expression(GEE~(mu*mol~CO[2]~m^-2~s^-1))) + xlab(expression(PPFD~(mu*mol~m^-2~s^-1)))+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(r=0.2,unit = "cm")),
        legend.title=element_blank())+
    scale_color_manual(values = c("darkgoldenrod", "grey23"))

##################### VPD_f + phen_s ###########################
GEE_df_VPD_f_phen_s <- data.frame(Effect(c("phen_s", "VPD_f"), lme_GEE))
GEE_df_VPD_f_phen_s$phen_stage <- factor(GEE_df_VPD_f_phen_s$phen_s, levels = c("late vegetative", "reproductive", "harvest", "post-harvest1", "post-harvest2",
                                                                              "post-harvest3", "post-harvest4", "post-harvest5", "post-harvest6", "post-harvest7", "post-harvest8", 
                                                                              "early vegetative",  "dormant"))
lsm <- lsmeans(lme_GEE, ~ phen_s:VPD_f, at = list(VPD_f = c(0, 1, 2, 3)))
 pairs(lsm, by = "VPD_f", adjust = "mvt")
 pairs(lsm, by = "phen_s", adjust = "mvt")

p_VPD_phen_s_GEE <- ggplot(GEE_df_VPD_f_phen_s, aes(x=VPD_f, y=fit, group=phen_s)) +
  geom_errorbar(data = GEE_df_VPD_f_phen_s, aes(ymin = fit-se, ymax = fit+se), size=0.3, width=0.3, position = position_dodge(width = 0.01)) +
  geom_line(data = GEE_df_VPD_f_phen_s, mapping=aes(x=VPD_f, y=fit, colour=phen_s), size=1,position = position_dodge(width = 0.01))+
  geom_point(data = GEE_df_VPD_f_phen_s, mapping=aes(x=VPD_f, y=fit, colour=phen_s), size=2, position = position_dodge(width = 0.01))+
  ylab(expression(GEE~(mu*mol~CO[2]~m^-2~s^-1))) + xlab(expression("VPD (kPa)")) +
  theme(axis.ticks.length=unit(-0.1, "cm"),
        legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(r=0.2,unit = "cm")),
        legend.title=element_blank())+
    scale_color_manual(values = c("#208eb7", "#9ecbf4", "#214a65", "#78d8bc", "#1c875c", "#74ee65", "#59a20c", "#bde267",
                                   "#7c8869", "#fe8f06", "#994440", "#c98d9f", "#ea1349", "#f1c039", "#4b481f"))

##################### PPFD_f + Type ###########################
GEE_df_VPD_f_Type <- data.frame(Effect(c("Type", "VPD_f"), lme_GEE))

lsm <- lsmeans(lme_GEE, ~ Type:VPD_f, at = list(VPD_f = c(0, 0.1, 1, 2, 3)))
 pairs(lsm, by = "VPD_f", adjust = "mvt")
 pairs(lsm, by = "Type", adjust = "mvt")

p_VPD_f_Type_GEE <- ggplot(GEE_df_VPD_f_Type, aes(x=VPD_f, y=fit, group=Type)) +
  geom_errorbar(data = GEE_df_VPD_f_Type, aes(ymin = fit-se, ymax = fit+se), size=0.3, width=0.3, position = position_dodge(width = 0.10)) +
  geom_line(data = GEE_df_VPD_f_Type, mapping=aes(x=VPD_f, y=fit, colour=Type), size=1, position = position_dodge(width = 0.10))+
  geom_point(data = GEE_df_VPD_f_Type, mapping=aes(x=VPD_f, y=fit, colour=Type), size=2, position = position_dodge(width = 0.10))+
  ylab(expression(GEE~(mu*mol~CO[2]~m^-2~s^-1))) + xlab(expression(VPD~(kPa)))+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(r=0.2,unit = "cm")),
        legend.title=element_blank())+
    scale_color_manual(values = c("darkgoldenrod", "grey23"))


#########################################################################################
###################### ++++++++++++++ GEE model ++++++++++++++++ ########################
#########################################################################################
lme_Reco <- lme(Reco_uStar ~  phen_s*Type+
                  phen_s*VPD_f+
                  phen_s*PPFD_m+
                  PPFD_m*VPD_f+
                  VPD_f+
                  PPFD_m*Type+
                  VPD_f*Type+
                  ToD , random = (~1 | Year/DoY),
                method="ML", control=lmeControl(msMaxIter=100, returnObject=TRUE),
                data=K_long, na.action = na.omit)

plot(lme_Reco) 
car::Anova(lme_Reco)
AIC(lme_Reco)
model_performance(lme_Reco)

############################# Plots Reco ##############################################
Reco_df_phen_stage_Type <- data.frame(Effect(c("phen_s", "Type"), lme_Reco))
Reco_df_phen_stage_Type$phen_stage <- factor(Reco_df_phen_stage_Type$phen_s, levels = c("late vegetative", "reproductive", "harvest", "post-harvest1", "post-harvest2",
                                                                              "post-harvest3", "post-harvest4", "post-harvest5", "post-harvest6", "post-harvest7",
                                                                              "post-harvest8","post-harvest9",  "post-harvest10", 
                                                                              "early vegetative",  "dormant"),
                                          labels = c("lv", "r", "h", "ph1", "ph2", "ph3", "ph4", "ph5", "ph6", "ph7",
                                                                              "ph8","ph9",  "ph10", "ev",  "d"))

lsm <- emmeans(lme_Reco, ~ Type:phen_s)
 pairs(lsm, by = "phen_s", adjust = "mvt")
 pairs(lsm, by = "Type", adjust = "mvt")

p_phen_stage_Type_Reco <- ggplot(Reco_df_phen_stage_Type, aes(x=phen_stage, y=fit, group=Type)) +
  geom_errorbar(data = Reco_df_phen_stage_Type, aes(ymin = fit-se, ymax = fit+se), size=0.3, width=0.3, position = position_dodge(width = 0.01)) +
  geom_point(data = Reco_df_phen_stage_Type, mapping=aes(x=phen_stage, y=fit, colour=Type), size=3, position = position_dodge(width = 0.01))+
  geom_line(data = Reco_df_phen_stage_Type, mapping=aes(x=phen_stage, y=fit, colour=Type), size=1, position = position_dodge(width = 0.01))+
  ylab(expression(R[eco]~(mu*mol~CO[2]~m^-2~s^-1))) + xlab("") +
  theme(axis.ticks.length=unit(-0.1, "cm"),
        legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(r=0.2,unit = "cm")),
        legend.title=element_blank())+
    scale_color_manual(values = c("darkgoldenrod", "grey23"))


##################### PPFD_f + Type ###########################
Reco_df_PPFD_VPD_f <- data.frame(Effect(c("PPFD_m", "VPD_f"), lme_Reco, xlevels = list(VPD_f = c(0, 1, 2, 3))))
Reco_df_PPFD_VPD_f$VPD_f <- factor(Reco_df_PPFD_VPD_f$VPD_f) # , levels = c("0", "0.5", "1", "1.5", "2", "2.5", "3")

lsm <- emmeans(lme_Reco, ~ VPD_f:PPFD_m, at = list(VPD_f = c(0, 1, 2, 3), PPFD_m = c(0, 500, 1000, 1500, 2000)))
 pairs(lsm, by = "PPFD_m", adjust = "mvt")
 pairs(lsm, by = "VPD_f", adjust = "mvt")

p_PPFD_f_VPD_f_Reco <- ggplot(Reco_df_PPFD_VPD_f, aes(x=PPFD_m, y=fit, group=VPD_f)) +
  geom_errorbar(data = Reco_df_PPFD_VPD_f, aes(ymin = fit-se, ymax = fit+se), size=0.3, width=0.4, position = position_dodge(width = 0.01)) +
  geom_point(data = Reco_df_PPFD_VPD_f, mapping=aes(x=PPFD_m, y=fit, colour=VPD_f),size=2, position = position_dodge(width = 0.01))+
  geom_line(data = Reco_df_PPFD_VPD_f, mapping=aes(x=PPFD_m, y=fit, colour=VPD_f), size=1, position = position_dodge(width = 0.01))+
  ylab(expression(R[eco]~(mu*mol~CO[2]~m^-2~s^-1))) + xlab(expression(PPFD~(mu*mol~m^-2~s^-1)))+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(r=0.2,unit = "cm")),
        legend.title=element_blank())+
    scale_color_manual(values = c("grey23", "#66a182", "#edae49", "lightblue"))


##################### PPFD_f + phen_s ###########################
Reco_df_PPFD_f_phen_s <- data.frame(Effect(c("phen_s", "PPFD_m"), lme_Reco))
Reco_df_PPFD_f_phen_s$phen_stage <- factor(Reco_df_PPFD_f_phen_s$phen_s, levels = c("late vegetative", "reproductive", "harvest", "post-harvest1", "post-harvest2",
                                                                              "post-harvest3", "post-harvest4", "post-harvest5", "post-harvest6", "post-harvest7", "post-harvest8", 
                                                                              "early vegetative",  "dormant"))

lsm <- lsmeans(lme_Reco, ~ phen_s:PPFD_m, at = list(PPFD_m = c(0, 500, 1000, 1500, 2000)))
 pairs(lsm, by = "PPFD_m", adjust = "mvt")
 pairs(lsm, by = "phen_s", adjust = "mvt")

p_PPFD_f_phen_s_Reco <- ggplot(Reco_df_PPFD_f_phen_s, aes(x=PPFD_m, y=fit, group=phen_s)) +
  geom_errorbar(data = Reco_df_PPFD_f_phen_s, aes(ymin = fit-se, ymax = fit+se), size=0.3, width=300, position = position_dodge(width = 10)) +
  geom_line(data = Reco_df_PPFD_f_phen_s, mapping=aes(x=PPFD_m, y=fit, colour=phen_s), size=1, position = position_dodge(width = 10))+
  geom_point(data = Reco_df_PPFD_f_phen_s, mapping=aes(x=PPFD_m, y=fit, colour=phen_s), size=2, position = position_dodge(width = 10))+
  ylab(expression(R[eco]~(mu*mol~CO[2]~m^-2~s^-1))) + xlab(expression(PPFD~(mu*mol~m^-2~s^-1)))+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(r=0.2,unit = "cm")),
        legend.title=element_blank())+
    scale_color_manual(values = c("#208eb7", "#9ecbf4", "#214a65", "#78d8bc", "#1c875c", "#74ee65", "#59a20c", "#bde267",
                                   "#7c8869", "#fe8f06", "#994440", "#c98d9f", "#ea1349", "#f1c039", "#4b481f"))


##################### PPFD_f + Type ###########################
Reco_df_PPFD_f_Type <- data.frame(Effect(c("Type", "PPFD_m"), lme_Reco))

lsm <- lsmeans(lme_Reco, ~ Type:PPFD_m , at = list(PPFD_m = c(0, 500, 1000, 1500, 2000)))
 pairs(lsm, by = "PPFD_m", adjust = "mvt")
 pairs(lsm, by = "Type", adjust = "mvt")

p_PPFD_f_Type_Reco <- ggplot(Reco_df_PPFD_f_Type, aes(x=PPFD_m, y=fit, group=Type)) +
  geom_errorbar(data = Reco_df_PPFD_f_Type, aes(ymin = fit-se, ymax = fit+se), size=0.3, width=300, position = position_dodge(width = 10)) +
  geom_line(data = Reco_df_PPFD_f_Type, mapping=aes(x=PPFD_m, y=fit, colour=Type), size=1, position = position_dodge(width = 10))+
  geom_point(data = Reco_df_PPFD_f_Type, mapping=aes(x=PPFD_m, y=fit, colour=Type), size=2, position = position_dodge(width = 10))+
  ylab(expression(R[eco]~(mu*mol~CO[2]~m^-2~s^-1))) +xlab(expression(PPFD~(mu*mol~m^-2~s^-1)))+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(r=0.2,unit = "cm")),
        legend.title=element_blank())+
    scale_color_manual(values = c("darkgoldenrod", "grey23"))


##################### VPD_f + phen_s ###########################
Reco_df_VPD_f_phen_s <- data.frame(Effect(c("phen_s", "VPD_f"), lme_Reco))
Reco_df_VPD_f_phen_s$phen_stage <- factor(Reco_df_VPD_f_phen_s$phen_s, levels = c("late vegetative", "reproductive", "harvest", "post-harvest1", "post-harvest2",
                                                                              "post-harvest3", "post-harvest4", "post-harvest5", "post-harvest6", "post-harvest7", "post-harvest8", 
                                                                              "early vegetative",  "dormant"))
lsm <- lsmeans(lme_Reco, ~ phen_s:VPD_f, at = list(VPD_f = c(0, 1, 2, 3)))
 pairs(lsm, by = "VPD_f", adjust = "mvt")
 pairs(lsm, by = "phen_s", adjust = "mvt")

p_VPD_phen_s_Reco <- ggplot(Reco_df_VPD_f_phen_s, aes(x=VPD_f, y=fit, group=phen_s)) +
  geom_errorbar(data = Reco_df_VPD_f_phen_s, aes(ymin = fit-se, ymax = fit+se), size=0.3, width=0.3, position = position_dodge(width = 0.01)) +
  geom_line(data = Reco_df_VPD_f_phen_s, mapping=aes(x=VPD_f, y=fit, colour=phen_s), size=1,position = position_dodge(width = 0.01))+
  geom_point(data = Reco_df_VPD_f_phen_s, mapping=aes(x=VPD_f, y=fit, colour=phen_s), size=2, position = position_dodge(width = 0.01))+
  ylab(expression(R[eco]~(mu*mol~CO[2]~m^-2~s^-1))) + xlab(expression("VPD (kPa)")) +
  theme(axis.ticks.length=unit(-0.1, "cm"),
        legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(r=0.2,unit = "cm")),
        legend.title=element_blank())+
    scale_color_manual(values = c("#208eb7", "#9ecbf4", "#214a65", "#78d8bc", "#1c875c", "#74ee65", "#59a20c", "#bde267",
                                   "#7c8869", "#fe8f06", "#994440", "#c98d9f", "#ea1349", "#f1c039", "#4b481f"))

##################### VPD_f + Type ###########################
Reco_df_VPD_f_Type <- data.frame(Effect(c("Type", "VPD_f"), lme_Reco))

lsm <- lsmeans(lme_Reco, ~ Type:VPD_f, at = list(VPD_f = c(0, 1, 2, 3)))
 pairs(lsm, by = "VPD_f", adjust = "mvt")
 pairs(lsm, by = "Type", adjust = "mvt")

p_VPD_Type_Reco <- ggplot(Reco_df_VPD_f_Type, aes(x=VPD_f, y=fit, group=Type)) +
  geom_errorbar(data = Reco_df_VPD_f_Type, aes(ymin = fit-se, ymax = fit+se), size=0.3, width=0.3, position = position_dodge(width = 0.01)) +
  geom_line(data = Reco_df_VPD_f_Type, mapping=aes(x=VPD_f, y=fit, colour=Type), size=1,position = position_dodge(width = 0.01))+
  geom_point(data = Reco_df_VPD_f_Type, mapping=aes(x=VPD_f, y=fit, colour=Type), size=2, position = position_dodge(width = 0.01))+
  ylab(expression(R[eco]~(mu*mol~CO[2]~m^-2~s^-1))) + xlab(expression("VPD (kPa)")) +
  theme(axis.ticks.length=unit(-0.1, "cm"),
        legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        axis.text.x = element_text(margin = margin(t = 0.2, unit = "cm")),
        axis.text.y = element_text(margin = margin(r=0.2,unit = "cm")),
        legend.title=element_blank())+
    scale_color_manual(values = c("darkgoldenrod", "grey23"))


########### main text ########
f <- ggplot_build(p_PPFD_f_phen_s_NEE)
f$plot$layers[[1]]$geom_params$width  <- 0.1/4*15 * diff(f$layout$panel_params[[1]]$x.range)
p_PPFD_f_phen_s_NEE <- f$plot

f <- ggplot_build(p_VPD_phen_s_NEE)
f$plot$layers[[1]]$geom_params$width  <- 0.1/4*15 * diff(f$layout$panel_params[[1]]$x.range)
p_VPD_phen_s_NEE <- f$plot

f <- ggplot_build(p_PPFD_f_phen_s_GEE)
f$plot$layers[[1]]$geom_params$width  <- 0.1/4*15 * diff(f$layout$panel_params[[1]]$x.range)
p_PPFD_f_phen_s_GEE <- f$plot

f <- ggplot_build(p_VPD_phen_s_GEE)
f$plot$layers[[1]]$geom_params$width  <- 0.1/4*15 * diff(f$layout$panel_params[[1]]$x.range)
p_VPD_phen_s_GEE <- f$plot

f <- ggplot_build(p_PPFD_f_phen_s_Reco)
f$plot$layers[[1]]$geom_params$width  <- 0.1/4*15 * diff(f$layout$panel_params[[1]]$x.range)
p_PPFD_f_phen_s_Reco <- f$plot

f <- ggplot_build(p_VPD_phen_s_Reco)
f$plot$layers[[1]]$geom_params$width  <- 0.1/4*15 * diff(f$layout$panel_params[[1]]$x.range)
p_VPD_phen_s_Reco <- f$plot

Figure7 <- cowplot::plot_grid(p_PPFD_f_phen_s_NEE, p_VPD_phen_s_NEE,
                        p_PPFD_f_phen_s_GEE, p_VPD_phen_s_GEE,
                        p_PPFD_f_phen_s_Reco, p_VPD_phen_s_Reco,
  labels = c("a", "b", "c", "d", "e", "f"), align = "hv",
  rel_widths = c(1, 1),
  nrow = 3, ncol = 2
)

#cowplot::save_plot(path.expand("/Users/susannewiesner/Documents/AFM/plots/Figure7.pdf"), 
#                   Figure7, base_height = 11, base_width = 8)

############################
f <- ggplot_build(p_PPFD_f_VPD_f_NEE)
f$plot$layers[[1]]$geom_params$width  <- 0.1 * diff(f$layout$panel_params[[1]]$x.range)
p_PPFD_f_VPD_f_NEE <- f$plot

f <- ggplot_build(p_PPFD_f_VPD_f_GEE)
f$plot$layers[[1]]$geom_params$width  <- 0.1 * diff(f$layout$panel_params[[1]]$x.range)
p_PPFD_f_VPD_f_GEE <- f$plot

f <- ggplot_build(p_PPFD_f_VPD_f_Reco)
f$plot$layers[[1]]$geom_params$width  <- 0.1 * diff(f$layout$panel_params[[1]]$x.range)
p_PPFD_f_VPD_f_Reco <- f$plot

f <- ggplot_build(p_PPFD_f_Type_NEE)
f$plot$layers[[1]]$geom_params$width  <- 0.05 * diff(f$layout$panel_params[[1]]$x.range)
p_PPFD_f_Type_NEE <- f$plot

f <- ggplot_build(p_PPFD_f_Type_GEE)
f$plot$layers[[1]]$geom_params$width  <- 0.05 * diff(f$layout$panel_params[[1]]$x.range)
p_PPFD_f_Type_GEE <- f$plot

f <- ggplot_build(p_PPFD_f_Type_Reco)
f$plot$layers[[1]]$geom_params$width  <- 0.05 * diff(f$layout$panel_params[[1]]$x.range)
p_PPFD_f_Type_Reco <- f$plot

f <- ggplot_build(p_VPD_Type_Reco)
f$plot$layers[[1]]$geom_params$width  <- 0.05 * diff(f$layout$panel_params[[1]]$x.range)
p_VPD_Type_Reco <- f$plot

f <- ggplot_build(p_VPD_f_Type_GEE)
f$plot$layers[[1]]$geom_params$width  <- 0.05 * diff(f$layout$panel_params[[1]]$x.range)
p_VPD_f_Type_GEE <- f$plot

f <- ggplot_build(p_VPD_f_Type_NEE)
f$plot$layers[[1]]$geom_params$width  <- 0.05 * diff(f$layout$panel_params[[1]]$x.range)
p_VPD_f_Type_NEE <- f$plot

Figure8 <- cowplot::plot_grid(p_PPFD_f_VPD_f_NEE,  p_PPFD_f_VPD_f_GEE, p_PPFD_f_VPD_f_Reco, p_PPFD_f_Type_NEE,
                        p_PPFD_f_Type_GEE, p_PPFD_f_Type_Reco, p_VPD_f_Type_NEE, p_VPD_f_Type_GEE, p_VPD_Type_Reco,
  labels = c("a", "b", "c", "d","e", "f", "g", "h", "i"), align = "hv",
  rel_widths = c(1, 1, 1),
  nrow = 3, ncol =3
)

#cowplot::save_plot(path.expand("/Users/susannewiesner/Documents/AFM/plots/Figure8.pdf"), 
#                   Figure8, base_height = 11, base_width = 12)


############################
Figure6 <- cowplot::plot_grid(p_phen_stage_Type_NEE,  p_phen_stage_Type_GEE, p_phen_stage_Type_Reco,
                        labels = c("a", "b", "c"), align = "hv", nrow = 3, ncol = 1)

#cowplot::save_plot(path.expand("/Users/susannewiesner/Documents/AFM/plots/Figure6.pdf"), 
#                   Figure6, base_height = 12, base_width = 9)


##################### Weather Stations ######################3
lodi_w <- read.csv("/Users/susannewiesner/Downloads/AFM/Weatherstations/Lodi_weather_precip.csv", skip = 7, header = F, na.strings = "M")
lodi_w_H <- fread("/Users/susannewiesner/Downloads/AFM/Weatherstations/Lodi_weather_precip.csv", skip = 5, header = T, nrows = 1, stringsAsFactors = F)
colnames(lodi_w) <- colnames(lodi_w_H)
lodi_w$date <- as.Date(lodi_w$Date, "%m/%d/%y")
lodi_w$temp_max <- as.numeric(lodi_w$temp_max)

baraboo_w <- read.csv("/Users/susannewiesner/Downloads/AFM/Weatherstations/Baraboo_precip.csv", skip = 7, header = F, na.strings = "M")
baraboo_w_H <- fread("/Users/susannewiesner/Downloads/AFM/Weatherstations/Baraboo_precip.csv", skip = 5, header = T, nrows = 1, stringsAsFactors = F)
colnames(baraboo_w) <- colnames(baraboo_w_H)
baraboo_w$date <- as.Date(baraboo_w$Date, "%m/%d/%y")
baraboo_w$precipitation <- as.numeric(baraboo_w$precipitation)
baraboo_w$temp_max <- as.numeric(baraboo_w$max_temp)

Sauk_w <- read.csv("/Users/susannewiesner/Downloads/AFM/Weatherstations/Sauk_city_precip.csv", skip = 7, header = F, na.strings = "M")
Sauk_w_H <- fread("/Users/susannewiesner/Downloads/AFM/Weatherstations/Sauk_city_precip.csv", skip = 5, header = T, nrows = 1, stringsAsFactors = F)
colnames(Sauk_w) <- colnames(Sauk_w_H)
Sauk_w$date <- as.Date(Sauk_w$Date, "%m/%d/%y")
Sauk_w$precipitation <- as.numeric(Sauk_w$precipitation)
Sauk_w$temp_max <- as.numeric(Sauk_w$max_temp)

#################### summarize daily data #######################
K_long_d <- K_long %>%
  dplyr::group_by(Year, DoY, Type, phen_s) %>%
  dplyr::summarize(Tair_m = mean(Tair_f, na.rm = T),
                   TairNA = mean(Tair_side, na.rm = T),
                   VPD_m = mean(VPD_f, na.rm = T),
                   VPDNA = mean(VPD_side, na.rm = T),
                   SWin_D = mean(SWIN_1_1_1, na.rm = T),
                   SWout_D = mean(SWOUT_1_1_1, na.rm = T),
                   LWout_D = mean(LWOUT_1_1_1, na.rm = T),
                   LWin_D = mean(LWIN_1_1_1, na.rm = T),
                   Rn_D = mean(Rn_1_1_1, na.rm = T),
                   Rn_D_s = sum(Rn_1_1_1, na.rm = T),
                   Rn_m = mean(Rn_Fill, na.rm = T),
                   Rn_s = sum(Rn_Fill, na.rm = T),
                   PPFD_m_NA = mean(PPFD_1_1_1, na.rm = T),
                   PPFD_m = mean(PPFD_m, na.rm = T),
                   PPFD_s = sum(PPFD_m, na.rm = T),
                   rain_s = sum(P_RAIN_1_1_1, na.rm = F),
                   P_m = mean(air_pressure, na.rm = T),
                   NEE_m = mean(NEE_uStar_f, na.rm = T),
                   NEE_s = sum(bigleaf::umolCO2.to.gC(NEE_uStar_f)/24/2, na.rm = T), #g C m-2
                   NEE_s05 = sum(bigleaf::umolCO2.to.gC(NEE_U05_f)/24/2, na.rm = T),
                   NEE_s50 = sum(bigleaf::umolCO2.to.gC(NEE_U50_f)/24/2, na.rm = T),
                   NEE_s95 = sum(bigleaf::umolCO2.to.gC(NEE_U95_f)/24/2, na.rm = T),
                   GEE_s = sum(bigleaf::umolCO2.to.gC(GPP_uStar_f)/24/2, na.rm = T),
                   GEE_s05 = sum(bigleaf::umolCO2.to.gC(GPP_U05_f)/24/2, na.rm = T),
                   GEE_s50 = sum(bigleaf::umolCO2.to.gC(GPP_U50_f)/24/2, na.rm = T),
                   GEE_s95 = sum(bigleaf::umolCO2.to.gC(GPP_U95_f)/24/2, na.rm = T),
                   Reco_s = sum(bigleaf::umolCO2.to.gC(Reco_uStar)/24/2, na.rm = T),
                   Reco_s05 = sum(bigleaf::umolCO2.to.gC(Reco_U05)/24/2, na.rm = T), # half hourly data added to daily
                   Reco_s50 = sum(bigleaf::umolCO2.to.gC(Reco_U50)/24/2, na.rm = T),
                   Reco_s95 = sum(bigleaf::umolCO2.to.gC(Reco_U95)/24/2, na.rm = T),
                   GEE_m = mean(bigleaf::umolCO2.to.gC(GPP_uStar_f)/24/2, na.rm = T),
                   GEE_LUEm = mean(GPP_uStar_f, na.rm = T), # umol CO2 m-2 s-2
                   GEE_LUEs = sum(GPP_uStar_f, na.rm = T),
                   Reco_m = mean(bigleaf::umolCO2.to.gC(Reco_uStar)/24/2, na.rm = T), # concerted to g C s-1
                   LE_s = sum(LE_F, na.rm = T),
                   ET_s = sum(ET, na.rm = T),
                   H_s = sum(H_F, na.rm = T),
                   LE_m = mean(LE_F, na.rm = T),
                   ET_m = mean(ET, na.rm = T),
                   H_m = mean(H_F, na.rm = T),
                   G_m = mean(G_f, na.rm = T),
                   G_s = sum(G_f, na.rm = T)) 

K_long_d = K_long_d %>% group_by(Type) %>% arrange(Year, DoY) %>% mutate(cs = cumsum(NEE_s))
K_long_d <- as.data.frame(K_long_d)
K_long_d$date <- as.Date(paste0(K_long_d$Year, K_long_d$DoY), format="%Y%j")
K_long_d$date2 <- as.POSIXct(K_long_d$date)
K_long_d <- tidyr::separate(K_long_d, date, into = c("Year","Month", "Day"), sep = "-", extra = "merge", remove=F)
K_long_d$Month <- factor(K_long_d$Month)
K_long_d$Year <- factor(K_long_d$Year)

K_long_d <- join(K_long_d, lodi_w[, c("date", "precipitation", "temp_max")], by = "date", type = "left", match = "first")
K_long_d$precipitation <- as.numeric(K_long_d$precipitation)
colnames(K_long_d)[52:53] <- c("precip_lodi", "mtemp_l")

K_long_d <- join(K_long_d, Sauk_w[, c("date", "precipitation", "temp_max")], by = "date", type = "left", match = "first")
K_long_d$precipitation <- as.numeric(K_long_d$precipitation)
colnames(K_long_d)[54:55] <- c("precip_sauk", "mtemp_sauk")

K_long_d$T_lodi <- ifelse(is.na(K_long_d$Tair_m), lead(K_long_d$mtemp_l, 2)*0.8643-4.2, K_long_d$Tair_m)
K_long_d$rain_sauk <- ifelse(is.na(K_long_d$rain_s), lead(K_long_d$precip_sauk, 2)*0.0254, K_long_d$rain_s)

K_long_d$ETs <- LE.to.ET(K_long_d$LE_s, K_long_d$Tair_m)*60*30 # converted to mm/day kg/m2/s is equivalent to mm/s
K_long_d$ETm <- LE.to.ET(K_long_d$LE_m, K_long_d$Tair_m)*60*30 # converted to mm/day

K_long_d$WUE <- K_long_d$GEE_s/K_long_d$ETs
K_long_d$WUEm <- K_long_d$GEE_s/K_long_d$ET_s

K_long_d$WUE <- ifelse(K_long_d$WUE>10, NA, K_long_d$WUE)
K_long_d$WUE <- ifelse(K_long_d$WUE<(-3), NA, K_long_d$WUE)
K_long_d$WUE <- ifelse(!is.finite(K_long_d$WUE), NA, K_long_d$WUE)

K_long_d$WUEm <- ifelse(K_long_d$WUEm>2, NA, K_long_d$WUEm)
K_long_d$WUEm <- ifelse(K_long_d$WUEm<(-1), NA, K_long_d$WUEm)
K_long_d$WUEm <- ifelse(!is.finite(K_long_d$WUEm), NA, K_long_d$WUEm)

K_long_d = K_long_d %>% group_by(Type) %>% arrange(Year, DoY) %>% mutate(cs = cumsum(NEE_s))
K_long_d = K_long_d %>% group_by(Type) %>% arrange(Year, DoY) %>% mutate(csET = cumsum(ETs))
K_long_d = K_long_d %>% group_by(Type) %>% arrange(Year, DoY) %>% mutate(csGEE = cumsum(GEE_s))
K_long_d = K_long_d %>% group_by(Type) %>% arrange(Year, DoY) %>% mutate(csReco = cumsum(Reco_s))

K_long_d$LUE <- ((K_long_d$GEE_s/K_long_d$PPFD_s*86400/1000000))*100  # GPP = g C m-2 d-1 / PPFD mol m-2 s-1 

constants = bigleaf::bigleaf.constants()
Tair <- K_long_d$Tair_m
P <- K_long_d$P_m/1000
VPD <- K_long_d$VPD_m
ET <- K_long_d$ETs

K_long_d$y <- bigleaf::psychrometric.constant(Tair, P)
K_long_d$d <- bigleaf::air.density(Tair, P)
gs <- bigleaf::surface.conductance(K_long_d, K_long_d$Tair_m, K_long_d$P_m/1000, VPD=K_long_d$VPD_m, LE=K_long_d$LE_m, formulation = "Flux-Gradient")
K_long_d$gs <- gs$Gs_mol

##############################
K_long_a <- K_long_d %>%
  dplyr::group_by(Year, Type) %>%
  dplyr::summarize(NEE_a = sum(NEE_s, na.rm = T),
                   NEE_a05 = sum(NEE_s05, na.rm = T),
                   NEE_a50 = sum(NEE_s50, na.rm = T),
                   NEE_a95 = sum(NEE_s95, na.rm = T),
                   GEE_a = sum(GEE_s, na.rm = T),
                   LE_a = sum(LE_m, na.rm = T),
                   H_a = sum(H_m, na.rm = T),
                   Reco_a = sum(Reco_s, na.rm = T),
                   ET_r = sum(ET_m, na.rm = T),
                   LUE_a= mean(LUE, na.rm = T),
                   WUE_a = mean(WUE, na.rm = T),
                   Rain_sauk = sum(rain_sauk, na.rm = T)*1000)

########### Plots of carbon fluxes ############
K_long_m <- K_long_d %>%
  dplyr::group_by(Year, Month, Type) %>%
  dplyr::summarize(NEE_a = sum(NEE_s, na.rm = T),
                   ET_a = sum(ETs, na.rm = T),
                   ET_a = mean(ETs, na.rm = T),
                   NEE_mean = mean(NEE_s, na.rm = T),
                   NEE_med = median(NEE_s, na.rm = T),
                   NEE_sd = sd(NEE_s, na.rm = T),
                   GEE_a = sum(GEE_s, na.rm = T),
                   GEE_mean = mean(GEE_s, na.rm = T),
                   GEE_sd = sd(GEE_s, na.rm = T),
                   Reco_a = sum(Reco_s, na.rm = T),
                   Reco_mean = mean(Reco_s, na.rm = T),
                   Reco_sd = sd(Reco_s, na.rm = T),
                   Rn_a = mean(Rn_m, na.rm = T),
                   Rn_mean = mean(Rn_s, na.rm = T),
                   Rn_s = sum(Rn_s, na.rm = T),
                   Rn_sd = sd(Rn_s, na.rm = T),
                   LE_a = mean(LE_m, na.rm = T),
                   LE_mean = mean(LE_s, na.rm = T),
                   LE_s = sum(LE_m, na.rm = T),
                   LE_sd = sd(LE_m, na.rm = T),
                   G_a = mean(G_m, na.rm = T),
                   G_s = sum(G_m, na.rm = T),
                   G_mean = mean(G_s, na.rm = T),
                   G_sd = sd(G_m, na.rm = T),
                   H_a = mean(H_m, na.rm = T),
                   H_mean = mean(H_s, na.rm = T),
                   H_s = sum(H_m, na.rm = T),
                   H_sd = sd(H_m, na.rm = T),
                   Tair_a = mean(Tair_m, na.rm = T),
                   Tair_sd = sd(Tair_m, na.rm = T),
                   VPD_a = mean(VPD_m, na.rm = T),
                   VPD_sd = sd(VPD_m, na.rm = T),
                   PPFD_a = mean(PPFD_m, na.rm = T),
                   PPFD_sd = sd(PPFD_m, na.rm = T),
                   Rain = sum(rain_sauk, na.rm = T),
                   VPD = mean(VPD_m, na.rm = T),
                   WUE = mean(WUE, na.rm = T),
                   LUE_m = mean(LUE, na.rm = T),
                   gs_m = mean(gs, na.rm = T))

K_long_C <- gather(K_long_m, Flux, C_flux, NEE_a, GEE_a, Reco_a,  factor_key=TRUE)
K_long_C$Flux <- factor(K_long_C$Flux, levels = c("Reco_a", "GEE_a", "NEE_a"))
K_long_C$Date <- as.POSIXct(as.Date(paste0(K_long_C$Year, "-", K_long_C$Month, "-1"), "%Y-%m-%d"))
K_long_C$Type <- factor(K_long_C$Type, levels = c("IWGm", "IWGb"))

K_p_C <- ggplot(subset(K_long_C, Year%in%c(2019, 2020)), aes(x=Date, y=C_flux, group=Flux)) + facet_grid(Type~.)+
   geom_col(data=subset(K_long_C, Year%in%c(2019, 2020)), mapping=aes(x=Date, y=C_flux, colour=NA, fill=Flux),
           position = position_dodge2(width = 0.9, preserve = "single")) +
  ylab(expression(C~Flux~(g~C~m^-2~s^-1))) + xlab("Month") + 
  theme(legend.position="top",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 10),
        axis.title=element_text(size=16),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=15, colour="black"),
        legend.title=element_blank())+
  scale_color_manual(values = c("grey23", "#66a182", "#edae49")) +
  scale_fill_manual(values = c("grey23", "#66a182", "#edae49")) +
  scale_x_datetime(date_labels = "%b-%Y", date_breaks = "4 month") 

########### Plots of energy fluxes ############
K_long_E <- gather(K_long_m, Flux, E_flux, Rn_a, LE_a, H_a, G_a, factor_key=TRUE)
K_long_E$Flux <- factor(K_long_E$Flux, levels = c("Rn_a", "LE_a", "H_a", "G_a"))
K_long_E$Type <- factor(K_long_E$Type, levels = c("IWGm", "IWGb"))
K_long_E$Date <- as.POSIXct(as.Date(paste0(K_long_E$Year, "-", K_long_E$Month, "-1"), "%Y-%m-%d"))

K_p_E <- ggplot(subset(K_long_E, Year%in%c(2019, 2020)), aes(x=Date, y=E_flux, group=Flux)) + facet_grid(Type~.)+
  geom_col(data=subset(K_long_E, Year%in%c(2019, 2020)), mapping=aes(x=Date, y=E_flux, colour=NA, fill=Flux),
           position = position_dodge2(width = 0.9, preserve = "single")) +
  ylab(expression(E~Flux~(W~m^-2))) + xlab("Month") + 
  theme(legend.position="top",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 10),
        axis.title=element_text(size=16),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=15, colour="black"),
        legend.title=element_blank())+
   scale_color_manual(values = c("grey40",  "steelblue3", "tomato1","#DFCC74")) +
   scale_fill_manual(values = c("grey40",  "steelblue3", "tomato1","#DFCC74")) +
  scale_x_datetime(date_labels = "%b-%Y", date_breaks = "4 month") 

Figure3 <- cowplot::plot_grid(K_p_C,  K_p_E, ncol = 1, 
                           nrow = 2, 
                           labels = c("a", "b", "c", "d", "e"), align = "hv")
#print
#cowplot::save_plot(path.expand("/Users/susannewiesner/Documents/Jonescenter/Energy_density_disturbance/plots/Figure3.pdf"), 
#                   Figure3, base_height = 12, base_width = 8)

######### energy balance closure #################
K_long_E <- subset(K_long, Hour>9)
K_long_E <- subset(K_long_E, Hour<17)
K_long_E$Rn_G <- rowSums(cbind(K_long_E$Rn_Fill, K_long_E$G), na.rm=TRUE)
K_long_E$Type <- factor(K_long_E$Type, levels = c("IWGm", "IWGb"))
K_long_E$Year <- as.factor(K_long_E$Year)

K_p_E2 <- ggplot(subset(K_long_E, Year%in%c(2019, 2020)), aes(x=H_F+LE_F, y=Rn_Fill, group=Year)) + facet_grid(.~Type)+
  geom_point(data=subset(K_long_E, Year%in%c(2019, 2020)), mapping=aes(x=H_F+LE_F, y=Rn_Fill, colour=Year, fill=Year))+
  ylab(expression(R[n]-G~(W~m^-2))) + xlab(expression(H+LE~(W~m^-2))) + #xlim(-30, 130)+
  theme(legend.position="top",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 10),
        axis.title=element_text(size=16),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=15, colour="black"),
        legend.title=element_blank())+
  scale_color_manual(values = c("black", "skyblue2")) +
  scale_fill_manual(values = c("black", "skyblue2")) 

#print
#cowplot::save_plot(path.expand("/Users/susannewiesner/Documents/Jonescenter/Energy_density_disturbance/plots/FigureA2.pdf"), 
#                   K_p_E2, base_height = 5, base_width = 8)

########### Plots of indices ############
K_long_m2 <- K_long_d %>%
  dplyr::group_by(Year, Month, Type) %>%
  dplyr::summarize(WUE_m = mean(WUE, na.rm = T),
                   WUE_sd = sd(WUE, na.rm = T),
                   LUE_m = mean(LUE, na.rm = T),
                   LUE_sd = sd(LUE, na.rm = T),
                   gs_m = mean(gs, na.rm = T),
                   gs_sd = sd(gs, na.rm = T))
                   
K_long_m2$date <- as.Date(paste0(K_long_m2$Year, "/", K_long_m2$Month, "/1"), "%Y/%m/%d")
K_long_m2 <- subset(K_long_m2, Year %in% c("2019", "2020"))
K_long_m2 <- as.data.frame(K_long_m2)

K_p_WUE <- ggplot(K_long_m2, aes(x=date, y=WUE_m, group = Type)) +
  geom_errorbar(data = K_long_m2, aes(ymin = WUE_m-WUE_sd, ymax = WUE_m+WUE_sd), size=0.3, width=20, position = position_dodge(width = 10)) +
  geom_point(data=K_long_m2, mapping=aes(x=date, y=WUE_m, colour = Type), size=2, position = position_dodge(width = 10))+
  ylab(expression(WUE~(g~C~kg^-1~H2O))) + xlab(" ") +
   theme(legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        legend.title=element_blank()) +
    scale_color_manual(values = c("darkgoldenrod", "grey23"))

K_p_LUE <- ggplot(K_long_m2, aes(x=date, y=LUE_m, group = Type)) +
  geom_errorbar(data = K_long_m2, aes(ymin = LUE_m-LUE_sd, ymax = LUE_m+LUE_sd), size=0.3, width=20, position = position_dodge(width = 10)) +
  geom_point(data=K_long_m2, mapping=aes(x=date, y=LUE_m, colour = Type), size=2, position = position_dodge(width = 10))+
  ylab(expression(LUE~(g~C~mol^-1~PPFD))) + xlab(" ") +
  theme(legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        legend.title=element_blank()) +
    scale_color_manual(values = c("darkgoldenrod", "grey23"))


K_p_gs <- ggplot(K_long_m2, aes(x=date, y=gs_m, group = Type)) +
  geom_errorbar(data = K_long_m2, aes(ymin = gs_m-gs_sd, ymax = gs_m+gs_sd), size=0.3, width=20, position = position_dodge(width = 10)) +
  geom_point(data=K_long_m2, mapping=aes(x=date, y=gs_m, colour = Type), size=2, position = position_dodge(width = 10))+
  ylab(expression(gs~(mol~m^-2~s^-1))) + xlab(" ") +
  theme(legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        legend.title=element_blank()) +
    scale_color_manual(values = c("darkgoldenrod", "grey23"))

Figure5 <- cowplot::plot_grid(K_p_WUE,  K_p_LUE, K_p_gs, ncol = 3, 
                           nrow = 1, 
                           labels = c("a", "b", "c", "d", "e"), align = "hv")
#print
#cowplot::save_plot(path.expand("/Users/susannewiesner/Documents/Jonescenter/Energy_density_disturbance/plots/Figure5.pdf"), 
#                   Figure5, base_height = 5, base_width = 16)

########### Monthly environmental variables ############
K_long_plot <- K_long_d %>%
  dplyr::group_by(Year, Month, Type) %>%
  dplyr::summarize(Tair_a = mean(Tair_m, na.rm = T),
                   Tair_sd = sd(Tair_m, na.rm = T),
                   VPD_a = mean(VPD_m, na.rm = T),
                   VPD_sd = sd(VPD_m, na.rm = T),
                   Rain_sauk = sum(rain_sauk, na.rm = T)*1000,
                   P_sauk = sum(precip_sauk*0.0254, na.rm = T)*1000,
                   PPFD_a = mean(PPFD_m, na.rm = T),
                   PPFD_sd = sd(PPFD_m, na.rm = T))

K_long_plot$date <- as.Date(paste0(K_long_plot$Year, "/", K_long_plot$Month, "/1"), "%Y/%m/%d")
K_long_plot <- subset(K_long_plot, Type == "IWGm")
K_long_plot <- subset(K_long_plot, Year %in% c("2019", "2020"))
K_long_plot <- as.data.frame(K_long_plot)

########### Plots of environmental variables ############
K_p_Tair <- ggplot(K_long_plot, aes(x=date, y=Tair_a)) +
  geom_errorbar(data = K_long_plot, aes(ymin = Tair_a-Tair_sd, ymax = Tair_a+Tair_sd), size=.3, width=12) +
  geom_point(data=K_long_plot, mapping=aes(x=date, y=Tair_a), size=2,color = "brown4")+
  ylab(expression(T[air]~(degree*C))) + xlab(" ") +
  theme(legend.position="top",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 10),
        axis.title=element_text(size=16),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=15, colour="black"),
        legend.title=element_blank())

K_p_VPD <- ggplot(K_long_plot, aes(x=date, y=VPD_a)) +
  geom_errorbar(data = K_long_plot, aes(ymin = VPD_a-VPD_sd, ymax = VPD_a+VPD_sd), size=0.3, width=12) +
  geom_point(data=K_long_plot, mapping=aes(x=date, y=VPD_a), size=2, colour = "azure4")+
  ylab(expression(VPD~(kPa))) + xlab(" ") +
  theme(legend.position="top",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 10),
        axis.title=element_text(size=16),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=15, colour="black"),
        legend.title=element_blank())

K_p_PPFD <- ggplot(K_long_plot, aes(x=date, y=PPFD_a)) +
  geom_errorbar(data = K_long_plot, aes(ymin = PPFD_a-PPFD_sd, ymax = PPFD_a+PPFD_sd), size=0.3, width=12) +
  geom_point(data=K_long_plot, mapping=aes(x=date, y=PPFD_a), size=2, colour = "darkseagreen4")+
  ylab(expression(PPFD~(mu*mol~m^-2~s^-1))) + xlab(" ") + 
  theme(legend.position="top",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 10),
        axis.title=element_text(size=16),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=15, colour="black"),
        legend.title=element_blank())

K_p_Rain <- ggplot(K_long_plot, aes(x=date, y=Rain_sauk)) +
  geom_bar(data=K_long_plot, mapping=aes(x=date, y=Rain_sauk), stat = "identity", size=1, fill = "deepskyblue4")+
  ylab(expression(Rain~(mm))) + xlab(" ") +
  theme(legend.position="top",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 10),
        axis.title=element_text(size=16),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=15, colour="black"),
        legend.title=element_blank())

#print
Figure2 <- cowplot::plot_grid(K_p_Rain, K_p_Tair, K_p_VPD, K_p_PPFD, ncol = 2, 
                           nrow = 2, 
                           labels = c("a", "b", "c", "d"), align = "hv")

#cowplot::save_plot(path.expand("/Users/susannewiesner/Downloads/AFM/plots/final_plots/Figure1.pdf"), 
#                   Figure2, base_height = 10, base_width = 10)

######### Plots of cumulative variables ############
K_p_csNEE <- ggplot(K_long_d, aes(x=date, y=cs, group=Type)) +
  geom_point(data=K_long_d, mapping=aes(x=date, y=cs, colour=Type), size=1)+
  ylab(expression(NEE~(g~C~m^-2))) + xlab(" ") +
  theme(legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 10),
        axis.title=element_text(size=24),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        legend.title=element_blank())+
  scale_color_manual(values = c("darkgoldenrod", "grey23"))

K_p_csET <- ggplot(K_long_d, aes(x=date, y=csET, group=Type)) +
  geom_point(data=K_long_d, mapping=aes(x=date, y=csET, colour=Type), size=1)+
  ylab(expression(ET~(kg~H2O~m^-2))) + xlab(" ") +
  theme(legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 10),
        axis.title=element_text(size=24),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        legend.title=element_blank())+
  scale_color_manual(values = c("darkgoldenrod", "grey23"))

K_p_csGEE <- ggplot(K_long_d, aes(x=date, y=csGEE, group=Type)) +
  geom_point(data=K_long_d, mapping=aes(x=date, y=csGEE,  colour=Type), size=1)+
  ylab(expression(GPP~(g~C~m^-2))) + xlab(" ") + ylim(-10, 4000) +
  theme(legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 10),
        axis.title=element_text(size=24),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        legend.title=element_blank())+
  scale_color_manual(values = c("darkgoldenrod", "grey23"))

K_p_csReco <- ggplot(K_long_d, aes(x=date, y=csReco, group=Type)) +
  geom_point(data=K_long_d, mapping=aes(x=date, y=csReco, colour=Type), size=1)+
  ylab(expression(R[eco]~(g~C~m^-2))) + xlab(" ") + ylim(-10, 4000) +
  theme(legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=24),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        legend.title=element_blank())+
  scale_color_manual(values = c("darkgoldenrod", "grey23"))

#print
Figure4 <- cowplot::plot_grid(K_p_csNEE, K_p_csGEE, K_p_csReco, K_p_csET, ncol = 2, 
                           nrow = 2, 
                           labels = c("a", "b", "c", "d", "e"), align = "hv")

#cowplot::save_plot(path.expand("/Users/susannewiesner/Downloads/AFM/plots/Figure4.pdf"), 
#                   Figure4, base_height = 11, base_width = 12)

########### ET ############
K_long_ET <- K_long_d %>%
  dplyr::group_by(Year, Month, Type) %>%
  dplyr::summarize(ET_mean = mean(ETs, na.rm = T),
                   ET_sd = sd(ETs, na.rm = T),
                   ET_r = mean(ET_m, na.rm = T),
                   ET_r_sd = sd(ET_m, na.rm = T),
                   ET_r2 = mean(ET_s, na.rm = T),
                   ET_r2_sd = sd(ET_s, na.rm = T))
                   
K_long_ET$date <- as.Date(paste0(K_long_ET$Year, "/", K_long_ET$Month, "/1"), "%Y/%m/%d")
K_long_ET <- subset(K_long_ET, Year %in% c("2019", "2020"))
K_long_ET <- as.data.frame(K_long_ET)

K_p_ET <- ggplot(K_long_ET, aes(x=date, y=ET_mean, group = Type)) + 
  geom_errorbar(data = K_long_ET, aes(ymin = ET_mean-ET_sd, ymax = ET_mean+ET_sd), size=0.3, width=20, position = position_dodge(width = 10)) +
  geom_point(data=K_long_ET, mapping=aes(x=date, y=ET_mean, colour = Type), size=2, position = position_dodge(width = 10))+
  geom_line(data=K_long_ET, mapping=aes(x=date, y=ET_mean, colour = Type), size=1, position = position_dodge(width = 10))+
  ylab(expression(ET~(g~C~kg^-1~H2O))) + xlab(" ") +
   theme(legend.position="top",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        legend.title=element_blank()) +
    scale_color_manual(values = c("darkgoldenrod", "grey23"))

#cowplot::save_plot(path.expand("/Users/susannewiesner/Documents/AFM/plots/FigureA3.pdf"), 
#                   K_p_ET, base_height = 6, base_width = 8)


########################## +++++++++++ monthly Fluxes +++++++++++ #####################
############ Tair climatologie #############
K_long_T <- subset(K_long_d, Year == c("2019"))
K_long_T <- subset(K_long_d, Year == c("2020"))

lme_T <- lm(TairNA ~ Month*Type, data=K_long_T, na.action = na.omit) # by year separately

lme_T <- lm(TairNA ~ Month*Type+Type*Year, data=K_long_d, na.action = na.omit) # by year and month

car::Anova(lme_T) 
anova(lme_T) 
AIC(lme_T)
model_performance(lme_T)

lsm <- lsmeans(lme_T, ~ Type:Month)
 pairs(lsm, by = "Month", adjust = "mvt")

plot_model(lme_T, type = "pred", terms = c("Month", "Type"))


############ VPD climatologie #############
K_long_V <- subset(K_long_d, Year == c("2019"))
K_long_V <- subset(K_long_d, Year == c("2020"))

lme_V <- lm(VPDNA ~ Month*Type, data=K_long_V, na.action = na.omit) # by year separately

lme_V <- lm(VPDNA ~ Month*Type+Type*Year, data=K_long_d, na.action = na.omit) # by year and month

car::Anova(lme_V)
anova(lme_V)
AIC(lme_V)
model_performance(lme_V)

lsm <- lsmeans(lme_V, ~ Type:Month)
 pairs(lsm, by = "Month", adjust = "Tukey")

 plot_model(lme_V, type = "pred", terms = c("Month", "Type"))

############ Tair #############
lm_T <- lm(Tair_m ~ Month*Year, data=K_long_d, na.action = na.omit)
car::Anova(lm_T) 
AIC(lm_T)
model_performance(lm_T)

lsm <- lsmeans(lm_T, ~ Year:Month)
 pairs(lsm, by = "Month", adjust = "mvt")
 pairs(lsm, by = "Year", adjust = "mvt")
 
plot_model(lm_T, type = "pred", terms = c("Month", "Year"))

############ VPD #############
lm_V <- lm(VPD_m ~ Month*Year, data=K_long_d, na.action = na.omit)
car::Anova(lm_V) 
AIC(lm_V)
model_performance(lm_V)

lsm <- lsmeans(lm_V, ~ Year:Month)
 pairs(lsm, by = "Month", adjust = "mvt")
 pairs(lsm, by = "Year", adjust = "mvt")
 
plot_model(lm_V, type = "pred", terms = c("Month", "Year"))

############ PPFD #############
lm_P <- lm(PPFD_m ~ Month*Year, data=K_long_d, na.action = na.omit)
car::Anova(lm_P)
AIC(lm_P)
model_performance(lm_P)

lsm <- lsmeans(lm_P, ~ Year:Month)
 pairs(lsm, by = "Month", adjust = "mvt")
 pairs(lsm, by = "Year", adjust = "mvt")
 
plot_model(lm_P, type = "pred", terms = c("Month", "Year"))

############ NEE #############
K_long_N <- subset(K_long_d, Year == c("2019"))
K_long_N <- subset(K_long_d, Year == c("2020"))

lme_NEE <- lm(NEE_s ~ Month*Type, data=K_long_N, na.action = na.omit)
car::Anova(lme_NEE)
AIC(lme_NEE)
model_performance(lme_NEE)

lsm <- lsmeans(lme_NEE, ~ Type:Month)
 pairs(lsm, by = "Month", adjust = "mvt")
plot_model(lme_NEE, type = "pred", terms = c("Month", "Type"))

############ GEE #############
K_long_G <- subset(K_long_d, Year == c("2019"))
K_long_G <- subset(K_long_d, Year == c("2020"))

lme_GEE <- lm(GEE_s ~ Month*Type, data=K_long_G, na.action = na.omit)
car::Anova(lme_GEE)
AIC(lme_GEE)
model_performance(lme_GEE)

lsm <- lsmeans(lme_GEE, ~ Type:Month)
 pairs(lsm, by = "Month", adjust = "mvt")

plot_model(lme_GEE, type = "pred", terms = c("Month", "Type"))

############ Reco #############
K_long_R <- subset(K_long_d, Year == c("2019"))
K_long_R <- subset(K_long_d, Year == c("2020"))

lme_Reco <- lm(Reco_s ~ Month*Type, data=K_long_R, na.action = na.omit)
car::Anova(lme_Reco)
AIC(lme_Reco)
model_performance(lme_Reco)

lsm <- lsmeans(lme_Reco, ~ Type:Month)
 pairs(lsm, by = "Month", adjust = "mvt")

plot_model(lme_Reco, type = "pred", terms = c("Month", "Type"))

############ LE #############
K_long_LE <- subset(K_long_d, Year == c("2019"))
K_long_LE <- subset(K_long_d, Year == c("2020"))

lme_LE <- lm(LE_m ~ Month*Type, data=K_long_LE, na.action = na.omit) 
lme_LE <- lme(LE_m ~ Month*Type, data=K_long_LE, random = (~1 | DoY), na.action = na.omit) 

car::Anova(lme_LE) 
AIC(lme_LE)
model_performance(lme_LE)

lsm <- lsmeans(lme_LE, ~ Type:Month)
 pairs(lsm, by = "Month", adjust = "mvt")

plot_model(lme_LE, type = "pred", terms = c("Month", "Type"))

############ H #############
K_long_H <- subset(K_long_d, Year == c("2019"))
K_long_H <- subset(K_long_d, Year == c("2020"))

lme_H <- lme(H_m ~ Month*Type, data=K_long_H, random = (~1 | DoY), na.action = na.omit)
car::Anova(lme_H)
AIC(lme_H)
model_performance(lme_H)

lsm <- lsmeans(lme_H, ~ Type:Month)
 pairs(lsm, by = "Month", adjust = "mvt")

plot_model(lme_H, type = "pred", terms = c("Month", "Type"))

############ G #############
K_long_G <- subset(K_long_d, Year == c("2019"))
K_long_G <- subset(K_long_d, Year == c("2020"))

lme_G <- lm(G_m ~ Month*Type, data=K_long_H, random = (~1 | DoY), na.action = na.omit)
car::Anova(lme_G)
AIC(lme_G)
model_performance(lme_G)

lsm <- lsmeans(lme_G, ~ Type:Month)
 pairs(lsm, by = "Month", adjust = "mvt")

plot_model(lme_G, type = "pred", terms = c("Month", "Type"))

############ ET #############
K_long_ET <- subset(K_long_d, Year == c("2019"))
K_long_ET <- subset(K_long_d, Year == c("2020"))

lme_ET <- lm(ET_s ~ Month*Type, data=K_long_ET, na.action = na.omit) 
car::Anova(lme_ET)
AIC(lme_ET)
model_performance(lme_ET)

lsm <- lsmeans(lme_ET, ~ Type:Month)
 pairs(lsm, by = "Month", adjust = "Type")

plot_model(lme_ET, type = "pred", terms = c("Month", "Type"))

############ WUE #############
K_long_WUE <- subset(K_long_d, Year == c("2019"))
K_long_WUE <- subset(K_long_d, Year == c("2020"))

lme_WUE <- lm(WUE ~ Month*Type, data=K_long_WUE, na.action = na.omit)
car::Anova(lme_WUE)
AIC(lme_WUE)
model_performance(lme_WUE)

lsm <- lsmeans(lme_WUE, ~ Type:Month)
 pairs(lsm, by = "Month", adjust = "mvt")

plot_model(lme_WUE, type = "pred", terms = c("Month", "Type"))

############ LUE #############
K_long_LUE <- subset(K_long_d, Year == c("2020"))

lme_LUE <- lm(LUE ~ Month*Type, data=K_long_LUE, na.action = na.omit) 
car::Anova(lme_LUE) 
AIC(lme_LUE)
model_performance(lme_LUE)

lsm <- lsmeans(lme_LUE, ~ Type:Month)
 pairs(lsm, by = "Month", adjust = "mvt")

plot_model(lme_LUE, type = "pred", terms = c("Month", "Type"))

############ gs #############
K_long_gs <- subset(K_long_d, Year == c("2019"))

lme_gs <- lm(gs ~ Month*Type, data=K_long_gs, na.action = na.omit)
car::Anova(lme_gs)
AIC(lme_gs)
model_performance(lme_gs)

lsm <- lsmeans(lme_gs, ~ Type:Month)
 pairs(lsm, by = "Month", adjust = "mvt")

plot_model(lme_gs, type = "pred", terms = c("Month", "Type"))


######################## ++++++++++++++++++ Uncertainty estimate ++++++++++++++++++ ##########################
############### Kernza monoculture IWGb #####################
KC_long$Time <- paste(floor(KC_long$Hour), round((KC_long$Hour-floor(KC_long$Hour))*60), sep=":")
KC_long$DateTime <- as.POSIXct(paste0(KC_long$Year, KC_long$DoY, " ", KC_long$Time), format="%Y%j %H:%M", tz="UTC")

KC_long %>% group_by(Year) %>% dplyr::filter(NEE_uStar_fqc == 0) %>% dplyr::summarise(nRec = sum(is.finite(NEE_uStar_f))
  , varSum = sum(NEE_uStar_fsd^2, na.rm = TRUE)
  , seMean = sqrt(varSum) / nRec
  , seMeanApprox = mean(NEE_uStar_fsd, na.rma = TRUE) / sqrt(nRec)
  ) %>% select(nRec, seMean, seMeanApprox)

results <- KC_long %>% 
  mutate(resid = ifelse(NEE_uStar_fqc == 0, NEE_uStar_orig - NEE_uStar_fall, NA ))

autoCorr <- computeEffectiveAutoCorr(results$resid)
nEff <- computeEffectiveNumObs(results$resid, na.rm = TRUE)
c( nEff = nEff, nObs = sum(is.finite(results$resid)))

acf(results$resid, na.action = na.pass, main = "")

results <- results %>% mutate(
  DateTime = KC_long$DateTime
  , DoY = as.POSIXlt(DateTime - 15*60)$yday # midnight belongs to the previous
)

results %>% filter(NEE_uStar_fqc == 0) %>% summarise(
  nRec = sum(is.finite(NEE_uStar_fsd))
  , varMean = sum(NEE_uStar_fsd^2, na.rm = TRUE) / nRec / (!!nEff - 1)
  , seMean = sqrt(varMean) 
  , seMeanApprox = mean(NEE_uStar_fsd, na.rm = TRUE) / sqrt(!!nEff - 1)
  ) %>% select(seMean, seMeanApprox)

results %>% group_by(Year, DoY) %>% filter(NEE_uStar_fqc == 0) %>% summarise(
  nRec = sum(is.finite(NEE_uStar_fsd))
  , varMean = sum(NEE_uStar_fsd^2, na.rm = TRUE) / nRec / (!!nEff - 1)
  , seMean = sqrt(varMean) 
  , seMeanApprox = mean(NEE_uStar_fsd, na.rm = TRUE) / sqrt(!!nEff - 1)
  ) %>% select(seMean, seMeanApprox)

plot(NEE_uStar_fsd ~ NEE_uStar_fall, slice(results, sample.int(nrow(results),400)))

aggDay <- results %>% group_by(DoY, Year) %>% 
  summarise(
    DateTime = first(DateTime)
    , nRec = sum( NEE_uStar_fqc == 0, na.rm = TRUE)
    , nEff = computeEffectiveNumObs(
       resid, effAcf = !!autoCorr, na.rm = TRUE)
    , NEE = mean(NEE_uStar_f, na.rm = TRUE)
    , sdNEE = if (nEff <= 1) NA_real_ else sqrt(
      mean(NEE_uStar_fsd^2, na.rm = TRUE) / (nEff - 1)) 
    , sdNEEuncorr = if (nRec == 0) NA_real_ else sqrt(
       mean(NEE_uStar_fsd^2, na.rm = TRUE) / (nRec - 1))
  )

aggMonth <- results %>% group_by(Year, Month) %>% 
  summarise(
    DateTime = first(DateTime)
    , nRec = sum( NEE_uStar_fqc == 0, na.rm = TRUE)
    , nEff = computeEffectiveNumObs(
       resid, effAcf = !!autoCorr, na.rm = TRUE)
    , NEE = mean(NEE_uStar_f, na.rm = TRUE)
    , sdNEE = if (nEff <= 1) NA_real_ else sqrt(
      mean(NEE_uStar_fsd^2, na.rm = TRUE) / (nEff - 1)) 
    , sdNEEuncorr = if (nRec == 0) NA_real_ else sqrt(
       mean(NEE_uStar_fsd^2, na.rm = TRUE) / (nRec - 1))
  )

aggYear <- results %>% group_by(Year) %>% 
  summarise(
    DateTime = first(DateTime)
    , nRec = sum( NEE_uStar_fqc == 0, na.rm = TRUE)
    , nEff = computeEffectiveNumObs(
       resid, effAcf = !!autoCorr, na.rm = TRUE)
    , NEE = mean(NEE_uStar_f, na.rm = TRUE)
    , sdNEE = if (nEff <= 1) NA_real_ else sqrt(
      mean(NEE_uStar_fsd^2, na.rm = TRUE) / (nEff - 1)) 
    , sdNEEuncorr = if (nRec == 0) NA_real_ else sqrt(
       mean(NEE_uStar_fsd^2, na.rm = TRUE) / (nRec - 1))
  )

aggDay$DateTime <- as.POSIXct(aggDay$DateTime, "%Y-%m-%d %H:%M:%S")
aggDay$Type <- "IWGb"

aggMonth$DateTime <- as.POSIXct(aggMonth$DateTime, "%Y-%m-%d %H:%M:%S")
aggMonth$Type <- "IWGb"

############### Kernza monoculture IWGm #####################
KO_long$Time <- paste(floor(KO_long$Hour), round((KO_long$Hour-floor(KO_long$Hour))*60), sep=":")
KO_long$DateTime <- as.POSIXct(paste0(KO_long$Year, KO_long$DoY, " ", KO_long$Time), format="%Y%j %H:%M", tz="UTC")

KO_long %>% filter(NEE_uStar_fqc == 0) %>% summarise(
  nRec = sum(is.finite(NEE_uStar_f))
  , varSum = sum(NEE_uStar_fsd^2, na.rm = TRUE)
  , seMean = sqrt(varSum) / nRec
  , seMeanApprox = mean(NEE_uStar_fsd, na.rma = TRUE) / sqrt(nRec)
  ) %>% select(nRec, seMean, seMeanApprox)

results2 <- KO_long %>% 
  mutate(
    resid = ifelse(NEE_uStar_fqc == 0, NEE_uStar_orig - NEE_uStar_fall, NA )
  )

results2 <- results2 %>% mutate(
  DateTime = KO_long$DateTime
  , DoY = as.POSIXlt(DateTime - 15*60)$yday # midnight belongs to the previous
)

autoCorr2 <- computeEffectiveAutoCorr(results2$resid)
nEff2 <- computeEffectiveNumObs(results2$resid, na.rm = TRUE)
c( nEff2 = nEff2, nObs = sum(is.finite(results2$resid)))

### uncertainty for effective number of observations 
results2 %>% filter(NEE_uStar_fqc == 0) %>% summarise(
  nRec = sum(is.finite(NEE_uStar_fsd))
  , varMean = sum(NEE_uStar_fsd^2, na.rm = TRUE) / nRec / (!!nEff - 1)
  , seMean = sqrt(varMean) 
  , seMeanApprox = mean(NEE_uStar_fsd, na.rm = TRUE) / sqrt(!!nEff - 1)
  ) %>% select(seMean, seMeanApprox)

results2 %>% group_by(Year) %>% filter(NEE_uStar_fqc == 0) %>% summarise(
  nRec = sum(is.finite(NEE_uStar_fsd))
  , varMean = sum(NEE_uStar_fsd^2, na.rm = TRUE) / nRec / (!!nEff - 1)
  , seMean = sqrt(varMean) 
  , seMeanApprox = mean(NEE_uStar_fsd, na.rm = TRUE) / sqrt(!!nEff - 1)
  ) %>% select(seMean, seMeanApprox)

aggDay2 <- results2 %>% group_by(DoY, Year) %>% 
  summarise(
    DateTime = first(DateTime)
    , nRec = sum( NEE_uStar_fqc == 0, na.rm = TRUE)
    , nEff = computeEffectiveNumObs(
       resid, effAcf = !!autoCorr, na.rm = TRUE)
    , NEE = mean(NEE_uStar_f, na.rm = TRUE)
    , sdNEE = if (nEff <= 1) NA_real_ else sqrt(
      mean(NEE_uStar_fsd^2, na.rm = TRUE) / (nEff - 1)) 
    , sdNEEuncorr = if (nRec == 0) NA_real_ else sqrt(
       mean(NEE_uStar_fsd^2, na.rm = TRUE) / (nRec - 1))
  )

aggMonth2 <- results2 %>% group_by(Year, Month) %>% 
  summarise(
    DateTime = first(DateTime)
    , nRec = sum( NEE_uStar_fqc == 0, na.rm = TRUE)
    , nEff = computeEffectiveNumObs(
       resid, effAcf = !!autoCorr, na.rm = TRUE)
    , NEE = mean(NEE_uStar_f, na.rm = TRUE)
    , sdNEE = if (nEff <= 1) NA_real_ else sqrt(
      mean(NEE_uStar_fsd^2, na.rm = TRUE) / (nEff - 1)) 
    , sdNEEuncorr = if (nRec == 0) NA_real_ else sqrt(
       mean(NEE_uStar_fsd^2, na.rm = TRUE) / (nRec - 1))
  )

aggYear2 <- results2 %>% group_by(Year) %>% 
  summarise(
    DateTime = first(DateTime)
    , nRec = sum( NEE_uStar_fqc == 0, na.rm = TRUE)
    , nEff = computeEffectiveNumObs(
       resid, effAcf = !!autoCorr, na.rm = TRUE)
    , NEE = mean(NEE_uStar_f, na.rm = TRUE)
    , sdNEE = if (nEff <= 1) NA_real_ else sqrt(
      mean(NEE_uStar_fsd^2, na.rm = TRUE) / (nEff - 1)) 
    , sdNEEuncorr = if (nRec == 0) NA_real_ else sqrt(
       mean(NEE_uStar_fsd^2, na.rm = TRUE) / (nRec - 1))
  )

aggDay2$DateTime <- as.POSIXct(aggDay2$DateTime, "%Y-%m-%d %H:%M:%S")
aggDay2$Type <- "IWGm"

aggMonth2$DateTime <- as.POSIXct(aggMonth2$DateTime, "%Y-%m-%d %H:%M:%S")
aggMonth2$Type <- "IWGm"

aggD <- rbind(aggDay, aggDay2)
aggD$NEE_perc_er <- aggD$sdNEE/aggD$NEE*100

aggM <- rbind(aggMonth, aggMonth2)

ggplot(aggD, aes(x=DateTime, y=NEE, group = Type)) +
  geom_point(data=aggD, mapping=aes(x=DateTime, y=NEE, colour = Type), size=1)+
  geom_errorbar(data = aggD, aes(ymin = NEE-sdNEE, ymax = NEE+sdNEE, colour =Type), size=0.3, width=0.4) +
  ylab(expression(NEE~(g~CO[2]~m^-2~s^-1))) + xlab("Date") + ylim (-20, 20) +
  theme(legend.position="top",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 10),
        axis.title=element_text(size=16),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=15, colour="black"),
        legend.title=element_blank())+
  scale_color_manual(values = c("darkgoldenrod", "grey46", "lightgray"))+
  scale_shape_manual(values = c(15, 16, 17))+
  scale_x_datetime(date_labels = "%b %d", date_breaks = "2 month")

themeTw <- theme_bw(base_size = 10) + theme(axis.title = element_text(size = 9))

ggplot(aggD, aes(DateTime, NEE, group = Type)) +
  geom_line(aes(colour = aggD$Type)) + ylim(-20, 20)+
  geom_ribbon(aes(ymin = NEE - 1.96*sdNEE, ymax = NEE + 1.96*sdNEE, fill = Type), alpha = 0.2) +
  scale_fill_discrete(labels = c('IWGb','IWGm')) +
  themeTw + theme(legend.position = c(0.95,0.95), legend.justification = c(1,1))

ggplot(aggM, aes(DateTime, NEE, group = Type)) +
  geom_line(aes(colour = Type)) + ylim(-20, 20)+
  geom_ribbon(aes(ymin = NEE - 1.96*sdNEE, ymax = NEE + 1.96*sdNEE, fill = Type), alpha = 0.2) +
  scale_fill_discrete(labels = c('IWGb','IWGm')) +
  themeTw + theme(legend.position = c(0.95,0.95), legend.justification = c(1,1))


NEE_agg <- read.csv("/Users/susannewiesner/Downloads/AFM/NEE_uncertainty_KC_KO.csv")
NEE_agg$DateTime <- as.Date(paste0(NEE_agg$Year, "-", NEE_agg$Month, "-1"), "%Y-%m-%d")
NEE_agg$lwr <- NEE_agg$NEE_f - NEE_agg$error
NEE_agg$upr <- NEE_agg$NEE_f + NEE_agg$error

p_uM <- ggplot(NEE_agg, aes(as.Date(DateTime), NEE_f, group = Type))+
    geom_point(aes(colour = Type))+
    geom_line(data=NEE_agg, aes(colour = Type))+
    geom_ribbon(data=NEE_agg, aes(ymin=lwr, ymax=upr, fill = Type), alpha=0.2, colour = NA)+
    ylab(expression(NEE~(mu*mol~CO["2"]~m^-2~s^-1))) + xlab(" ") +
   theme(legend.position="none",panel.grid.minor = element_blank(),
        panel.grid.major =  element_blank(),
        legend.text = element_text(size = 15),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = 'white', colour="black"),
        legend.key = element_rect(fill = NA, colour = NA),
        axis.text = element_text(size=17, colour="black"),
        legend.title=element_blank()) +
    scale_color_manual(values = c("darkgoldenrod", "grey23"))+
    scale_fill_manual(values = c("darkgoldenrod", "grey23"))

#cowplot::save_plot("/Users/susannewiesner/Downloads/AFM/FigureA1.pdf",  p_uM, base_height = 7, base_width = 12)


################# calculate Ustar threshold uncertainty ################
results <- rbind(KC_long, KO_long)

results <- results %>% 
  mutate(
    resid = ifelse(NEE_uStar_fqc == 0, NEE_uStar_orig - NEE_uStar_fall, NA )
  )

autoCorr <- computeEffectiveAutoCorr(results$resid)
nEff <- computeEffectiveNumObs(results$resid, na.rm = TRUE)
c( nEff = nEff, nObs = sum(is.finite(results$resid)))

results %>% filter(NEE_uStar_fqc == 0) %>% summarise(
  nRec = sum(is.finite(NEE_uStar_fsd))
  , varMean = sum(NEE_uStar_fsd^2, na.rm = TRUE) / nRec / (!!nEff - 1)
  , seMean = sqrt(varMean) 
  , seMeanApprox = mean(NEE_uStar_fsd, na.rm = TRUE) / sqrt(!!nEff - 1)
  ) %>% select(seMean, seMeanApprox)

#################################################################################
results <- results %>% mutate(
  DateTime = results$DateTime
  , DoY = as.POSIXlt(DateTime - 15*60)$yday # midnight belongs to the previous
)

results$DoY <- as.character(results$DoY)
results <- tidyr::separate(results, DateTime, into = c("Date","Time"), sep = " ", extra = "merge", remove=F)
results <- tidyr::separate(results, Date, into = c("Year","Month", "Day"), sep = "-", extra = "merge", remove=F)

########################### Systematic Uncertainty ############################3
uStarSuffixes <- c("uStar", "U05", "U50", "U95")

NEEAggCO2 <- results %>% dplyr::group_by(Year, Type) %>%
  dplyr::summarise(
    NEE_f = mean(NEE_uStar_f),
    NEE_f05 = mean(NEE_U05_f),
    NEE_f50 = mean(NEE_U50_f),
    NEE_f95 = mean(NEE_U95_f)
  )

NEEAggCO2$NEE_m <- rowMeans(cbind(NEEAggCO2$NEE_f, NEEAggCO2$NEE_f05, NEEAggCO2$NEE_f50, NEEAggCO2$NEE_f95), na.rm=TRUE) 
NEEAggCO2$NEE_sd <- rowSds(cbind(NEEAggCO2$NEE_f, NEEAggCO2$NEE_f05, NEEAggCO2$NEE_f50, NEEAggCO2$NEE_f95), na.rm=TRUE) 
NEEAggCO2$NEE_var <- rowVars(cbind(NEEAggCO2$NEE_f, NEEAggCO2$NEE_f05, NEEAggCO2$NEE_f50, NEEAggCO2$NEE_f95), na.rm=TRUE)

#write.csv(NEEAggCO2, "/Users/susannewiesner/Downloads/AFM/NEE_annual_uncertainty_IWGm_IWGb.csv")

################## +++++++++++++++++ biomass model +++++++++++++++++ ###################
K_bio <- read.csv("/Users/susannewiesner/Documents/Kernza/data/Kernza_biomass.csv")
K_bio$Type <- as.factor(K_bio$Type)

lm_b <- lm(Biomass.Dry.Weight..g.~Type, K_bio)
Anova(lm_b)
lsm <- emmeans(lm_b, ~ Type)
b.emm.s <- emmeans(lm_b, "Type")
pairs(b.emm.s)
contrast(b.emm.s, "poly")

lm_g <- lm(grain..g.~Type, K_bio)
Anova(lm_g)

lm_s <- lm(spikes..g.~Type, K_bio)
Anova(lm_s)

lm_h <- lm(HI~Type, K_bio)
Anova(lm_h)

lm_m <- lm(X..Moisture.Content~Type, K_bio)
Anova(lm_m)     
     
K_bio_a <- K_bio %>%
  dplyr::group_by(Type) %>%
  dplyr::summarize(Bio_a = mean(Biomass.Dry.Weight..g., na.rm = T)*4,
                   Grain_a = mean(grain..g., na.rm = T)*4,
                   Spikes_a = mean(spikes..g., na.rm = T)*4,
                   HI_a = mean(HI, na.rm = T),
                   Moist_a= mean(X..Moisture.Content, na.rm = T),
                   Bio_sd = sd(Biomass.Dry.Weight..g.*4, na.rm = T),
                   Grain_sd = sd(grain..g.*4, na.rm = T),
                   Spikes_sd = sd(spikes..g.*4, na.rm = T),
                   HI_sd = sd(HI, na.rm = T),
                   Moist_sd= sd(X..Moisture.Content, na.rm = T))

################## +++++++++++++++++ quadrat cover model +++++++++++++++++ ###################
K_quadrat <- read.csv("/Users/susannewiesner/Documents/Kernza/data/KC_KO_quadrat_final_results.csv")
K_quadrat$Op <- as.factor(K_quadrat$Op)
K_quadrat$Type <- as.factor(K_quadrat$Type)

lm_b <- lme(Percent.green~0+Type, random = (~1|Op), K_quadrat)
plot(lm_b)
summary(lm_b)
Anova(lm_b)

lsm <- emmeans(lm_b, ~ Type)
b.emm.s <- emmeans(lm_b, "Type")
pairs(b.emm.s)
contrast(b.emm.s, "poly")

lm_K <- lme(Kernza_P~0+Type, random = (~1|Op), K_quadrat)
plot(lm_K)
summary(lm_K)
Anova(lm_K)
b.emm.s <- emmeans(lm_K, "Type")
pairs(b.emm.s)

lm_D <- lme(Non_Leg_p~0+Type, random = (~1|Op), K_quadrat)
plot(lm_D)
summary(lm_D)
Anova(lm_D)
b.emm.s <- emmeans(lm_D, "Type")
pairs(b.emm.s)

lm_C <- lme(Legume_p~0+Type, random = (~1|Op), K_quadrat)
plot(lm_C)
summary(lm_C)
Anova(lm_C)
b.emm.s <- emmeans(lm_C, "Type")
pairs(b.emm.s)

lm_CNL <- lme(Legume_p+Non_Leg_p~0+Type, random = (~1|Op), K_quadrat)
plot(lm_CNL)
summary(lm_CNL)
Anova(lm_CNL)
b.emm.s <- emmeans(lm_CNL, "Type")
pairs(b.emm.s)

################################## get C:N ratio from Holly et al. 2014 data ###################################
#https://data.nal.usda.gov/dataset/carbon-dioxide-methane-nitrous-oxide-and-ammonia-emissions-digested-and-separated-dairy-manure-during-storage-and-land-application
# DOI 10.15482/USDA.ADC/1332491
Holly2014 <- read.csv("/Users/susannewiesner/Downloads/AFM/Manure_application_2020.csv")

lme_cn <- lme(C.N~Total_N_., random = (~1|Date/Replicate), Holly2014)
plot(lme_cn)
summary(lme_cn)
intervals(lme_cn, "sigma")
emmeans(lme_cn, "Total_N_.") # 0.717
r.rmse <- sqrt(mean(lme_cn$residuals^2))

predictSE(lme_cn, newdata = data.frame(Total_N_. = c(0.13))) ####### predict C:N ratio for N percent of 0.13%