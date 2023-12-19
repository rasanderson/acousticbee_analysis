# Selby statistical modelling approach
# 
rm(list = ls())

library(ggplot2)
library(dplyr)
library(lubridate)
library(DHARMa)
library(glmmTMB)

rawd <- read.csv("data/All_Data_Master_Shortened.csv")
# Three rows with missing values causing trouble need removing
rawd <- rawd[!is.na(rawd$ACI),]
rawd$timestamp    <- ymd_hms(rawd$timestamp)
rawd$date         <- ymd(rawd$date)
rawd$time         <- hms(rawd$time)
rawd$colony       <- as.factor(rawd$colony)
rawd$cbpv_status  <- as.factor(rawd$cbpv_status)
rawd$cbpv_status1 <- as.factor(rawd$cbpv_status1)
# CBPV coded None, Unsure, Low, Medium, Very High
rawd$cbpv_status  <- ordered(rawd$cbpv_status, levels = c("N", "US", "L", "M", "VH"))
rawd$cbpv_status1  <- ordered(rawd$cbpv_status1, levels = c("N", "US", "L", "M", "VH"))

selby_dat <- filter(rawd, date <= "2022-08-31" & site == "Selby")
selby_dat <- filter(selby_dat, !is.na(rain_mm_h_mean))
selby_dat <- filter(selby_dat, !is.na(temperature))


# Create PCA of acoustics ----
library(vegan)
tmp2 <- selby_dat[, 14:41]
tmp2 <- decostand(tmp2, method = "hellinger")
acoustic_pca <- rda(tmp2)
plot(acoustic_pca, display="sites")
acoustic_sco <- data.frame(scores(acoustic_pca, display="sites"))
detach("package:vegan", unload = TRUE)

# Acoustic patterns with varroa
# And by disease severity with acoustics ----
ggplot(selby_dat, aes(y = varroa_per_300_bees1, x = acoustic_sco$PC1, colour=colony)) +
  xlab("Acoustic PCA1 (66.7%)") +
  geom_smooth() 
ggplot(selby_dat, aes(y = varroa_per_300_bees1, x = acoustic_sco$PC2, colour=colony)) +
  xlab("Acoustic PCA2 (27.0%)") +
  geom_smooth()


# Assemble data for statistical models
sin_hr <- sin(2*pi*selby_dat$hour/24)
cos_hr <- cos(2*pi*selby_dat$hour/24)
sin_day <- sin(2*pi*selby_dat$day/365)
cos_day <- cos(2*pi*selby_dat$day/365)
selby_dat <- data.frame(cbind(selby_dat), sin_hr, cos_hr,
                        sin_day, cos_day, acoustic_sco)

# Basic lm / glm ignoring colony ----
selby_glm1 <- glmmTMB(varroa_per_300_bees1 ~ sin_hr + cos_hr + sin_day + cos_day,
                      data = selby_dat)
summary(selby_glm1)
simulateResiduals(selby_glm1, plot=TRUE)
selby_glm2 <- glmmTMB(log(varroa_per_300_bees1+1) ~ sin_hr + cos_hr + sin_day + cos_day,
                      data = selby_dat)
summary(selby_glm2)
simulateResiduals(selby_glm2, plot=TRUE)
selby_glm3 <- glmmTMB(log(varroa_per_300_bees1+1) ~ sin_day + cos_day + 
                       PC1 + PC2,
                     data = selby_dat)
summary(selby_glm3)
simulateResiduals(selby_glm3, plot=TRUE)
selby_glm4 <- glmmTMB(log(varroa_per_300_bees1+1) ~ sin_day + cos_day + 
                        PC1 + PC2 + rain_mm_h_mean + temperature + windsp,
                      data = selby_dat)
summary(selby_glm4)
simulateResiduals(selby_glm4, plot=TRUE)
library(nlme)
selby_glm4 <- gls(log(varroa_per_300_bees1+1) ~ sin_day + cos_day + 
                        PC1 + PC2 + rain_mm_h_mean + temperature + windsp,
                      data = selby_dat[!is.na(selby_dat$rain_mm_h_mean),])
summary(selby_glm4)
plot(ACF(selby_glm4))

# Mixed effect
selby_glm5 <- glmmTMB(log(varroa_per_300_bees1+1) ~ sin_day + cos_day + 
                        PC1 + PC2 + rain_mm_h_mean + temperature + windsp +
                        (1 | colony),
                      data = selby_dat)
summary(selby_glm5)
simulateResiduals(selby_glm5, plot=TRUE)
