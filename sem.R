# SEM approach
# After chat with Steve suggests piecewiseSEM given high intercorrelations
# Also include day_no (maybe sin_day, cos_day can be omitted) with log varroa
# Check residuals for ACF

# Initial setup ----
rm(list = ls())

library(ggplot2)
library(tidyr)
library(dplyr)
library(lubridate)
library(DHARMa)
library(glmmTMB)
library(gridExtra)
library(MuMIn)
library(Metrics)
library(piecewiseSEM)
library(solrad)
library(chillR)
source("night_vs_day.R")

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

# Selby before 31st August. Colony 10 has no records in dataset until 22nd June
latitude_selby <- 53.475
selby_dat <- filter(rawd, date <= "2022-08-31" & site == "Selby" & colony != 10)
selby_dat <- filter(selby_dat, !is.na(rain_mm_h_mean))
selby_dat <- filter(selby_dat, !is.na(temperature))
selby_dat <- mutate(selby_dat, day_no = yday(date))
selby_dat <- mutate(selby_dat, day_length = DayLength(day_no, latitude_selby))
daynight <- light_dark(day_no = selby_dat$day_no,
                       hour = selby_dat$hour,
                       latitude = latitude_selby)
selby_dat <- mutate(selby_dat, daynight = daynight)
# Create PCA of acoustics ----
library(vegan)
tmp2 <- selby_dat[, 14:41]
tmp2 <- decostand(tmp2, method = "hellinger")
acoustic_pca <- rda(tmp2)
plot(acoustic_pca, display="sites")
plot(acoustic_pca, display="species")
acoustic_sco <- data.frame(scores(acoustic_pca, display="sites", ch=1:3))
detach("package:vegan", unload = TRUE)
## PCA plots for colonies, acoustics, varroa ----
ggplot(cbind(selby_dat, acoustic_sco), aes(x = colony, y = PC1)) +
  geom_boxplot()
ggplot(cbind(selby_dat, acoustic_sco), aes(x = colony, y = PC2)) +
  geom_boxplot()
ggplot(cbind(selby_dat, acoustic_sco), aes(x = date, y = PC1, colour = colony)) +
  geom_point() +
  geom_smooth()
ggplot(cbind(selby_dat, acoustic_sco), aes(x = date, y = PC2, colour = colony)) +
  geom_point() +
  geom_smooth()
ggplot(selby_dat, aes(x = colony, y = RMS)) +
  geom_boxplot()
ggplot(selby_dat, aes(x = colony, y = freq_mode)) +
  geom_boxplot()
ggplot(selby_dat, aes(x = colony, y = IQR)) +
  geom_boxplot()
ggplot(selby_dat, aes(x = date, y = RMS, colour = colony)) +
  geom_point() +
  geom_smooth()
ggplot(selby_dat, aes(x = date, y = freq_mode, colour = colony)) +
  geom_point() +
  geom_smooth()
ggplot(selby_dat, aes(x = date, y = IQR, colour = colony)) +
  geom_point() +
  geom_smooth()
ggplot(cbind(selby_dat, acoustic_sco), aes(x = PC1, y = varroa_per_300_bees1, colour = colony)) +
  geom_point() +
  geom_smooth()
ggplot(cbind(selby_dat, acoustic_sco), aes(x = PC2, y = varroa_per_300_bees1, colour = colony)) +
  geom_point() +
  geom_smooth()
 
# Plot changes of key variables with time ----
ggplot(selby_dat, aes(x = date, y = log(varroa_per_300_bees1+1), colour = colony)) +
  geom_point() +
  geom_smooth()
ggplot(selby_dat, aes(x = date, y = RMS, colour = colony)) +
  geom_point() +
  geom_smooth() 
ggplot(selby_dat, aes(x = date, y = freq_mode, colour = colony)) +
  geom_point() +
  geom_smooth()
ggplot(selby_dat, aes(x = date, y = IQR, colour = colony)) +
  geom_point() +
  geom_smooth()

