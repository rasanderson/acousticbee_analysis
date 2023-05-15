# Try a simple lme-style approach

# Just operate a Selby initially (can later try Hexham for first part) as there
# seems to be a steadily increasing amount of Varroa.

# For simplicity will use weather plus acoustic PC1 and PC2 as predictors, hive
# a random effect and log numbers of varroa as response

# Data pre-process
rm(list = ls())

library(ggplot2)
library(dplyr)
library(lubridate)
library(tidyr)
library(scales)
library(nlme)

# Read in, clean and pre-process data ####
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
data_subset <- filter(rawd, date <= "2022-08-31")
data_subset$day_of_year <- yday(data_subset$date)
data_subset$cos_hour <- cos((2 * pi * data_subset$hour ) / 24)
data_subset$sin_hour <- sin((2 * pi * data_subset$hour ) / 24)
data_subset$grid_easting2 <- data_subset$grid_easting^2
data_subset$grid_northing2 <- data_subset$grid_northing^2
# Note there are 175 NA for station_easting/northing, but 320 for grid_easting/n
data_subset <- filter(data_subset, site == "Selby") %>% 
  drop_na(grid_easting) %>% 
  drop_na(grid_northing) %>% 
  drop_na(rain_mm_h_mean) %>% 
  drop_na(windsp) %>% 
  drop_na(temperature) %>% 
  filter(colony != 10) # Colony 10 at Selby has no data for first 100 days?

# Setup data for "meta-" variables ####
meta_variables <- c("varroa_per_300_bees1", "day_of_year", "cos_hour", "sin_hour",
                    "grid_easting", "grid_northing", "rain_mm_h_mean", "windsp",
                    "temperature", "grid_easting2", "grid_northing2")
acoustic_variables <- c("varroa_per_300_bees1", "ACI", "RMS", "H", "mean_dfreq",
                        "M", "Ht", "freq_mean", "sd", "freq_median", "sem",
                        "freq_mode", "Q25", "Q75", "IQR", "cent", "skewness",
                        "kurtosis", "sfm", "sh", "prec", "time.P1", "time.M",
                        "time.P2", "time.IPR", "freq.P1", "freq.M", "freq.P2",
                        "freq.IPR")
meta_data     <- select(data_subset, all_of(meta_variables))
acoustic_data <- select(data_subset, all_of(acoustic_variables))

# PCA of acoustic data
library(vegan)
tmp <- acoustic_data[, -1]
tmp <- decostand(tmp, method = "hellinger")
acoustic_pca <- rda(tmp)
plot(acoustic_pca, display = "species")
plot(acoustic_pca, display = "sites")
acoustic_sco <- data.frame(scores(acoustic_pca, display="sites"))

ggplot(data_subset, aes(x=acoustic_sco$PC1, y=acoustic_sco$PC2, colour=colony)) +
  geom_point()

ggplot(data_subset, aes(x=date, y=acoustic_sco$PC1, colour=colony)) +
  geom_point()
ggplot(data_subset, aes(x=date, y=acoustic_sco$PC2, colour=colony)) +
  geom_point()
detach("package:vegan", unload = TRUE)

# Simple lme
varroa_dat <- cbind(data_subset, acoustic_PC1 = acoustic_sco$PC1, acoustic_PC2 = acoustic_sco$PC2)
# Overview
ggplot(varroa_dat, aes(x = day_of_year, y = varroa_per_300_bees1)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_continuous(trans = log1p_trans()) +
  facet_wrap(~colony)

ggplot(varroa_dat, aes(x = day_of_year, y = varroa_per_300_bees1)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_continuous(trans = log1p_trans())


varroa_lme1 <- lme(log(varroa_per_300_bees1+1) ~ cos_hour +
                     sin_hour + rain_mm_h_mean + windsp + temperature,
                   random = ~day_of_year|colony, data = varroa_dat)
summary(varroa_lme1)

varroa_lm1 <- lm(log(varroa_per_300_bees1+1) ~ cos_hour + sin_hour + 
                  rain_mm_h_mean * windsp * temperature * day_of_year,
                data = varroa_dat)
summary(varroa_lm1)
varroa_lm2 <- update(varroa_lm1, . ~ . - rain_mm_h_mean:windsp:temperature:day_of_year)
summary(varroa_lm2)
varroa_lm3 <- update(varroa_lm2, . ~ . - rain_mm_h_mean:temperature:day_of_year)
summary(varroa_lm3)
varroa_lm4 <- update(varroa_lm3, . ~ . - rain_mm_h_mean:windsp:day_of_year)
summary(varroa_lm4)
varroa_lm5 <- update(varroa_lm4, . ~ . - rain_mm_h_mean:windsp:temperature)
summary(varroa_lm5)
varroa_lm6 <- update(varroa_lm5, . ~ . - rain_mm_h_mean:windsp)
summary(varroa_lm6)
varroa_lm7 <- update(varroa_lm6, . ~ . - temperature:day_of_year)
summary(varroa_lm7)
varroa_lm8 <- update(varroa_lm7, . ~ . - windsp:temperature)
summary(varroa_lm8)
varroa_lm9 <- update(varroa_lm8, . ~ . - windsp:temperature:day_of_year)
summary(varroa_lm9)
varroa_lm10 <- update(varroa_lm9, . ~ . - rain_mm_h_mean)
summary(varroa_lm10)

varroa_resid <- residuals(varroa_lm10)
acoustic_lm1 <- lm(varroa_resid ~ acoustic_PC1 + acoustic_PC2, data = varroa_dat)
summary(acoustic_lm1)

ggplot(varroa_dat, aes(x = acoustic_PC1, y = varroa_resid, colour = colony)) + 
  geom_point() + 
  geom_smooth()
ggplot(varroa_dat, aes(x = acoustic_PC2, y = varroa_resid, colour = colony)) + 
  geom_point() + 
  geom_smooth()
ggplot(varroa_dat, aes(x = data_subset$IQR, y = varroa_resid, colour = colony)) + 
  geom_point() + 
  geom_smooth()
