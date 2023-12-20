# Selby statistical modelling approach
# 
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
selby_glm4a <- glm(log(varroa_per_300_bees1+1) ~ sin_day + cos_day + 
                        PC1 + PC2 + rain_mm_h_mean + temperature + windsp,
                      data = selby_dat)
summary(selby_glm4a)

library(nlme)
selby_glm4b <- gls(log(varroa_per_300_bees1+1) ~ sin_day + cos_day + 
                        PC1 + PC2 + rain_mm_h_mean + temperature + windsp,
                      data = selby_dat[!is.na(selby_dat$rain_mm_h_mean),])
summary(selby_glm4b)
plot(ACF(selby_glm4b))

# Mixed effect
selby_glm5 <- glmmTMB(log(varroa_per_300_bees1+1) ~ sin_day + cos_day + 
                        PC1 + PC2 + rain_mm_h_mean + temperature + windsp +
                        (1 | colony),
                      data = selby_dat)
summary(selby_glm5)
simulateResiduals(selby_glm5, plot=TRUE)


# Try testing Selby-derived model glm4a vs Hexham ----
hexham_dat <- filter(rawd, (date >= "2022-05-25" & date <= "2022-06-21") & site == "Hexham")
hexham_dat <- filter(hexham_dat, !is.na(rain_mm_h_mean))
hexham_dat <- filter(hexham_dat, !is.na(temperature))

library(vegan)
tmp2 <- hexham_dat[, 14:41]
tmp2 <- decostand(tmp2, method = "hellinger")
hexham_sco <- data.frame(predict(acoustic_pca, newdata = tmp2, type = "wa"))
sin_hr <- sin(2*pi*hexham_dat$hour/24)
cos_hr <- cos(2*pi*hexham_dat$hour/24)
sin_day <- sin(2*pi*hexham_dat$day/365)
cos_day <- cos(2*pi*hexham_dat$day/365)
hexham_dat <- data.frame(cbind(hexham_dat), sin_hr, cos_hr,
                        sin_day, cos_day, 
                        PC1 = hexham_sco$PC1, PC2 = hexham_sco$PC2)

hexham_pred <- predict(selby_glm4a, newdata = hexham_dat, type = "response")

hexham_dat <- cbind(hexham_dat, hexham_pred)
p1 <- ggplot(hexham_dat, aes(x = date, y = varroa_per_300_bees1, colour = colony)) +
  geom_smooth() 
p2 <- ggplot(hexham_dat, aes(x = date, y = exp(hexham_pred)-1, colour = colony)) +
  geom_smooth() 
grid.arrange(p1, p2, nrow = 1)
hexham_obs_pred_lm <- lm(log(varroa_per_300_bees1 + 1) ~ I(exp(hexham_pred)-1), data = hexham_dat)
summary(hexham_obs_pred_lm)
hexham_obs_pred_lme <- lme(log(varroa_per_300_bees1 + 1) ~ I(exp(hexham_pred)-1),
                           random = ~1|colony, data = hexham_dat)
summary(hexham_obs_pred_lme)
r.squaredLR(hexham_obs_pred_lme)

# Try creating model for 3 out of 4 Selby nests ----
selby_dat <- mutate(selby_dat, PC1 = NULL, PC2 = NULL)
## Omit Selby 1, setup model ----
selby_omit <- filter(selby_dat, colony != 1)
selby_test <- filter(selby_dat, colony == 1)
selby_acoustic_omit <- selby_omit[, 14:41]
selby_acoustic_omit <- decostand(selby_acoustic_omit, method = "hellinger")
selby_omit_pca <- rda(selby_acoustic_omit)
plot(selby_omit_pca, display="sites")
selby_omit_acoustic_sco <- data.frame(scores(selby_omit_pca, display="sites"))
selby_omit <- data.frame(cbind(selby_omit, selby_omit_acoustic_sco))
selby_omit_glm <- glm(log(varroa_per_300_bees1+1) ~ sin_day + cos_day + 
                     PC1 + PC2 + rain_mm_h_mean + temperature + windsp,
                   data = selby_omit)
summary(selby_omit_glm)
## Test Selby1 model ----
selby_test_acoustic <- selby_test[, 14:41]
selby_test_acoustic <- decostand(selby_test_acoustic, method = "hellinger")
selby_test_sco <- data.frame(predict(selby_omit_pca, newdata = selby_test_acoustic, type = "wa"))
selby_test <- data.frame(cbind(selby_test), 
                         PC1 = selby_test_sco$PC1, PC2 = selby_test_sco$PC2)
selby_pred <- predict(selby_omit_glm, newdata = selby_test, type = "response")

selby_test <- cbind(selby_test, selby_pred)
p1 <- ggplot(selby_test, aes(x = date, y = varroa_per_300_bees1, colour = colony)) +
  geom_smooth() 
p2 <- ggplot(selby_test, aes(x = date, y = exp(selby_pred)-1, colour = colony)) +
  geom_smooth() 
grid.arrange(p1, p2, nrow = 1)
selby_obs_pred_lm <- lm(log(varroa_per_300_bees1 + 1) ~ I(exp(selby_pred)-1), data = selby_test)
summary(selby_obs_pred_lm)
## Omit Selby 4 ----
selby_omit <- filter(selby_dat, colony != 4)
selby_test <- filter(selby_dat, colony == 4)
selby_acoustic_omit <- selby_omit[, 14:41]
selby_acoustic_omit <- decostand(selby_acoustic_omit, method = "hellinger")
selby_omit_pca <- rda(selby_acoustic_omit)
plot(selby_omit_pca, display="sites")
selby_omit_acoustic_sco <- data.frame(scores(selby_omit_pca, display="sites"))
selby_omit <- data.frame(cbind(selby_omit, selby_omit_acoustic_sco))
selby_omit_glm <- glm(log(varroa_per_300_bees1+1) ~ sin_day + cos_day + 
                        PC1 + PC2 + rain_mm_h_mean + temperature + windsp,
                      data = selby_omit)
summary(selby_omit_glm)
## Test Selby4 model ----
selby_test_acoustic <- selby_test[, 14:41]
selby_test_acoustic <- decostand(selby_test_acoustic, method = "hellinger")
selby_test_sco <- data.frame(predict(selby_omit_pca, newdata = selby_test_acoustic, type = "wa"))
selby_test <- data.frame(cbind(selby_test), 
                         PC1 = selby_test_sco$PC1, PC2 = selby_test_sco$PC2)
selby_pred <- predict(selby_omit_glm, newdata = selby_test, type = "response")

selby_test <- cbind(selby_test, selby_pred)
p1 <- ggplot(selby_test, aes(x = date, y = varroa_per_300_bees1, colour = colony)) +
  geom_smooth() 
p2 <- ggplot(selby_test, aes(x = date, y = exp(selby_pred)-1, colour = colony)) +
  geom_smooth() 
grid.arrange(p1, p2, nrow = 1)
selby_obs_pred_lm <- lm(log(varroa_per_300_bees1 + 1) ~ I(exp(selby_pred)-1), data = selby_test)
summary(selby_obs_pred_lm)
## Omit Selby 6 ----
selby_omit <- filter(selby_dat, colony != 6)
selby_test <- filter(selby_dat, colony == 6)
selby_acoustic_omit <- selby_omit[, 14:41]
selby_acoustic_omit <- decostand(selby_acoustic_omit, method = "hellinger")
selby_omit_pca <- rda(selby_acoustic_omit)
plot(selby_omit_pca, display="sites")
selby_omit_acoustic_sco <- data.frame(scores(selby_omit_pca, display="sites"))
selby_omit <- data.frame(cbind(selby_omit, selby_omit_acoustic_sco))
selby_omit_glm <- glm(log(varroa_per_300_bees1+1) ~ sin_day + cos_day + 
                        PC1 + PC2 + rain_mm_h_mean + temperature + windsp,
                      data = selby_omit)
summary(selby_omit_glm)
## Test Selby6 model ----
selby_test_acoustic <- selby_test[, 14:41]
selby_test_acoustic <- decostand(selby_test_acoustic, method = "hellinger")
selby_test_sco <- data.frame(predict(selby_omit_pca, newdata = selby_test_acoustic, type = "wa"))
selby_test <- data.frame(cbind(selby_test), 
                         PC1 = selby_test_sco$PC1, PC2 = selby_test_sco$PC2)
selby_pred <- predict(selby_omit_glm, newdata = selby_test, type = "response")

selby_test <- cbind(selby_test, selby_pred)
p1 <- ggplot(selby_test, aes(x = date, y = varroa_per_300_bees1, colour = colony)) +
  geom_smooth() 
p2 <- ggplot(selby_test, aes(x = date, y = exp(selby_pred)-1, colour = colony)) +
  geom_smooth() 
grid.arrange(p1, p2, nrow = 1)
selby_obs_pred_lm <- lm(log(varroa_per_300_bees1 + 1) ~ I(exp(selby_pred)-1), data = selby_test)
summary(selby_obs_pred_lm)
## Omit Selby 8
selby_omit <- filter(selby_dat, colony != 8)
selby_test <- filter(selby_dat, colony == 8)
selby_acoustic_omit <- selby_omit[, 14:41]
selby_acoustic_omit <- decostand(selby_acoustic_omit, method = "hellinger")
selby_omit_pca <- rda(selby_acoustic_omit)
plot(selby_omit_pca, display="sites")
selby_omit_acoustic_sco <- data.frame(scores(selby_omit_pca, display="sites"))
selby_omit <- data.frame(cbind(selby_omit, selby_omit_acoustic_sco))
selby_omit_glm <- glm(log(varroa_per_300_bees1+1) ~ sin_day + cos_day + 
                        PC1 + PC2 + rain_mm_h_mean + temperature + windsp,
                      data = selby_omit)
summary(selby_omit_glm)
## Test Selby8 model ----
selby_test_acoustic <- selby_test[, 14:41]
selby_test_acoustic <- decostand(selby_test_acoustic, method = "hellinger")
selby_test_sco <- data.frame(predict(selby_omit_pca, newdata = selby_test_acoustic, type = "wa"))
selby_test <- data.frame(cbind(selby_test), 
                         PC1 = selby_test_sco$PC1, PC2 = selby_test_sco$PC2)
selby_pred <- predict(selby_omit_glm, newdata = selby_test, type = "response")

selby_test <- cbind(selby_test, selby_pred)
p1 <- ggplot(selby_test, aes(x = date, y = varroa_per_300_bees1, colour = colony)) +
  geom_smooth() 
p2 <- ggplot(selby_test, aes(x = date, y = exp(selby_pred)-1, colour = colony)) +
  geom_smooth() 
grid.arrange(p1, p2, nrow = 1)
selby_obs_pred_lm <- lm(log(varroa_per_300_bees1 + 1) ~ I(exp(selby_pred)-1), data = selby_test)
summary(selby_obs_pred_lm)


# Try a lag model for all 4 Selby colonies ----
## Create lagged datasets
selby_lead1 <- selby_dat %>% 
  filter(colony == 1) %>% 
  mutate(lead_varroa = lead(varroa_per_300_bees1, n = 20)) %>% 
  drop_na(lead_varroa)
selby_lead4 <- selby_dat %>% 
  filter(colony == 4) %>% 
  mutate(lead_varroa = lead(varroa_per_300_bees1, n = 20)) %>% 
  drop_na(lead_varroa)
selby_lead6 <- selby_dat %>% 
  filter(colony == 6) %>% 
  mutate(lead_varroa = lead(varroa_per_300_bees1, n = 20)) %>% 
  drop_na(lead_varroa)
selby_lead8 <- selby_dat %>% 
  filter(colony == 8) %>% 
  mutate(lead_varroa = lead(varroa_per_300_bees1, n = 20)) %>% 
  drop_na(lead_varroa)
selby_lead <- rbind(selby_lead1, selby_lead4, selby_lead6, selby_lead8)
# Now create selby model with 1 day lag varroa (20 records) ----
#selby_omit <- filter(selby_dat, colony != 1)
#selby_test <- filter(selby_dat, colony == 1)
selby_omit <- selby_lead
selby_acoustic_omit <- selby_omit[, 14:41]
selby_acoustic_omit <- decostand(selby_acoustic_omit, method = "hellinger")
selby_omit_pca <- rda(selby_acoustic_omit)
plot(selby_omit_pca, display="sites")
selby_omit_acoustic_sco <- data.frame(scores(selby_omit_pca, display="sites"))
selby_omit <- data.frame(cbind(selby_omit, selby_omit_acoustic_sco))
selby_omit_lm <- lm(log(varroa_per_300_bees1+1) ~ sin_day + cos_day + 
                        PC1 + PC2 + rain_mm_h_mean + temperature + windsp +
                        lead_varroa,
                      data = selby_omit)
summary(selby_omit_lm)
# Test lagged (lead) model for each Selby colony ----
## Lagged Selby1 ----
selby_omit <- filter(selby_lead, colony != 1)
selby_test <- filter(selby_lead, colony == 1)
selby_acoustic_omit <- selby_omit[, 14:41]
selby_acoustic_omit <- decostand(selby_acoustic_omit, method = "hellinger")
selby_omit_pca <- rda(selby_acoustic_omit)
plot(selby_omit_pca, display="sites")
selby_omit_acoustic_sco <- data.frame(scores(selby_omit_pca, display="sites"))
selby_omit <- data.frame(cbind(selby_omit, selby_omit_acoustic_sco))
selby_omit_lm <- glm(log(varroa_per_300_bees1+1) ~ sin_day + cos_day + 
                      PC1 + PC2 + rain_mm_h_mean + temperature + windsp +
                      lead_varroa,
                    data = selby_omit)
summary(selby_omit_glm)
selby_test_acoustic <- selby_test[, 14:41]
selby_test_acoustic <- decostand(selby_test_acoustic, method = "hellinger")
selby_test_sco <- data.frame(predict(selby_omit_pca, newdata = selby_test_acoustic, type = "wa"))
selby_test <- data.frame(cbind(selby_test), 
                         PC1 = selby_test_sco$PC1, PC2 = selby_test_sco$PC2)
selby_pred <- predict(selby_omit_glm, newdata = selby_test, type = "response")
selby_test <- cbind(selby_test, selby_pred)
p1 <- ggplot(selby_test, aes(x = date, y = varroa_per_300_bees1, colour = colony)) +
  geom_smooth() 
p2 <- ggplot(selby_test, aes(x = date, y = exp(selby_pred)-1, colour = colony)) +
  geom_smooth() 
grid.arrange(p1, p2, nrow = 1)
selby_obs_pred_lm <- lm(log(varroa_per_300_bees1 + 1) ~ I(exp(selby_pred)-1), data = selby_test)
summary(selby_obs_pred_lm)
## Lagged Selby4 ----
selby_omit <- filter(selby_lead, colony != 4)
selby_test <- filter(selby_lead, colony == 4)
selby_acoustic_omit <- selby_omit[, 14:41]
selby_acoustic_omit <- decostand(selby_acoustic_omit, method = "hellinger")
selby_omit_pca <- rda(selby_acoustic_omit)
plot(selby_omit_pca, display="sites")
selby_omit_acoustic_sco <- data.frame(scores(selby_omit_pca, display="sites"))
selby_omit <- data.frame(cbind(selby_omit, selby_omit_acoustic_sco))
selby_omit_lm <- glm(log(varroa_per_300_bees1+1) ~ sin_day + cos_day + 
                       PC1 + PC2 + rain_mm_h_mean + temperature + windsp +
                       lead_varroa,
                     data = selby_omit)
summary(selby_omit_glm)
selby_test_acoustic <- selby_test[, 14:41]
selby_test_acoustic <- decostand(selby_test_acoustic, method = "hellinger")
selby_test_sco <- data.frame(predict(selby_omit_pca, newdata = selby_test_acoustic, type = "wa"))
selby_test <- data.frame(cbind(selby_test), 
                         PC1 = selby_test_sco$PC1, PC2 = selby_test_sco$PC2)
selby_pred <- predict(selby_omit_glm, newdata = selby_test, type = "response")
selby_test <- cbind(selby_test, selby_pred)
p1 <- ggplot(selby_test, aes(x = date, y = varroa_per_300_bees1, colour = colony)) +
  geom_smooth() 
p2 <- ggplot(selby_test, aes(x = date, y = exp(selby_pred)-1, colour = colony)) +
  geom_smooth() 
grid.arrange(p1, p2, nrow = 1)
selby_obs_pred_lm <- lm(log(varroa_per_300_bees1 + 1) ~ I(exp(selby_pred)-1), data = selby_test)
summary(selby_obs_pred_lm)
## Lagged Selby4 ----
selby_omit <- filter(selby_lead, colony != 6)
selby_test <- filter(selby_lead, colony == 6)
selby_acoustic_omit <- selby_omit[, 14:41]
selby_acoustic_omit <- decostand(selby_acoustic_omit, method = "hellinger")
selby_omit_pca <- rda(selby_acoustic_omit)
plot(selby_omit_pca, display="sites")
selby_omit_acoustic_sco <- data.frame(scores(selby_omit_pca, display="sites"))
selby_omit <- data.frame(cbind(selby_omit, selby_omit_acoustic_sco))
selby_omit_lm <- glm(log(varroa_per_300_bees1+1) ~ sin_day + cos_day + 
                       PC1 + PC2 + rain_mm_h_mean + temperature + windsp +
                       lead_varroa,
                     data = selby_omit)
summary(selby_omit_glm)
selby_test_acoustic <- selby_test[, 14:41]
selby_test_acoustic <- decostand(selby_test_acoustic, method = "hellinger")
selby_test_sco <- data.frame(predict(selby_omit_pca, newdata = selby_test_acoustic, type = "wa"))
selby_test <- data.frame(cbind(selby_test), 
                         PC1 = selby_test_sco$PC1, PC2 = selby_test_sco$PC2)
selby_pred <- predict(selby_omit_glm, newdata = selby_test, type = "response")
selby_test <- cbind(selby_test, selby_pred)
p1 <- ggplot(selby_test, aes(x = date, y = varroa_per_300_bees1, colour = colony)) +
  geom_smooth() 
p2 <- ggplot(selby_test, aes(x = date, y = exp(selby_pred)-1, colour = colony)) +
  geom_smooth() 
grid.arrange(p1, p2, nrow = 1)
selby_obs_pred_lm <- lm(log(varroa_per_300_bees1 + 1) ~ I(exp(selby_pred)-1), data = selby_test)
summary(selby_obs_pred_lm)
## Lagged Selby4 ----
selby_omit <- filter(selby_lead, colony != 8)
selby_test <- filter(selby_lead, colony == 8)
selby_acoustic_omit <- selby_omit[, 14:41]
selby_acoustic_omit <- decostand(selby_acoustic_omit, method = "hellinger")
selby_omit_pca <- rda(selby_acoustic_omit)
plot(selby_omit_pca, display="sites")
selby_omit_acoustic_sco <- data.frame(scores(selby_omit_pca, display="sites"))
selby_omit <- data.frame(cbind(selby_omit, selby_omit_acoustic_sco))
selby_omit_lm <- glm(log(varroa_per_300_bees1+1) ~ sin_day + cos_day + 
                       PC1 + PC2 + rain_mm_h_mean + temperature + windsp +
                       lead_varroa,
                     data = selby_omit)
summary(selby_omit_glm)
selby_test_acoustic <- selby_test[, 14:41]
selby_test_acoustic <- decostand(selby_test_acoustic, method = "hellinger")
selby_test_sco <- data.frame(predict(selby_omit_pca, newdata = selby_test_acoustic, type = "wa"))
selby_test <- data.frame(cbind(selby_test), 
                         PC1 = selby_test_sco$PC1, PC2 = selby_test_sco$PC2)
selby_pred <- predict(selby_omit_glm, newdata = selby_test, type = "response")
selby_test <- cbind(selby_test, selby_pred)
p1 <- ggplot(selby_test, aes(x = date, y = varroa_per_300_bees1, colour = colony)) +
  geom_smooth() 
p2 <- ggplot(selby_test, aes(x = date, y = exp(selby_pred)-1, colour = colony)) +
  geom_smooth() 
grid.arrange(p1, p2, nrow = 1)
selby_obs_pred_lm <- lm(log(varroa_per_300_bees1 + 1) ~ I(exp(selby_pred)-1), data = selby_test)
summary(selby_obs_pred_lm)
