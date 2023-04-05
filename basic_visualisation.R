# Basic visualisation

rm(list = ls())

library(ggplot2)
library(dplyr)
library(lubridate)

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

tmp <- filter(rawd, date <= "2022-08-31")

ggplot(tmp, aes(x=date, y=ACI, colour=colony)) + # Acoustic complexity index
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=RMS, colour=colony)) + # Root mean square
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=H, colour=colony)) + # Shannon (?)
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=mean_dfreq, colour=colony)) + # Mean dominant frequency
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=M, colour=colony)) + # Median of amplitude envelope
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=Ht, colour=colony)) + # Temporal entropy (noisy = 1, quiet = 0)
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=freq_mean, colour=colony)) + # Mean frequency
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=sd, colour=colony)) + # Standard deviation of mean frequency
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=freq_median, colour=colony)) + # Median frequency
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=sem, colour=colony)) + # Spectrum standard error of the mean frequency
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=freq_mode, colour=colony)) + # Freq mode (dominant frequency)
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=Q25, colour=colony)) + # Spectral first quartile
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=Q75, colour=colony)) + # Spectral third quartile
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=IQR, colour=colony)) + # Spectral inter-quartile range
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=cent, colour=colony)) + # Spectral centroid
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=skewness, colour=colony)) + # Skewness (assymmetry)
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=kurtosis, colour=colony)) + # Kurtosis (peakedness)
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=sfm, colour=colony)) + # Spectral flatness measure
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=sh, colour=colony)) + # Shannon and Reyni spectral entropy
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=prec, colour=colony)) + # Frequency precision of the spectrum
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=time.P1, colour=colony)) + # STFT time initial percentile
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=time.M, colour=colony)) + # STFT time median
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=time.P2, colour=colony)) + # STFT time terminal percentile
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=time.IPR, colour=colony)) + # STFT time inter-percentile range
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=freq.P1, colour=colony)) + # STFT freq initial percentile
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=freq.M, colour=colony)) + # STFT frequency median
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=freq.P2, colour=colony)) + # STFT frequency terminal percentile
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=freq.IPR, colour=colony)) + # STFT frequency inter-percentile range
  geom_point() +
  facet_wrap(~site)

# Check acoustic measures in a PCA out of curiosity
library(vegan)
tmp2 <- tmp[, 14:41]
tmp2 <- decostand(tmp2, method = "hellinger")
acoustic_pca <- rda(tmp2)
plot(acoustic_pca, display = "species")
plot(acoustic_pca, display = "sites")
acoustic_sco <- data.frame(scores(acoustic_pca, display="sites"))

ggplot(tmp, aes(x=date, y=acoustic_sco$PC1, colour=colony)) +
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=acoustic_sco$PC2, colour=colony)) +
  geom_point() +
  facet_wrap(~site)
detach("package:vegan", unload = TRUE)

# Compare some of disease levels over time
ggplot(tmp, aes(x=date, y=varroa_per_300_bees1, colour=colony)) +
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=cbpv_status1, colour=colony)) +
  geom_point() +
  facet_wrap(~site)
# And by disease severity with acoustics
ggplot(tmp, aes(x = varroa_per_300_bees1, y = acoustic_sco$PC1, colour=colony)) +
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x = varroa_per_300_bees1, y = acoustic_sco$PC2, colour=colony)) +
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x = cbpv_status1, y = acoustic_sco$PC1, colour=colony)) +
  geom_boxplot() +
  facet_wrap(~site)
ggplot(tmp, aes(x = cbpv_status1, y = acoustic_sco$PC2, colour=colony)) +
  geom_boxplot() +
  facet_wrap(~site)

# A few correlations of acoustic stuff
tmp3 <- cbind(tmp2, acoustic_sco)
library(ggcorrplot)
corr <- cor(tmp3)
ggcorrplot(corr, hc.order = TRUE)


# Some comparisons of acoustics with weather
ggplot(tmp, aes(x = temperature, y = ACI, colour = colony)) +
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x = windsp, y = ACI, colour = colony)) +
  geom_point() +
  facet_wrap(~site)
ggplot(tmp4, aes(x = rain_mm_h_mean, y = mean_dfreq, colour = colony)) +
  geom_point() +
  facet_wrap(~site)
