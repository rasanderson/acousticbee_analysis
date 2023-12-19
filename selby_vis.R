# Visualise Selby data prior to analysis
# 
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

tmp <- filter(rawd, date <= "2022-08-31" & site == "Selby")
# Acoustic indices with date ----
ggplot(tmp, aes(x=date, y=ACI, colour=as.factor(hour))) + # Acoustic complexity index
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=RMS, colour=as.factor(hour))) + # Root mean square
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=H, colour=as.factor(hour))) + # Shannon (?)
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=mean_dfreq, colour=as.factor(hour))) + # Mean dominant frequency
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=M, colour=as.factor(hour))) + # Median of amplitude envelope
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=Ht, colour=as.factor(hour))) + # Temporal entropy (noisy = 1, quiet = 0)
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=freq_mean, colour=as.factor(hour))) + # Mean frequency
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=sd, colour=as.factor(hour))) + # Standard deviation of mean frequency
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=freq_median, colour=as.factor(hour))) + # Median frequency
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=sem, colour=as.factor(hour))) + # Spectrum standard error of the mean frequency
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=freq_mode, colour=as.factor(hour))) + # Freq mode (dominant frequency)
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=Q25, colour=as.factor(hour))) + # Spectral first quartile
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=Q75, colour=as.factor(hour))) + # Spectral third quartile
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=IQR, colour=as.factor(hour))) + # Spectral inter-quartile range
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=cent, colour=as.factor(hour))) + # Spectral centroid
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=skewness, colour=as.factor(hour))) + # Skewness (assymmetry)
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=kurtosis, colour=as.factor(hour))) + # Kurtosis (peakedness)
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=sfm, colour=as.factor(hour))) + # Spectral flatness measure
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=sh, colour=as.factor(hour))) + # Shannon and Reyni spectral entropy
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=prec, colour=as.factor(hour))) + # Frequency precision of the spectrum
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=time.P1, colour=as.factor(hour))) + # STFT time initial percentile
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=time.M, colour=as.factor(hour))) + # STFT time median
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=time.P2, colour=as.factor(hour))) + # STFT time terminal percentile
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=time.IPR, colour=as.factor(hour))) + # STFT time inter-percentile range
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=freq.P1, colour=as.factor(hour))) + # STFT freq initial percentile
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=freq.M, colour=as.factor(hour))) + # STFT frequency median
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=freq.P2, colour=as.factor(hour))) + # STFT frequency terminal percentile
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=freq.IPR, colour=as.factor(hour))) + # STFT frequency inter-percentile range
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)

# Check acoustic measures in a PCA out of curiosity ----
library(vegan)
tmp2 <- tmp[, 14:41]
tmp2 <- decostand(tmp2, method = "hellinger")
acoustic_pca <- rda(tmp2)
plot(acoustic_pca, display = "species")
plot(acoustic_pca, display = "sites")
acoustic_sco <- data.frame(scores(acoustic_pca, display="sites"))

ggplot(tmp, aes(x=acoustic_sco$PC1, y=acoustic_sco$PC2, colour=colony)) +
  geom_point() +
  facet_wrap(~as.factor(hour))

ggplot(tmp, aes(x=date, y=acoustic_sco$PC1, colour=as.factor(hour))) +
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
ggplot(tmp, aes(x=date, y=acoustic_sco$PC2, colour=as.factor(hour))) +
  geom_smooth(se=FALSE) +
  facet_wrap(~colony)
detach("package:vegan", unload = TRUE)

ggplot(tmp, aes(x=date, y=acoustic_sco$PC1, colour=colony)) +
  geom_point() +
  facet_wrap(~site)
ggplot(tmp, aes(x=date, y=acoustic_sco$PC2, colour=colony)) +
  geom_point() +
  facet_wrap(~site)
detach("package:vegan", unload = TRUE)

# Compare some of disease levels over time ----
ggplot(tmp, aes(x=date, y=varroa_per_300_bees1, colour=colony)) +
  #geom_point() +
  geom_smooth(se=FALSE) 
# Not enough CBPV to work with at Selby
ggplot(tmp, aes(x=date, y=cbpv_status1, colour=colony)) +
  geom_point() +
  facet_wrap(~site)

# And by disease severity with acoustics ----
ggplot(tmp, aes(y = varroa_per_300_bees1, x = acoustic_sco$PC1, colour=colony)) +
  xlab("Acoustic PCA1 (66.7%)") +
  geom_smooth() 
ggplot(tmp, aes(y = varroa_per_300_bees1, x = acoustic_sco$PC2, colour=colony)) +
  xlab("Acoustic PCA2 (27.0%)") +
  geom_smooth()

ggplot(tmp, aes(y = varroa_per_300_bees1, x = acoustic_sco$PC1, colour=colony)) +
  xlab("Acoustic PCA1 (66.7%)") +
  geom_smooth() +
  facet_wrap(~as.factor(hour))
ggplot(tmp, aes(y = varroa_per_300_bees1, x = acoustic_sco$PC2, colour=colony)) +
  xlab("Acoustic PCA2 (27.0%)") +
  geom_smooth() +
  facet_wrap(~as.factor(hour))

