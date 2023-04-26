# Single site simple_regr
#
# Deep learning models that are regression based, but just pick one location
rm(list = ls())
location <- "Selby"

library(ggplot2)
library(dplyr)
library(lubridate)
library(tidyr)

# Read in and pre-process data ####
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

# Begin with just the "meta-data" about the environment ####
# Work on the Varroa as more variable response
# However, no records at Thorganby as not assessed. Remove missing values.
# Note there are 175 NA for station_easting/northing, but 320 for grid_easting/n
data_subset <- filter(data_subset, site != "Thorganby") %>% 
  filter(site == location) %>% 
  drop_na(grid_easting) %>% 
  drop_na(grid_northing) %>% 
  drop_na(rain_mm_h_mean) %>% 
  drop_na(windsp) %>% 
  drop_na(temperature)

meta_variables <- c("varroa_per_300_bees1", "day_of_year", "cos_hour", "sin_hour",
                    "grid_easting", "grid_northing", "rain_mm_h_mean", "windsp",
                    "temperature")
meta_subset <- select(data_subset, all_of(meta_variables))
library(keras)
library(rsample)
split <- initial_split(meta_subset, 0.8)
train_dataset <- training(split)
test_dataset <- testing(split)
# # Smoothed freq histograms and scatterplots (slow)
# train_dataset %>%
#   select(day_of_year, cos_hour, sin_hour, temperature) %>%
#   GGally::ggpairs()
# # Range etc.
# skimr::skim(meta_subset)

# Now split features from labels
train_features <- train_dataset %>% select(-varroa_per_300_bees1)
test_features <- test_dataset %>% select(-varroa_per_300_bees1)

train_labels <- train_dataset %>% select(varroa_per_300_bees1)
test_labels <- test_dataset %>% select(varroa_per_300_bees1)

# Normalise the data
normalizer <- layer_normalization(axis = -1L)
# The adapt() function fits the state of the normalized data and fixes it so
# it doesn't change subsequently
normalizer %>% adapt(as.matrix(train_features))
print(normalizer$mean)
first <- as.matrix(train_features[1,])
cat('First example:', first)
cat('Normalized:', as.matrix(normalizer(first)))

# As relatively simple, use a single function to define and compile model
build_and_compile_model <- function(norm) {
  model <- keras_model_sequential() %>%
    norm() %>%
    layer_dense(64, activation = 'relu') %>%
    layer_dense(64, activation = 'relu') %>%
    layer_dense(64, activation = 'relu') %>%
    layer_dense(1)
  
  model %>% compile(
    metrics = list("accuracy"),
    loss = 'mean_absolute_error',
    optimizer = optimizer_adam(0.001)
  )
  
  model
}

# Apply to normalized data
dnn_model <- build_and_compile_model(normalizer)
summary(dnn_model)
# Now the slow bit. Overtraining after about 20 epochs
history <- dnn_model %>% fit(
  as.matrix(train_features),
  as.matrix(train_labels),
  validation_split = 0.2,
  verbose = 1,
  epochs = 50
)
plot(history)
# Results on test dataset
test_results <- list()
test_results[['dnn_model']] <- dnn_model %>% evaluate(
  as.matrix(test_features),
  as.matrix(test_labels),
  verbose = 0
)
sapply(test_results, function(x) x)
# Make some predictions
test_predictions <- predict(dnn_model, as.matrix(test_features))
ggplot(data.frame(pred = as.numeric(test_predictions), varroa = test_labels$varroa_per_300_bees1)) +
  geom_point(aes(x = pred, y = varroa)) +
  geom_abline(intercept = 0, slope = 1, color = "blue") +
  geom_smooth(aes(x = pred, y = varroa), method = "lm", color = "red", se = FALSE)
# Error distribution
qplot(test_predictions - test_labels$varroa_per_300_bees1, geom = "density")



# Deep learning model for acoustic indices ####
acoustic_variables <- c("varroa_per_300_bees1", "ACI", "RMS", "H", "mean_dfreq",
                        "M", "Ht", "freq_mean", "sd", "freq_median", "sem",
                        "freq_mode", "Q25", "Q75", "IQR", "cent", "skewness",
                        "kurtosis", "sfm", "sh", "prec", "time.P1", "time.M",
                        "time.P2", "time.IPR", "freq.P1", "freq.M", "freq.P2",
                        "freq.IPR")
acoustic_subset <- select(data_subset, all_of(acoustic_variables))
library(keras)
library(rsample)
split <- initial_split(acoustic_subset, 0.8)
train_dataset <- training(split)
test_dataset <- testing(split)
# # Smoothed freq histograms and scatterplots (slow)
# train_dataset %>%
#   select(day_of_year, cos_hour, sin_hour, temperature) %>%
#   GGally::ggpairs()
# # Range etc.
# skimr::skim(meta_subset)

# Now split features from labels
train_features <- train_dataset %>% select(-varroa_per_300_bees1)
test_features <- test_dataset %>% select(-varroa_per_300_bees1)

train_labels <- train_dataset %>% select(varroa_per_300_bees1)
test_labels <- test_dataset %>% select(varroa_per_300_bees1)

# Normalise the data
normalizer <- layer_normalization(axis = -1L)
# The adapt() function fits the state of the normalized data and fixes it so
# it doesn't change subsequently
normalizer %>% adapt(as.matrix(train_features))
print(normalizer$mean)
first <- as.matrix(train_features[1,])
cat('First example:', first)
cat('Normalized:', as.matrix(normalizer(first)))

# As relatively simple, use a single function to define and compile model
build_and_compile_model <- function(norm) {
  model <- keras_model_sequential() %>%
    norm() %>%
    layer_dense(64, activation = 'relu') %>%
    layer_dense(64, activation = 'relu') %>%
    layer_dense(64, activation = 'relu') %>%
    layer_dense(1)
  
  model %>% compile(
    metrics = list("accuracy"),
    loss = 'mean_absolute_error',
    optimizer = optimizer_adam(0.001)
  )
  
  model
}

# Apply to normalized data
dnn_model <- build_and_compile_model(normalizer)
summary(dnn_model)
# Now the slow bit. Overtraining after about 20 epochs
history <- dnn_model %>% fit(
  as.matrix(train_features),
  as.matrix(train_labels),
  validation_split = 0.2,
  verbose = 1,
  epochs = 25
)
plot(history)
# Results on test dataset
test_results <- list()
test_results[['dnn_model']] <- dnn_model %>% evaluate(
  as.matrix(test_features),
  as.matrix(test_labels),
  verbose = 0
)
sapply(test_results, function(x) x)
# Make some predictions
test_predictions <- predict(dnn_model, as.matrix(test_features))
ggplot(data.frame(pred = as.numeric(test_predictions), varroa = test_labels$varroa_per_300_bees1)) +
  geom_point(aes(x = pred, y = varroa)) +
  geom_abline(intercept = 0, slope = 1, color = "blue") +
  geom_smooth(aes(x = pred, y = varroa), method = "lm", color = "red", se = FALSE)
# Error distribution
qplot(test_predictions - test_labels$varroa_per_300_bees1, geom = "density")
