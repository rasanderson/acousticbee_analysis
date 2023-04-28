# Simple regression-based DL model
#
# Deep learning models that are regression based
rm(list = ls())


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
data_subset$grid_easting2 <- data_subset$grid_easting^2
data_subset$grid_northing2 <- data_subset$grid_northing^2

# Begin with just the "meta-data" about the environment ####
# Work on the Varroa as more variable response
# However, no records at Thorganby as not assessed. Remove missing values.
# Note there are 175 NA for station_easting/northing, but 320 for grid_easting/n
data_subset <- filter(data_subset, site != "Thorganby") %>% 
  drop_na(grid_easting) %>% 
  drop_na(grid_northing) %>% 
  drop_na(rain_mm_h_mean) %>% 
  drop_na(windsp) %>% 
  drop_na(temperature) 

# Split data first, so meta- and audio- using same data points
library(keras)
library(rsample)
split <- initial_split(data_subset, 0.8)
split_train_dataset <- training(split)
split_test_dataset <- testing(split)

meta_variables <- c("varroa_per_300_bees1", "day_of_year", "cos_hour", "sin_hour",
                    "grid_easting", "grid_northing", "rain_mm_h_mean", "windsp",
                    "temperature", "grid_easting2", "grid_northing2")
meta_train_dataset <- select(split_train_dataset, all_of(meta_variables))
meta_test_dataset <- select(split_test_dataset, all_of(meta_variables))

# Now split features from labels
meta_train_features <- meta_train_dataset %>% select(-varroa_per_300_bees1)
meta_test_features <- meta_test_dataset %>% select(-varroa_per_300_bees1)

meta_train_labels <- meta_train_dataset %>% select(varroa_per_300_bees1)
meta_test_labels <- meta_test_dataset %>% select(varroa_per_300_bees1)

# Normalise the data
meta_normalizer <- layer_normalization(axis = -1L)
# The adapt() function fits the state of the normalized data and fixes it so
# it doesn't change subsequently
meta_normalizer %>% adapt(as.matrix(meta_train_features))
print(meta_normalizer$mean)
first <- as.matrix(meta_train_features[1,])
cat('First example:', first)
cat('Normalized:', as.matrix(meta_normalizer(first)))



# Deep learning model for acoustic indices ####
acoustic_variables <- c("varroa_per_300_bees1", "ACI", "RMS", "H", "mean_dfreq",
                        "M", "Ht", "freq_mean", "sd", "freq_median", "sem",
                        "freq_mode", "Q25", "Q75", "IQR", "cent", "skewness",
                        "kurtosis", "sfm", "sh", "prec", "time.P1", "time.M",
                        "time.P2", "time.IPR", "freq.P1", "freq.M", "freq.P2",
                        "freq.IPR")
acoustic_train_dataset <- select(split_train_dataset, all_of(acoustic_variables))
acoustic_test_dataset <- select(split_test_dataset, all_of(acoustic_variables))
# # Smoothed freq histograms and scatterplots (slow)
# train_dataset %>%
#   select(day_of_year, cos_hour, sin_hour, temperature) %>%
#   GGally::ggpairs()
# # Range etc.
# skimr::skim(meta_subset)

# Now split features from labels
acoustic_train_features <- acoustic_train_dataset %>% select(-varroa_per_300_bees1)
acoustic_test_features <- acoustic_test_dataset %>% select(-varroa_per_300_bees1)

acoustic_train_labels <- acoustic_train_dataset %>% select(varroa_per_300_bees1)
acoustic_test_labels <- acoustic_test_dataset %>% select(varroa_per_300_bees1)

# Normalise the data
acoustic_normalizer <- layer_normalization(axis = -1L)
# The adapt() function fits the state of the normalized data and fixes it so
# it doesn't change subsequently
acoustic_normalizer %>% adapt(as.matrix(acoustic_train_features))
print(acoustic_normalizer$mean)
first <- as.matrix(acoustic_train_features[1,])
cat('First example:', first)
cat('Normalized:', as.matrix(acoustic_normalizer(first)))

# For simplicity, just scale the features
acoustic_train_features_scale <- scale(acoustic_train_features)
meta_train_features_scale     <- scale(meta_train_features)

# As relatively simple, use a single function to define and compile model
meta_dnn_model <- layer_input(shape = 10, name = "meta_dnn") %>% 
  layer_dense(64, activation = 'relu', name = "meta1") %>% 
  layer_dropout(0.25) %>% 
  layer_dense(10, activation = 'relu', name = "meta1") 

# As relatively simple, use a single function to define and compile model
acoustic_dnn_model <- layer_input(shape = 28, name = "acoustic_dnn") %>%
    layer_dense(64, activation = 'relu', name = "acoustic1") %>% 
    layer_dropout(0.25) %>% 
    layer_dense(28, activation = 'relu', name = "acoustic1") 

merge_inputs <-
  layer_concatenate(list(meta_dnn_model, acoustic_dnn_model), name = "conc") %>% 
  layer_dropout(0.25) %>% 
  layer_dense(38, activation = 'relu', name = "conc_dense")

merge_outputs <- merge_inputs %>% 
  layer_dense(1, name = "output")

combined_dnn <- keras_model(
  inputs = list(meta_dnn_model, acoustic_dnn_model),
  outputs = merge_outputs
)

combined_dnn %>% compile(
  loss = 'mean_absolute_error',
  optimizer = optimizer_adam(0.001),
)

summary(combined_dnn)
plot(combined_dnn, show_shapes = TRUE, show_layer_names = TRUE, expand_nested = TRUE)

# Outputs for target are same for both meta- and acoustic- datasets
output_targets <- acoustic_train_labels

# Now fit the combined model
combined_dnn_history <- combined_dnn %>% fit(
  x = list(as.matrix(meta_train_features_scale), as.matrix(acoustic_train_features_scale)),
  y = list(output_targets),
  validation_split = 0.2,
  verbose = 1,
  epochs = 50
)
plot(combined_dnn_history)


# Results on test dataset
# labels are identical for meta- and acoustic test datasets
combined_test_labels <- acoustic_test_labels

combined_test_results <- list()
combined_test_results[['dnn_model']] <- combined_dnn %>% evaluate(
  x = list(as.matrix(scale(meta_test_features)), as.matrix(scale(acoustic_test_features))),
  y = as.matrix(combined_test_labels),
  verbose = 1
)


sapply(combined_test_results, function(x) x)
# Make some predictions
combined_test_predictions <- predict(combined_dnn, list(as.matrix(scale(meta_test_features)), as.matrix(scale(acoustic_test_features))))
test_obs_preds <- data.frame(pred = as.numeric(combined_test_predictions), varroa = combined_test_labels$varroa_per_300_bees1)
ggplot(test_obs_preds) +
  geom_point(aes(x = pred, y = varroa)) +
  geom_abline(intercept = 0, slope = 1, color = "blue") +
  geom_smooth(aes(x = pred, y = varroa), method = "lm", color = "red", se = FALSE)
# Error distribution
ggplot(test_obs_preds, aes(x = pred - varroa)) +
  geom_density()
