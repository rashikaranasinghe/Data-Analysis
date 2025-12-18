###############################
######## SARIMA Modeling ######
###############################

# SARIMA Model contains 2 components 
#   1. ARIMA component
#   2. seasonal component 
# ARIMA(p,d,q), SARIMA(P,D,Q,S)
# p - auto regressive lags
# d - order of intergations
# q - moving average lags

# P - seasonal AR lags
# D - seasonal differences
# Q - searsonal MA lags
# S - length of seasonal cycle

# use SARIMA model when data contains seasonal fluctuations
# SARIMA take into account the seasonality in the data and make more accurate predictions for each month

# Observe ACF and PACF of the data to determine whether to use SARIMA model or not. If data shows a cyclic pattern then SARIMA model is the most appropriate one.


# set the workind directory
setwd("/Users/rashikaranasinghe/Documents/My_website/SARIMA Modeling")

# Install packages if not already installed
# remotes::install_github("hrbrmstr/cdcfluview")

# Load libraries
library(cdcfluview)
library(tidyverse)
library(lubridate)
library(forecast)
library(ggplot2)

###### Download publicly available respiratory infection data
# CDC FluView provides weekly ILI (Influenza-Like Illness) data for the US.
ili_data <- ilinet(region = "national", years = 2015:2024)
dim(ili_data)# [1] 521  16
head(ili_data)

###### Prepare time-series data
ili_ts_data <- ili_data %>%
  select(week_start, weighted_ili) %>%
  arrange(week_start)
head(ili_ts_data)
dim(ili_ts_data)

# Convert to time series object
ili_ts <- ts(
  ili_ts_data$weighted_ili,
  start = c(year(min(ili_ts_data$week_start)), 
            week(min(ili_ts_data$week_start))),
  frequency = 52)



####### Exploratory visualization
quartz()
plot(ili_ts, lwd=2, col=2)

## OR ###

autoplot(ili_ts) +
  labs(
    title = "Weekly Respiratory Infection Trends (ILI Proxy)",
    x = "Year",
    y = "ILI (%)"
  ) +
  theme_minimal()


## Test whether the assumptions are met for the SARIMA modle
# SARIMA models assume the data are stationary (or can be made stationary)

# Do Augmented Dickey–Fuller (ADF) test for that:
# Null hypothesis (H₀): The series has a unit root → it is non-stationary (data wander over time and don’t return to a stable level)
# Alternative hypothesis (H₁): The series does not have a unit root → it is stationary (data stay around a stable level and don’t drift over time)

adf.test(ili_ts, alternative = "stationary")
###########################
# Augmented Dickey-Fuller Test
# data:  ili_ts
# Dickey-Fuller = -5.9415, Lag order = 8, p-value = 0.01
# alternative hypothesis: stationary
## p-value < 0.05, reject the null hypothesis, so the data are stationary and stable enough for the model to work properly.
############################

### Decomposing a time series into its constituent components
# breaking the data into simpler parts to understand its structure.
# These components are:
#   Trend – the long-term upward or downward movement
#   Seasonality – repeating patterns (e.g., yearly or weekly, daily cycles)
#   Irregular / noise/ residual – random fluctuations or outliars
# 
# By separating these pieces, you can see what part of the data is trend, what part is seasonal, and what is random. SARIMA then models the remaining structure (after accounting for trend and seasonality) to make better forecasts.

# This is done by using 
#   Seasonal-Trend decomposition - LOESS(STL) [LOESS - smooth local regression method]
#   seasonal decomposition of time series (STL)

dec <- decompose(ili_ts, type = "additive")
plot(dec, col = 2, lwd = 3)
acf(ili_ts, main="ACF", lwd = 2, col=2)
pacf(ili_ts, main="PACF", lwd = 3, col=3)
# see a cycling pattern - SARIMA model can be used

####### Stationarity check
ndiffs(ili_ts)     # Non-seasonal differencing
nsdiffs(ili_ts)   # Seasonal differencing

###### Fit SARIMA model
# SARIMA models were used because they effectively capture the strong seasonal patterns, temporal autocorrelation, and non-stationarity inherent in respiratory infection surveillance data, while providing interpretable forecasts with quantified uncertainty.s
sarima_model <- auto.arima(
  ili_ts,
  seasonal = TRUE,
  trace = TRUE,
  stepwise = FALSE,
  approximation = FALSE, # Run all the combinations
  ic="aic", test = "kpss")

sarima_model <- auto.arima(
  ili_ts,
  trace = TRUE,
  approximation = FALSE, # Run all the combinations
  ic="aic", test = "kpss")

summary(sarima_model)
################
# Series: ili_ts 
# ARIMA(2,0,0)(0,1,1)[52] 
# (2,0,0) is the non-seasonal componendt (ARIMA(2,0,0))
# (0,1,1) represents the seasonal component (SARIMA(2,0,0)(0,1,1)[52])
# 
# Coefficients:
#   ar1     ar2     sma1
# 1.5058  -0.565  -0.5796
# s.e.  0.0380   0.038   0.0480
# 
# sigma^2 = 0.07864:  log likelihood = -79.99
# AIC=167.98   AICc=168.07   BIC=184.59
# 
# Training set error measures:
#   ME      RMSE      MAE          MPE
# Training set 0.01423669 0.2652079 0.146277 -0.002439618
# MAPE      MASE       ACF1
# Training set 6.003585 0.1684528 0.03111112
##################

arimaorder(sarima_model)
###########
# p         d         q         P         D         Q Frequency 
# 2         0         0         2         1         0        52 
##########

###### Model diagnostics 
# to check wether there is any leftover autocorrelation patterns in the data
#######  plot residuals #########
res <- residuals(sarima_model)
res
acf(res) # no spike pass the dashed line means the model is good
pacf(res) # same
# The dashed lines are significance limits.
# If no spikes cross the dashed lines, it means there is no significant autocorrelation left in the residuals.

####### Test it statistically ######
checkresiduals(sarima_model)
############
# Ljung-Box test
# Null hypothesis: Residuals are random (no remaining autocorrelation)
# data:  Residuals from ARIMA(2,0,0)(2,1,0)[52]
# Q* = 67.987, df = 100, p-value = 0.994
# Model df: 4.   Total lags used: 104
# p-value > 0.05, fail to reject the H0, so there is no evidence of leftover patterns in the residuals.The model has captured the trend and seasonality well.
############


####### Forecast future respiratory infection patterns
forecast_12_months <- forecast(
  sarima_model,
  h = 52 # 104 for 2 more years
)
forecast_12_months
plot(forecast_12_months)

####### Visualize forecasts
quartz()
autoplot(forecast_12_months) +
  labs(
    title = "SARIMA Forecast of Respiratory Infections",
    x = "Year",
    y = "Predicted ILI (%)"
  ) +
  theme_minimal()







