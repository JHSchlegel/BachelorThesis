#TODO: 1) OLS regression

#TODO: 2) bootstrap OLS residuals

#TODO: 3) estimate distribution F of factor returns via Copula GARCH

#TODO: 4) simulate rt from F& epsilont from Fe

#TODO: 5) calculate stock returns w/FFC

#TODO:6) calculate portfolio return

#TODO: VaR_t^p=-percentile(rpf,100p)

# TODO: DCC which distribution? how to backtransform? 


library(rugarch)
library(rmgarch)
library(tidyverse)
library(xts)
library(parallel)



# Import Data and Create xts Objects --------------------------------------


FFCFactors_df <- read.csv("./Data/FFCFactors.csv", header = TRUE)
stocks_plret_df <- read.csv("./Data/StockPlrets.csv", header = TRUE)
portfolio_plret_df <- data.frame(Date = stocks_plret_df$Date, Portfolio = apply(stocks_plret_df[,-1], 1, mean))

## Create Dataframe with all stocks and market factors
joined_df <- stocks_plret_df %>% 
  full_join(FFCFactors_df)

## Convert to time series
FFCFactors_ts <- xts(FFCFactors_df[,-1], order.by = as.Date(FFCFactors_df$Date))
stocks_plret_ts <- xts(stocks_plret_df[,-1], order.by = as.Date(stocks_plret_df$Date))
portfolio_plret_ts <- xts(portfolio_plret_df[,-1], order.by = as.Date(portfolio_plret_df$Date))

## Coefficients of Carhart Model
