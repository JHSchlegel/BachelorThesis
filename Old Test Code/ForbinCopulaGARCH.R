#TODO: 1) OLS regression

#TODO: 2) bootstrap OLS residuals

#TODO: 3) estimate distribution F of factor returns via Copula GARCH

#TODO: 4) simulate rt from F& errort from Fe

#TODO: 5) calculate stock returns w/FFC

#TODO:6) calculate portfolio return

#TODO: VaR_t^p=-percentile(rpf,100p)

# TODO: DCC which distribution? how to backtransform? 


library(rugarch)
library(rmgarch)
library(tidyverse)
library(xts)
library(parallel)


#--------------------------------------------------------#
########### Import Data and Create xts Objects ###########
#--------------------------------------------------------#

## Import data from csv files created in EDA.R
FFCFactors_df <- read.csv("./Data/FFCFactors.csv", header = TRUE)
stocks_plret_df <- read.csv("./Data/StockPlrets.csv", header = TRUE)
portfolio_plret_df <- data.frame(Date = stocks_plret_df$Date, Portfolio = apply(stocks_plret_df[,-1], 1, mean))

## Create dataframe with all stocks and market factors
joined_df <- stocks_plret_df %>% 
  full_join(FFCFactors_df)

head(joined_df)

## Convert to time series
FFCFactors_ts <- xts(FFCFactors_df[,-1], order.by = as.Date(FFCFactors_df$Date))
stocks_plret_ts <- xts(stocks_plret_df[,-1], order.by = as.Date(stocks_plret_df$Date))
portfolio_plret_ts <- xts(portfolio_plret_df[,-1], order.by = as.Date(portfolio_plret_df$Date))


#-----------------------------------------------------------------------------#
########### Coefficients and Residuals of Carhart Four-Factor Model ###########
#-----------------------------------------------------------------------------#

n_dates <- dim(FFCFactors_df)[1]
coefs_mat <- matrix(0L, nrow = 10, ncol = 5) # empty matrix to store coefficients
error_mat <- matrix(0L, nrow = n_dates, ncol = 10) # empty matrix to store residuals
for (i in 1:10){
  # columns 2-11 are the shares
  fit <- lm((joined_df[,i+1]-RF) ~ Mkt.RF + SMB+ HML + Mom,data = joined_df) 
  coefs_mat[i,] <- coef(fit)
  error_mat[,i] <- resid(fit)
}
coefs_df <- data.frame(coefs_mat); error_df <- data.frame(Date = joined_df$Date, error_mat)
rownames(coefs_df) <- c("KO", "XOM", "GE", "IBM", "CVX", "RTX", "PG", "CAT", "BA", "MRK")
colnames(error_df) <- c("Date", "KO", "XOM", "GE", "IBM", "CVX", "RTX", "PG", "CAT", "BA", "MRK")
colnames(coefs_df) <- c("Intercept", "Market", "Size", "Value", "Momentum")
head(coefs_df); head(error_df)

# Each row in error_vec_resampled is one draw of a vector of error terms
# to make sure dependencies are kept, all elements in a row are from the same time t
set.seed(42)
N_boot <- 200000
bootind <- sample.int(n_dates, size = N_boot, replace = TRUE)
head(error_df[bootind,])
error_vec_resampled <- error_df[bootind,-1] # now without date

# Plot bootstrapped error distributions
error_vec_resampled %>% 
  gather() %>% 
  ggplot(aes(value)) +
  geom_histogram(bins = 1000, color = "white")+
  facet_wrap(~key, scale = "free_x")
