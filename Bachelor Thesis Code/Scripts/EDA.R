#=============================================================================#
########################## Exploratory Data Analysis ##########################
#=============================================================================#

library(tidyverse)
library(yfR)
library(zoo)
library(psych)
library(tseries)

theme_set(theme_light())

# TODO: check for multivariate normality
# TODO: tsdisplay to check for stationarity


#---------------------------------------------------------------#
########### Load Factors and Create Factors Dataframe ###########
#---------------------------------------------------------------#

FFFactors <- read.csv("Data/FamaFrenchFactorsDaily.csv", header = TRUE)
View(FFFactors)
MomFactor <- read.csv("Data/MomentumFactorDaily.csv", header = TRUE)
View(MomFactor)


## Check whether the dates are the same for the period 2nd January 2001 to 30th December 2011
Date_fff_all <- as.Date(FFFactors$X, format = "%Y%m%d") 
Date_fff <- Date_fff_all[(Date_fff_all >= "2001-01-02") & (Date_fff_all <= "2011-12-30")]
Date_mom_all <- as.Date(MomFactor$X, format = "%Y%m%d")
Date_mom <- Date_mom_all[(Date_mom_all >= "2001-01-02") & (Date_mom_all <= "2011-12-30")]
all.equal(Date_fff, Date_mom)
# All dates are the same for this period

## Construct Fama-French factors return dataframe for the period in question
FFFactors_df <- FFFactors %>% 
  mutate(Date = X) %>% 
  select(-X) %>% 
  relocate(Date, .before = "Mkt.RF") %>% 
  relocate(RF, .before = "Mkt.RF") %>% 
  mutate(Date = as.Date(Date, format = "%Y%m%d")) %>%
  filter((Date >= "2001-01-02") & (Date <= "2011-12-30"))
View(FFFactors_df)

## Construct momentum factor return dataframe for the period in question
MomFactor_df <- MomFactor %>% 
  mutate(Date = X) %>% 
  select(-X) %>% 
  relocate(Date, .before = "Mom") %>% 
  mutate(Date = as.Date(Date, format = "%Y%m%d")) %>% 
  filter((Date >= "2001-01-02") & (Date <= "2011-12-30"))
View(MomFactor_df)

## one dataframe for all factors
FFCFactors_df <- data.frame(FFFactors_df, Mom = MomFactor_df$Mom)
View(FFCFactors_df)
## Save newly created dataframe as csv file to allow for easy importing
write.csv(FFCFactors_df, "Data\\FFCFactors.csv", row.names = FALSE)



#-------------------------------------------------------#
########### Load Stock Percentage Log Returns ###########
#-------------------------------------------------------#


stocks_plret_df <- read.csv("Data/StockPlrets.csv", header = TRUE)
View(stocks_plret_df)

all.equal(FFCFactors_df$Date, stocks_plret_df$Date)
# all dates are the same for this period


## The following (commented out) section shows how this csv file has been created


#------------------------------------------------------------------------------------------------#
########### Creation of Stock Percentage Log Returns Dataframe with Yahoo Finance Data ###########
#------------------------------------------------------------------------------------------------#

# # For numerical reasons of stability, the GARCH Models will be fitted using percentage log returns
# # import log returns from yahoo;
# # UTC now is called RTX due to fusion in April of 2020
# tickers <- c("KO", "XOM", "GE", "IBM", "CVX", "RTX", "PG", "CAT", "BA", "MRK")
# stock_data_log <- yf_get(tickers, first_date = "2000-12-29", last_date = "2011-12-31",
#                      freq_data = "daily",  type_return = "log")
# head(stock_data_log)
# 
# ## Only the first log return for each of the stocks is NA
# ## i.e. our actual log returns begin on the 2nd of January 2001
# sum(is.na(stock_data_log$ret_adjusted_prices))
# 
# ## Bind and assign the stock name, date and log return of adjusted prices to the stock name
# for (ticker_nr in 1:10){
#   assign(tickers[ticker_nr], stock_data_log %>% filter(ticker == tickers[ticker_nr]) %>%
#            select(ticker, ref_date, ret_adjusted_prices) %>% na.omit())
# }
# ## View PG as an example
# View(PG)
# 
# 
# ## Create dataframe with the log returns of all stocks
# stocks_lret_df <- data.frame(Date = KO$ref_date, KO = KO$ret_adjusted_prices, XOM = XOM$ret_adjusted_prices,
#                              GE = GE$ret_adjusted_prices, IBM = IBM$ret_adjusted_prices, CVX = CVX$ret_adjusted_prices,
#                              UTC = RTX$ret_adjusted_prices, PG = PG$ret_adjusted_prices, CAT = CAT$ret_adjusted_prices,
#                              BA = BA$ret_adjusted_prices, MRK = MRK$ret_adjusted_prices)
# View(stocks_lret_df)
# 
# 
# 
# ## Multiply all columns except Date by 100 to get percentage log returns
# stocks_plret_df <- data.frame(Date = KO$ref_date, 100*stocks_lret_df[,-1])
# View(stocks_plret_df)
# 
# # ## Save the created dataframe as a csv file to allow for easy importing and guarantee reproducibility
# write.csv(stocks_plret_df, "Data\\StockPlrets.csv", row.names = FALSE)

#-----------------------------------------------------------------------------------------------#
########### Calculate Portfolio Percentage Log Returns of Equally Weighted Portfolio ###########
#-----------------------------------------------------------------------------------------------#


portfolio_weights <- rep(1/10, 10)
portfolio_plret_df <- data.frame(Date = stocks_plret_df$Date, 
                                 Portfolio = apply(stocks_plret_df[,-1], 1, function(x)x%*%portfolio_weights))
View(portfolio_plret_df)
## Alternatively: portfolio_plret_df <- data.frame(Date = stocks_plret_df$Date, Portfolio = apply(stocks_plret_df[,-1], 1, mean))

## Visualise portfolio percentage log returns
par(mfrow=c(1,2), oma = c(0.2, 0.2, 0.2, 0.2), mar = c(5.1, 4.1, 4.1, 2.1),
    cex.main = 0.9)
plot(as.Date(portfolio_plret_df$Date), portfolio_plret_df$Portfolio, type = "h", 
     main = "", xlab = "Date",
     ylab = "Percentage Log Returns")
title("Panel A: Daily Percentage Log Returns", line = 0.5)
abline(h = 0, col = "grey")
# periods of higher and lower volatility clearly visible
hist(portfolio_plret_df$Portfolio, breaks = 50,  xlab = "Percentage Log Returns", 
     main = "", col  ="grey95", xlim = c(-10, 12))
title("Panel B: Histogram of Daily Percentage Log Returns", line = 0.5)
# non-normal with very extreme values on either side

max_return_ind <- which.max(portfolio_plret_df$Portfolio)
portfolio_plret_df[max_return_ind, 1] #2008-10-28

min_return_ind <- which.min(portfolio_plret_df$Portfolio)
portfolio_plret_df[min_return_ind, 1] # 2008-10-15
# both the highest and lowest return were observed in October of 2008                

#--------------------------------------------------#
########### Calculate Summary Statistics ###########
#--------------------------------------------------#


#' Summary Statistics
#' @param dataframe Dataframe for which we want to calculate the summary statistics
#' @param multiple.rets Boolean whether there are multiple return columns or not
summary_statistics <- function(dataframe, multiple.rets =  TRUE){
  descr_stats <- psych::describe(dataframe[,-1])
  ## Conduct Jarque-Bera-Test
  ifelse(multiple.rets == TRUE,
         JB_test <- apply(dataframe[,-1], 2, FUN = jarque.bera.test),
         JB_test <- jarque.bera.test(dataframe[,-1]))
  ## Extract Jarque-Bera test statistic
  ifelse(multiple.rets == TRUE, 
         JB_test_stat <- matrix(unlist(lapply(JB_test, function(x)return(x[1])))),
         JB_test_stat <- matrix(unlist(JB_test[1])))
  tab <- data.frame(Mean = descr_stats[,3], Median = descr_stats[,5], SD = descr_stats[,4], MAD = descr_stats[,7],
                    Min = descr_stats[,8], Max = descr_stats[,9], Skew = descr_stats[,11], Kurt = descr_stats[,12], JB = JB_test_stat)
  tab <- data.frame(Stock = colnames(dataframe)[-1], tab %>% round(3))
  return(tab)
}
summary_statistics(stocks_plret_df)
summary_statistics(FFCFactors_df[,-2])
summary_statistics(portfolio_plret_df, multiple.rets = FALSE)



#---------------------------------------------------------------------------------------#
########### Examining Relationship and Correlation Between Factors/ Shares ##############
#---------------------------------------------------------------------------------------#
## Scatterplot Matrix
if (!require(GGally)) install.packages("GGally")

# function to add smoother to scatterplot matrix:
add_smooth <- function(dataframe, mapping, method = "loess", ...){
  ggplot(data = dataframe, mapping = mapping)+
  geom_point(alpha = 0.4)+
  geom_smooth(method = method, ...)
}
FFCFactors_df %>% 
  select(-Date, -RF) %>% 
  GGally::ggpairs(lower = list(continuous = add_smooth))

stocks_plret_df %>% 
  select(-Date) %>% 
  GGally::ggpairs(lower = list(continuous = add_smooth))

## Correlation Heatmap
if (!require(pheatmap)) install.packages("pheatmap")
FFCFactors_df %>% 
  select(-Date, -RF) %>%
  cor() %>% 
  pheatmap(display_numbers = TRUE)
  # round(digits = 3) %>% 
  # reshape2::melt() %>% 
  # ggplot(aes(Var1, Var2, fill = value))+
  # geom_tile()+
  # theme(panel.grid.major = element_blank(), panel.border = element_blank(),
  #       axis.ticks = element_blank())+
  # scale_fill_gradient2(high = "blue", low = "red", mid = "white", midpoint = 0,
  #                      name = "Correlation", limit = c(-1,1))+
  # geom_text(aes(Var1, Var2, label = value), color = "black")+
  # labs(x = "", y = "", title  ="Correlation Heatmap of Factors")

# momentum and market show strongest correlation w/ -0.465, others barely correlated

stocks_plret_df %>% 
  select(-Date) %>% 
  cor() %>% 
  pheatmap(display_numbers = TRUE)
  # round(digits = 3) %>% 
  # reshape2::melt() %>%
  # ggplot(aes(Var1, Var2, fill = value))+
  # geom_tile()+
  # theme(panel.grid.major = element_blank(), panel.border = element_blank(),
  #       axis.ticks = element_blank())+
  # scale_fill_gradient2(high = "blue", low = "red", mid = "white", midpoint = 0,
  #                      name = "Correlation", limit = c(-1,1))+
  # geom_text(aes(Var1, Var2, label = value), color = "black")+
  # labs(x = "", y = "", title  ="Correlation Heatmap of Shares")
  
# all positively correlated; correlations between shares higher than correlations between factors


## Network Plot of Correlation
# highly correlated variables are clustered together in network plots:
if (!require(corrr)) install.packages("corrr")
FFCFactors_df %>% 
  select(-Date, -RF) %>% 
  correlate() %>% 
  rearrange() %>% 
  network_plot(min_cor = 0, colors = c("blue", "white", "red"))
# market and momentum show strongest linear relationship
# SMB almost uncorrelated to other factors
# very slight correlation between Mom & HML and Mkt.RF & HML

stocks_plret_df %>% 
  select(-Date) %>% 
  correlate() %>% 
  rearrange() %>% 
  network_plot(min_cor = 0, colors = c("blue", "white", "red"))
# all shares positively correlated with each other
# can identify some clusters based on correlation:
# GE, UTC, IBM, BA and CAT tend to be highly correlated; all technology
# CVX and XOM are clustered; both energy
# PG and KO are clustered; both consumption goods
# MRK does not belong to any of the other clusters; MRK is pharma
# returns show strong linear dependencies within a certain sector and slightly lower
# correlation between sectors


#------------------------------------------#
########### ACF Plots Factors ##############
#------------------------------------------#

par(mfrow = c(2,2),  cex = 0.8, oma = c(0.2, 0.2, 0.2, 0.2),
    mar = c(2, 4.1, 4.1, 2.1))
names <- c("Panel A: Market Factor", "Panel B: Size Factor", 
           "Panel C: Value Factor", "Panel D: Momentum Factor")
for (i in 1:4){
  factors_acf <- acf(FFCFactors_df[,i+2], plot = FALSE)
  abs_factors_acf <- acf(abs(FFCFactors_df[,i+2]), plot = FALSE)
  plot(factors_acf, type = "l", ci.col = "black", lty = 1, lwd = 1.3,
       xlab = "Lag", main = "", ylab = "ACF")
  title(names[i], line = 0.5)
  lines(x = 0:(length(abs_factors_acf$acf)-1), abs_factors_acf$acf, lty =3,
        lwd = 1.3)
  legend("topright", legend = c("returns", "absolute returns"), lty = c(1,4),
         bty = "n", lwd = c(1.3, 1.3))
}


acf(FFCFactors_df[,3:6], ci.col = "black") # little to no auto-& crosscorrelation between returns
acf(abs(FFCFactors_df[,3:6]), ci.col = "black") # high auto-& crosscorrelation between absolute returns


#-------------------------------------------------------#
########### ACF Plots Portfolio and Shares ##############
#-------------------------------------------------------#


par(mfrow = c(4,3), cex = 0.8, oma = c(0.2, 0.2, 0.2, 0.2),
    mar = c(1, 4.1, 4.1, 1))

pf_acf <- acf(portfolio_plret_df$Portfolio, plot = FALSE)
abs_pf_acf <- acf(abs(portfolio_plret_df$Portfolio), plot = FALSE)
plot(pf_acf, type = "l", ci.col = "black", lty = 1, lwd = 1.3,
     xlab = "Lag", main = "", ylab = "ACF")
title("Panel A: Portfolio", line = 0.5)
lines(x = 0:(length(abs_pf_acf$acf)-1), abs_pf_acf$acf, lty =3,
      lwd = 1.3)
legend("topright", legend = c("returns", "absolute returns"), lty = c(1,4),
       bty = "n", lwd = c(1.3, 1.3))

for (i in 1:10){
  stocks_acf <- acf(stocks_plret_df[,i+1], plot = FALSE)
  abs_stocks_acf <- acf(abs(stocks_plret_df[,i+1]), plot = FALSE)
  plot(stocks_acf, type = "l", ci.col = "black", lty = 1, lwd = 1.3,
       xlab = "Lag", main = "", ylab = "ACF")
  title(paste0("Panel ", LETTERS[i+1], ": ", colnames(stocks_plret_df)[i+1]), line = 0.5)
  lines(x = 0:(length(abs_stocks_acf$acf)-1), abs_stocks_acf$acf, lty =3,
        lwd = 1.3)
  legend("topright", legend = c("returns", "absolute returns"), lty = c(1,4),
         bty = "n", lwd = c(1.3, 1.3))
}

acf(stocks_plret_df[,-1], ci.col = "black") # little to no auto-& crosscorrelation between returns
acf(abs(stocks_plret_df[,-1]), ci.col = "black") # high auto-& crosscorrelation between absolute returns







