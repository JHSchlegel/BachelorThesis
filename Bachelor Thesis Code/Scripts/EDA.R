#=================================#
#### Exploratory Data Analysis ####
#=================================#

library(tidyverse)
library(yfR)
library(zoo)
library(psych)
library(tseries)
library(broom)

# TODO: check whether adjusted prices or not
# TODO: tsplot; especially acf plot w/ lines for absolute and "normal"returns
# TODO: check for multivariate normality



# Load Factors and Create Factors Dataframe -------------------------------


FFFactors <- read.csv("./Data/FamaFrenchFactorsDaily.csv", header = TRUE)
View(FFFactors)
MomFactor <- read.csv("./Data/MomentumFactorDaily.csv", header = TRUE)
View(MomFactor)

# TODO: include this section or not?
# Check whether the dates are the same for the period 2nd January 2001 to 30th December 2011
Date.fff.all <- as.Date(FFFactors.df$X, format = "%Y%m%d") 
Date.fff <- Date.fff.all[(Date.fff.all >= "2001-01-02") & (Date.fff.all <= "2011-12-30")]
Date.mom.all <- as.Date(MomFactor.df$X, format = "%Y%m%d")
Date.mom <- Date.mom.all[(Date.mom.all >= "2001-01-02") & (Date.mom.all <= "2011-12-30")]
all.equal(Date.fff, Date.mom)
# All dates are the same for this period

## Construct factor return dataframe for the period in question
FFFactors.df <- FFFactors %>% 
  mutate(Date = X) %>% 
  select(-X) %>% 
  relocate(Date, .before = "Mkt.RF") %>% 
  relocate(RF, .before = "Mkt.RF") %>% 
  mutate(Date = as.Date(Date, format = "%Y%m%d")) %>%
  filter((Date >= "2001-01-02") & (Date <= "2011-12-30"))
View(FFFactors.df)

MomFactor.df <- MomFactor %>% 
  mutate(Date = X) %>% 
  select(-X) %>% 
  relocate(Date, .before = "Mom") %>% 
  mutate(Date = as.Date(Date, format = "%Y%m%d")) %>% 
  filter((Date >= "2001-01-02") & (Date <= "2011-12-30"))
View(MomFactor.df)

## one dataframe for all factors
FFCFactors.df <- data.frame(FFFactors.df, Mom = MomFactor.df$Mom)
View(FFCFactors.df)
## Save newly created dataframe as csv file to allow for easy importing
write.csv(FFCFactors.df, "Data\\FFCFactors.csv", row.names = FALSE)


# Load Stock Percentage Log Returns ---------------------------------------


stocks.plret.df <- read.csv("./Data/StockPlrets.csv", header = TRUE)
View(stocks.plret.df)

## The following (commented out) section shows how this csv file has been created


# Create Stock Percentage Log Returns Dataframe with Yahoo Finance Data ---------------------------


# import log returns from yahoo;
# UTC now is called RTX due to fusion in April of 2020
# tickers <- c("KO", "XOM", "GE", "IBM", "CVX", "RTX", "PG", "CAT", "BA", "MRK")
# stock_data <- yf_get(tickers, first_date = "2000-12-29", last_date = "2011-12-30", 
#                      freq_data = "daily",  type_return = "log")
# head(stock_data)
# ## Only the first log return for each of the stocks is NA
# ## i.e. our actual log returns begin on the 2nd of January 2001
# sum(is.na(stock_data$ret_adjusted_prices))
# 
# ## Bind and assign the stock name, date and log return of adjusted prices to the stock name
# for (ticker_nr in 1:10){
#   assign(tickers[ticker_nr], stock_data %>% filter(ticker == tickers[ticker_nr]) %>% 
#            select(ticker, ref_date, ret_adjusted_prices) %>% na.omit())
# }
# ## View PG as an example
# View(PG)
# 
# ## Create dataframe with the log returns of all stocks
# stocks.lret.df <- data.frame(Date = KO$ref_date, KO = KO$ret_adjusted_prices, XOM = XOM$ret_adjusted_prices, 
#                              GE = GE$ret_adjusted_prices, IBM = IBM$ret_adjusted_prices, CVX = CVX$ret_adjusted_prices, 
#                              UTC = RTX$ret_adjusted_prices, PG = PG$ret_adjusted_prices, CAT = CAT$ret_adjusted_prices,
#                              BA = BA$ret_adjusted_prices, MRK = MRK$ret_adjusted_prices)
# View(stocks.lret.df)
# 
# ## Multiply all columns except Date by 100 to get percentage log returns
# stocks.plret.df <- data.frame(Date = KO$ref_date, apply(stocks.lret.df[,-1], 2, function(x)100*x))
# View(stocks.plret.df)
# 
# ## Save the created dataframe as a csv file to allow for easy importing and guarantee reproducibility
# write.csv(stocks.plret.df, "Data\\StockPlrets.csv", row.names = FALSE)


# Calculate Portfolio Percentage Log Returns of Equally Weighted Portfolio ---------------


portfolio.weights <- rep(1/10, 10)
portfolio.plret.df <- data.frame(Date = stocks.plret.df$Date, 
                                 Portfolio = apply(stocks.plret.df[,-1], 1, function(x)x%*%portfolio.weights))
View(portfolio.plret.df)
## Alternatively: portfolio.plret.df <- data.frame(Date = stocks.plret.df$Date, Portfolio = apply(stocks.plret.df[,-1], 1, mean))

## Visualise portfolio percentage log-returns
par(mfrow=c(1,2))
plot(as.Date(portfolio.plret.df$Date), portfolio.plret.df$Portfolio, type = "h", 
     main = "Daily Percentage Log-Returns of Portfolio", xlab = "Date", ylab = "Percentage Log-Returns")
abline(h = 0, col = "grey")
# periods of higher and lower volatility clearly visible
hist(portfolio.plret.df$Portfolio, breaks = 50, prob = TRUE, xlab = "Percentage Log-Returns", 
     main = "Histogram of Daily Percentage Log-Returns", col  ="grey95")
lines(density(portfolio.plret.df$Portfolio), col = "black", lty = 1)
# non-normal with very extreme values on either side


# Calculate Summary Statistics --------------------------------------------


#' Summary Statistics
#' @param dataframe Dataframe for which we want to calculate the summary statistics
#' @param multiple.rets Boolean whether there are multiple return columns or not
summary.statistics <- function(dataframe, multiple.rets =  TRUE){
  descr.stats <- psych::describe(dataframe[,-1])
  ## Conduct Jarque-Bera-Test
  ifelse(multiple.rets == TRUE,
         JB.test <- apply(dataframe[,-1], 2, FUN = jarque.bera.test),
         JB.test <- jarque.bera.test(dataframe[,-1]))
  ## Extract Jarque-Bera test statistic
  ifelse(multiple.rets == TRUE, 
         JB.test.stat <- matrix(lapply(JB.test, function(x)return(x[1]))),
         JB.test.stat <- matrix(JB.test[1]))
  tab <- data.frame(Mean = descr.stats[,3], Median = descr.stats[,5], SD = descr.stats[,4], MAD = descr.stats[,7],
                    Min = descr.stats[,8], Max = descr.stats[,9], Skew = descr.stats[,11], Kurt = descr.stats[,12], JB = JB.test.stat)
  tab <- data.frame(Stock = colnames(dataframe)[-1], tab)
  return(tab)
}
summary.statistics(stocks.plret.df)
summary.statistics(FFCFactors.df)
summary.statistics(portfolio.plret.df, multiple.rets = FALSE)




# ACF Plots Factors -------------------------------------------------------

par(mfrow = c(2,2),  cex = 0.8, oma = c(0.2, 0.2, 0.2, 0.2))
names <- c("Market", "Size", "Value", "Momentum")
for (i in 1:4){
  factors_acf <- acf(FFCFactors.df[,i+2], plot = FALSE)
  abs_factors_acf <- acf(abs(FFCFactors.df[,i+2]), plot = FALSE)
  plot(factors_acf, type = "l", ci.col = "black", lty = 1, lwd = 1.3,
       xlab = "Lag", main = "", ylab = "ACF")
  title(names[i], line = 0.5)
  lines(x = 0:(length(abs_factors_acf$acf)-1), abs_factors_acf$acf, lty =3,
        lwd = 1.3)
  legend("topright", legend = c("returns", "absolute returns"), lty = c(1,4),
         bty = "n", lwd = c(1.3, 1.3))
}


acf(FFCFactors.df[,3:6], ci.col = "black") # little to no auto-& crosscorrelation between returns
acf(abs(FFCFactors.df[,3:6]), ci.col = "black") # high auto-& crosscorrelation between absolute returns



# ACF Plots Shares --------------------------------------------------------

par(mfrow = c(4,3))
for (i in 1:10){
  stocks_acf <- acf(stocks.plret.df[,i+1], plot = FALSE)
  abs_stocks_acf <- acf(abs(stocks.plret.df[,i+1]), plot = FALSE)
  plot(stocks_acf, type = "l", ci.col = "black", lty = 1, lwd = 1.3,
       xlab = "Lag", main = "", ylab = "ACF")
  title(colnames(stocks.plret.df)[i+1], line = 0.5)
  lines(x = 0:(length(abs_stocks_acf$acf)-1), abs_stocks_acf$acf, lty =3,
        lwd = 1.3)
  legend("topright", legend = c("returns", "absolute returns"), lty = c(1,4),
         bty = "n", lwd = c(1.3, 1.3))
}

acf(stocks.plret.df[,-1], ci.col = "black") # little to no auto-& crosscorrelation between returns
acf(abs(stocks.plret.df[,-1]), ci.col = "black") # high auto-& crosscorrelation between absolute returns


