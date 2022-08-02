#################################
### Exploratory Data Analysis ###
#################################


# Don't forget to download newest version!
# check whether dates are the same for stocks and factor returns
library(tidyverse)
library(magrittr)
library(yfR)
library(zoo)
library(psych)
library(tseries)
library(broom)
library(xts)

# Load Factors
FFFactors.df <- read.csv("./Data/FamaFrenchFactorsDaily.csv", header = TRUE)
View(FFFactors.df)
MomFactor.df <- read.csv("./Data/MomentumFactorDaily.csv", header = TRUE)
View(MomFactor)
dim(FFFactors.df)
dim(MomFactor.df)

# Check whether the dates are the same for the period 2nd January 2001 to 30th December 2011
Date.fff.all <- as.Date(FFFactors.df$X, format = "%Y%m%d") 
Date.fff <- Date.fff.all[(Date.fff.all >= "2001-01-02") & (Date.fff.all <= "2011-12-30")]
Date.mom.all <- as.Date(MomFactor.df$X, format = "%Y%m%d")
Date.mom <- Date.mom.all[(Date.mom.all >= "2001-01-02") & (Date.mom.all <= "2011-12-30")]
all.equal(Date.fff, Date.mom)
# All dates are the same for the period in question

FFFactors.df <- FFFactors.df %<>% 
  rename(Date = X, MRP = Mkt.RF, SMB=SMB, HML=HML, RF=RF) %>% 
  mutate(Date = as.Date(Date, format = "%Y%m%d")) %>%
  filter((Date >= "2001-01-02") & (Date <= "2011-12-30"))
View(FFFactors.df)

MomFactor.df <- MomFactor.df %<>% 
  rename(Date = X, Mom = Mom) %>% 
  mutate(Date = as.Date(Date, format = "%Y%m%d")) %>% 
  filter((Date >= "2001-01-02") & (Date <= "2011-12-30"))
View(MomFactor.df)

FFCFactors.df <- data.frame(FFFactors.df, Mom = MomFactor.df$Mom)
View(FFCFactors.df)

########### RTX instead of UTC due to fusion in April of 2020
# Return from 2001-01-02 inclusive?
tickers <- c("KO", "XOM", "GE", "IBM", "CVX", "RTX", "PG", "CAT", "BA", "MRK")
stock_data <- yf_get(tickers, first_date = "2000-12-29", last_date = "2011-12-30", 
                     freq_data = "daily",  type_return = "log")
stock_data

# Only the first log return for each of the stocks is NA
# i.e. our actual log returns begin on the 1st of January 2001
sum(is.na(stock_data$ret_adjusted_prices))

for (ticker_nr in 1:10){
  assign(tickers[ticker_nr], stock_data %>% filter(ticker == tickers[ticker_nr]) %>% 
           select(ticker, ref_date, ret_adjusted_prices) %>% na.omit())
}
View(PG)

stocks.lret.df <- data.frame(Date = KO$ref_date, KO = KO$ret_adjusted_prices, XOM = XOM$ret_adjusted_prices, GE = GE$ret_adjusted_prices,
                        IBM = IBM$ret_adjusted_prices, CVX = CVX$ret_adjusted_prices, RTX = RTX$ret_adjusted_prices, PG = PG$ret_adjusted_prices,
                        CAT = CAT$ret_adjusted_prices, BA = BA$ret_adjusted_prices, MRK = MRK$ret_adjusted_prices)

stocks.lret.df

stocks.perclret.df <- data.frame(Date = KO$ref_date, apply(stocks.lret.df[,-1], 2, function(x)100*x))
stocks.perclret.df


summary.statistics <- function(dataframe, Datecol = TRUE){
  if (Datecol == TRUE){
    descr.stats <- psych::describe(dataframe[,-1])
    JB.test <- apply(dataframe[,-1], 2, FUN = jarque.bera.test)
  }
  else{
    descr.stats <- psych::describe(dataframe)
    JB.test <- apply(dataframe, 2, FUN = jarque.bera.test)
  }
  JB.test.stat <- matrix(lapply(JB.test, function(x)return(x[1])))
  tab <- data.frame(Mean = descr.stats[,3], Median = descr.stats[,5], SD = descr.stats[,4], MAD = descr.stats[,7],
                    Min = descr.stats[,8], Max = descr.stats[,9], Skew = descr.stats[,11], Kurt = descr.stats[,12], JB = JB.test.stat)
  rownames(tab) <- rownames(descr.stats)
  return(tab)
}
summary.statistics(stocks.perclret.df)
summary.statistics(FFCFactors.df)

# Calculate portfolio returns of equally weighted portfolio
portfolio.weights <- rep(1/10, 10)
portfolio.plret.df <- data.frame(Date = stocks.perclret.df$Date, plret = apply(stocks.perclret.df[,-1], 1, function(x)x%*%portfolio.weights))
# Alternatively: portfolio.perclret <- apply(stocks.perclret.df[,-1], 1, mean)

plot(portfolio.plret.df$Date, portfolio.plret.df$plret, type = "h", main = "Daily Percentage Log-Returns of Portfolio", xlab = "Date", ylab = "Percentage Log-Returns")
abline(h = 0, col = "grey")

# Convert to time series
portfolio.plret.ts <- xts(portfolio.plret.df$plret, order.by = portfolio.plret.df$Date)
plot(portfolio.plret.ts)

stocks.plret.ts <- xts(stocks.perclret.df[,-1], order.by = stocks.perclret.df$Date)
stocks.ts