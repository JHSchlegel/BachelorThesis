#=================================#
#### Forecasting Portfolio VaR ####
#=================================#

portfolio.plret.df <- data.frame(Date = stocks.plret.df$Date, Portfolio = apply(stocks.plret.df[,-1], 1, mean))

# Convert to time series
portfolio.plret.ts <- xts(portfolio.plret.df$plret, order.by = as.Date(portfolio.plret.df$Date))

stocks.plret.ts <- xts(stocks.plret.df[,-1], order.by = as.Date(stocks.plret.df$Date))
stocks.ts

