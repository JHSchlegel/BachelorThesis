#=================================#
#### Forecasting Portfolio VaR ####
#=================================#


library(rugarch)
library(rmgarch)
library(xts)
library(parallel)



# Import Data and Create xts Objects --------------------------------------


FFCFactors.df <- read.csv("./Data/FFCFactors.csv", header = TRUE)
stocks.plret.df <- read.csv("./Data/StockPlrets.csv", header = TRUE)
portfolio.plret.df <- data.frame(Date = stocks.plret.df$Date, Portfolio = apply(stocks.plret.df[,-1], 1, mean))

## Convert to time series
FFCFactors.ts <- xts(FFCFactors.df[,-1], order.by = as.Date(FFCFactors.df$Date))
stocks.plret.ts <- xts(stocks.plret.df[,-1], order.by = as.Date(stocks.plret.df$Date))
portfolio.plret.ts <- xts(portfolio.plret.df[,-1], order.by = as.Date(portfolio.plret.df$Date))

## Detect how many cores are available for parallelizing
numcores <- detectCores()


# Univariate Normal GARCH(1,1) --------------------------------------------


cl <- makePSOCKcluster(numcores)
uni_normal.spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), mean.model = list(armaOrder=c(0,0), include.mean=TRUE),
                              distribution.model = "norm")
uni_normal_GARCH <- ugarchroll(uni_normal.spec, portfolio.plret.ts, n.start = 1000, refit.every = 1,
                       refit.window = "moving", solver = "hybrid",
                       calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05),
                       cluster = cl, keep.coef = FALSE)

stopCluster(cl)

uni_normal_GARCH.VaR <- data.frame(Date = portfolio.plret.df$Date[-c(1:1000)], alpha.1perc = uni_normal_GARCH@forecast$VaR[,1],
                                   alpha.5perc = uni_normal_GARCH@forecast$VaR[,2])
head(uni_normal_GARCH.VaR)
write.csv(uni_normal_GARCH.VaR, "Data\\Uni_Normal_GARCH.VaR.csv", row.names = FALSE)



# Univariate EWMA ---------------------------------------------------------


lambda <- 0.94
cl <- makePSOCKcluster(numcores)
ewma.spec <- ugarchspec(variance.model=list(model="iGARCH", garchOrder=c(1,1)), 
                        mean.model=list(armaOrder=c(0,0), include.mean=TRUE),  
                        distribution.model="norm", fixed.pars=list(omega=0, alpha = 1-lambda))
ewma_GARCH <- ugarchroll(uni_normal.spec, portfolio.plret.ts, n.start = 1000, refit.every = 1,
                               refit.window = "moving", solver = "hybrid",
                               calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05),
                               cluster = cl, keep.coef = FALSE)
stopCluster(cl)
ewma_GARCH.VaR <- data.frame(Date = portfolio.plret.df$Date[-c(1:1000)], alpha.1perc = ewma_GARCH@forecast$VaR[,1],
                             alpha.5perc = ewma_GARCH@forecast$VaR[,2])
head(ewma_GARCH.VaR)
write.csv(ewma_GARCH.VaR, "Data\\EWMA_GARCH.VaR.csv", row.names = FALSE)


# Univariate t-GJR-GARCH(1,1) ---------------------------------------------



# Multivariate Normal DCC-GARCH(1,1) --------------------------------------


