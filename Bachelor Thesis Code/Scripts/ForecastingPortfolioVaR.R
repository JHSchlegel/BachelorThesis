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

uni_normal.meanModel <- list(armaOrder=c(0,0), include.mean=TRUE)
uni_normal.varModel <- list(model = "sGARCH", garchOrder = c(1,1))
uni_normal.spec <- ugarchspec(variance.model = uni_normal.varModel, mean.model = uni_normal.meanModel,
                              distribution.model = "norm")
cl <- makePSOCKcluster(numcores)
uni_normal_GARCH.roll <- ugarchroll(uni_normal.spec, portfolio.plret.ts, n.start = 1000, refit.every = 1,
                       refit.window = "moving", solver = "hybrid",
                       calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05),
                       cluster = cl, keep.coef = FALSE)

stopCluster(cl)

uni_normal_GARCH_VaR <- data.frame(Date = portfolio.plret.df$Date[-c(1:1000)], alpha.1perc = uni_normal_GARCH.roll@forecast$VaR[,1],
                                   alpha.5perc = uni_normal_GARCH.roll@forecast$VaR[,2])
head(uni_normal_GARCH_VaR)
write.csv(uni_normal_GARCH_VaR, "Data\\VaR\\Uni_Normal_GARCH_VaR.csv", row.names = FALSE)



# Univariate EWMA ---------------------------------------------------------


lambda <- 0.94
ewma.meanModel <- list(armaOrder=c(0,0), include.mean=TRUE)
ewma.varModel <- list(model="iGARCH", garchOrder=c(1,1))
ewma.spec <- ugarchspec(variance.model= ewma.varModel, mean.model= ewma.meanModel,  
                        distribution.model="norm", fixed.pars=list(omega=0, alpha = 1-lambda))
cl <- makePSOCKcluster(numcores)
ewma_GARCH.roll <- ugarchroll(uni_normal.spec, portfolio.plret.ts, n.start = 1000, refit.every = 1,
                               refit.window = "moving", solver = "hybrid",
                               calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05),
                               cluster = cl, keep.coef = FALSE)
stopCluster(cl)
ewma_GARCH_VaR <- data.frame(Date = portfolio.plret.df$Date[-c(1:1000)], alpha.1perc = ewma_GARCH.roll@forecast$VaR[,1],
                             alpha.5perc = ewma_GARCH.roll@forecast$VaR[,2])
head(ewma_GARCH_VaR)
write.csv(ewma_GARCH_VaR, "Data\\VaR\\EWMA_GARCH_VaR.csv", row.names = FALSE)


# Univariate t-GJR-GARCH(1,1) ---------------------------------------------

uni_t_gjr.meanModel <- list(armaOrder = c(0,0), include.mean = TRUE)
uni_t_gjr.varModel <- list(model = "gjrGARCH", garchOrder = c(1,1))
uni_t_gjr.spec <- ugarchspec(variance.model = uni_t_gjr.varModel, mean.model = uni_t_gjr.meanModel, distribution.model = "std")
cl <- makePSOCKcluster(numcores)
uni_t_gjr_GARCH.roll <- ugarchroll(uni_normal.spec, portfolio.plret.ts, n.start = 1000, refit.every = 1,
                              refit.window = "moving", solver = "hybrid",
                              calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05),
                              cluster = cl, keep.coef = FALSE)
stopCluster(cl)
uni_t_gjr_GARCH_VaR <- data.frame(Date = portfolio.plret.df$Date[-c(1:1000)], alpha.1perc = uni_t_gjr_GARCH.roll@forecast$VaR[,1],
                             alpha.5perc = uni_t_gjr_GARCH.roll@forecast$VaR[,2])
head(uni_t_gjr_GARCH_VaR)
write.csv(uni_t_gjr_GARCH_VaR, "Data\\VaR\\Uni_t_GJR_GARCH.csv", row.names = FALSE)


# Univariate skew-t-NGARCH(1,1) -------------------------------------------


uni_skewt_NGARCH.meanModel <- list(armaOrder = c(0, 0), include.mean = TRUE)
uni_skewt_NGARCH.varModel <- list(model = "fGARCH", submodel = "NGARCH", garchOrder = c(1,1))
uni_skewt_NGARCH.spec <- ugarchspec(variance.model = uni_skewt_NGARCH.varModel, mean.model = uni_skewt_NGARCH.meanModel, 
                                    distribution.model = "sstd")
cl <- makePSOCKcluster(numcores)
uni_skewt_NGARCH.roll <- ugarchroll(uni_normal.spec, portfolio.plret.ts, n.start = 1000, refit.every = 1,
                              refit.window = "moving", solver = "hybrid",
                              calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05),
                              cluster = cl, keep.coef = FALSE)
stopCluster(cl)
uni_skewt_NGARCH_VaR <- data.frame(Date = portfolio.plret.df$Date[-c(1:1000)], alpha.1perc = uni_skewt_NGARCH.roll@forecast$VaR[,1],
                                  alpha.5perc = uni_skewt_NGARCH.roll@forecast$VaR[,2])
head(uni_skewt_NGARCH_VaR)
write.csv(uni_skewt_NGARCH_VaR, "Data\\VaR\\Uni_Skewt_NGARCH.csv", row.names = FALSE)


###################### Don't forget univariate mix normal garch
######################function for dataframe and head()





# Multivariate Normal DCC-GARCH(1,1) --------------------------------------


