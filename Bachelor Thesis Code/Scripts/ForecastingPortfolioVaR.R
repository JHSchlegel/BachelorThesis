#=============================================================================#
########################## Forecasting Portfolio VaR ##########################
#=============================================================================#

if (!require(parallel)) install.packages("parallel")
if (!require(rugarch)) install.packages("rugarch")
if (!require(rmgarch)) install.packages("rmgarch")
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(xts)) install.packages("xts")
if (!require(MSGARCH)) install.packages("MSGARCH")


#--------------------------------------------------------#
########### Import Data and Create xts Objects ###########
#--------------------------------------------------------#

FFCFactors_df <- read.csv("./Data/FFCFactors.csv", header = TRUE)
stocks_plret_df <- read.csv("./Data/StockPlrets.csv", header = TRUE)
portfolio_plret_df <- data.frame(Date = stocks_plret_df$Date, Portfolio =rowMeans(stocks_plret_df[,-1]))


## Convert to time series
FFCFactors_ts <- xts(FFCFactors_df[,-1], order.by = as.Date(FFCFactors_df$Date))
stocks_plret_ts <- xts(stocks_plret_df[,-1], order.by = as.Date(stocks_plret_df$Date))
portfolio_plret_ts <- xts(portfolio_plret_df[,-1], order.by = as.Date(portfolio_plret_df$Date))

## Detect how many cores are available for parallelizing
# subtract one so that not all cores are used
numcores <- detectCores()-1


#--------------------------------------------------#
########### Univariate Normal GARCH(1,1) ###########
#--------------------------------------------------#

uni_normal_meanModel <- list(armaOrder=c(0,0), include.mean=TRUE)
uni_normal_varModel <- list(model = "sGARCH", garchOrder = c(1,1))
uni_normal_spec <- ugarchspec(variance.model = uni_normal_varModel, mean.model = uni_normal_meanModel,
                              distribution.model = "norm")
cl <- makePSOCKcluster(numcores)
uni_normal_GARCH_roll <- ugarchroll(uni_normal_spec, portfolio_plret_ts, refit.every = 1,
                       refit.window = "moving", n.start = 1000, solver = "hybrid",
                       calculate_var = TRUE, VaR.alpha = c(0.01, 0.05),
                       cluster = cl, keep.coef = FALSE)

stopCluster(cl)

uni_normal_GARCH_VaR <- data.frame(Date = portfolio_plret_df$Date[-c(1:1000)], alpha_1perc = uni_normal_GARCH_roll@forecast$VaR[,1],
                                   alpha_5perc = uni_normal_GARCH_roll@forecast$VaR[,2])

write.csv(uni_normal_GARCH_VaR, "Data\\VaR\\Uni_Normal_GARCH_VaR.csv", row.names = FALSE)




#------------------------------------------------------#
########### Univariate EWMA w/ lambda = 0.94 ###########
#------------------------------------------------------#


lambda <- 0.94
ewma_meanModel <- list(armaOrder=c(0,0), include.mean=TRUE)
ewma_varModel <- list(model="iGARCH", garchOrder=c(1,1))
uni_ewma_spec <- ugarchspec(variance.model= ewma_varModel, mean.model= ewma_meanModel,  
                        distribution.model="norm", fixed.pars=list(omega=0, alpha = 1-lambda))
cl <- makePSOCKcluster(numcores)
uni_ewma_roll <- ugarchroll(uni_ewma_spec, portfolio_plret_ts, n.start = 1000, refit.every = 1,
                               refit.window = "moving", solver = "hybrid",
                               calculate_var = TRUE, VaR.alpha = c(0.01, 0.05),
                               cluster = cl, keep.coef = FALSE)
stopCluster(cl)
uni_ewma_VaR <- data.frame(Date = portfolio_plret_df$Date[-c(1:1000)], alpha_1perc = uni_ewma_roll@forecast$VaR[,1],
                             alpha_5perc = uni_ewma_roll@forecast$VaR[,2])
head(uni_ewma_VaR)
write.csv(uni_ewma_VaR, "Data\\VaR\\Uni_EWMA_VaR.csv", row.names = FALSE)


#-------------------------------------------------#
########### Univariate t-GJR-GARCH(1,1) ###########
#-------------------------------------------------#

uni_t_gjr_meanModel <- list(armaOrder = c(0,0), include.mean = TRUE)
uni_t_gjr_varModel <- list(model = "gjrGARCH", garchOrder = c(1,1))
uni_t_gjr_spec <- ugarchspec(variance.model = uni_t_gjr_varModel, mean.model = uni_t_gjr_meanModel, distribution.model = "std")
cl <- makePSOCKcluster(numcores)
uni_t_gjr_GARCH_roll <- ugarchroll(uni_t_gjr_spec, portfolio_plret_ts, n.start = 1000, refit.every = 1,
                              refit.window = "moving", solver = "hybrid",
                              calculate_var = TRUE, VaR.alpha = c(0.01, 0.05),
                              cluster = cl, keep.coef = FALSE)
stopCluster(cl)
uni_t_gjr_GARCH_VaR <- data.frame(Date = portfolio_plret_df$Date[-c(1:1000)], alpha_1perc = uni_t_gjr_GARCH_roll@forecast$VaR[,1],
                             alpha_5perc = uni_t_gjr_GARCH_roll@forecast$VaR[,2])
head(uni_t_gjr_GARCH_VaR)
write.csv(uni_t_gjr_GARCH_VaR, "Data\\VaR\\Uni_t_GJR_GARCH.csv", row.names = FALSE)


#-----------------------------------------------------#
########### Univariate skewt-GJR-GARCH(1,1) ###########
#-----------------------------------------------------#

uni_skewt_gjr_meanModel <- list(armaOrder = c(0,0), include.mean = TRUE)
uni_skewt_gjr_varModel <- list(model = "gjrGARCH", garchOrder = c(1,1))
uni_skewt_gjr_spec <- ugarchspec(variance.model = uni_skewt_gjr_varModel, mean.model = uni_skewt_gjr_meanModel, distribution.model = "sstd")
cl <- makePSOCKcluster(numcores)
uni_skewt_gjr_GARCH_roll <- ugarchroll(uni_skewt_gjr_spec, portfolio_plret_ts, n.start = 1000, refit.every = 1,
                                   refit.window = "moving", solver = "hybrid",
                                   calculate_var = TRUE, VaR.alpha = c(0.01, 0.05),
                                   cluster = cl, keep.coef = FALSE)
stopCluster(cl)
uni_skewt_gjr_GARCH_VaR <- data.frame(Date = portfolio_plret_df$Date[-c(1:1000)], alpha_1perc = uni_skewt_gjr_GARCH_roll@forecast$VaR[,1],
                                  alpha_5perc = uni_skewt_gjr_GARCH_roll@forecast$VaR[,2])
head(uni_skewt_gjr_GARCH_VaR)
write.csv(uni_skewt_gjr_GARCH_VaR, "Data\\VaR\\Uni_skewt_GJR_GARCH.csv", row.names = FALSE)


#--------------------------------------------------#
########### Univariate Skewt-NGARCH(1,1) ###########
#--------------------------------------------------#


uni_skewt_NGARCH_meanModel <- list(armaOrder = c(0, 0), include.mean = TRUE)
uni_skewt_NGARCH_varModel <- list(model = "fGARCH", submodel = "NGARCH", garchOrder = c(1,1))
uni_skewt_NGARCH_spec <- ugarchspec(variance.model = uni_skewt_NGARCH_varModel, mean.model = uni_skewt_NGARCH_meanModel, 
                                    distribution.model = "sstd")
cl <- makePSOCKcluster(numcores)
uni_skewt_NGARCH_roll <- ugarchroll(uni_skewt_NGARCH_spec, portfolio_plret_ts, n.start = 1000, refit.every = 1,
                              refit.window = "moving", solver = "hybrid",
                              calculate_var = TRUE, VaR.alpha = c(0.01, 0.05),
                              cluster = cl, keep.coef = FALSE)
stopCluster(cl)
uni_skewt_NGARCH_VaR <- data.frame(Date = portfolio_plret_df$Date[-c(1:1000)], alpha_1perc = uni_skewt_NGARCH_roll@forecast$VaR[,1],
                                  alpha_5perc = uni_skewt_NGARCH_roll@forecast$VaR[,2])
head(uni_skewt_NGARCH_VaR)
write.csv(uni_skewt_NGARCH_VaR, "Data\\VaR\\Uni_Skewt_NGARCH.csv", row.names = FALSE)


#--------------------------------------#
########### Mix-Normal GARCH ###########
#--------------------------------------#
window_size <- 1000
n_windows <- dim(portfolio_plret_df)[1]-window_size

ms2garch_VaR <- matrix(0L, nrow = n_windows, ncol = 2)
ms2garch_ES <- matrix(0L, nrow = n_windows, ncol = 2)


ms2garch_spec <- CreateSpec(variance.spec = list(model = "sGARCH"),
           distribution.spec = list(distribution = "norm"),
           switch.spec = list(K = 2))

for (i in 1:n_windows){
  fit_ml <- FitML(spec = ms2garch_spec, data = portfolio_plret_ts[i:(1000+i-1)])
  ms2garch_VaR[i,] <- Risk(fit_ml, alpha = c(0.01, 0.05), nahead = 1)$VaR
  ms2garch_ES[i,] <- Risk(fit_ml, alpha = c(0.01, 0.05), nahead = 1)$ES
  message("completed: ", i, " of ", n_windows)
}
ms2garch_VaR_df <- data.frame(Date = portfolio_plret_df[-c(1:1000),1],
                           alpha_1perc = ms2garch_VaR[,1], alpha_5perc = ms2garch_VaR[,2])
ms2garch_VaR_df %>% head()
rugarch::VaRTest(alpha = 0.01, portfolio_plret_df[-c(1:1000),2], ms2garch_VaR[,1])
ms2garch_ES

write.csv(ms2garch_VaR_df, "Data\\VaR\\MS2GARCH.csv", row.names = FALSE)



#--------------------------------------------------------#
########### Multivariate Normal DCC-GARCH(1,1) ###########
#--------------------------------------------------------#

portfolio.weights <- rep(1/10, 10)
dcc_meanModel <- list(armaOrder = c(0,0), include.mean = TRUE)
dcc_varModel <- list(model = "sGARCH", garchOrder = c(1,1))
dcc_uspec <- ugarchspec(variance.model = dcc_varModel, mean.model = dcc_meanModel, distribution.model = "norm")
dcc_mspec <- multispec(replicate(10, dcc_uspec))
dcc_spec <- dccspec(dcc_mspec, VAR = FALSE, model = "DCC", dccOrder = c(1,1), distribution =  "mvnorm")

cl <- makePSOCKcluster(numcores)
multi_dcc_roll <- dccroll(spec = dcc_spec, data = stocks_plret_ts, n.ahead = 1, n.start = 1000, refit.every = 1, refit.window = "moving", window.size = 1000,  cluster = cl)
stopCluster(cl)

n <- length(stocks_plret_df[,1])
n.window <- n-1000
dcc_sigma_pf <- numeric(n.window)
dcc_cov <- rcov(multi_dcc_roll) # rcov to extract forecasted Covariance matrix; rcor would extract forecasted Correlation matrix
for (i in 1:n.window) dcc_sigma_pf[i] <- portfolio.weights %*% dcc_cov[,,i] %*% portfolio.weights
multi_dcc_VaR <- data.frame(Date = portfolio_plret_df$Date[-c(1:1000)])
multi_dcc_VaR$alpha_1perc <- sapply(dcc_sigma_pf, function(x) x*qnorm(0.01, 0, 1, lower.tail = TRUE))
multi_dcc_VaR$alpha_5perc <-  sapply(dcc_sigma_pf, function(x) x*qnorm(0.05, 0, 1, lower.tail = TRUE))
head(multi_dcc_VaR)
write.csv(multi_dcc_VaR, "Data\\VaR\\Multi_DCC_GARCH.csv", row.names = FALSE)



dcc_NGARCH_meanModel <- list(armaOrder = c(0,0), include.mean = TRUE)
dcc_NGARCH_varModel <- list(model = "fGARCH", submodel = "NGARCH", garchOrder = c(1,1))
dcc_NGARCH_uspec <- ugarchspec(variance.model = dcc_NGARCH_varModel, mean.model = dcc_NGARCH_meanModel, distribution.model = "sstd")
dcc_NGARCH_mspec <- multispec(replicate(10, dcc_NGARCH_uspec))
dcc_NGARCH_spec <- dccspec(dcc_NGARCH_mspec, VAR = FALSE, model = "DCC", dccOrder = c(1,1), distribution =  "mvnorm")

cl <- makePSOCKcluster(numcores)
multi_dcc_NGARCH_roll <- dccroll(spec = dcc_NGARCH_spec, data = stocks_plret_ts, n.ahead = 1, n.start = 1000, refit.every = 1, refit.window = "moving", window.size = 1000,  cluster = cl)
stopCluster(cl)

n <- length(stocks_plret_df[,1])
n.window <- n-1000
dcc_NGARCH_sigma_pf <- numeric(n.window)
dcc_NGARCH_cov <- rcov(multi_dcc_NGARCH_roll) # rcov to extract forecasted Covariance matrix; rcor would extract forecasted Correlation matrix
for (i in 1:n.window) dcc_NGARCH_sigma_pf[i] <- portfolio.weights %*% dcc_NGARCH_cov[,,i] %*% portfolio.weights
multi_dcc_NGARCH_VaR <- data.frame(Date = portfolio_plret_df$Date[-c(1:1000)])
multi_dcc_NGARCH_VaR$alpha_1perc <- sapply(dcc_NGARCH_sigma_pf, function(x) x*qnorm(0.01, 0, 1, lower.tail = TRUE))
multi_dcc_NGARCH_VaR$alpha_5perc <-  sapply(dcc_NGARCH_sigma_pf, function(x) x*qnorm(0.05, 0, 1, lower.tail = TRUE))
head(multi_dcc_NGARCH_VaR)
write.csv(multi_dcc_NGARCH_VaR, "Data\\VaR\\Multi_dcc_NGARCH_GARCH.csv", row.names = FALSE)
