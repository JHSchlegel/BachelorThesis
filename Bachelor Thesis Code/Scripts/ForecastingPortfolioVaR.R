#=============================================================================#
########################## Forecasting Portfolio VaR ##########################
#=============================================================================#

if (!require(parallel)) install.packages("parallel") # parallel computing
if (!require(rugarch)) install.packages("rugarch") # univariate GARCH models
if (!require(rmgarch)) install.packages("rmgarch") # dcc garch
if (!require(tidyverse)) install.packages("tidyverse") # data manipulation
if (!require(xts)) install.packages("xts") # for multivariate time series
if (!require(MSGARCH)) install.packages("MSGARCH") # for MN GARCH


#--------------------------------------------------------#
########### Import Data and Create xts Objects ###########
#--------------------------------------------------------#

FFCFactors_df <- read.csv("./Data/FFCFactors.csv", header = TRUE)
stocks_plret_df <- read.csv("./Data/StockPlrets.csv", header = TRUE)
portfolio_plret_df <- read.csv("./Data/PortfolioPlrets.csv", header = TRUE)



## Convert to time series
FFCFactors_ts <- xts(FFCFactors_df[,-1], order.by = as.Date(FFCFactors_df$Date))
stocks_plret_ts <- xts(stocks_plret_df[,-1],
                       order.by = as.Date(stocks_plret_df$Date))
portfolio_plret_ts <- xts(portfolio_plret_df[,-1],
                          order.by = as.Date(portfolio_plret_df$Date))

## Detect how many cores are available for parallelizing
# subtract one so that not all cores are used
numcores <- detectCores()-1


#--------------------------------------------------#
########### Univariate Normal GARCH(1,1) ###########
#--------------------------------------------------#

# no ARFIMA, just use unconditional mean
uni_normal_meanModel <- list(armaOrder=c(0,0), include.mean=TRUE) 
uni_normal_varModel <- list(model = "sGARCH", garchOrder = c(1,1))
uni_normal_spec <- ugarchspec(
  variance.model = uni_normal_varModel, mean.model = uni_normal_meanModel, 
  distribution.model = "norm"
  )
cl <- makePSOCKcluster(numcores)

## Rolling window fitting
uni_normal_GARCH_roll <- ugarchroll(
  uni_normal_spec, portfolio_plret_ts, refit.every = 1, refit.window = "moving", 
  n.start = 1000, solver = "hybrid", calculate_var = TRUE, 
  VaR.alpha = c(0.01, 0.05), cluster = cl, keep.coef = FALSE
  )

stopCluster(cl)

uni_normal_GARCH_VaR <- data.frame(
  Date = portfolio_plret_df$Date[-c(1:1000)], 
  alpha_1perc = uni_normal_GARCH_roll@forecast$VaR[,1],
  alpha_5perc = uni_normal_GARCH_roll@forecast$VaR[,2]
  )

head(uni_normal_GARCH_VaR)
write.csv(uni_normal_GARCH_VaR, "Data\\VaR\\Uni_Normal_GARCH_VaR.csv", 
          row.names = FALSE)




#------------------------------------------------------#
########### Univariate EWMA w/ lambda = 0.94 ###########
#------------------------------------------------------#


lambda <- 0.94 # as in RiskMetrics technical document
# no ARFIMA, just use unconditional mean
ewma_meanModel <- list(armaOrder=c(0,0), include.mean=TRUE)
ewma_varModel <- list(model="iGARCH", garchOrder=c(1,1))
uni_ewma_spec <- ugarchspec(
  variance.model= ewma_varModel, mean.model= ewma_meanModel,  
  distribution.model="norm", fixed.pars=list(omega=0, alpha = 1-lambda)
  )
cl <- makePSOCKcluster(numcores)
uni_ewma_roll <- ugarchroll(
  uni_ewma_spec, portfolio_plret_ts, n.start = 1000, refit.every = 1,
  refit.window = "moving", solver = "hybrid", calculate_var = TRUE, 
  VaR.alpha = c(0.01, 0.05), cluster = cl, keep.coef = FALSE
  )
stopCluster(cl)
uni_ewma_VaR <- data.frame(
  Date = portfolio_plret_df$Date[-c(1:1000)], 
  alpha_1perc = uni_ewma_roll@forecast$VaR[,1],
  alpha_5perc = uni_ewma_roll@forecast$VaR[,2]
  )


head(uni_ewma_VaR)
write.csv(uni_ewma_VaR, "Data\\VaR\\Uni_EWMA_VaR.csv", row.names = FALSE)


#-------------------------------------------------#
########### Univariate t-GJR-GARCH(1,1) ###########
#-------------------------------------------------#

# no ARFIMA, just use unconditional mean
uni_t_gjr_meanModel <- list(armaOrder = c(0,0), include.mean = TRUE)
uni_t_gjr_varModel <- list(model = "gjrGARCH", garchOrder = c(1,1))
uni_t_gjr_spec <- ugarchspec(
  variance.model = uni_t_gjr_varModel, mean.model = uni_t_gjr_meanModel, 
  distribution.model = "std"
  )
cl <- makePSOCKcluster(numcores)
uni_t_gjr_GARCH_roll <- ugarchroll(
  uni_t_gjr_spec, portfolio_plret_ts, n.start = 1000, refit.every = 1,
  refit.window = "moving", solver = "hybrid",
  calculate_var = TRUE, VaR.alpha = c(0.01, 0.05),
  cluster = cl, keep.coef = FALSE
  )
stopCluster(cl)
uni_t_gjr_GARCH_VaR <- data.frame(
  Date = portfolio_plret_df$Date[-c(1:1000)], 
  alpha_1perc = uni_t_gjr_GARCH_roll@forecast$VaR[,1],
  alpha_5perc = uni_t_gjr_GARCH_roll@forecast$VaR[,2]
  )


head(uni_t_gjr_GARCH_VaR)
write.csv(uni_t_gjr_GARCH_VaR, "Data\\VaR\\Uni_t_GJR_GARCH.csv",
          row.names = FALSE)


#-----------------------------------------------------#
########### Univariate skewt-GJR-GARCH(1,1) ###########
#-----------------------------------------------------#

# no ARFIMA, just use unconditional mean
uni_skewt_gjr_meanModel <- list(armaOrder = c(0,0), include.mean = TRUE)
uni_skewt_gjr_varModel <- list(model = "gjrGARCH", garchOrder = c(1,1))
uni_skewt_gjr_spec <- ugarchspec(
  variance.model = uni_skewt_gjr_varModel, mean.model = uni_skewt_gjr_meanModel, 
  distribution.model = "sstd"
  )
cl <- makePSOCKcluster(numcores)
uni_skewt_gjr_GARCH_roll <- ugarchroll(
  uni_skewt_gjr_spec, portfolio_plret_ts, n.start = 1000, refit.every = 1,
  refit.window = "moving", solver = "hybrid", calculate_var = TRUE, 
  VaR.alpha = c(0.01, 0.05), cluster = cl, keep.coef = FALSE
  )
stopCluster(cl)
uni_skewt_gjr_GARCH_VaR <- data.frame(
  Date = portfolio_plret_df$Date[-c(1:1000)],
  alpha_1perc = uni_skewt_gjr_GARCH_roll@forecast$VaR[,1],
  alpha_5perc = uni_skewt_gjr_GARCH_roll@forecast$VaR[,2]
  )
head(uni_skewt_gjr_GARCH_VaR)
write.csv(uni_skewt_gjr_GARCH_VaR, "Data\\VaR\\Uni_skewt_GJR_GARCH.csv", 
          row.names = FALSE)


#--------------------------------------------------#
########### Univariate Skewt-NGARCH(1,1) ###########
#--------------------------------------------------#

# no ARFIMA, just use unconditional mean
uni_skewt_NGARCH_meanModel <- list(armaOrder = c(0, 0), include.mean = TRUE)
uni_skewt_NGARCH_varModel <- list(model = "fGARCH", submodel = "NGARCH", 
                                  garchOrder = c(1,1))
uni_skewt_NGARCH_spec <- ugarchspec(
  variance.model = uni_skewt_NGARCH_varModel, 
  mean.model = uni_skewt_NGARCH_meanModel, distribution.model = "sstd"
  )
cl <- makePSOCKcluster(numcores)
uni_skewt_NGARCH_roll <- ugarchroll(
  uni_skewt_NGARCH_spec, portfolio_plret_ts, n.start = 1000, refit.every = 1,
  refit.window = "moving", solver = "hybrid",calculate_var = TRUE, 
  VaR.alpha = c(0.01, 0.05), cluster = cl, keep.coef = FALSE)
stopCluster(cl)

uni_skewt_NGARCH_VaR <- data.frame(
  Date = portfolio_plret_df$Date[-c(1:1000)],
  alpha_1perc = uni_skewt_NGARCH_roll@forecast$VaR[,1],
  alpha_5perc = uni_skewt_NGARCH_roll@forecast$VaR[,2]
  )
head(uni_skewt_NGARCH_VaR)
write.csv(uni_skewt_NGARCH_VaR, "Data\\VaR\\Uni_Skewt_NGARCH.csv", 
          row.names = FALSE)


#--------------------------------------#
########### Mix-Normal GARCH ###########
#--------------------------------------#

# MN(k, g) as in Mixed Normal Conditional Heteroskedasticity by Haas, Mittnik
# and Paolella (2004). doi: https://doi.org/10.1093/JJFINEC/NBH009

## create cutom optimizer function for msgarch::fitML
# default method in MSGARCH is "BFGS" which sometimes caused errors
f_custom_optim <- function(vPw, f_nll, spec, data, do.plm){
  out <- stats::optim(vPw, f_nll, spec = spec, data = data,
                      do.plm = do.plm, method = "L-BFGS-B")
  return(out)
}

#' Diagonal Mix-Normal GARCH
#'
#' Implements MN(k,g) GARCH models as specified in Haas, Mittnik and 
#' Paolella (2004)
#'
#' @param k how many component densities
#' @param g how many of component densities follow sGARCH process. The remaining
#' k-g components are restricted to be constant
#'
#' @return dataframe of VaR with date in first column, 1% VaR in second column
#' and 5% VaR in third column
mn_k_g_garch <- function(k, g){
  window_size <- 1000
  n_windows <- dim(portfolio_plret_df)[1]-window_size
  mn_k_g_garch_VaR <- matrix(0L, nrow = n_windows, ncol = 2)
  # mn_k_g_garch_ES <- matrix(0L, nrow = n_windows, ncol = 2)
  
  
  mn_k_g_garch_spec <- CreateSpec(
    # g GARCH terms, k-g constant
    variance.spec = list(model = rep("sGARCH", g)), 
    # mixture of k normal distributions
    distribution.spec = list(distribution = rep("norm", k)),
    switch.spec = list(do.mix = TRUE) # do.mix for mixed normal GARCH (diagonal)
    )
  
  message("===================================================================")
  message("                                MN(",k,",",g,")                    ")
  message("===================================================================")
  for (i in 316:n_windows){
    # rolling window fitting
    fit_ml <- tryCatch(
      {MSGARCH::FitML(
      spec = mn_k_g_garch_spec, 
      data = portfolio_plret_ts[i:(1000+i-1)],
      ctr = list(OptimFUN = f_custom_optim)
      )},
      error = function(cond) {
        MSGARCH::FitML(
          spec = mn_k_g_garch_spec, 
          data = portfolio_plret_ts[i:(1000+i-1)]
        ) # use BFGS instead of "L-BFGS-B"
      })
    mn_k_g_garch_VaR[i,] <- MSGARCH::Risk(fit_ml, alpha = c(0.01, 0.05), 
                                          nahead = 1)$VaR
    # mn_k_g_garch_ES[i,] <- Risk(fit_ml, alpha = c(0.01, 0.05), nahead = 1)$ES
    message("completed: ", i, " of ", n_windows)
  }
  mn_k_g_garch_VaR_df <- data.frame(Date = portfolio_plret_df[-c(1:1000),1],
                                    alpha_1perc = mn_k_g_garch_VaR[,1], 
                                    alpha_5perc = mn_k_g_garch_VaR[,2])
  invisible(mn_k_g_garch_VaR_df)
}

## MN(2,2)
mn_2_2_garch <- mn_k_g_garch(2,2)
write.csv(mn_2_2_garch, "Data\\VaR\\Uni_MN_2_2_GARCH.csv", row.names = FALSE)

## MN(3,2)
mn_3_2_garch <- mn_k_g_garch(3,2)
write.csv(mn_3_2_garch, "Data\\VaR\\Uni_MN_3_2_GARCH.csv", row.names = FALSE)

## MN(3,3)
mn_3_3_garch <- mn_k_g_garch(3,3) #316
write.csv(mn_3_3_garch, "Data\\VaR\\Uni_MN_3_2_GARCH.csv", row.names = FALSE)

rugarch::VaRTest(alpha = 0.01, portfolio_plret_df[-c(1:1000),2], 
                 mn_2_2_garch[,2])




#--------------------------------------------------------#
########### Multivariate Normal DCC-GARCH(1,1) ###########
#--------------------------------------------------------#

portfolio.weights <- rep(1/10, 10) # weight vector for 1/N portfolio
# no ARFIMA, just use unconditional mean
dcc_meanModel <- list(armaOrder = c(0,0), include.mean = TRUE) 
dcc_varModel <- list(model = "sGARCH", garchOrder = c(1,1)) # variance model
dcc_uspec <- ugarchspec(
  variance.model = dcc_varModel, 
  mean.model = dcc_meanModel, 
  distribution.model = "norm"
  )
# model all stocks with same mean model and variance model
dcc_mspec <- multispec(replicate(10, dcc_uspec))
# normal DCC GARCH(1,1)
dcc_spec <- dccspec(dcc_mspec, VAR = FALSE, model = "DCC", dccOrder = c(1,1), 
                    distribution =  "mvnorm")

cl <- makePSOCKcluster(numcores)
multi_dcc_roll <- dccroll(
  spec = dcc_spec, data = stocks_plret_ts[c(1:1005),], n.ahead = 1, n.start = 1000, 
  refit.every = 1, refit.window = "moving", window.size = 1000,  cluster = cl
  )
stopCluster(cl)

n <- length(stocks_plret_df[,1])
n_windows <- n-1000
dcc_sigma_pf <- numeric(n_windows)
dcc_mu_stocks <- matrix(NA, nrow = n_windows, ncol = 10)
# rcov to extract forecasted Covariance matrix; rcor would extract forecasted
# Correlation matrix
dcc_cov <- rcov(multi_dcc_roll) 

## calculating portfolio variance and extracting means of stocks
for (i in 1:5){#n_windows) {
  dcc_sigma_pf[i] <- portfolio.weights %*% dcc_cov[,,i] %*% portfolio.weights
  dcc_mu_stocks[i,] <- multi_dcc_roll@mforecast[[i]]@mforecast$mu
}
## Calculating portfolio mean
dcc_mu_pf <- apply(dcc_mu_stocks, 1, function(x) 100*log(mean(exp(x/100)-1)+1))

## Calculating VaR and saving it in data frame  
multi_dcc_VaR <- data.frame(Date = portfolio_plret_df$Date[-c(1:1000)])
multi_dcc_VaR$alpha_1perc <- dcc_mu_pf + sapply(dcc_sigma_pf, function(x) 
  x*qnorm(0.01, 0, 1, lower.tail = TRUE))
multi_dcc_VaR$alpha_5perc <-  dcc_mu_pf +sapply(dcc_sigma_pf, function(x) 
  x*qnorm(0.05, 0, 1, lower.tail = TRUE))


head(multi_dcc_VaR)
write.csv(multi_dcc_VaR, "Data\\VaR\\Multi_DCC_GARCH.csv", row.names = FALSE)

