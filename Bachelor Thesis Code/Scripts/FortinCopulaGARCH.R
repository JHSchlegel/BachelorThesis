#==========================================================================#
############################### Copula GARCH ###############################
#==========================================================================#

# R version 4.2.2 (2022-10-31 ucrt)

# !!! Important Notice !!! #
# the idea to use rugarch::pdist and rugarch::qdist for PIT/ quantile 
# was inspired by the GitHub of the rmgarch package. The backtransformation
# to returns (they used sqrtm instead of Cholesky though) as well the general
# structure of fitting a DCC Copula GARCH was taken from this video:
# https://www.youtube.com/watch?v=OMjnDnGJOqY&list=LL&index=130

if (!require(tidyverse)) install.packages("tidyverse")
if (!require(rugarch)) install.packages("rugarch") # for univariate GARCH models
# rugarch allows for parallel programming using parallel package
if (!require(parallel)) install.packages("parallel") 
if (!require(rmgarch)) install.packages("rmgarch") # for DCC GARCH models
if (!require(copula)) install.packages("copula") # for copulas
if (!require(Rcpp)) install.packages("Rcpp") # to integrate C++ code in R
theme_set(theme_light())


sourceCpp("Scripts/Columnwise_Sum.cpp")
source("Scripts/Skew_t_Copula.R")

#--------------------------------------------------------#
########### Import Data and Create xts Objects ###########
#--------------------------------------------------------#

## Import data from csv files created in EDA.R
FFCFactors_df <- read.csv("./Data/FFCFactors.csv", header = TRUE)
stocks_plret_df <- read.csv("./Data/StockPlrets.csv", header = TRUE)
portfolio_plret_df <- read.csv("./Data/PortfolioPlrets.csv", header = TRUE)


## Convert to Matrix
# Matrices only include numeric values; hence especially indexing by row much 
# faster than indexing by row in dataframe
FFCFactors_mat <- as.matrix(FFCFactors_df[,-1])
stocks_plret_mat <- as.matrix(stocks_plret_df[,-1])
portfolio_plret_mat <- as.matrix(portfolio_plret_df[,-1])


## Create dataframe with all stocks and market factors for regression
joined_df <- stocks_plret_df %>% 
  full_join(FFCFactors_df, by = "Date")

head(joined_df)
dim(joined_df)
joined_df$Date

## Convert to time series
if (!require(xts)) install.packages("xts") # handling of multivariate time series
Factors_ts <- xts(FFCFactors_df[,-c(1,2)], 
                  order.by = as.Date(FFCFactors_df$Date))
stocks_plret_ts <- xts(stocks_plret_df[,-1], 
                       order.by = as.Date(stocks_plret_df$Date))
portfolio_plret_ts <- xts(portfolio_plret_df[,-1], 
                          order.by = as.Date(portfolio_plret_df$Date))


#---------------------------------------#
########### Microbenchmarking ###########
#---------------------------------------#

# Have to calculate simulated portfolio returns from simulated factor returns
# check which way it is the fastest using microbenchmark package
# for completeness sake, I included all the benchmarking I did even though
# some parts might not be directly relevant for the final function I wrote

if (!require(microbenchmark)) install.packages("microbenchmark")


mbm_idx <- microbenchmark(
  stocks_plret_ts[50:2000,2:9],
  stocks_plret_df[50:2000,2:9],
  stocks_plret_mat[50:2000,1:8],
  times = 1000
)
mbm_idx
mbm_idx %>% autoplot()
# xts and matrix way faster when trying to access subset of rows

cppFunction("double DotProdCpp(NumericVector x, NumericVector y){
  double result = 0;
  int n = x.length();
  for (int i=0; i<n; i++){
    result+=x[i]*y[i];
  }
  return result;
}")

set.seed(42)
a <- 3 # e.g. intercept or risk free rate
x <- rnorm(4)
y <- rnorm(4)
mbm_dot_small <- microbenchmark(
  base = a+x%*%y,
  crossprod = a+crossprod(x,y),
  cpp_dot_prd = a+DotProdCpp(x,y),
  times = 100
)
mbm_dot_small
mbm_dot_small %>% autoplot()
# Use %*% for dot product of simulated factor returns and OLS factor coefs


x <- rnorm(100000)
y <- rnorm(100000)
mbm_dot_large <- microbenchmark(
  base = a+x%*%y,
  crossprod = a+crossprod(x,y),
  cpp_dot_prd = a+DotProdCpp(x,y),
  times = 100
)
mbm_dot_large
mbm_dot_large %>% autoplot()
# use C++ function to calculate dot product of large vectors

numcores <- detectCores()-1

parapply_test <- function(cl, test_rets, test_coefs){
  cl <- makePSOCKcluster(numcores)
  clusterExport(cl, list("test_rets", "test_coefs"))
  sim_test <- matrix(0L, nrow = 2e5, ncol = 10)
  for (i in 1:10) sim_test[,i] <- parApply(cl, test_rets, 1, 
                function(x) test_coefs[i,1] + test_coefs[i,-1]%*%x)
  stopCluster(cl)
}


apply_test <- function(test_rets, test_coefs){
  sim_test <- matrix(0L, nrow = 2e5, ncol = 10)
  for (i in 1:10) sim_test[,i] <- apply(test_rets, 1,
                              function(x)test_coefs[i,1] + test_coefs[i,-1]%*%x)
  sim_test
}

columnwise_sum <- function(test_rets, test_coefs){
  sim_test <- matrix(0L, nrow = 2e5, ncol = 10)
  for (i in 1:10) sim_test[,i] <- test_coefs[i,1] + 
      test_coefs[i,2] * test_rets[,1] + test_coefs[i,3] * test_rets[,2]  + 
      test_coefs[i,4] * test_rets[,3] + test_coefs[i,5] * test_rets[,4]
  sim_test
}




test_rets <- matrix(rnorm(4*2e5), nrow = 2e5, ncol = 4)
test_errors <- matrix(rnorm(10*2e5), nrow = 2e5, ncol = 10)
test_coefs <- matrix(rnorm(5*10), nrow = 10, ncol = 5)



mbm_apply <- microbenchmark(
  parapply = parapply_test(cl, test_rets, test_coefs),
  apply = apply_test(test_rets, test_coefs),
  columnwise = columnwise_sum(test_rets, test_coefs),
  times = 20
)
mbm_apply
mbm_apply %>% autoplot()
# apply faster than parapply in this case, probably due to communication time
# between threads when parallelizing
# columnwise summing almost a hundred times faster


in_loop <- function(constant, test_rets, test_coefs, test_errors){
  sim_test <- matrix(0L, nrow = 2e5, ncol = 10)
  for (i in 1:10) sim_test[,i] <- constant + test_coefs[i,1] +
      test_coefs[i,2] * test_rets[,1] +
      test_coefs[i,3] * test_rets[,2] + test_coefs[i,4]*test_rets[,3] + 
      test_coefs[i,5] * test_rets[,4] + test_errors[,i]
  sim_test
}
out_loop <- function(constant, test_rets, test_coefs, test_errors){
  sim_test <- matrix(0L, nrow = 2e5, ncol = 10)
  for (i in 1:10) sim_test[,i] <- constant + test_coefs[i,1] + 
      test_coefs[i,2] * test_rets[,1] + test_coefs[i,3] * test_rets[,2] + 
      test_coefs[i,4] * test_rets[,3] + test_coefs[i,5] * test_rets[,4]
  sim_test + test_errors
}




x <- columnwise_sum_cpp(rep(a, 2e5), test_rets, test_coefs, test_errors, 2e5)
y <- in_loop(a, test_rets, test_coefs, test_errors)
z <- out_loop(a, test_rets, test_coefs, test_errors)
all.equal(x,y)
all.equal(x, z)
# all functions yield the same result

mbm_loop <- microbenchmark(
  inloop = in_loop(a, test_rets, test_coefs, test_errors),
  outloop = out_loop(a, test_rets, test_coefs, test_errors),
  cpp = columnwise_sum_cpp(rep(a, 2e5), test_rets, test_coefs, test_errors, 2e5),
  times = 50
)
mbm_loop
mbm_loop %>% autoplot()
# outside of loop is faster than inside the loop
# c++ 7x faster than columnwise implementation w/ summation of error matrix 
# outside the loop

sim_test <- columnwise_sum_cpp(rep(a, 2e5), test_rets, test_coefs, 
                               test_errors, 2e5)

cl <- makePSOCKcluster(numcores)
clusterExport(cl, list("sim_test"))
mbm_pf <- microbenchmark(
  apply = apply(sim_test, 1, mean),
  parallel = parApply(cl, sim_test, 1, mean),
  rowmeans = rowMeans(sim_test),
  times = 50
)
stopCluster(cl)
mbm_pf
mbm_pf %>% autoplot()
# rowMeans the fastest; when testing it was even faster than a quick C++ 
# implementation I made

# Hence, for calculating the returns of the stocks, the C++ function 
# columnwise_sum_cpp will be used
# To get the portfolio returns from the stock returns, rowMeans will be used



#-----------------------------------------------------------------------------#
########### Coefficients and Residuals of Carhart Four-Factor Model ###########
#-----------------------------------------------------------------------------#

if (!require(tseries)) install.packages(tseries) # for JB test
if (!require(sandwich)) install.packages(sandwich)
n_dates <- dim(FFCFactors_df)[1]
coefs_mat <- matrix(0L, nrow = 10, ncol = 5) # empty matrix for coefficients
error_mat <- matrix(0L, nrow = n_dates, ncol = 10) # empty matrix for residuals
adj_R2 <- numeric(10)
JB <- numeric(10)
for (i in 1:10){
  # columns 2-11 are the shares
  fit <- lm((joined_df[,i+1]-RF) ~ Mkt.RF + SMB+ HML + Mom,data = joined_df)
  par(mfrow=c(2,2))
  plot(fit)
  coefs_mat[i,] <- coef(fit)
  error_mat[,i] <- resid(fit)
  adj_R2[i] <- summary(fit)$adj.r.squared
  JB[i] <- tseries::jarque.bera.test(resid(fit))$statistic
}
# residuals vs fitted looks fine
# qq-plots show leptokurtic behavior
# scale location plots show heteroscedastic behaviour
# some points with high leverage but still w/in 0.5 boundary for Cook's distance

coefs_df <- data.frame(coefs_mat); error_df <- data.frame(Date = joined_df$Date,
                                                          error_mat)
rownames(coefs_df) <- c("KO", "XOM", "GE", "IBM", "CVX", "UTC", "PG", "CAT", 
                        "BA", "MRK")
colnames(error_df) <- c("Date", "KO", "XOM", "GE", "IBM", "CVX", "UTC", "PG", 
                        "CAT", "BA", "MRK")
colnames(coefs_df) <- c("Intercept", "Market", "Size", "Value", "Momentum")
head(coefs_df); head(error_df)

OLS_table <- data.frame(coefs_df, adj_R2 = adj_R2, JB = JB) %>% round(3)
OLS_table

source("Scripts/EDA.R") # for summary_statistics function
summary_statistics(error_df) #highly non-normal OLS residuals

## Plot and investigate error distribution
error_df[,-1] %>% 
  gather() %>% 
  ggplot(aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 80, fill = "grey69", 
                 color = "white")+
  geom_density()+
  facet_wrap(~key, scale = "free_x")
# very long left tail; unimodal w/ mode around zero

apply(error_mat, 2, min)
apply(error_mat, 2, max)
# minimal values are a lot smaller in absolute value than maximal values for 
# all shares
# can observe that the bootstrapped error distributions have a long left tail

if (!require(pheatmap)) install.packages("pheatmap")
error_df[,-1] %>% 
  cor() %>% 
  pheatmap(display_numbers = TRUE)
# most residuals only barely correlated
# XOM& CVX residuals show strongest correlation, followed by RTX& BA and KO&PG

if (!require(corrr)) install.packages("corrr")
error_df[,-1] %>% 
  correlate() %>% 
  rearrange() %>% 
  network_plot(min_cor = 0, colors = c("blue", "white", "red"))
# 3 or 4 clusters
# CVX and XOM together; both oil companies
# Boeing and UTC together; both aircraft/technology
# KO and PG; both consumption goods
# others are less clear
# can clearly see that error terms show dependencies based on the sector

## Bootstrap OLS Residuals Distribution
set.seed(42)
N_boot <- 200000
bootind <- sample.int(n_dates, size = N_boot, replace = TRUE)
head(error_df[bootind,])
error_vec_resampled <- error_df[bootind,-1] # now without date
head(error_vec_resampled)

## Plot and investigate bootstrapped error distribution
error_vec_resampled %>%
  gather() %>% 
  ggplot(aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 80, fill = "grey69", 
                 color = "white")+
  geom_density()+
  facet_wrap(~key, scale = "free_x")
# distribution very similar to original just slightly less smooth 

error_vec_resampled %>% 
  cor() %>% 
  pheatmap(display_numbers = TRUE)
# correlations almost the same as in original OLS residuals


error_df %>% 
  correlate() %>% 
  rearrange() %>% 
  network_plot(min_cor = 0, colors = c("blue", "white", "red"))
# can identify similar clusters as in original OLS residuals

# Hence bootstrapping accurately preserves crosscorrelation between residuals

#-------------------------------------------------------------------------------#
########### Bootstrap Resample Error Distributions of OLS Regressions ###########
#-------------------------------------------------------------------------------#

# Each row in error_vec_resampled is one draw of a vector of error terms
# to make sure dependencies are kept, all elements in a row are from the same
# time t

## Hence we will be using the statements below in the final function to preserve
## relations between different stocks
set.seed(42)
N_boot <- 200000
bootind <- sample.int(n_dates, size = N_boot, replace = TRUE)
error_vec_resampled <- error_mat[bootind,] # now without date
head(error_vec_resampled)


#-------------------------------------------------------------#
########### Copula GARCH as in Fortin et al. (2022) ###########
#-------------------------------------------------------------#

#' Copula GARCH with Factor Returns
#' 
#' The function uses some global variables specified above. Hence, quickly run 
#' all sections but microbenchmarking quickly before using this function
#'
#' @param DCC_corr_mat boolean: whether DCC-GARCH corr matrix is used for Copula 
#' fitting or not
#' @param pseudo_obs boolean: whether pseudo observations or PIT should be used
#' @param ugarch_dist choose between "norm" and "sstd" for distribution of 
#' GARCH marginals
#' @param ugarch_model choose between "GARCH" and "NGARCH" for conditional 
#' variance model of marginals
#' @param copula_dist copula distribution: choose from "norm", "t" and "skewt" 
#' distribution
#' @param error_mat Tx10 Matrix with all OLS errors of Carhart four-factor model
#' By default use error mat calculated above in section about Carhart four-
#' factor model
#' @return dataframe with dates in first column, 0.01 VaR in second column and
#' 0.05 VaR in third column
fortin_cgarch_VaR <- function(DCC_corr_mat, pseudo_obs = TRUE, 
                              ugarch_dist = c("norm", "sstd"),
                              ugarch_model = c("GARCH", "NGARCH"),
                              copula_dist = c("norm", "t", "skewt"),
                              error_mat = error_mat){
  N_sim <- 200000 # nr. of random draws from copula
  N_boot <- 200000 # nr. of bootstrap resamples of OLS residuals
  n_ahead <- 1 # 1 step ahead forecast
  n_window <- length(FFCFactors_df[,1])-1000 # how many rolling windows
  numcores <- detectCores() - 1 # use all but 1 core
  
  ugarch_dist <- match.arg(ugarch_dist)
  
  ugarch_model = match.arg(ugarch_model)
  
  copula_dist <- match.arg(copula_dist)
  
  ## DCC specification
  meanModel <- list(armaOrder = c(0, 0), include.mean = TRUE)
  varModel <- switch(ugarch_model,
                     "NGARCH" = list(model = "fGARCH", submodel = "NGARCH",
                                     garchOrder = c(1,1)),
                     "GARCH" = list(model = "sGARCH", garchOrder = c(1,1))
  )
  uspec <- ugarchspec(varModel, mean.model = meanModel, 
                      distribution.model = ugarch_dist)
  mspec <- multispec(replicate(4, uspec))
  dcc_spec <- dccspec(mspec, VAR = FALSE, model = "DCC", dccOrder = c(1,1), 
                      distribution =  "mvnorm")
  
  
  VaR_cop <- matrix(0L, nrow = n_window, ncol = 2) # empty matrix to store VaR
  
  
  ## rolling window VaR forecast
  for (i in 1:n_window){
    cl <- makePSOCKcluster(numcores)
    dcc_fit <- dccfit(dcc_spec, data = Factors_ts[i:(1000+i-1),],cluster = cl)
    dcc_fcst <- dccforecast(dcc_fit, cluster = cl)
    stopCluster(cl)
    
    ## Extract standardized residuals
    garch_dcc_res <- dcc_fit@mfit$stdresid
    
    
    ## Extract distribution of marginals (will be used later on 
    ## for quantile transformation)
    marginal_dist <- dcc_fit@model$umodel$modeldesc$distribution
    if (pseudo_obs == TRUE){
      for (k in 1:4){
        # get skewness parameter and store as skew_k
        skew <- paste("skew", k, sep = "_")
        assign(skew, ifelse(dcc_fit@model$umodel$modelinc["skew", k]>0,
                            dcc_fit@model$mpars["skew", k],
                            0))
        
        # get shape parameter and store as shape_k
        shape <- paste("shape", k, sep = "_")
        assign(shape, ifelse(dcc_fit@model$umodel$modelinc["shape", k]>0,
                             dcc_fit@model$mpars["shape", k],
                             0))
      }
    }
    
    else {
      u_res <- matrix(0L, nrow = dim(garch_dcc_res)[1], ncol = 4)
      for (k in 1:4){
        # get skewness parameter and store as skew_k
        skew <-  ifelse(dcc_fit@model$umodel$modelinc["skew", k]>0, 
                        dcc_fit@model$mpars["skew", k], 
                        0)
        skew_k <- paste("skew", k, sep = "_")
        assign(skew_k, skew) # name skew_1, ..., skew_4 to reuse for qdist later
        
        # get shape parameter and store as shape_k
        shape <- ifelse(dcc_fit@model$umodel$modelinc["shape", k]>0, 
                        dcc_fit@model$mpars["shape", k], 
                        0)
        shape_k <- paste("shape", k, sep = "_")
        assign(shape_k, shape)  # name shape_1, ..., shape_4 to reuse for qdist 
        
        ## PIT to get u's that lie within unit cube
        u_res[,k] <- rugarch::pdist(distribution = marginal_dist[k], 
                                    q  = garch_dcc_res[,k], 
                                    mu = 0,
                                    sigma = 1,
                                    shape = shape,
                                    skew = skew)
      }
    }
    
    
    
    ## Extract forecasted cov matrix and mean
    dcc_fcst_cov <- matrix(dcc_fcst@mforecast$H[[1]], nrow = 4) 
    dcc_fcst_cor <- cov2cor(dcc_fcst_cov)
    dcc_fcst_mean <- matrix(dcc_fcst@mforecast$mu, nrow = 1)
    
    ## Create pseudoobservations from standardized residuals that lie within
    ## unit cube
    pobs_res <- apply(garch_dcc_res, 2, function(x) copula::pobs(x))
    
    ## Fit copula
    
    # use pobs or PIT values for fitting depending on user input
    ifelse(pseudo_obs == TRUE, data <- pobs_res, data <- u_res)
    # start values for skew t
    start <- list(rho = numeric(6), delta=numeric(4), nu = 6)
    
    if (DCC_corr_mat == FALSE){
      if (copula_dist == "skewt"){
        if (i-1%%20==0){
          cop_fit <- tryCatch(
            {stcop.mle(data, start=start, 
                       control = list(reltol=1e-4, maxit = 10000))},
            error = function(cond){
              stcop_MLE <- stcop.mle(data, start=start,
                                     control = list(reltol=1e-4, maxit = 10000),
                                     solver = "Nelder-Mead")
              message("Had to use Nelder-Mead")
              return(stcop_MLE)
            }
          )
        }
      }
      else {
        cop_fit <- switch(copula_dist,
                          "norm" = fitCopula(
                            normalCopula(dim = 4), data = data
                          ),
                          "t" = fitCopula(
                            tCopula(dim = 4, df.fixed = FALSE), data = data
                          )
                          
        )}
    }
    
    
    if (DCC_corr_mat == TRUE){
      cop_fit <- switch(copula_dist,
          "norm" = fitCopula(
            normalCopula(dim = 4, param = dcc_fcst_cor[lower.tri(dcc_fcst_cor)], 
                       dispstr = "un"), data = data
      ),
          "t" = fitCopula(
            tCopula(dim = 4, param = dcc_fcst_cor[lower.tri(dcc_fcst_cor)],
                    df.fixed = FALSE, dispstr = "un"), data = data)
                          
    )
    }

    
    ## Simulation from copula:
    if (copula_dist == "skewt"){
      rho <- cop_fit$dp$rho
      delta <- cop_fit$dp$delta
      df <- cop_fit$dp$nu
      
      set.seed(i)
      cop_sim <- rstcop(N_sim, rho, delta, df)$u
    }
    else{
      set.seed(i)
      cop_sim <- rCopula(N_sim, cop_fit@copula)
    }
    
    
    ## Quantile transform
    sim_residuals <- cbind(
      rugarch::qdist(distribution = marginal_dist[1], cop_sim[,1], 
                     mu = 0, sigma = 1, skew = skew_1, shape = shape_1),
      rugarch::qdist(distribution = marginal_dist[2], cop_sim[,2], 
                     mu = 0, sigma = 1, skew = skew_2, shape = shape_2), 
      rugarch::qdist(distribution = marginal_dist[3], cop_sim[,3], 
                     mu = 0, sigma = 1, skew = skew_3, shape = shape_3), 
      rugarch::qdist(distribution = marginal_dist[4], cop_sim[,4], 
                     mu = 0, sigma = 1, skew = skew_4, shape = shape_4)
    )
    
    
    # returns are X_t = mu_t+sigma_t*epsilon_t
    # where sigma_t is calculated using Cholesky decomposition of forecasted
    # sigma matrix
    factor_rets <- matrix(0L, nrow = N_sim, ncol=4)
    chol_cov <- chol(dcc_fcst_cov)
    factor_rets <- matrix(rep(dcc_fcst_mean, each = N_sim), ncol = 4) + 
      t(chol_cov%*%t(sim_residuals))
    
    set.seed(i)
    bootind <- sample(i:(1000+i-1), size = N_boot, replace = TRUE)
    error_vec_resampled <- error_mat[bootind,] 
    
    ## sim_returns = alpha+beta%*%factor returns + bootstrapped OLS residuals
    # for more information about function see cpp file with function
    sim_returns <- columnwise_sum_cpp(
      rep(FFCFactors_mat[1000+i-1,1], 2e5), 
      factor_rets, 
      coefs_mat, 
      error_vec_resampled, 
      2e5
      )
    
    # calculate portfolio percentage log returns for equally weighted portfolio
    
    fractional_arithmetic_returns <- apply(sim_returns, 2,
                                           function(x)exp((x/100))-1)
    sim_pf_returns <- rowMeans(fractional_arithmetic_returns)
    sim_pf_plreturns <- sapply(sim_pf_returns, function(x) 100*log(x+1))
    
    VaR_cop[i,] <- quantile(sim_pf_plreturns, c(0.01, 0.05))
    message("Completed: ", i, " of ", n_window)
    message(VaR_cop[i,])
  }
  VaR_cop_df <- data.frame(Date = portfolio_plret_df$Date[-c(1:1000)],
                           alpha_0.01 = VaR_cop[,1], alpha_0.05 = VaR_cop[,2])
  invisible(VaR_cop_df)
}


## NGARCH normal copula
VaR_cop_norm_NGARCH_df <- fortin_cgarch_VaR(
  DCC_corr_mat = FALSE, pseudo_obs = TRUE, ugarch_dist = "sstd",
  ugarch_model = "NGARCH", copula_dist = "norm", error_mat = error_mat
)

write.csv(VaR_cop_norm_NGARCH_df, 
          "Data\\VaR\\Fortin_cop_norm_NGARCH.csv", row.names = FALSE)

## NGARCH t copula
VaR_cop_t_NGARCH_df <- fortin_cgarch_VaR(
  DCC_corr_mat = FALSE, pseudo_obs = TRUE, ugarch_dist = "sstd",
  ugarch_model = "NGARCH", copula_dist = "t", error_mat = error_mat
)

write.csv(VaR_cop_t_NGARCH_df, 
          "Data\\VaR\\Fortin_cop_t.csv", row.names = FALSE)





## sGARCH normal Copula
VaR_cop_norm_sGARCH_df <- fortin_cgarch_VaR(
  DCC_corr_mat = FALSE, pseudo_obs = TRUE, ugarch_dist = "norm",
  ugarch_model = "GARCH", copula_dist = "norm", error_mat = error_mat
  )

write.csv(VaR_cop_norm_sGARCH_df, 
          "Data\\VaR\\Fortin_cop_norm_sGARCH.csv", row.names = FALSE)


## sGARCH t copula
VaR_cop_t_sGARCH_df <- fortin_cgarch_VaR(
  DCC_corr_mat = FALSE, pseudo_obs = TRUE, ugarch_dist = "norm",
  ugarch_model = "GARCH", copula_dist = "t", error_mat = error_mat
  )
write.csv(VaR_cop_t_sGARCH_df, 
          "Data\\VaR\\Fortin_cop_t_sGARCH.csv", row.names = FALSE)

## sGARCH skewt copula
VaR_cop_skewt_sGARCH_df <- fortin_cgarch_VaR(
  DCC_corr_mat = FALSE, pseudo_obs = F, ugarch_dist = "norm",
  ugarch_model = "GARCH", copula_dist = "skewt", error_mat = error_mat
  )
write.csv(VaR_cop_skewt_sGARCH_df, 
          "Data\\VaR\\Fortin_cop_skewt_sGARCH.csv", row.names = FALSE)


