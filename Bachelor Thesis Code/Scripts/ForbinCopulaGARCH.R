#TODO: 1) OLS regression

#TODO: 2) bootstrap OLS residuals

#TODO: 3) estimate distribution F of factor returns via Copula GARCH

#TODO: 4) simulate rt from F& errort from Fe

#TODO: 5) calculate stock returns w/FFC

#TODO:6) calculate portfolio return

#TODO: VaR_t^p=-percentile(rpf,100p)

# TODO: DCC which distribution? how to backtransform? 

# TODO: set.seed(i) for resampling

if (!require(tidyverse)) install.packages("tidyverse")
if (!require(rugarch)) install.packages("rugarch")
if (!require(rmgarch)) install.packages("rmgarch")
if (!require(parallel)) install.packages("parallel")
if (!require(copula)) install.packages("copula")
if (!require(Rcpp)) install.packages("Rcpp")
theme_set(theme_light())


sourceCpp("Scripts/CppFunctions.cpp")
#--------------------------------------------------------#
########### Import Data and Create xts Objects ###########
#--------------------------------------------------------#

## Import data from csv files created in EDA.R
FFCFactors_df <- read.csv("./Data/FFCFactors.csv", header = TRUE)
stocks_plret_df <- read.csv("./Data/StockPlrets.csv", header = TRUE)
portfolio_plret_df <- data.frame(Date = stocks_plret_df$Date, Portfolio = apply(stocks_plret_df[,-1], 1, mean))

## Convert to Matrix
# Matrices only include numeric values; hence especially indexing by row much faster than indexing by row in dataframe
FFCFactors_mat <- as.matrix(FFCFactors_df[,-1])
stocks_plret_mat <- as.matrix(stocks_plret_df[,-1])
portfolio_plret_mat <- as.matrix(portfolio_plret_df[,-1])


## Create dataframe with all stocks and market factors
joined_df <- stocks_plret_df %>% 
  full_join(FFCFactors_df)

head(joined_df)

## Convert to time series
if (!require(xts)) install.packages("xts") # allows handling of multivariate time series
Factors_ts <- xts(FFCFactors_df[,-c(1,2)], order.by = as.Date(FFCFactors_df$Date))
stocks_plret_ts <- xts(stocks_plret_df[,-1], order.by = as.Date(stocks_plret_df$Date))
portfolio_plret_ts <- xts(portfolio_plret_df[,-1], order.by = as.Date(portfolio_plret_df$Date))


#---------------------------------------#
########### Microbenchmarking ###########
#---------------------------------------#

# Have to calculate simulated portfolio returns from simulated factor returns
# check which way it is the fastest using microbenchmarking

if (!require(microbenchmark)) install.packages("microbenchmark")


mbm_idx <- microbenchmark(
  stocks_plret_ts[50:2000,2:9],
  stocks_plret_df[50:2000,2:9],
  stocks_plret_mat[50:2000,1:8],
  times = 1000
)
mbm_idx
plot(mbm_idx)
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
plot(mbm_dot_small)
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
plot(mbm_dot_large)
# use C++ function to calculate dot product of large vectors

numcores <- detectCores()-1

parapply_test <- function(cl, test_rets, test_coefs){
  cl <- makePSOCKcluster(numcores)
  clusterExport(cl, list("test_rets", "test_coefs"))
  sim_test <- matrix(0L, nrow = 2e5, ncol = 10)
  for (i in 1:10) sim_test[,i] <- parApply(cl, test_rets, 1, function(x)test_coefs[i,1]+test_coefs[i,-1]%*%x)
  stopCluster(cl)
}


apply_test <- function(constant, test_rets, test_coefs){
  sim_test <- matrix(0L, nrow = 2e5, ncol = 10)
  for (i in 1:10) sim_test[,i] <- constant + apply(test_rets, 1, function(x)test_coefs[i,1]+test_coefs[i,-1]%*%x)
  sim_test
}

columnwise_sum <- function(constant, test_rets, test_coefs){
  sim_test <- matrix(0L, nrow = 2e5, ncol = 10)
  for (i in 1:10) sim_test[,i] <- constant +test_coefs[i,1]+test_coefs[i,2]*test_rets[,1]+
      test_coefs[i,3]*test_rets[,2]+test_coefs[i,4]*test_rets[,3]+test_coefs[i,5]*test_rets[,4]
  sim_test
}




a <- 3
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
plot(mbm_apply)
# apply faster than parapply in this case, probably due to communication time between threads when parallelizingmbm
# columnwise summing almost a hundred times faster


in_loop <- function(constant, test_rets, test_coefs, test_errors){
  sim_test <- matrix(0L, nrow = 2e5, ncol = 10)
  for (i in 1:10) sim_test[,i] <- constant+ test_coefs[i,1]+test_coefs[i,2]*test_rets[,1]+
      test_coefs[i,3]*test_rets[,2]+test_coefs[i,4]*test_rets[,3]+test_coefs[i,5]*test_rets[,4]+test_errors[,i]
  sim_test
}
out_loop <- function(constant, test_rets, test_coefs, test_errors){
  sim_test <- matrix(0L, nrow = 2e5, ncol = 10)
  for (i in 1:10) sim_test[,i] <- constant+ test_coefs[i,1]+test_coefs[i,2]*test_rets[,1]+
      test_coefs[i,3]*test_rets[,2]+test_coefs[i,4]*test_rets[,3]+test_coefs[i,5]*test_rets[,4]
  sim_test+test_errors
}

# cppFunction("NumericMatrix columnwise_sum_cpp(NumericVector constant, NumericMatrix rets, NumericMatrix coefs, NumericMatrix sim_errors, int N_Sim){
#   NumericMatrix sim_ret(N_Sim, 10);
#   int n = 10;
#   NumericVector unit(N_Sim, 1.0);
#   for (int i=0; i<n; i++){
#     sim_ret(_,i) = constant + coefs(i,0)*unit+ coefs(i,1)*rets(_,0)+coefs(i,2)*rets(_,1)+coefs(i,3)*rets(_,2)+coefs(i,4)*rets(_,3)+sim_errors(_,i);
#   }
#   return sim_ret;
# }")

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
plot(mbm_loop)
# outside of loop is faster than inside the loop
# c++ 7x faster than columnwise implementation w/ summation of error matrix outside the loop

sim_test <- columnwise_sum_cpp(rep(a, 2e5), test_rets, test_coefs, test_errors, 2e5)

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
plot(mbm_pf)
# rowMeans the fastest; when testing it was even faster than a quick C++ implementation I did

# Hence, for calculating the returns of the stocks, the C++ function columnwise_sum_cpp will be used
# To get the portfolio returns from the stock returns, rowMeans will be used

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
rownames(coefs_df) <- c("KO", "XOM", "GE", "IBM", "CVX", "UTC", "PG", "CAT", "BA", "MRK")
colnames(error_df) <- c("Date", "KO", "XOM", "GE", "IBM", "CVX", "UTC", "PG", "CAT", "BA", "MRK")
colnames(coefs_df) <- c("Intercept", "Market", "Size", "Value", "Momentum")
head(coefs_df); head(error_df)

## Plot and investigate error distribution
error_df %>% 
  select(-Date) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  geom_histogram(aes(y = ..density..), bins = 80, fill = "grey69", color = "white")+
  geom_density()+
  facet_wrap(~key, scale = "free_x")
# very long left tail; unimodal w/ mode around zero

apply(error_mat, 2, min)
apply(error_mat, 2, max)
# minimal values are a lot smaller in absolute value than maximal values for all shares
# can observe that the bootstrapped error distributions have a long left tail

if (!require(pheatmap)) install.packages("pheatmap")
error_df[,-1] %>% 
  cor() %>% 
  pheatmap(display_numbers = TRUE)
# most residuals only barely correlated
# XOM& CVX residuals show strongest correlation, followed by RTX& BA and KO&PG

if (!require(corrr)) install.packages("corrr")
error_df %>% 
  select(-Date) %>% 
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
  geom_histogram(aes(y = ..density..), bins = 80, fill = "grey69", color = "white")+
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
# to make sure dependencies are kept, all elements in a row are from the same time t
set.seed(42)
N_boot <- 200000
bootind <- sample.int(n_dates, size = N_boot, replace = TRUE)
error_vec_resampled <- error_mat[bootind,] # now without date
head(error_vec_resampled)

# TODO check ARMA model for multivariate case

numcores <- detectCores()-1 # all cores available except one

N_sim <- 200000
n_ahead <- 1
meanModel <- list(armaOrder = c(3, 0))
varModel <- list(model = "fGARCH", submodel = "NGARCH", garchOrder = c(1,1))
uspec <- ugarchspec(varModel, mean.model = meanModel, distribution.model = "sstd")
mspec <- multispec(replicate(4, uspec))
dcc_spec <- dccspec(mspec, VAR = FALSE, model = "DCC", dccOrder = c(1,1), distribution =  "mvnorm")

library(expm)

n_window <- length(FFCFactors_df[,1])-1000

VaR_cop_norm <- matrix(0L, nrow = n_window, ncol = 2)

for (i in 1:n_window){
  cl <- makePSOCKcluster(numcores)
  garch_dcc_fit <- dccfit(dcc_spec, data = Factors_ts[i:(1000+i-1),],cluster = cl)
  garch_dcc_fcst <- dccforecast(garch_dcc_fit, cluster = cl)
  stopCluster(cl)
  
  dcc_fcst_cov <- matrix(garch_dcc_fcst@mforecast$H[[1]], nrow = 4) # or rcov(garch_dcc_fcst)  
  dcc_fcst_mu <- matrix(garch_dcc_fcst@mforecast$mu, nrow = 1)
  garch_dcc_res <- garch_dcc_fit@mfit$stdresid
  res1 <- garch_dcc_res[,1]
  res2 <- garch_dcc_res[,2]
  res3 <- garch_dcc_res[,3]
  res4 <- garch_dcc_res[,4]
  pobs_res <- apply(garch_dcc_res, 2, function(x) copula::pobs(x))
  # par(mfrow = c(2,2))
  # for (i in 1:4)plot(garch_dcc_res[,i], type = "l")
  cop_n <- fitCopula(normalCopula(dim = 4), data = pobs_res)
  
  # simulation:
  cop_sim <- rCopula(N_sim, cop_n@copula)
  cop_sim_df <- data.frame(cop_sim)
  
  res_sim <- cbind(qnorm(cop_sim[,1]), qnorm(cop_sim[,2]), qnorm(cop_sim[,3]), qnorm(cop_sim[,4]))
  res_sim_df <- as.matrix(res_sim)
  #head(res_sim)
  # par(mfrow = c(2,2))
  # for (i in 1:4) plot(res_sim[1:100,i], type = "l")
  # for qt: df = length(resx)-1
  
  
  # use Chol if we want X'*X=A or X*X'=A; sqrtm if we want X*X = A
  # returns are X_t = mu_t+sigma_t*epsilon_t
  logret <- matrix(0L, nrow = N_sim, ncol=4)
  sqrt_h <- sqrtm(dcc_fcst_cov)
  #sqrt_h
  chol_h <- chol(dcc_fcst_cov)
  #logret <- matrix(rep(dcc_fcst_mu, each = N_sim), ncol = 4)+t(sqrt_h%*%t(res_sim))
  logret <- matrix(rep(dcc_fcst_mu, each = N_sim), ncol = 4)+t(chol_h%*%t(res_sim))
  
  # hist(logret[,1], breaks = 50)
  # dim(logret)
  
  # TODO: cholesky vs sqrt.m
  # cholesky smaller VaR than sqrt.m
  
  #### TODO: insert time for RF
  
  set.seed(i)
  bootind <- sample.int(n_dates, size = N_boot, replace = TRUE)
  error_vec_resampled <- error_mat[bootind,] 
  
  sim_rets <- columnwise_sum_cpp(rep(FFCFactors_mat[1000+i-1,1], 2e5), 
                                 logret, coefs_mat, error_vec_resampled, 2e5)
  # calculate portfolio log returns for equally weighted portfolio
  sim_plrets <- rowMeans(sim_rets)
  
  VaR_cop_norm[i,] <- quantile(sim_plrets, c(0.01, 0.05))
  message("completed: ", i, " of ", n_window)
}
VaR_cop_norm_df <- data.frame(Date = portfolio_plret_df$Date[-c(1:1000)], alpha_0.01 = VaR_cop_norm[,1], alpha_0.05 = VaR_cop_norm[,2])
write.csv(VaR_cop_norm, "Data\\VaR\\Multi_cop_norm_VaR.csv", row.names = FALSE)

VaRplot(0.05, portfolio_plret_ts[-c(1:1000)], VaR_cop_norm[,2])
VaRTest(0.05, portfolio_plret_ts[-c(1:1000)], VaR_cop_norm[,2])

VaRplot(0.01, portfolio_plret_ts[-c(1:1000)], VaR_cop_norm[,1])
VaRTest(0.01, portfolio_plret_ts[-c(1:1000)], VaR_cop_norm[,1])

length(VaR_cop_norm)
length(portfolio_plret_ts[-c(1:1000)])



VaR_cop_t <- matrix(0L, nrow = n_window, ncol = 2)
dcc_spec_t <- dccspec(mspec, VAR = FALSE, model = "DCC", dccOrder = c(1,1), distribution =  "mvt")

for (i in 1:n_window){
  cl <- makePSOCKcluster(numcores)
  garch_dcc_fit <- dccfit(dcc_spec_t, data = Factors_ts[i:(1000+i-1),], cluster = cl)
  garch_dcc_fcst <- dccforecast(garch_dcc_fit, cluster = cl)
  stopCluster(cl)
  
  dcc_fcst_cov <- matrix(garch_dcc_fcst@mforecast$H[[1]], nrow = 4) # or rcov(garch_dcc_fcst)  
  dcc_fcst_mu <- matrix(garch_dcc_fcst@mforecast$mu, nrow = 1)
  garch_dcc_res <- garch_dcc_fit@mfit$stdresid
  res1 <- garch_dcc_res[,1]
  res2 <- garch_dcc_res[,2]
  res3 <- garch_dcc_res[,3]
  res4 <- garch_dcc_res[,4]
  pobs_res <- apply(garch_dcc_res, 2, function(x) copula::pobs(x))
  # par(mfrow = c(2,2))
  # for (i in 1:4)plot(garch_dcc_res[,i], type = "l")
  cop_t <- fitCopula(tCopula(dim = 4, df.fixed = TRUE), data = pobs_res)
  
  # simulation:
  cop_sim <- rCopula(N_sim, cop_t@copula)
  cop_sim_df <- data.frame(cop_sim)
  
  res_sim <- cbind(qt(cop_sim[,1], df = 4), qt(cop_sim[,2], df = 4), qt(cop_sim[,3], df = 4), qt(cop_sim[,4], df = 4))
  res_sim_df <- as.matrix(res_sim)
  #head(res_sim)
  # par(mfrow = c(2,2))
  # for (i in 1:4) plot(res_sim[1:100,i], type = "l")
  # for qt: df = length(resx)-1
  
  
  # use Chol if we want X'*X=A or X*X'=A; sqrtm if we want X*X = A
  # returns are X_t = mu_t+sigma_t*epsilon_t
  logret <- matrix(0L, nrow = N_sim, ncol=4)
  sqrt_h <- sqrtm(dcc_fcst_cov)
  #sqrt_h
  chol_h <- chol(dcc_fcst_cov)
  #logret <- matrix(rep(dcc_fcst_mu, each = N_sim), ncol = 4)+t(sqrt_h%*%t(res_sim))
  logret <- matrix(rep(dcc_fcst_mu, each = N_sim), ncol = 4)+t(chol_h%*%t(res_sim))
  
  # hist(logret[,1], breaks = 50)
  # dim(logret)
  
  # TODO: cholesky vs sqrt.m
  # cholesky smaller VaR than sqrt.m
  
  #### TODO: insert time for RF
  
  set.seed(i)
  bootind <- sample.int(n_dates, size = N_boot, replace = TRUE)
  error_vec_resampled <- error_mat[bootind,] 
  
  sim_rets <- columnwise_sum_cpp(rep(FFCFactors_mat[1000+i-1,1], 2e5), 
                                 logret, coefs_mat, error_vec_resampled, 2e5)
  # calculate portfolio log returns for equally weighted portfolio
  sim_plrets <- rowMeans(sim_rets)
  
  VaR_cop_t[i,] <- quantile(sim_plrets, c(0.01, 0.05))
  message("completed: ", i, " of ", n_window)
}
VaR_cop_t_df <- data.frame(Date = portfolio_plret_df$Date[-c(1:1000)], alpha_0.01 = VaR_cop_t[,1], alpha_0.05 = VaR_cop_t[,2])
write.csv(VaR_cop_norm, "Data\\VaR\\Multi_cop_t_VaR.csv", row.names = FALSE)
