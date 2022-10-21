#TODO: 1) OLS regression

#TODO: 2) bootstrap OLS residuals

#TODO: 3) estimate distribution F of factor returns via Copula GARCH

#TODO: 4) simulate rt from F& errort from Fe

#TODO: 5) calculate stock returns w/FFC

#TODO:6) calculate portfolio return

#TODO: VaR_t^p=-percentile(rpf,100p)

# TODO: DCC which distribution? how to backtransform? 

if (!require(tidyverse)) install.packages("tidyverse")
if (!require(rugarch)) install.packages("rugarch")
if (!require(rmgarch)) install.packages("rmgarch")
if (!require(parallel)) install.packages("parallel")
if (!require(copula)) install.packages("copula")
theme_set(theme_light())

#--------------------------------------------------------#
########### Import Data and Create xts Objects ###########
#--------------------------------------------------------#

## Import data from csv files created in EDA.R
FFCFactors_df <- read.csv("./Data/FFCFactors.csv", header = TRUE)
stocks_plret_df <- read.csv("./Data/StockPlrets.csv", header = TRUE)
portfolio_plret_df <- data.frame(Date = stocks_plret_df$Date, Portfolio = apply(stocks_plret_df[,-1], 1, mean))

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
if (!require(Rcpp)) install.packages("Rcpp")

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
# use c++ function to calculate dot product of large vectors

numcores <- detectCores()-1

parapply_test <- function(cl, test_rets, test_coefs){
  sim_test <- matrix(0L, nrow = N_sim, ncol = 10)
  for (i in 1:10) sim_test[,i] <- parApply(cl, test_rets, 1, function(x)test_coefs[i,1]+test_coefs[i,-1]%*%x)
}

apply_test <- function(test_rets, test_coefs){
  sim_test <- matrix(0L, nrow = N_sim, ncol = 10)
  for (i in 1:10) sim_test[,i] <- apply(test_rets, 1, function(x)test_coefs[i,1]+test_coefs[i,-1]%*%x)
}

test_rets <- matrix(rnorm(4*200000), nrow = 200000, ncol = 4)
test_errors <- matrix(rnorm(10*200000), nrow = 200000, ncol = 10)
test_coefs <- matrix(rnorm(5*10), nrow = 10, ncol = 5)


cl <- makePSOCKcluster(numcores)
clusterExport(cl, list("test_rets", "test_coefs", "i"))
mbm_apply <- microbenchmark(
  parapply = parapply_test(cl, test_rets, test_coefs),
  apply = apply_test(test_rets, test_coefs),
  times = 10
)
stopCluster(cl)
mbm_apply
# apply faster than parapply in this case

apply_test_inloop <- function(constant, test_rets, test_coefs, test_errors){
  sim_test <- matrix(0L, nrow = N_sim, ncol = 10)
  for (i in 1:10) sim_test[,i] <- constant + apply(test_rets, 1, function(x)test_coefs[i,1]+test_coefs[i,-1]%*%x)+test_errors[,i]
}

apply_test_outloop <- function(constant, test_rets, test_coefs, test_errors){
  sim_test <- matrix(0L, nrow = N_sim, ncol = 10)
  for (i in 1:10) sim_test[,i] <- constant + apply(test_rets, 1, function(x)test_coefs[i,1]+test_coefs[i,-1]%*%x)
  sim_test+test_errors
}

mbm_loop <- microbenchmark(
  inloop = apply_test_inloop(a, test_rets, test_coefs, test_errors),
  outloop = apply_test_outloop(a, test_rets, test_coefs, test_errors),
  times = 10
)
mbm_loop
# out of loop faster

sim_test <- matrix(0L, nrow = N_sim, ncol = 10)
for (i in 1:10) sim_test[,i] <- constant + apply(test_rets, 1, function(x)test_coefs[i,1]+test_coefs[i,-1]%*%x)
  sim_test+test_errors

  

cl <- makePSOCKcluster(numcores)
clusterExport(cl, list("sim_test"))
mbm_pf <- microbenchmark(
  apply = apply(sim_test, 1, mean),
  parallel = parApply(cl, sim_test, 1, mean)
)
stopCluster(cl)
# parallel faster for mean

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
rownames(coefs_df) <- c("KO", "XOM", "GE", "IBM", "CVX", "RTX", "PG", "CAT", "BA", "MRK")
colnames(error_df) <- c("Date", "KO", "XOM", "GE", "IBM", "CVX", "RTX", "PG", "CAT", "BA", "MRK")
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

#-------------------------------------------------------------------------------#
########### Bootstrap Resample Error Distributions of OLS Regressions ###########
#-------------------------------------------------------------------------------#

# Each row in error_vec_resampled is one draw of a vector of error terms
# to make sure dependencies are kept, all elements in a row are from the same time t
set.seed(42)
N_boot <- 200000
bootind <- sample.int(n_dates, size = N_boot, replace = TRUE)
head(error_df[bootind,])
error_vec_resampled <- error_df[bootind,-1] # now without date
head(error_vec_resampled)

# TODO check ARMA model for multivariate case

numcores <- detectCores()-1 # all cores available except one

N_sim <- 200000
n_ahead <- 1
n <- length(stocks_plret_df[,-1])-n_ahead+1
meanModel <- list(armaOrder = c(3, 0))
varModel <- list(model = "fGARCH", submodel = "NGARCH", garchOrder = c(1,1))
uspec <- ugarchspec(varModel, mean.model = meanModel, distribution.model = "sstd")
mspec <- multispec(replicate(4, uspec))
dcc_spec <- dccspec(mspec, VAR = FALSE, model = "DCC", dccOrder = c(1,1), distribution =  "mvnorm")
cl <- makePSOCKcluster(numcores)
garch_dcc_fit <- dccfit(dcc_spec, data = Factors_ts ,cluster = cl)


loglik_dcc <- garch_dcc_fit@mfit$log.likelihoods
AICdcc <- infocriteria(garch_dcc_fit)[1]


garch_dcc_fcst <- dccforecast(garch_dcc_fit, cluster = cl)
stopCluster(cl)
dcc_fcst_cov <- matrix(garch_dcc_fcst@mforecast$H[[1]], nrow = 4) # or rcov(garch_dcc_fcst)  
dcc_fcst_mu <- matrix(garch_dcc_fcst@mforecast$mu, nrow = 1)
garch_dcc_res <- garch_dcc_fit@mfit$stdresid
head(garch_dcc_res)
res1 <- garch_dcc_res[,1]
res2 <- garch_dcc_res[,2]
res3 <- garch_dcc_res[,3]
res4 <- garch_dcc_res[,4]
pobs_res <- apply(garch_dcc_res, 2, function(x) copula::pobs(x))
par(mfrow = c(2,2))
for (i in 1:4)plot(garch_dcc_res[,i], type = "l")
cop_n <- fitCopula(normalCopula(dim = 4), data = pobs_res)

# simulation:
cop_sim <- rCopula(N_sim, cop_n@copula)
cop_sim_df <- data.frame(cop_sim)

res_sim <- cbind(qnorm(cop_sim[,1]), qnorm(cop_sim[,2]), qnorm(cop_sim[,3]), qnorm(cop_sim[,4]))
res_sim_df <- as.data.frame(res_sim)
head(res_sim)
par(mfrow = c(2,2))
for (i in 1:4) plot(res_sim[1:100,i], type = "l")
# for qt: df = length(resx)-1

library(expm)
# use Chol if we want X'*X=A or X*X'=A; sqrtm if we want X*X = A
# returns are X_t = mu_t+sigma_t*epsilon_t
logret <- matrix(0L, nrow = N_sim, ncol=4)
sqrt_h <- sqrtm(dcc_fcst_cov)
sqrt_h
logret <- matrix(rep(dcc_fcst_mu, each = N_sim), ncol = 4)+t(sqrt_h%*%t(res_sim))

hist(logret[,1], breaks = 50)
dim(logret)


#### TODO: insert time for RF



sim_rets <- matrix(0L, nrow = N_sim, ncol = 10)
for (i in 1:10) sim_rets[,i] <- apply(logret, 1, function(x)coefs_mat[i,1]+coefs_mat[i,-1]%*%x)
sim_rets <- FFCFactors_df$RF[time]+sim_rets+error_vec_resampled
cl <- makePSOCKcluster(numcores)
clusterExport(cl, "logret")
sim_plrets <- parApply(cl, sim_rets, 1, mean)
stopCluster(cl)
# crossprod slightly faster than x%*%y

quantile(sim_plrets, c(0.01, 0.05))
