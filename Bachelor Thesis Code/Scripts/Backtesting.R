#=========================================================================#
############################### Backtesting ###############################
#=========================================================================#

# TODO: ranking of mean loss
# TODO: CPA test matrix
# TODO: VaR list implementation for LR tests
# TODO: check p_value corrections
# TODO: check whether p or 1-p for coverage etc.

#------------------------------------#
########### Importing Data ###########
#------------------------------------#

## Portfolio Plrets
stocks_plret_df <- read.csv("./Data/StockPlrets.csv", header = TRUE)
portfolio_plret_df <- data.frame(Date = stocks_plret_df$Date, 
                                 Portfolio = rowMeans(stocks_plret_df[,-1]))

## VaR
Uni_Normal_GARCH_VaR <- read.csv("./Data/VaR/Uni_Normal_GARCH_VaR.csv", 
                                 header = TRUE)
Uni_EWMA_VaR <- read.csv("./Data/VaR/Uni_EWMA_VaR.csv", 
                         header = TRUE)
Uni_t_GJR_GARCH_VaR <- read.csv("./Data/VaR/Uni_t_GJR_GARCH.csv", 
                                header = TRUE)


#-----------------------------------------------------------------------------#
###### Tests for Independence and Conditional and Unconditional Coverage ######
#-----------------------------------------------------------------------------#

#' Test of unconditional coverage
#'
#' @param p VaR percentile e.g. 1%
#' @param VaR VaR forecasts of a model (only the column with p% VaR values)
#' @param plrets portfolio returns dataframe with dates in first column& returns
#'in second column
#'
#' @return test statistic of likelihood ratio test of unconditional coverage
LR_uc <- function(p, VaR, plrets = portfolio_plret_df[-c(1:1000),2]){
  indicator <- ifelse(plrets-VaR<0, 1, 0)
  n1 <- sum(indicator)
  n0 <- length(VaR) - n1
  
  lik_p <- (1 - p)^n0 * p^n1
  
  pi_mle <- n1 / (n0 + n1)
  lik_pi_mle <- (1 - pi_mle)^n0 * pi_mle^n1
  
  LR <- -2 * log(lik_p / lik_pi_mle)
  return(LR)
}

uc <- LR_uc(0.01, Uni_t_GJR_GARCH_VaR[,2])
ugarch_uc <- VaRTest(alpha = 0.01, portfolio_plret_df[-c(1:1000),2], 
                     Uni_t_GJR_GARCH_VaR[,2], conf.level = 0.95)$uc.LRstat
uc==ugarch_uc
# get same value as rugarch implementation


#' Test of independence
#'
#' @param p VaR percentile e.g. 1%
#' @param VaR VaR forecasts of a model (only the column with p% VaR values)
#' @param plrets portfolio returns dataframe with dates in first column& returns
#'in second column
#'
#' @return test statistic of likelihood ratio test of independence
LR_ind <- function(p, VaR, plrets = portfolio_plret_df[-c(1:1000),2]){
  indicator <- as.numeric(ifelse(plrets-VaR<0, 1, 0))
  tab <- table(indicator[-length(indicator)], indicator[-1])
  n00 <- tab[1,1]
  n01 <- tab[1,2]
  n10 <- tab[2,1]
  n11 <- tab[2,2]
  
  pi_MLE_01 <- n01/(n00+n01)
  pi_MLE_11 <- n11/(n10+n11)
  lik_Pi1_MLE <- (1 - pi_MLE_01)^n00 * pi_MLE_01^n01 * (1 - pi_MLE_11)^n10 * pi_MLE_11^n11
  
  pi2_MLE <- (n01 + n11) / sum(tab)
  lik_Pi2_MLE <- (1 - pi2_MLE)^(n00 + n10) * pi2_MLE^(n01 + n11)
  
  LR <- -2 * log(lik_Pi2_MLE / lik_Pi1_MLE)
  return(LR)
}


#' Test of conditional coverage
#'
#' @param p VaR percentile e.g. 1%
#' @param VaR VaR forecasts of a model (only the column with p% VaR values)
#' @param plrets portfolio returns dataframe with dates in first column& returns
#'in second column
#'
#' @return test statistic of likelihood ratio test of coonditional coverage
LR_cc <- function(p, VaR, plrets = portfolio_plret_df[-c(1:1000),2]){
  uc <- LR_uc(p, VaR)
  ind <- LR_ind(p, VaR)
  LR <- uc + ind
  return(LR)
}
cc <- LR_cc(0.01, Uni_t_GJR_GARCH_VaR[, 2])
cc==VaRTest(alpha = 0.01, portfolio_plret_df[-c(1:1000),2], 
            Uni_t_GJR_GARCH_VaR[,2], conf.level = 0.95)$cc.LRstat
# get same value as rugarch implementation


## Create class to return separate list for each test
setClass(Class="LR_tests",
         representation(
           cc  = "list",
           ind = "list",
           uc  = "list"
         )
)


#' Test of unconditional coverage
#'
#' Implements backtesting as described in Christoffersen (1998) i.e. implements
#' the LR test of unconditional coverage, the LR test of independence and the LR
#' test of conditional coverage.
#'
#' @param p VaR percentile e.g. 1%
#' @param VaR VaR forecasts of a model (only the column with p% VaR values)
#' @param plrets portfolio returns dataframe with dates in first column& returns
#'in second column
#' @param conf_level the confidence level of the test
#'
#' @return returns instance of class "LR_tests" i.e. a list for each of the three
#' tests that includes the critical value, the test statistic, the p-value and
#' the decision i.e. reject or not
VaR_LR_tests <- function(p, VaR, plrets = portfolio_plret_df[-c(1:1000),2], conf_level = 0.95){
  LR_uc <- LR_uc(p, VaR)
  LR_ind <- LR_ind(p, VaR)
  LR_cc <- LR_cc(p, VaR)
  
  crit_val_uc <- crit_val_ind <- qchisq(conf_level, df = 1)
  crit_val_cc <- qchisq(conf_level, df = 2)
  
  p_val_uc <- 1 - pchisq(LR_uc, df = 1)
  p_val_ind <- 1 - pchisq(LR_ind, df = 1)
  p_val_cc <- 1 - pchisq(LR_cc, df = 2)
  
  reject_uc <- ifelse(p_val_uc < 1 - conf_level, TRUE, FALSE)
  reject_ind <- ifelse(p_val_ind < 1 - conf_level, TRUE, FALSE)
  reject_cc <- ifelse(p_val_cc < 1 - conf_level, TRUE, FALSE)
  
  return(new("LR_tests",
             cc  = list(crit_val_cc = crit_val_cc, LR_cc = LR_cc, p_val_cc = p_val_cc, reject_cc = reject_cc),
             ind = list(crit_val_ind = crit_val_ind, LR_ind = LR_ind, p_val_ind = p_val_ind, reject_ind = reject_ind),
             uc  = list(crit_val_uc = crit_val_uc, LR_uc = LR_uc, p_val_uc = p_val_uc, reject_uc = reject_uc)))
}

VaR_LR_tests(0.01, Uni_t_GJR_GARCH_VaR[, 2])


#----------------------------------------------------------------#
########### Calculate Exceedances and Nominal Coverage ###########
#----------------------------------------------------------------#

#' Nominal Coverage
#' 
#' Calculates and returns the nominal coverage
#'
#' @param VaR Value at risk forecasts of a model
#' @param plrets portfolio returns dataframe with dates in first column& returns
#' in second column
#'
#' @return nominal coverage
nominal_coverage <- function(VaR, plrets = portfolio_plret_df[-c(1:1000),2]){
  indicator <- ifelse(plrets-VaR<0, 1, 0)
  coverage <- sum(indicator)/length(VaR)
  return(nominal_coverage = coverage)
}

nominal_coverage(Uni_EWMA_VaR[,2])

if (!require(tidyverse)) install.packages("tidyverse")
if (!require(lubridate)) install.packages("lubridate") # to extract year from Date

#' Exceedances
#' 
#' Calculates and returns the total number of exceedances as well as the exceedances
#' per year
#'
#' @param VaR Value at risk forecasts of a model
#' @param plrets portfolio returns dataframe with dates in first column& returns
#'in second column
#'
#' @return list with total number of exceedances and exceedences per year
exceedances <- function(VaR, plrets = portfolio_plret_df[-c(1:1000),]){
  indicator <- ifelse(plrets[,2]-VaR<0, 1, 0)
  indicator_df <- data.frame(Date = plrets[,1], Exceedance = as.factor(indicator))
  exc_per_year  <- indicator_df %>% 
    mutate(year = year(Date)) %>% 
    select(Exceedance, year) %>% 
    count(year, Exceedance) %>% 
    mutate(n = ifelse(Exceedance==1, n, 0)) %>% 
    select(-Exceedance) %>% 
    group_by(year) %>% 
    summarise(n = sum(n))
  return(list(total_exc = sum(indicator), exc_per_year = exc_per_year))
}

exceedances(Uni_EWMA_VaR[,2])




#' Exceedances Table
#'
#' Create tables for the 1% VaR and the 5% VaR that consist of the total number
#' of exceedances and the exceedences per year 
#'
#' @param VaR_list list of VaR dataframes with date in first column, 1% VaR in
#' second column and 5% VaR in third column
#' @param plrets portfolio returns dataframe with dates in first column& returns
#'in second column
#'
#' @return list of exceedance tables for 1% VaR and 5% VaR
exceedances_table <- function(VaR_list, plrets = portfolio_plret_df[-c(1:1000),]){
  n <- length(VaR_list)
  matrix_99 <- matrix(0L, nrow = n, ncol = 9)
  for (i in 1:n){
    matrix_99[i, 1] <- unlist(exceedances(VaR_list[[i]][,2])$total_exc)
    matrix_99[i, 2:9] <- unlist(exceedances(VaR_list[[i]][,2])$exc_per_year[,2])
  }
  table_99 <- data.frame(matrix_99)
  colnames(table_99) <- c("Total_Exc", "2004", "2005", "2006", "2007",
                          "2008", "2009", "2010", "2011")
  rownames(table_99) <- names(VaR_list)
  
  matrix_95 <- matrix(0L, nrow = n, ncol = 9)
  for (i in 1:n){
    matrix_95[i, 1] <- unlist(exceedances(VaR_list[[i]][,3])$total_exc)
    matrix_95[i, 2:9] <- unlist(exceedances(VaR_list[[i]][,3])$exc_per_year[,2])
  }
  table_95 <- data.frame(matrix_95)
  colnames(table_95) <- c("Total_Exc", "2004", "2005", "2006", "2007",
                          "2008", "2009", "2010", "2011")
  rownames(table_95) <- names(VaR_list)
  
  return(list(table_99 = table_99, table_95 = table_95))
}

test_VaR_list <- list(EWMA = Uni_EWMA_VaR, Normal_GARCH = Uni_Normal_GARCH_VaR,
                      t_GJR = Uni_t_GJR_GARCH_VaR)
exceedances_table(test_VaR_list)$table_99
exceedances_table(test_VaR_list)$table_95


#---------------------------------------------------------------#
########### Table for Exceedances and LR Test Pvalues ###########
#---------------------------------------------------------------#

#' Performance table
#' 
#' Create a summary table for backtesting that includes nominal coverage,
#' exceedances (over the years) and LR tests of Christoffersen (1998). 
#'
#' @param VaR_list list of VaR dataframes with date in first column, 1% VaR in
#' second column and 5% VaR in third column
#' @param plrets portfolio returns dataframe with dates in first column& returns
#'in second column
#' @param conf_level confidence level for LR tests of Christoffersen (1998). By
#' default 95%
#'
#' @return list with performance table for 1% VaR and for 5% VaR
performance_table <- function(VaR_list, plrets = portfolio_plret_df[-c(1:1000),],
                              conf_level = 0.95){
  n <- length(VaR_list)
  coverage_99 <- matrix(0L, nrow = n, ncol = 1)
  tests_99 <- matrix(0L, nrow = n, ncol = 3)
  for (i in 1:n){
    coverage_99[i, 1] <- nominal_coverage(VaR_list[[i]][,2], plrets[,2])
    tests_99[i, 1] <- unlist(VaR_LR_tests(0.01, VaR_list[[i]][,2], plrets[,2])@uc$p_val_uc)
    tests_99[i, 2] <- unlist(VaR_LR_tests(0.01, VaR_list[[i]][,2], plrets[,2])@ind$p_val_ind)
    tests_99[i, 3] <- unlist(VaR_LR_tests(0.01, VaR_list[[i]][,2], plrets[,2])@cc$p_val_cc)
  }
  exceed_99 <- exceedances_table(VaR_list, plrets)$table_99
  performance_table_99 <- data.frame(coverage_99, tests_99, exceed_99)
  colnames(performance_table_99) <- c("coverage_99", "uc", "ind", "cc", 
                                      "Total_Exc", "2004", "2005", "2006", 
                                      "2007", "2008", "2009", "2010", "2011")
  rownames(performance_table_99) <- names(VaR_list)
  
  coverage_95 <- matrix(0L, nrow = n, ncol = 1)
  tests_95 <- matrix(0L, nrow = n, ncol = 3)
  for (i in 1:n){
    coverage_95[i, 1] <- nominal_coverage(VaR_list[[i]][,3], plrets[,2])
    tests_95[i, 1] <- unlist(VaR_LR_tests(0.05, VaR_list[[i]][,3], plrets[,2])@uc$p_val_uc)
    tests_95[i, 2] <- unlist(VaR_LR_tests(0.05, VaR_list[[i]][,3], plrets[,2])@ind$p_val_ind)
    tests_95[i, 3] <- unlist(VaR_LR_tests(0.05, VaR_list[[i]][,3], plrets[,2])@cc$p_val_cc)
  }
  exceed_95 <- exceedances_table(VaR_list, plrets)$table_95
  performance_table_95 <- data.frame(coverage_95, tests_95, exceed_95)
  colnames(performance_table_95) <- c("coverage_95", "uc", "ind", "cc", 
                                      "Total_Exc", "2004", "2005", "2006", 
                                      "2007", "2008", "2009", "2010", "2011")
  rownames(performance_table_95) <- names(VaR_list)
  
  return(list(performance_table_99 = performance_table_99 %>% round(3), 
              performance_table_95 = performance_table_95 %>% round(3)))
}


performance_table(test_VaR_list)$performance_table_99 
performance_table(test_VaR_list)$performance_table_95



#---------------------------------------------------------------#
########### CPA Test as in Giacomini and White (2006) ###########
#---------------------------------------------------------------#


#' Loss VaR
#' 
#' Calculates loss function for VaR i.e. tick loss function for quantile 
#' regression
#' 
#' @param VaR dataframe w/ dates in first column, 99% VaR in second column and 
#' 95% VaR in third column
#' @param plrets portfolio returns; by default portfolio returns from t=1001 
#' until t=T
#' @param percentile 1- VaR percentiles; by default c(0.99, 0.95)
#'
#' @return list w/ losses for the percentiles as first two elements and 
#' mean losses for the percentiles as last two elements
loss_VaR <- function(VaR, plrets = portfolio_plret_df[-c(1:1000),2], 
                     percentile = c(0.99, 0.95)){
  indicator_99 <- ifelse(plrets-VaR[2]<0, 1, 0)
  indicator_95 <- ifelse(plrets-VaR[3]<0, 1, 0)
  loss_99 <- (plrets-VaR[2])*(percentile[1]-indicator_99)
  loss_95 <- (plrets-VaR[3])*(percentile[2]-indicator_95)
  return(list(
    loss_99=loss_99, 
    loss_95 = loss_95, 
    mean_loss_99 = colMeans(loss_99), 
    mean_loss_95 = colMeans(loss_95)
    ))
}

Uni_t_GJR_loss <- loss_VaR(Uni_t_GJR_GARCH_VaR)
Uni_EWMA_loss <- loss_VaR(Uni_EWMA_VaR)
Uni_EWMA_loss$mean_loss_99
Uni_t_GJR_loss$mean_loss_99


## Compare w/ VaRloss function from rugarch package for e.g. GJR GARCH:
if (!require(rugarch)) install.packages("rugarch")
all.equal(VaRloss(0.95, portfolio_plret_df[-c(1:1000),2], Uni_t_GJR_GARCH_VaR[,3]), 
          as.numeric(as.matrix(100*Uni_t_GJR_loss$loss_95)))
# rugarch VaRloss is 100* the loss calculated above


#' Ranking VaR forecasts
#' 
#' Ranking VaR forecasts in ascending order based on their average VaR/ tick loss.
#'
#' @param VaR_list list of VaR dataframes with date in first column, 1% VaR in
#' second column and 5% VaR in third column
#' @param plrets portfolio returns; by default portfolio returns from t=1001 
#' until t=T
#' @param percentile 1- VaR percentiles; by default c(0.99, 0.95)
#'
#' @return list with tables for 1% and 5% VaR ranking
VaR_loss_ranking <- function(VaR_list, plrets = portfolio_plret_df[-c(1:1000),],
                             percentile = c(0.99, 0.95)){
  n <- length(VaR_list)
  matrix_99 <- matrix(0L, nrow = n, ncol = 1)
  matrix_95 <- matrix(0L, nrow = n, ncol = 1)
  for (i in 1:n){
    matrix_99[i, 1] <- unlist(loss_VaR(VaR_list[[i]])$mean_loss_99)
    matrix_95[i, 1] <- unlist(loss_VaR(VaR_list[[i]])$mean_loss_95)
  }
  table_99 <- data.frame(matrix_99)
  colnames(table_99) <- c("mean_VaR_loss")
  rownames(table_99) <- names(VaR_list)
  table_99 <- table_99 %>% arrange(mean_VaR_loss) # arange in ascending order
  
  
  table_95 <- data.frame(matrix_95)
  colnames(table_95) <- c("mean_VaR_loss")
  rownames(table_95) <- names(VaR_list)
  table_95 <- table_95 %>% arrange(mean_VaR_loss) 
  return(list(table_99 = table_99, table_95 = table_95))
}

VaR_loss_ranking(test_VaR_list)$table_99
VaR_loss_ranking(test_VaR_list)$table_95


## Create class to return two lists in CPA_test function
setClass(Class="CPA",
         representation(
           VaR_99="list",
           VaR_95="list"
         )
)

#' Conditional Predictive Ability Test by Giacomini and White (2006)
#' 
#' Implements CPA test to allow for binary model comparisons in predictive ability
#'
#' @param VaR1 VaR forecasts of model i (whole dataframe)
#' @param VaR2 VaR forecasts of model j (whole dataframe)
#' @param plrets portfolio returns; by default portfolio returns from t=1001 
#' until t=T
#' @param percentile 1- VaR percentiles (by default: 99% and 95%)
#'
#' @return instance of class "CPA" i.e. for the 1% and the 5% VaR it returns
#' a list of the test statistic, pvalue, critical value and the test decision
#' (i.e. which model has higher predictive ability)
CPA_test <- function(VaR1, VaR2, plrets = portfolio_plret_df[-c(1:1000),2], 
                                          percentile = c(0.99, 0.95)){
  
  loss1 <- loss_VaR(VaR1)
  loss2 <- loss_VaR(VaR2)
  
  ## 99%
  d_ij_99 <- as.matrix(loss1$loss_99)-as.matrix(loss2$loss_99)
  T <- length(d_ij_99)
  
  h_tminus1_99 <- matrix(0L, ncol = T-1, nrow = 2)
  for (i in 1:(T-1)) h_tminus1_99[,i] <- rbind(1,d_ij_99[i])
  
  Z_99 <- matrix(0L, ncol =T-1, nrow = 2)
  for (i in 1:(T-1)) Z_99[,i] <- h_tminus1_99[,i]*d_ij_99[i+1]
  
  Z_bar_99 <- matrix(0L, nrow = 2, ncol = 1)
  for (i in 1:2) Z_bar_99[i] <- 1/T*sum(Z_99[i,])
  
  outer_Z_99 <- array(0L, dim = c(2,2,(T-1))) # dim: rows, cols, time
  for (i in 1:(T-1)) outer_Z_99[,,i] <- Z_99[,i]%*%t(Z_99[,i])
  Omega_inv_99 <- matrix(0L, nrow = 2, ncol =2)
  for (i in 1:2)for (j in 1:2) Omega_inv_99[i,j] <- 1/T*sum(outer_Z_99[i,j,])
  
  GW_ij_99 <- T * t(Z_bar_99) %*% Omega_inv_99 %*% Z_bar_99
  
  
  ## 95%
  d_ij_95 <- as.matrix(loss1$loss_95)-as.matrix(loss2$loss_95)
  T <- length(d_ij_95)
  
  h_tminus1_95 <- matrix(0L, ncol = T-1, nrow = 2)
  for (i in 1:(T-1)) h_tminus1_95[,i] <- rbind(1,d_ij_95[i])
  
  Z_95 <- matrix(0L, ncol =T-1, nrow = 2)
  for (i in 1:(T-1)) Z_95[,i] <- h_tminus1_95[,i]*d_ij_95[i+1]
  
  Z_bar_95 <- matrix(0L, nrow = 2, ncol = 1)
  for (i in 1:2) Z_bar_95[i] <- 1/T*sum(Z_95[i,])
  
  outer_Z_95 <- array(0L, dim = c(2,2,(T-1))) # dim: rows, cols, time
  for (i in 1:(T-1)) outer_Z_95[,,i] <- Z_95[,i]%*%t(Z_95[,i])
  Omega_inv_95 <- matrix(0L, nrow = 2, ncol =2)
  for (i in 1:2)for (j in 1:2) Omega_inv_95[i,j] <- 1/T*sum(outer_Z_95[i,j,])
  
  GW_ij_95 <- T * t(Z_bar_95) %*% Omega_inv_95 %*% Z_bar_95
  
  ## P values
  p_val_99 <- 1-pchisq(GW_ij_99,2)
  p_val_95 <- 1-pchisq(GW_ij_95,2)
  
  better_99 <- NULL; worse_99 <- NULL
  better_95 <- NULL; worse_95 <- NULL
  
  c <- 0 # threshold as in Giacomini and White (2006)
  
  if (p_val_99<=0.05){
    signif_99 <- TRUE
    fit_99 <- lm(d_ij_99[-1]~t(h_tminus1_99)[,-1])
    delta_99 <- matrix(coef(fit_99), nrow = 2, ncol = 1)
    indicator_seq_99 <- numeric(0)
    for (i in 1:(T-1)) indicator_seq_99[i] <- ifelse(t(delta_99)%*%h_tminus1_99[,i]>c, 1, 0)
    seq_mean_99 <- 1/T*sum(indicator_seq_99)
    ifelse(seq_mean_99>0.5, 
           c(better_99 <- deparse(substitute(VaR1)), worse_99 <-  deparse(substitute(VaR2))),
           c(worse_99 <- deparse(substitute(VaR1)), better_99 <-  deparse(substitute(VaR2))))
    message(better_99, " significantly outperforms ", worse_99, " for the 99% VaR")
  }
  else{signif_99 <- FALSE}
  
  
  
  if (p_val_95<=0.05){
    signif_95 <- TRUE
    fit_95 <- lm(d_ij_95[-1]~t(h_tminus1_95)[,-1])
    delta_95 <- matrix(coef(fit_95), nrow = 2, ncol = 1)
    indicator_seq_95 <- numeric(0)
    for (i in 1:(T-1)) indicator_seq_95[i] <- ifelse(t(delta_95)%*%h_tminus1_95[,i]>c, 1, 0)
    seq_mean_95 <- 1/T*sum(indicator_seq_95)
    ifelse(seq_mean_95>0.5, 
           c(better_95 <- deparse(substitute(VaR1)), worse_95 <-  deparse(substitute(VaR2))),
           c(worse_95 <- deparse(substitute(VaR1)), better_95 <-  deparse(substitute(VaR2))))
    message(better_95, " significantly outperforms ", worse_95, " for the 95% VaR")
  }
  else{signif_95 <- FALSE}
  
  return(new("CPA", 
             VaR_99 = list(GW_ij_99 = GW_ij_99, p_val_99 = p_val_99, 
                           signif_99 = signif_99, better_99 = better_99, worse_99 = worse_99),
             VaR_95 = list(GW_ij_95 = GW_ij_95, p_val_95 = p_val_95, 
                       signif_95 = signif_95, better_95 = better_95, worse_95 = worse_95)
         ))
}
CPA_test(Uni_t_GJR_GARCH_VaR, Uni_EWMA_VaR)
CPA_test(Uni_t_GJR_GARCH_VaR, Uni_Normal_GARCH_VaR)
CPA_test(Uni_Normal_GARCH_VaR, Uni_EWMA_VaR)


CPA_table <- function()
