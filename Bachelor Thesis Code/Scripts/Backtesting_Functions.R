#==========================================================================#
########################## Backtesting Functions  ##########################
#==========================================================================#

# R version 4.2.2 (2022-10-31 ucrt)

# !!!!! Important Note !!!!! #
# The code for the Conditional Predictive Ability Test by Giacomini and White
# (2006) is based on the Matlab code provided by the authors. I translated it
# to R such that it covers my specific case. Other forecasting horizons etc.
# are not implemented. The Matlab code of the Giacomini and White can be found
# here:
# http://www.runmycode.org/companion/view/88

# I compared my implementation with theirs for different models as well as
# different loss functions and always got the same result

# For the LR tests of Christoffersen I have consulted the GitHub page for the
# rugarch package to compare my own implementation with the rugarch one. I
# wanted to implement it on my own since implementing sth on my own in R
# generally helps me understand concepts better. After comparing my LR tests
# to the rugarch GitHub, I made slight adjustments inspired by the rugarch
# GitHub: 
# 1. using the table function for the LR test of independence
# 2. for numerical reasons: calculating the test statistic of the LR
# test of conditional coverage as the sum of LR_uc and LR_ind instead of 
# calculating it using the corresponding likelihood functions


#------------------------------------#
########### Importing Data ###########
#------------------------------------#

## Portfolio Plrets
stocks_plret_df <- read.csv("./Data/StockPlrets.csv", header = TRUE)
portfolio_plret_df <- read.csv("./Data/PortfolioPlrets.csv", header = TRUE)


## Don't run when importing
if (sys.nframe() == 0) {
  ## VaR
  # Univariate:
  Uni_Normal_GARCH_VaR <- read.csv("./Data/VaR/Uni_Normal_GARCH_VaR.csv", 
                                   header = TRUE)
  Uni_EWMA_VaR <- read.csv("./Data/VaR/Uni_EWMA_VaR.csv", 
                           header = TRUE)
  Uni_t_GJR_GARCH_VaR <- read.csv("./Data/VaR/Uni_t_GJR_GARCH.csv", 
                                  header = TRUE)
  Uni_Skewt_GJR_GARCH_VaR <- read.csv("./Data/VaR/Uni_Skewt_GJR_GARCH.csv", 
                                      header = TRUE)
  Uni_Skewt_NGARCH_VaR <- read.csv("./Data/VaR/Uni_Skewt_NGARCH.csv", 
                                   header = TRUE)
  
  # Multivariate
  Multi_DCC_GARCH_VaR <- read.csv("./Data/VaR/Multi_DCC_GARCH.csv",
                                  header = TRUE)
  
  Multi_Fortin_Normal_VaR <- read.csv("./Data/VaR/Multi_cop_norm_VaR.csv", 
                                      header = TRUE)
  
  
  Multi_Fortin_t_VaR <- read.csv("./Data/VaR/Multi_cop_t_VaR.csv", 
                                 header = TRUE)
}



#--------------------------------------------------#
########### VaR Exceedence Plot Function ###########
#--------------------------------------------------#
# before the more formal tests it is always a good idea to first have a look at 
# a graphical representation


# load lubridate for year() function to extract year from Date
if (!require(lubridate)) install.packages("lubridate")
if (!require(tidyverse)) install.packages("tidyverse")

#' VaR Exceedence Plot
#' 
#' Plots returns and line of VaR and marks every point where the return is 
#' below the VaR line as an exceedence and displays how many exceedences there 
#' are in each year
#' 
#' @param dataframe dataframe with VaR
#' @param VaR_in_col_nr integer indicating in which column of dataframe the VaR
#'  is
#' @param pf_plrets dataframe of portfolio percentage returns
#' @param VaR_percentile VaR level in percent w/o percentage sign
#' @param modelname Name of model as character
#' 
#' @return Exceedences plot which highlights exceedances and
#'  displays number of exceedances per year
VaR_exceed_plot <- function(dataframe, VaR_in_col_nr, pf_plrets, VaR_percentile, 
                            modelname){
  VaR_df <- data.frame(Date = as.Date(dataframe[,1]),
      VaR = dataframe[,VaR_in_col_nr],
      Exceedance = as.factor(pf_plrets[-c(1:1000),2]<dataframe[,VaR_in_col_nr])
    )
  exceedances_per_year  <- VaR_df %>% 
    mutate(year = year(Date)) %>% 
    select(Exceedance, year) %>% 
    count(year, Exceedance) %>% 
    mutate(n = ifelse(Exceedance==TRUE, n, 0)) %>% 
    select(-Exceedance) %>% 
    group_by(year) %>% 
    summarise(n = sum(n))
  
  ggplot(VaR_df, aes(x = Date, y = VaR))+
    geom_point(aes(x = Date, y = pf_plrets[-c(1:1000), 2], 
                   color = Exceedance, shape = Exceedance), size =1.5, 
               alpha = 2)+
    scale_shape_manual(values = c(20, 4), name="",
                       labels = c("Lower than VaR", "Greater than VaR"))+
    scale_color_manual(values = c("gray", "red"), name = "", 
                       labels = c("Lower than VaR", "Greater than VaR"))+
    geom_line(alpha = 0.7)+
    labs(y = "Daily Portfolio Returns", x = "Date", 
         title = paste0(VaR_percentile, "% VaR Exceedances Plot for ",modelname))+
    theme_light()+
    theme(legend.position = c(.15, .8), 
          legend.background = element_rect(color = NA), 
          legend.key = element_rect(color = "transparent"))+
    annotate("text", x = as.Date("2005-01-15"), y = -7, size = 3, hjust = 0,
             label = paste("Number of Exceedances per Year:", "\n2004:", 
                           exceedances_per_year$n[1], "\n2005:", 
                           exceedances_per_year$n[2],
                           "\n2006:", exceedances_per_year$n[3], 
                           "\n2007:", exceedances_per_year$n[4]))+
    annotate("text", x = as.Date("2005-10-15"), y = -7, size = 3, hjust = 0,
             label = paste(" ","\n2008:", exceedances_per_year$n[5],"\n2009:", 
                           exceedances_per_year$n[6],
                           "\n2010:", exceedances_per_year$n[7], "\n2011:", 
                           exceedances_per_year$n[8]))
}

## Don't run when importing
if (sys.nframe() == 0) {
  VaR_exceed_plot(Uni_Normal_GARCH_VaR, 3, portfolio_plret_df, 
                  VaR_percentile = 5,"Uni_Normal_GARCH")
  VaR_exceed_plot(Uni_Normal_GARCH_VaR, 2, portfolio_plret_df,
                  VaR_percentile = 1, "Uni_Normal_GARCH")
}


#-----------------------------------------------------------------------------#
###### Tests for Independence and Conditional and Unconditional Coverage ######
#-----------------------------------------------------------------------------#

## Load rugarch to compare tests to tests implemented in rugarch
if (!require(rugarch)) install.packages("rugarch")

#' Test of unconditional coverage
#'
#' @param p VaR quantile e.g. 1%
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

## Don't run when importing
if (sys.nframe() == 0) {
  uc1 <- LR_uc(0.01, Uni_t_GJR_GARCH_VaR[,2])
  ugarch_uc1 <- VaRTest(alpha = 0.01, portfolio_plret_df[-c(1:1000),2], 
                        Uni_t_GJR_GARCH_VaR[,2], conf.level = 0.95)$uc.LRstat
  uc1==ugarch_uc1
  
  
  uc2 <- LR_uc(0.05, Uni_Normal_GARCH_VaR[,3])
  ugarch_uc2 <- VaRTest(alpha = 0.05, portfolio_plret_df[-c(1:1000),2], 
                        Uni_Normal_GARCH_VaR[,3], conf.level = 0.95)$uc.LRstat
  uc2==ugarch_uc2
  # get same value as rugarch implementation
}



#' Test of independence
#'
#' @param p VaR quantile e.g. 1%
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
  lik_Pi1_MLE <- (1 - pi_MLE_01)^n00 * pi_MLE_01^n01 * 
    (1 - pi_MLE_11)^n10 * pi_MLE_11^n11
  
  pi2_MLE <- (n01 + n11) / sum(tab)
  lik_Pi2_MLE <- (1 - pi2_MLE)^(n00 + n10) * pi2_MLE^(n01 + n11)
  
  LR <- -2 * log(lik_Pi2_MLE / lik_Pi1_MLE)
  return(LR)
}


#' Test of conditional coverage
#'
#' As in rugarch, the cc test statistic is for numerical reasons calculated as
#' the sum of ind& uc test statistics.
#'
#' @param p VaR quantile e.g. 1%
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

## Don't run when importing
if (sys.nframe() == 0) {
  cc1 <- LR_cc(0.01, Uni_t_GJR_GARCH_VaR[, 2])
  cc1==VaRTest(alpha = 0.01, portfolio_plret_df[-c(1:1000),2], 
               Uni_t_GJR_GARCH_VaR[,2], conf.level = 0.95)$cc.LRstat
  
  cc2 <- LR_cc(0.05, Uni_Normal_GARCH_VaR[, 3])
  cc2==VaRTest(alpha = 0.05, portfolio_plret_df[-c(1:1000),2], 
               Uni_Normal_GARCH_VaR[,3], conf.level = 0.95)$cc.LRstat
}


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
#' @param p VaR quantile e.g. 1%
#' @param VaR VaR forecasts of a model (only the column with p% VaR values)
#' @param plrets portfolio returns dataframe with dates in first column& returns
#'in second column
#' @param conf_level the confidence level of the test
#'
#' @return returns instance of class "LR_tests" i.e. a list for each of the three
#' tests that includes the critical value, the test statistic, the p-value and
#' the decision i.e. reject or not
VaR_LR_tests <- function(p, VaR, plrets = portfolio_plret_df[-c(1:1000),2],
                         conf_level = 0.95){
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
             cc  = list(crit_val_cc = crit_val_cc, LR_cc = LR_cc, 
                        p_val_cc = p_val_cc, reject_cc = reject_cc),
             ind = list(crit_val_ind = crit_val_ind, LR_ind = LR_ind,
                        p_val_ind = p_val_ind, reject_ind = reject_ind),
             uc  = list(crit_val_uc = crit_val_uc, LR_uc = LR_uc,
                        p_val_uc = p_val_uc, reject_uc = reject_uc)))
}

## Don't run when importing
if (sys.nframe() == 0) {
  VaR_LR_tests(0.01, Uni_Normal_GARCH_VaR[, 2])
}


#----------------------------------------------------------------#
########### Calculate Exceedances and Nominal Coverage ###########
#----------------------------------------------------------------#

#' Empirical Coverage
#' 
#' Calculates and returns the nominal coverage
#'
#' @param VaR Value at risk forecasts of a model
#' @param plrets portfolio returns dataframe with dates in first column& returns
#' in second column
#'
#' @return nominal coverage
empirical_coverage <- function(VaR, plrets = portfolio_plret_df[-c(1:1000),2]){
  indicator <- ifelse(plrets-VaR<0, 1, 0)
  coverage <- sum(indicator)/length(VaR)
  return(empirical_coverage = coverage)
}

## Don't run when importing
if (sys.nframe() == 0) {
  empirical_coverage(Uni_Normal_GARCH_VaR[,2])
  empirical_coverage(Uni_t_GJR_GARCH_VaR[,2])
}


if (!require(tidyverse)) install.packages("tidyverse")
# lubridate to extract year from Date
if (!require(lubridate)) install.packages("lubridate")

#' Exceedances
#' 
#' Calculates and returns the total number of exceedances as well as the
#'  exceedances per year
#'
#' @param VaR Value at risk forecasts of a model
#' @param plrets portfolio returns dataframe with dates in first column& returns
#'in second column
#'
#' @return list with total number of exceedances and exceedences per year
exceedances <- function(VaR, plrets = portfolio_plret_df[-c(1:1000),]){
  indicator <- ifelse(plrets[,2]-VaR<0, 1, 0)
  indicator_df <- data.frame(Date = plrets[,1],
                             Exceedance = as.factor(indicator)
                             )
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

## Don't run when importing
if (sys.nframe() == 0) {
  exceedances(Uni_Normal_GARCH_VaR[,2])
  exceedances(Uni_t_GJR_GARCH_VaR[,2])
}





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
exceedances_table <- function(VaR_list, 
                              plrets = portfolio_plret_df[-c(1:1000),]){
  n <- length(VaR_list)
  matrix_01 <- matrix(0L, nrow = n, ncol = 9)
  for (i in 1:n){
    matrix_01[i, 1] <- unlist(exceedances(VaR_list[[i]][,2])$total_exc)
    matrix_01[i, 2:9] <- unlist(exceedances(VaR_list[[i]][,2])$exc_per_year[,2])
  }
  table_01 <- data.frame(matrix_01)
  colnames(table_01) <- c("Total_Exc", "2004", "2005", "2006", "2007",
                          "2008", "2009", "2010", "2011")
  rownames(table_01) <- names(VaR_list)
  
  matrix_05 <- matrix(0L, nrow = n, ncol = 9)
  for (i in 1:n){
    matrix_05[i, 1] <- unlist(exceedances(VaR_list[[i]][,3])$total_exc)
    matrix_05[i, 2:9] <- unlist(exceedances(VaR_list[[i]][,3])$exc_per_year[,2])
  }
  table_05 <- data.frame(matrix_05)
  colnames(table_05) <- c("Total_Exc", "2004", "2005", "2006", "2007",
                          "2008", "2009", "2010", "2011")
  rownames(table_05) <- names(VaR_list)
  
  return(list(table_01 = table_01, table_05 = table_05))
}

## Don't run when importing
if (sys.nframe() == 0) {
  test_VaR_list <- list(EWMA = Uni_EWMA_VaR, Normal_GARCH = Uni_Normal_GARCH_VaR,
                        t_GJR = Uni_t_GJR_GARCH_VaR, 
                        Skewt_GJR = Uni_Skewt_GJR_GARCH_VaR,
                        skewt_NGARCH = Uni_Skewt_NGARCH_VaR,
                        normal_DCC_GARCH = Multi_DCC_GARCH_VaR)
  exceedances_table(test_VaR_list)$table_01
  exceedances_table(test_VaR_list)$table_05
}


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
  coverage_01 <- matrix(0L, nrow = n, ncol = 1)
  tests_01 <- matrix(0L, nrow = n, ncol = 3)
  for (i in 1:n){
    coverage_01[i, 1] <- empirical_coverage(VaR_list[[i]][,2], plrets[,2])
    tests_01[i, 1] <- unlist(VaR_LR_tests(0.01, VaR_list[[i]][,2], 
                                          plrets[,2])@uc$p_val_uc)
    tests_01[i, 2] <- unlist(VaR_LR_tests(0.01, VaR_list[[i]][,2], 
                                          plrets[,2])@ind$p_val_ind)
    tests_01[i, 3] <- unlist(VaR_LR_tests(0.01, VaR_list[[i]][,2], 
                                          plrets[,2])@cc$p_val_cc)
  }
  exceed_01 <- exceedances_table(VaR_list, plrets)$table_01
  performance_table_01 <- data.frame(coverage_01, tests_01, exceed_01)
  colnames(performance_table_01) <- c("coverage_1%", "uc", "ind", "cc", 
                                      "Total_Exc", "2004", "2005", "2006", 
                                      "2007", "2008", "2009", "2010", "2011")
  rownames(performance_table_01) <- names(VaR_list)
  
  coverage_05 <- matrix(0L, nrow = n, ncol = 1)
  tests_05 <- matrix(0L, nrow = n, ncol = 3)
  for (i in 1:n){
    coverage_05[i, 1] <- empirical_coverage(VaR_list[[i]][,3], plrets[,2])
    tests_05[i, 1] <- unlist(VaR_LR_tests(0.05, VaR_list[[i]][,3], 
                                          plrets[,2])@uc$p_val_uc)
    tests_05[i, 2] <- unlist(VaR_LR_tests(0.05, VaR_list[[i]][,3],
                                          plrets[,2])@ind$p_val_ind)
    tests_05[i, 3] <- unlist(VaR_LR_tests(0.05, VaR_list[[i]][,3], 
                                          plrets[,2])@cc$p_val_cc)
  }
  exceed_05 <- exceedances_table(VaR_list, plrets)$table_05
  performance_table_05 <- data.frame(coverage_05, tests_05, exceed_05)
  colnames(performance_table_05) <- c("coverage_5%", "uc", "ind", "cc", 
                                      "Total_Exc", "2004", "2005", "2006", 
                                      "2007", "2008", "2009", "2010", "2011")
  rownames(performance_table_05) <- names(VaR_list)
  
  return(list(performance_table_01 = performance_table_01 %>% round(3), 
              performance_table_05 = performance_table_05 %>% round(3)))
}

## Don't run when importing
if (sys.nframe() == 0) {
  performance_table(test_VaR_list)$performance_table_01 
  performance_table(test_VaR_list)$performance_table_05
}




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
#' @param quantiles VaR quantiles; by default c(0.01, 0.05)
#'
#' @return list w/ losses for the percentiles as first two elements and 
#' mean losses for the percentiles as last two elements
loss_VaR <- function(VaR, plrets = portfolio_plret_df[-c(1:1000),2], 
                     quantiles = c(0.01, 0.05)){
  indicator_01 <- ifelse(plrets-VaR[,2]<0, 1, 0)
  indicator_05 <- ifelse(plrets-VaR[,3]<0, 1, 0)
  loss_01 <- (plrets-VaR[,2])*(quantiles[1]-indicator_01)
  loss_05 <- (plrets-VaR[,3])*(quantiles[2]-indicator_05)
  return(list(
    loss_01=loss_01, 
    loss_05 = loss_05, 
    mean_loss_01 = mean(loss_01), 
    mean_loss_05 = mean(loss_05)
  ))
}

## Don't run when importing
if (sys.nframe() == 0) {
  Uni_t_GJR_loss <- loss_VaR(Uni_t_GJR_GARCH_VaR)
  Uni_Normal_loss <- loss_VaR(Uni_Normal_GARCH_VaR)
  mean_loss_sgarch <- Uni_Normal_loss$mean_loss_01
  mean_loss_gjr <- Uni_t_GJR_loss$mean_loss_01
  
  ## Double checking with MCS:
  # Note: I only saw that MCS package had implemented tick loss function for VaR
  # towards the very end of my thesis. Hence, I have implemented my own function
  # for tick loss long before I discovered the LossVaR function in MCS
  library(MCS)
  mcs_mean_loss_sgarch <- mean(MCS::LossVaR(
    portfolio_plret_df[-c(1:1000),2], Uni_Normal_GARCH_VaR[,2], type = "normal", 
    tau = 0.01))
  mcs_mean_loss_gjr <-  mean(MCS::LossVaR(
    portfolio_plret_df[-c(1:1000),2], Uni_t_GJR_GARCH_VaR[,2], type = "normal", 
    tau = 0.01))
  mcs_mean_loss_sgarch == mean_loss_sgarch
  mcs_mean_loss_gjr == mean_loss_gjr
}



#' Ranking VaR Forecasts
#' 
#' Ranking VaR forecasts in ascending order based on their average VaR tick loss.
#'
#' @param VaR_list list of VaR dataframes with date in first column, 1% VaR in
#' second column and 5% VaR in third column
#' @param plrets portfolio returns; by default portfolio returns from t=1001 
#' until t=T
#' @param quantiles VaR quantiles; by default c(0.01, 0.05)
#'
#' @return list with tables for 1% and 5% VaR ranking
VaR_loss_ranking <- function(VaR_list, plrets = portfolio_plret_df[-c(1:1000),],
                             quantiles = c(0.01, 0.05)){
  n <- length(VaR_list)
  matrix_99 <- matrix(0L, nrow = n, ncol = 1)
  matrix_95 <- matrix(0L, nrow = n, ncol = 1)
  for (i in 1:n){
    matrix_99[i, 1] <- unlist(loss_VaR(VaR_list[[i]])$mean_loss_01)
    matrix_95[i, 1] <- unlist(loss_VaR(VaR_list[[i]])$mean_loss_05)
  }
  table_01 <- data.frame(matrix_99)
  colnames(table_01) <- c("mean_VaR_loss")
  rownames(table_01) <- names(VaR_list)
  table_01 <- table_01 %>% arrange(mean_VaR_loss) # arrange in ascending order
  
  
  table_05 <- data.frame(matrix_95)
  colnames(table_05) <- c("mean_VaR_loss")
  rownames(table_05) <- names(VaR_list)
  table_05 <- table_05 %>% arrange(mean_VaR_loss) 
  return(list(table_01 = table_01, table_05 = table_05))
}

## Don't run when importing
if (sys.nframe() == 0) {
  VaR_loss_ranking(test_VaR_list)$table_01
  VaR_loss_ranking(test_VaR_list)$table_05
}



## Create class to return two lists in CPA_test function
setClass(Class="CPA",
         representation(
           VaR_01="list",
           VaR_05="list"
         )
)

#' Conditional Predictive Ability Test by Giacomini and White (2006)
#' 
#' Implements CPA test to allow for binary model comparisons in predictive 
#' ability with a confidence level of 95%
#'
#' @param VaR1 VaR forecasts of model i (whole dataframe)
#' @param VaR2 VaR forecasts of model j (whole dataframe)
#' @param plrets portfolio returns; by default portfolio returns from t=1001 
#' until t=T
#'
#' @return instance of class "CPA" i.e. for the 1% and the 5% VaR it returns
#' a list of the test statistic, pvalue, critical value and the test decision
#' (i.e. which model has higher predictive ability and whether this difference
#' is significant or not)
CPA_test <- function(VaR1, VaR2, plrets = portfolio_plret_df[-c(1:1000),2]){
  loss1 <- loss_VaR(VaR1)
  loss2 <- loss_VaR(VaR2)
  
   
  ## 99%
  d_ij_01 <- as.matrix(loss1$loss_01)-as.matrix(loss2$loss_01)#loss differential
  T <- length(d_ij_01)
  h_tminus1_01 <- cbind(matrix(1, ncol = 1, nrow = T-1), 
                        matrix(d_ij_01[1:T-1], nrow = T-1, ncol = 1))
  
  
  
  lossdiff <- d_ij_01[2:T] # loss differential from 2nd observation onwards
  
  reg_01 <- matrix(1, nrow = nrow(h_tminus1_01), ncol = ncol(h_tminus1_01))
  for (jj in 1:2) reg_01[, jj] <- h_tminus1_01[, jj]*lossdiff
  
  # since forecasting horizon is 1, test stat can be calculated as n*R^2
  # from the regression of one on ht_minus1_01*lossdiff
  fit_01 <- lm(matrix(1, nrow = T-1, ncol =1)~0+reg_01)
  r2_01 <- summary(fit_01)$r.squared
  n <- T-1-1+1 # n=T-tau-m1+1 in paper
  test_stat_01 <- n*r2_01
  
  
  ## 95%
  d_ij_05 <- as.matrix(loss1$loss_05)-as.matrix(loss2$loss_05)#loss differential
  T <- length(d_ij_05)
  h_tminus1_05 <- cbind(matrix(1, ncol = 1, nrow = T-1), 
                        matrix(d_ij_05[1:T-1], nrow = T-1, ncol = 1))
  
  
  
  lossdiff <- d_ij_05[2:T] # loss differential from 2nd observation onwards
  
  reg_05 <- matrix(1, nrow = nrow(h_tminus1_05), ncol = ncol(h_tminus1_05))
  for (jj in 1:2) reg_05[, jj] <- h_tminus1_05[, jj]*lossdiff
  
  # since forecasting horizon is 1, test stat can be calculated as n*R^2
  # from the regression of one on ht_minus1_05*lossdiff
  fit_05 <- lm(matrix(1, nrow = T-1, ncol =1)~0+reg_05)
  r2_05 <- summary(fit_05)$r.squared
  n <- T-1-1+1 # n=T-tau-m1+1 in paper
  test_stat_05 <- n*r2_05
  
  ## Critical Values
  crit_val_01 <- qchisq(0.95, 2)
  crit_val_05 <- qchisq(0.95, 2)
  
  ## Calculate p-values
  p_val_01 <- 1-pchisq(abs(test_stat_01),2)
  p_val_05 <- 1-pchisq(abs(test_stat_05),2)
  
  mean_diff_loss_01 <- loss1$mean_loss_01-loss2$mean_loss_01
  signif_01 <- p_val_01<0.05
  ifelse(mean_diff_loss_01<0, c(better_01 <- "row"), c(better_01 <- "col"))
  
  mean_diff_loss_05 <- loss1$mean_loss_05-loss2$mean_loss_05
  signif_05 <- p_val_05<0.05
  ifelse(mean_diff_loss_05<0, c(better_05 <- "row"), c(better_05 <- "col"))
  

  
  if (signif_01 == TRUE){
    message(better_01, " performs significantly better")
  }
  
  if (signif_05 == TRUE){
    message(better_05, " performs significantly better")
  }
  
  return(new("CPA", 
             VaR_01 = list(test_stat_01 = test_stat_01, 
                           crit_val_01 = crit_val_01, p_val_01 = p_val_01, 
                           signif_01 = signif_01, better_01 = better_01),
             VaR_05 = list(test_stat_05 = test_stat_05, 
                           crit_val_05 = crit_val_05, p_val_05 = p_val_05, 
                           signif_05 = signif_05, better_05 = better_05)
  ))
}

## Don't run when importing
if (sys.nframe() == 0) {
  CPA_test(Uni_t_GJR_GARCH_VaR, Uni_Skewt_NGARCH_VaR)
  CPA_test(Uni_t_GJR_GARCH_VaR, Uni_EWMA_VaR)
}


#' CPA Test table
#'
#' @param VaR_list list of VaR dataframes with date in first column, 1% VaR in
#' second column and 5% VaR in third column
#' @param plrets portfolio returns; by default portfolio returns from t=1001 
#' until t=T
#' 
#' @return list of tables for the two VaR levels. Each entry in the table includes
#' the p-value and which model performed better (if one did)
CPA_table <- function(VaR_list, plrets = portfolio_plret_df[-c(1:1000),2]){
  CPA_matrix_01 <- matrix(nrow = length(VaR_list)-1, ncol = length(VaR_list)-1)
  CPA_matrix_05 <- matrix(nrow = length(VaR_list)-1, ncol = length(VaR_list)-1)
  rows <- VaR_list[-length(VaR_list)]
  cols <- VaR_list[-1]
  
  
  for (i in seq_along(rows)){
    for (j in seq_along(cols)){
      if (i<=j){
        p_val_01 <- CPA_test(rows[[i]], cols[[j]], 
                             plrets = plrets)@VaR_01$p_val_01
        better_01 <- CPA_test(rows[[i]], cols[[j]], 
                              plrets = plrets)@VaR_01$better_01
        
        p_val_05 <- CPA_test(rows[[i]], cols[[j]], 
                             plrets = plrets)@VaR_05$p_val_05
        better_05 <- CPA_test(rows[[i]], cols[[j]], 
                              plrets = plrets)@VaR_05$better_05
        
        CPA_matrix_01[i,j] <- paste(as.character(round(as.numeric(p_val_01),3)),
                                    better_01, sep  ="; ")
        CPA_matrix_05[i,j] <- paste(as.character(round(as.numeric(p_val_05),3)),
                                    better_05, sep  ="; ")
        }
    }
  }
  CPA_table_01 <- data.frame(CPA_matrix_01)
  colnames(CPA_table_01) <- names(cols)
  rownames(CPA_table_01) <- names(rows)
  
  CPA_table_05 <- data.frame(CPA_matrix_05)
  colnames(CPA_table_05) <- names(cols)
  rownames(CPA_table_05) <- names(rows)
  
  return(list(CPA_table_01 = CPA_table_01,
              CPA_table_05 = CPA_table_05))
}

## Don't run when importing
if (sys.nframe() == 0) {
  CPA_table(test_VaR_list)
}


