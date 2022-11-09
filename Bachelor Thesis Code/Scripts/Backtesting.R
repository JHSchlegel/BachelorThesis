#=========================================================================#
############################### Backtesting ###############################
#=========================================================================#

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
#' @param percentile VaR percentiles; by default c(99, 95)
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


## Compare w/ VaRloss function from rugarch package:
if (!require(rugarch)) install.packages("rugarch")
all.equal(VaRloss(0.95, portfolio_plret_df[-c(1:1000),2], Uni_t_GJR_GARCH[,3]), 
          as.numeric(as.matrix(100*Uni_t_GJR_loss$loss_95)))


#' Conditional Predictive Ability Test by Giacomini and White (2006)
#'
#' @param VaR1 
#' @param VaR2 
#' @param plrets 
#' @param percentile 
#'
#' @return
#' @export
#'
#' @examples
CPA_test <- function(VaR1, VaR2, plrets = portfolio_plret_df[-c(1:1000),2], 
                                          percentile = c(0.99, 0.95)){
  T <- length(plrets)
  loss1 <- loss_VaR(VaR1)
  loss2 <- loss_VaR(VaR2)
  d_ij_99 <- as.matrix(loss1$loss_99)-as.matrix(loss2$loss_99)
  d_ij_95 <- as.matrix(loss1$loss_95)-as.matrix(loss2$loss_95)
  h_tminus1_99 <- c(rep(1, length(d_ij_99)-1), d_ij_99[-length(d_ij_99)])
  h_tminus1_95 <- c(rep(1, length(d_ij_95)-1), d_ij_99[-length(d_ij_95)])
  Z_99 <- h_tminus1_99*d_ij_99[-1]
  Z_95 <- h_tminus1_95*d_ij_95[-1]
  Z_bar_99 <- mean(Z_99)
  Z_bar_95 <- mean(Z_95)
  Omega_inv_99 <- mean(Z_99%*%t(Z_99))
  Omega_inv_95 <- mean(Z_95%*%t(Z_95))
  GW_ij_99 <- T*Z_bar_99*Omega_inv_99*Z_bar_99
  GW_ij_95 <- T*Z_bar_95*Omega_inv_95*Z_bar_95
  c(pchisq(GW_ij_99, 2), pchisq(GW_ij_95, 2))
}

CPA_test(Uni_t_GJR_GARCH_VaR, Uni_EWMA_VaR)

plrets <- portfolio_plret_df[-c(1:1000),2]
VaR1 <- Uni_t_GJR_GARCH
percentile <- c(0.99, 0.95)


loss1 <- loss_VaR(Uni_t_GJR_GARCH_VaR)
loss2 <- loss_VaR(Uni_EWMA_VaR)
d_ij_99 <- as.matrix(loss1$loss_99)-as.matrix(loss2$loss_99)
d_ij_95 <- as.matrix(loss1$loss_95)-as.matrix(loss2$loss_95)
h_tminus1_99 <- matrix(cbind(rep(1, length(d_ij_99)-1), d_ij_99[-length(d_ij_99)]), ncol = 2)
h_tminus1_95 <- matrix(cbind(rep(1, length(d_ij_95)-1), d_ij_99[-length(d_ij_95)]), ncol = 2)
Z_99 <- h_tminus1_99*d_ij_99[-1]
Z_95 <- h_tminus1_95*d_ij_95[-1]


Z_bar_99 <- 1/T * matrix(apply(Z_99, 2, sum), nrow = 2, ncol = 1)
Z_bar_95 <- 1/T * matrix(apply(Z_95, 2, sum), nrow = 2, ncol = 1)

Omega_inv_99 <- 1/T * matrix(colSums(matrix(apply(Z_99, 1, function(x)x%*%t(x)), ncol = 4)), ncol = 2, nrow = 2)
Omega_inv_95 <- 1/T * matrix(colSums(matrix(apply(Z_95, 1, function(x)x%*%t(x)), ncol = 4)), ncol = 2, nrow = 2)


GW_ij_99 <- T * t(Z_bar_99) %*% Omega_inv_99 %*% Z_bar_99
GW_ij_95 <- T * t(Z_bar_95) %*% Omega_inv_95 %*% Z_bar_95

GW_ij_99
GW_ij_95

pchisq(GW_ij_99, 2)
pchisq(GW_ij_95, 2)
