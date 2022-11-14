#=========================================================================#
############################### Backtesting ###############################
#=========================================================================#

# TODO: ranking of mean loss

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


## Compare w/ VaRloss function from rugarch package for e.g. GJR GARCH:
if (!require(rugarch)) install.packages("rugarch")
all.equal(VaRloss(0.95, portfolio_plret_df[-c(1:1000),2], Uni_t_GJR_GARCH_VaR[,3]), 
          as.numeric(as.matrix(100*Uni_t_GJR_loss$loss_95)))


## Create class to return two lists in CPA_test function
setClass(Class="CPA",
         representation(
           VaR_99="list",
           VaR_95="list"
         )
)

#' Conditional Predictive Ability Test by Giacomini and White (2006)
#'
#' @param VaR1 
#' @param VaR2 
#' @param plrets 
#' @param percentile 
#'
#' @return
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
  
  ifelse(loss1$mean_loss_99<loss2$mean_loss_99, 
         c(better_99 <- deparse(substitute(VaR1)), worse_99 <-  deparse(substitute(VaR2))),
         c(worse_99 <- deparse(substitute(VaR1)), better_99 <-  deparse(substitute(VaR2)))
         )
  
  ifelse(loss1$mean_loss_95<loss2$mean_loss_95, 
         c(better_95 <- deparse(substitute(VaR1)), worse_95 <-  deparse(substitute(VaR2))),
         c(worse_95 <- deparse(substitute(VaR1)), better_95 <-  deparse(substitute(VaR2)))
  )
  
  if (p_val_99<=0.05){
    signif_99 <- TRUE
    message(better_99, " significantly outperforms ", worse_99, " for the 99% VaR")
  }
  else{signif_99 <- FALSE}
  
  if (p_val_95<=0.05){
    signif_95 <- TRUE
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
