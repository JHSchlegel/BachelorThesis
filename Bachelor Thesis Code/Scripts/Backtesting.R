#=========================================================================#
############################### Backtesting ###############################
#=========================================================================#

## Import my Backtesting Functions
source("Scripts/Backtesting_Functions.R")

#------------------------------------#
########### Importing Data ###########
#------------------------------------#

## Portfolio Plrets
stocks_plret_df <- read.csv("./Data/StockPlrets.csv", header = TRUE)
portfolio_plret_df <- read.csv("./Data/PortfolioPlrets.csv", header = TRUE)

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






all_VaR_list <- list(EWMA = Uni_EWMA_VaR, Normal_GARCH = Uni_Normal_GARCH_VaR,
                      t_GJR = Uni_t_GJR_GARCH_VaR, 
                      Skewt_GJR = Uni_Skewt_GJR_GARCH_VaR,
                      skewt_NGARCH = Uni_Skewt_NGARCH_VaR,
                      normal_DCC_GARCH = Multi_DCC_GARCH_VaR)
#------------------------------------------#
########### Graphical Inspection ###########
#------------------------------------------#

VaR_exceed_plot(Uni_Normal_GARCH_VaR, 3, portfolio_plret_df, alpha = 5,
                "Uni_Normal_GARCH")

#----------------------------------------------------------------#
########### Exceedances, Coverage and LR Tests Pvalues ###########
#----------------------------------------------------------------#
performance_table(all_VaR_list)$performance_table_99
performance_table(all_VaR_list)$performance_table_95


#----------------------------------------------------------------#
########### CPA Tests as in Giacomini and White (2006) ###########
#----------------------------------------------------------------#
CPA_table(passed_LRtests_VaR_list_99)$CPA_table_99
CPA_table(passed_LRtests_VaR_list_95)$CPA_table_95