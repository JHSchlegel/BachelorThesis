#=================================#
#### Forecasting Portfolio VaR ####
#=================================#

if (!require(parallel)) install.packages("parallel")
if (!require(rugarch)) install.packages("rugarch")
if (!require(rmgarch)) install.packages("rmgarch")
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(xts)) install.packages("xts")



# TODO: use AR(1)-GARCH(1,1) for univiariate models

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
########### VaR Exceedance Plot Function ###########
#--------------------------------------------------#

# load lubridate for year() function to extract year from Date
if (!require(lubridate)) install.packages("lubridate") 
#' @param dataframe dataframe with VaR
#' @param VaR_in_col_nr integer indicating in which column of dataframe the VaR is
#' @param pf_plrets dataframe of portfolio percentage returns
#' @param alpha VaR level in percent w/o percentage sign
#' @param modelname Name of model as character
#' 
#' @return Exceedences plot which highlights exceedances and
#'  displays number of exceedances per year
VaR_exceed_plot <- function(dataframe, VaR_in_col_nr, pf_plrets, alpha, modelname){
  VaR_df <- data.frame(Date = as.Date(dataframe[,1]), VaR = dataframe[,VaR_in_col_nr], 
                       Exceedance = as.factor(pf_plrets[-c(1:1000),2]<dataframe[,VaR_in_col_nr]))
  exceedances_per_year  <- VaR_df %>% 
    mutate(year = year(Date)) %>% 
    select(Exceedance, year) %>% 
    count(year, Exceedance) %>% 
    mutate(n = ifelse(Exceedance==TRUE, n, 0)) %>% 
    select(-Exceedance) %>% 
    group_by(year) %>% 
    summarise(n = sum(n))
  
  ggplot(VaR_df, aes(x = Date, y = VaR))+
    geom_point(aes(x = Date, y = pf_plrets[-c(1:1000), 2], color = Exceedance, shape = Exceedance), size =1.5, alpha = 2)+
    scale_shape_manual(values = c(20, 4), name="", labels = c("Lower than VaR", "Greater than VaR"))+
    scale_color_manual(values = c("gray", "red"), name = "", labels = c("Lower than VaR", "Greater than VaR"))+
    geom_line(alpha = 0.7)+
    labs(y = "Daily Portfolio Returns", x = "Date", 
         title = paste0(alpha, "% VaR Exceedances Plot for ",modelname))+
    theme_light()+
    theme(legend.position = c(.15, .8), legend.background = element_rect(color = NA), legend.key = element_rect(color = "transparent"))+
    annotate("text", x = as.Date("2005-01-15"), y = -7, size = 3, hjust = 0,
             label = paste("Number of Exceedances per Year:", "\n2004:", exceedances_per_year$n[1], "\n2005:", exceedances_per_year$n[2],
                           "\n2006:", exceedances_per_year$n[3], "\n2007:", exceedances_per_year$n[4]))+
    annotate("text", x = as.Date("2005-10-15"), y = -7, size = 3, hjust = 0,
             label = paste(" ","\n2008:", exceedances_per_year$n[5],"\n2009:", exceedances_per_year$n[6],
                           "\n2010:", exceedances_per_year$n[7], "\n2011:", exceedances_per_year$n[8]))
}


# Univariate Normal GARCH(1,1) --------------------------------------------

uni_normal_meanModel <- list(armaOrder=c(0,0), include_mean=TRUE)
uni_normal_varModel <- list(model = "sGARCH", garchOrder = c(1,1))
uni_normal_spec <- ugarchspec(variance.model = uni_normal_varModel, mean.model = uni_normal_meanModel,
                              distribution.model = "norm")
cl <- makePSOCKcluster(numcores)
uni_normal_GARCH_roll <- ugarchroll(uni_normal_spec, portfolio_plret_ts, n.start = 1000, refit.every = 1,
                       refit.window = "moving", solver = "hybrid",
                       calculate_var = TRUE, VaR.alpha = c(0.01, 0.05),
                       cluster = cl, keep.coef = FALSE)

stopCluster(cl)

uni_normal_GARCH_VaR <- data.frame(Date = portfolio_plret_df$Date[-c(1:1000)], alpha_1perc = uni_normal_GARCH_roll@forecast$VaR[,1],
                                   alpha_5perc = uni_normal_GARCH_roll@forecast$VaR[,2])
head(uni_normal_GARCH_VaR)
write.csv(uni_normal_GARCH_VaR, "Data\\VaR\\Uni_Normal_GARCH_VaR.csv", row.names = FALSE)

VaR_exceed_plot(uni_normal_GARCH_VaR, 3, portfolio_plret_df, alpha = 95, "Uni_Normal_GARCH")
VaR_exceed_plot(uni_normal_GARCH_VaR, 2, portfolio_plret_df, alpha = 99, "Uni_Normal_GARCH")

# Univariate EWMA ---------------------------------------------------------


lambda <- 0.94
uni_ewma_meanModel <- list(armaOrder=c(0,0), include_mean=TRUE)
uni_ewma_varModel <- list(model="iGARCH", garchOrder=c(1,1))
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


# Univariate t-GJR-GARCH(1,1) ---------------------------------------------

uni_t_gjr_meanModel <- list(armaOrder = c(0,0), include_mean = TRUE)
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


# Univariate skew-t-NGARCH(1,1) -------------------------------------------


uni_skewt_NGARCH_meanModel <- list(armaOrder = c(0, 0), include_mean = TRUE)
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

VaR_exceed_plot(uni_skewt_NGARCH_VaR, 3, portfolio_plret_df, alpha = 95, "Uni_Skewt_NGARCH")
VaR_exceed_plot(uni_skewt_NGARCH_VaR, 2, portfolio_plret_df, alpha = 99, "Uni_Skewt_NGARCH")


###################### Don't forget univariate mix normal garch
######################function for dataframe and head()


# Multivariate Normal DCC-GARCH(1,1) --------------------------------------


portfolio.weights <- rep(1/10, 10)
dcc_meanModel <- list(armaOrder = c(0,0), include_mean = TRUE)
dcc_varModel <- list(model = "sGARCH", garchOrder = c(1,1))
dcc.uspec <- ugarchspec(variance.model = dcc_varModel, mean.model = dcc_meanModel, distribution.model = "norm")
dcc_mspec <- multispec(replicate(10, dcc.uspec))
dcc_spec <- dccspec(dcc_mspec, VAR = FALSE, model = "DCC", dccOrder = c(1,1), distribution =  "mvnorm")

cl <- makeCluster(cl)
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
write.csv(multi_dcc_VaR, "Data\\VaR\\Multi_DCC_GARCH.csv")
