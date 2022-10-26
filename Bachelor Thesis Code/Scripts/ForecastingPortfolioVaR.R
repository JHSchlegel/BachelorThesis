#=================================#
#### Forecasting Portfolio VaR ####
#=================================#


if (!require(rugarch)) install.packages("rugarch")
ir (!require(rmgarch)) install.packages("rmgarch")
library(xts)
library(parallel)

#--------------------------------------------------#
########### VaR Exceedance Plot Function ###########
#--------------------------------------------------#

# load lubridate for year() function to extract year from Date
if (!require(lubridate)) install.packages("lubridate") 
#' @param dataframe dataframe with VaR
#' @param VaR_in_col_nr integer indicating in which column of dataframe the VaR is
#' @param pf_plrets dataframe of portfolio percentage returns
VaR_exceed_plot <- function(dataframe, VaR_in_col_nr, pf_plrets){
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
    labs(y = "Daily Portfolio Returns", x = "Date")+
    theme_light()+
    theme(legend.position = c(.15, .8), legend.background = element_rect(color = NA), legend.key = element_rect(color = "transparent"))+
    annotate("text", x = as.Date("2005-01-15"), y = -9, size = 3, hjust = 0,
             label = paste("Number of Exceedances per Year:\n2004:", exceedances_per_year$n[1],
                           "\n2005:", exceedances_per_year$n[2], "\n2006:", exceedances_per_year$n[3],
                           "\n2007:", exceedances_per_year$n[4], "\n2008:", exceedances_per_year$n[5],
                           "\n2009:", exceedances_per_year$n[6], "\n2010:", exceedances_per_year$n[7],
                           "\n2011:", exceedances_per_year$n[8]))
}

# TODO: use AR(1)-GARCH(1,1) for univiariate models

# Import Data and Create xts Objects --------------------------------------


FFCFactors.df <- read.csv("./Data/FFCFactors.csv", header = TRUE)
stocks.plret.df <- read.csv("./Data/StockPlrets.csv", header = TRUE)
portfolio.plret.df <- data.frame(Date = stocks.plret.df$Date, Portfolio = apply(stocks.plret.df[,-1], 1, mean))

## Convert to time series
FFCFactors.ts <- xts(FFCFactors.df[,-1], order.by = as.Date(FFCFactors.df$Date))
stocks.plret.ts <- xts(stocks.plret.df[,-1], order.by = as.Date(stocks.plret.df$Date))
portfolio.plret.ts <- xts(portfolio.plret.df[,-1], order.by = as.Date(portfolio.plret.df$Date))

## Detect how many cores are available for parallelizing
numcores <- detectCores()


# Univariate Normal GARCH(1,1) --------------------------------------------

uni_normal.meanModel <- list(armaOrder=c(0,0), include.mean=TRUE)
uni_normal.varModel <- list(model = "sGARCH", garchOrder = c(1,1))
uni_normal.spec <- ugarchspec(variance.model = uni_normal.varModel, mean.model = uni_normal.meanModel,
                              distribution.model = "norm")
cl <- makePSOCKcluster(numcores)
uni_normal_GARCH.roll <- ugarchroll(uni_normal.spec, portfolio.plret.ts, n.start = 1000, refit.every = 1,
                       refit.window = "moving", solver = "hybrid",
                       calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05),
                       cluster = cl, keep.coef = FALSE)

stopCluster(cl)

uni_normal_GARCH_VaR <- data.frame(Date = portfolio.plret.df$Date[-c(1:1000)], alpha.1perc = uni_normal_GARCH.roll@forecast$VaR[,1],
                                   alpha.5perc = uni_normal_GARCH.roll@forecast$VaR[,2])
head(uni_normal_GARCH_VaR)
write.csv(uni_normal_GARCH_VaR, "Data\\VaR\\Uni_Normal_GARCH_VaR.csv", row.names = FALSE)



# Univariate EWMA ---------------------------------------------------------


lambda <- 0.94
uni_ewma.meanModel <- list(armaOrder=c(0,0), include.mean=TRUE)
uni_ewma.varModel <- list(model="iGARCH", garchOrder=c(1,1))
uni_ewma.spec <- ugarchspec(variance.model= ewma.varModel, mean.model= ewma.meanModel,  
                        distribution.model="norm", fixed.pars=list(omega=0, alpha = 1-lambda))
cl <- makePSOCKcluster(numcores)
uni_ewma.roll <- ugarchroll(uni_ewma.spec, portfolio.plret.ts, n.start = 1000, refit.every = 1,
                               refit.window = "moving", solver = "hybrid",
                               calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05),
                               cluster = cl, keep.coef = FALSE)
stopCluster(cl)
uni_ewma_VaR <- data.frame(Date = portfolio.plret.df$Date[-c(1:1000)], alpha.1perc = uni_ewma.roll@forecast$VaR[,1],
                             alpha.5perc = uni_ewma.roll@forecast$VaR[,2])
head(uni_ewma_VaR)
write.csv(uni_ewma_VaR, "Data\\VaR\\Uni_EWMA_VaR.csv", row.names = FALSE)


# Univariate t-GJR-GARCH(1,1) ---------------------------------------------

uni_t_gjr.meanModel <- list(armaOrder = c(0,0), include.mean = TRUE)
uni_t_gjr.varModel <- list(model = "gjrGARCH", garchOrder = c(1,1))
uni_t_gjr.spec <- ugarchspec(variance.model = uni_t_gjr.varModel, mean.model = uni_t_gjr.meanModel, distribution.model = "std")
cl <- makePSOCKcluster(numcores)
uni_t_gjr_GARCH.roll <- ugarchroll(uni_t_gjr.spec, portfolio.plret.ts, n.start = 1000, refit.every = 1,
                              refit.window = "moving", solver = "hybrid",
                              calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05),
                              cluster = cl, keep.coef = FALSE)
stopCluster(cl)
uni_t_gjr_GARCH_VaR <- data.frame(Date = portfolio.plret.df$Date[-c(1:1000)], alpha.1perc = uni_t_gjr_GARCH.roll@forecast$VaR[,1],
                             alpha.5perc = uni_t_gjr_GARCH.roll@forecast$VaR[,2])
head(uni_t_gjr_GARCH_VaR)
write.csv(uni_t_gjr_GARCH_VaR, "Data\\VaR\\Uni_t_GJR_GARCH.csv", row.names = FALSE)


# Univariate skew-t-NGARCH(1,1) -------------------------------------------


uni_skewt_NGARCH.meanModel <- list(armaOrder = c(0, 0), include.mean = TRUE)
uni_skewt_NGARCH.varModel <- list(model = "fGARCH", submodel = "NGARCH", garchOrder = c(1,1))
uni_skewt_NGARCH.spec <- ugarchspec(variance.model = uni_skewt_NGARCH.varModel, mean.model = uni_skewt_NGARCH.meanModel, 
                                    distribution.model = "sstd")
cl <- makePSOCKcluster(numcores)
uni_skewt_NGARCH.roll <- ugarchroll(uni_skewt_NGARCH.spec, portfolio.plret.ts, n.start = 1000, refit.every = 1,
                              refit.window = "moving", solver = "hybrid",
                              calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05),
                              cluster = cl, keep.coef = FALSE)
stopCluster(cl)
uni_skewt_NGARCH_VaR <- data.frame(Date = portfolio.plret.df$Date[-c(1:1000)], alpha.1perc = uni_skewt_NGARCH.roll@forecast$VaR[,1],
                                  alpha.5perc = uni_skewt_NGARCH.roll@forecast$VaR[,2])
head(uni_skewt_NGARCH_VaR)
write.csv(uni_skewt_NGARCH_VaR, "Data\\VaR\\Uni_Skewt_NGARCH.csv", row.names = FALSE)


###################### Don't forget univariate mix normal garch
######################function for dataframe and head()


# Multivariate Normal DCC-GARCH(1,1) --------------------------------------


portfolio.weights <- rep(1/10, 10)
dcc.meanModel <- list(armaOrder = c(0,0), include.mean = TRUE)
dcc.varModel <- list(model = "sGARCH", garchOrder = c(1,1))
dcc.uspec <- ugarchspec(variance.model = dcc.varModel, mean.model = dcc.meanModel, distribution.model = "norm")
dcc.mspec <- multispec(replicate(10, dcc.uspec))
dcc.spec <- dccspec(dcc.mspec, VAR = FALSE, model = "DCC", dccOrder = c(1,1), distribution =  "mvnorm")

cl <- makeCluster(cl)
multi_dcc.roll <- dccroll(spec = dcc.spec, data = stocks.plret.ts, n.ahead = 1, n.start = 1000, refit.every = 1, refit.window = "moving", window.size = 1000,  cluster = cl)
stopCluster(cl)

n <- length(stocks.plret.df[,1])
n.window <- n-1000
dcc_sigma_pf <- numeric(n.window)
dcc.cov <- rcov(multi_dcc.roll) # rcov to extract forecasted Covariance matrix; rcor would extract forecasted Correlation matrix
for (i in 1:n.window) dcc_sigma_pf[i] <- portfolio.weights %*% dcc.cov[,,i] %*% portfolio.weights
multi_dcc_VaR <- data.frame(Date = portfolio.plret.df$Date[-c(1:1000)])
multi_dcc_VaR$alpha.1perc <- sapply(dcc_sigma_pf, function(x) x*qnorm(0.01, 0, 1, lower.tail = TRUE))
multi_dcc_VaR$alpha.5perc <-  sapply(dcc_sigma_pf, function(x) x*qnorm(0.05, 0, 1, lower.tail = TRUE))
head(multi_dcc_VaR)
write.csv(multi_dcc_VaR, "Data\\VaR\\Multi_DCC_GARCH.csv")
