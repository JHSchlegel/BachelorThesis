library(rugarch)
library(xdcclarge)
library(rmgarch)
library(xts)
#load data
data(us_stocks)
n<-3
Rtn<-log(us_stocks[-1,1:n]/us_stocks[-nrow(us_stocks),1:n])
# Step 1:GARCH Parameter Estimation with rugarch
spec = ugarchspec()
mspec = multispec( replicate(spec, n = n) )
fitlist = multifit(multispec = mspec, data = Rtn)
ht<-sigma(fitlist)^2
residuals<-residuals(fitlist)
# Step 2:cDCC-GARCH Parameter Estimation with xdcclarge
cDCC<-cdcc_estimation(ini.para=c(0.05,0.93) ,ht ,residuals)
#Time varying correlation matrix Rt at time t
(Rt<-matrix(cDCC$cdcc_Rt,n,n))

parameters <- cDCC$result$par

head(Rtn)
date <- as.Date(rownames(Rtn), format  ="%Y/%m/%d")
any(is.na(date))
Rtn <- xts(Rtn, order.by = date)
u.fit <- ugarchfit(spec, data = Rtn[,1])
u.fcst <- ugarchforecast(u.fit, n.ahead = 1, n.roll = F, out.sample = 0)
u.fcst@forecast$sigmaFor

dcc.spec <- dccspec(mspec)
dcc.fit <- dccfit(dcc.spec, data = Rtn)
dcc.fit@mfit$H

dcc.fcst <- dccforecast(dcc.fit, n.ahead = 1)
dcc.fcst@mforecast$H

dcc.fit@mfit$coef

