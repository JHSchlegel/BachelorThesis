library(forecast)
library(zoo)
spi <- read.csv("C:/Users/jansc/Desktop/Universitaet/Bachelor 6. Semester/Bachelor Thesis/Application/Literature Prediction/hspitr (5).csv", sep=";")
spi$DATE <- as.Date(spi$DATE, format="%d.%m.%Y")
spi <- data.frame(lr=spi$Log.Return)
colnames(spi) <- spi$DATE
spi <- ts(spi, start=1987-01-01, freq=frequency(spi))
plot.ts(spi)
spi_training <- spi[-c(1:5)]
spi_training2 <- spi[c(6:26)]
mod <- auto.arima(spi_training)
mod20 <- auto.arima(spi_training2)
prediction_5d <- predict(mod, 5)$pred[1:5]
prediction_20d <- predict(mod20, 5)$pred[1:5]
spi[1:5]-prediction_5d
spi[1:5]-prediction_20d
sum((spi[1:5]-prediction_5d)^2)
mean(abs((spi[1:5]-prediction_5d)/spi[1:5]))
mean(abs((spi[1:5]-prediction_20d)/spi[1:5]))



spi <- read.csv("C:/Users/jansc/Desktop/Universitaet/Bachelor 6. Semester/Bachelor Thesis/Application/Literature Prediction/hspitr (5).csv", sep=";")
spi$DATE <- as.Date(spi$DATE, format="%d.%m.%Y")
spi <- data.frame(l=spi$Return)
colnames(spi) <- spi$DATE
spi <- ts(spi, start=1987-01-01, freq=frequency(spi))
plot.ts(spi)
spi_training <- spi[-c(1:5)]
spi_training2 <- spi[c(6:26)]
mod <- auto.arima(spi_training)
mod20 <- auto.arima(spi_training2)
prediction_5d <- predict(mod, 5)$pred[1:5]
prediction_20d <- predict(mod20, 5)$pred[1:5]
spi[1:5]-prediction_5d
spi[1:5]-prediction_20d
sum((spi[1:5]-prediction_5d)^2)
mean(abs((spi[1:5]-prediction_5d)/spi[1:5]))
mean(abs((spi[1:5]-prediction_20d)/spi[1:5]))


spi <- read.csv("C:/Users/jansc/Desktop/Universitaet/Bachelor 6. Semester/Bachelor Thesis/Application/Literature Prediction/hspitr (5).csv", sep=";")
spi$DATE <- as.Date(spi$DATE, format="%d.%m.%Y")
spi <- data.frame(l=spi$Close)
colnames(spi) <- spi$DATE
spi <- ts(spi, start=1987-01-01, freq=frequency(spi))
plot.ts(spi)
spi_training <- spi[-c(1:5)]
spi_training2 <- spi[c(6:26)]
mod <- auto.arima(spi_training)
mod20 <- auto.arima(spi_training2)
prediction_5d <- predict(mod, 5)$pred[1:5]
prediction_20d <- predict(mod20, 5)$pred[1:5]
spi[1:5]-prediction_5d
spi[1:5]-prediction_20d
sum((spi[1:5]-prediction_5d)^2)
mean(abs((spi[1:5]-prediction_5d)/spi[1:5]))
mean(abs((spi[1:5]-prediction_20d)/spi[1:5]))
