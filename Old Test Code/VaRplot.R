library(ggplot2)
library(readr)
library(tidyverse)
library(rugarch)

VaR_cop_norm_df <- read.csv("Data/VaR/Multi_cop_norm_VaR.csv", header = TRUE)

cop_df <- cop_df %>% 
  mutate(Date = portfolio_plret_df[,1],
         alpha_0.01 = V1, alpha_0.05 = V2) %>% 
  select(-V1, -V2)

dataframe <- data.frame(Date = portfolio_plret_df$Date[-c(1:1000)], alpha_0.01 = VaR_cop_norm_df[,1], alpha_0.05 = VaR_cop_norm_df[,2])
head(dataframe)

VaR_exceed_plot <- function(dataframe, VaR_in_col_nr, plrets){
VaR_df <- data.frame(Date = as.Date(dataframe[,1]), VaR = dataframe[,VaR_in_col_nr], 
                     Exceedance = as.factor(plrets[-c(1:1000),2]<dataframe[,VaR_in_col_nr]))
ggplot(VaR_df, aes(x = Date, y = VaR))+
# TODO number of exceedances per year in plot
  geom_point(aes(x = Date, y = plrets[-c(1:1000), 2], color = Exceedance, shape = Exceedance), size =1.5, alpha = 2)+
  scale_shape_manual(values = c(20, 4), name="", labels = c("Lower than VaR", "Greater than VaR"))+
  scale_color_manual(values = c("gray", "red"), name = "", labels = c("Lower than VaR", "Greater than VaR"))+
  geom_line(alpha = 0.7)+
  scale_linetype_manual(values=c("solid", "solid"), name = "", labels = c("VaR", "VaR"))+
  labs(y = "Daily Portfolio Returns", x = "Date")+
  theme_light()+
  theme(legend.position = c(.15, .8), legend.background = element_rect(color = NA), legend.key = element_rect(color = "transparent"))
}
which(ifelse(plrets[-c(1:1000),2]<VaR_df$VaR, 3, 20)==20)
VaR_df$VaR  
VaRplot(0.01, portfolio_plret_ts[-c(1:1000)], VaR_df$VaR)
