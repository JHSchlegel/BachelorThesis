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


#' @param dataframe dataframe with VaR
#' @param VaR_in_col_nr integer indicating in which column of dataframe the VaR is
#' @param pf_plrets dataframe of portfolio percentage returns
VaR_exceed_plot <- function(dataframe, VaR_in_col_nr, pf_plrets){
  VaR_df <- data.frame(Date = as.Date(dataframe[,1]), VaR = dataframe[,VaR_in_col_nr], 
                       Exceedance = as.factor(pf_plrets[-c(1:1000),2]<dataframe[,VaR_in_col_nr]))
  library(lubridate) # for year() function to extract year from Date
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
    scale_linetype_manual(values=c("solid", "solid"), name = "", labels = c("VaR", "VaR"))+
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


which(ifelse(pf_plrets[-c(1:1000),2]<VaR_df$VaR, 3, 20)==20)
VaR_df$VaR  
VaRplot(0.01, portfolio_plret_ts[-c(1:1000)], VaR_df$VaR)


  
exceedances_per_year
exceedances_per_year$year
format(Date, "%Y")

