library("rvest")
idx <- "^SSMI"
getConstituents <- function(idx){
  url <- paste0("https://finance.yahoo.com/quote/",idx ,"/components?p=", idx)
  dt <- read_html(url)
  const <<- dt %>% html_nodes("body") %>% html_nodes("div") %>% .[[2]] %>% 
    html_nodes("div") %>% .[[1]] %>% html_nodes("section") %>% html_nodes("table") %>% 
    html_table() %>% as.data.frame()
  View(const)
  const_names <<- const[,1]
}
getConstituents(idx)
const_names
const_names_vec <- c(const_names)
class(const_names_vec)


library(BatchGetSymbols)
library(ggplot2)

first.date <- Sys.Date()-365
last.date <- Sys.Date()

#df.SP500 <- GetSP500Stocks()
#tickers <- df.SP500$Tickers
tickers <- const_names_vec

l.out <- BatchGetSymbols(tickers = tickers,
                         first.date = first.date,
                         last.date = last.date)
l.out$df.tickers

print(l.out$df.control)
print(l.out$df.tickers)

p <- ggplot(l.out$df.tickers, aes(x = ref.date, y = price.close))
p <- p + geom_line()
p <- p + facet_wrap(~ticker, scales = 'free_y') 
print(p)