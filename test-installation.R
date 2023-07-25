eval(parse("plot-tools.R", encoding="UTF-8"))
eval(parse("auto-fit.R", encoding="UTF-8"))
eval(parse("auto-select.R", encoding="UTF-8"))
eval(parse("forecasting.R", encoding="UTF-8"))
eval(parse('arima-simulation.R', encoding='UTF-8'))

library(fpp2)
library(tseries)
library(TSA)
library(seastests)
library(forecast)

library(plotly)
library(forecast)

library(prettydoc)
library(stringi)
library(stringr)
library(polynom)
library(parallel)


load("data/patatas.dat")
sales <- patatas[,1]   
xregs <- ts(matrix(patatas[,2]))
colnames(xregs) <- c('price')
chips_model <- drm.select(sales, xregs, show_info=F)

# point predictions
preds <- forecast_model(sales, xregs, chips_model, h=10, mode='bootstrap')
fig <- plot_forecast(preds)
save_image(fig, file='test-figure.pdf')
