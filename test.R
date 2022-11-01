source('auto-fit.R')
source('forecasting.R')
source('auto-select.R')

# Time series analysis
library(fpp2)
library(tseries)
library(TSA)
library(seastests)
library(forecast)

# Plot libraries
library(plotly)

# Auxiliars
library(prettydoc)
library(stringi)
library(stringr)
library(polynom)
library(parallel)
library(reticulate)
reticulate::py_run_string('import kaleido, sys')

cataluna <- read.csv("data/evolucion_gripe_covid.csv")
sarscov2 <- ts(cataluna$sarscov2)
xregs <- ts(cataluna[, c('vac1218', 'vac1845', 'vac4565', 'vac6580', 'vac80')])
adjust <- drm.select(sarscov2, xregs, 
                     ic='aicc', alpha=0.05, st_method='auto.arima', show_info=T)


preds <- forecast_model(sarscov2, xregs, adjust, h=10, 
                        mode='bootstrap', levels=c(50, 75, 90))
fig <- plot_forecast(preds)

# change x axis
time_vec <- window(cataluna$fecha, start=start(preds$x))
fig <- fig %>% layout(
    xaxis = list(
        ticktext = c('2021', '2022'),
        tickvals = time(time_vec)[1] + 
            c(
                which(str_detect(time_vec, '2021'))[1],
                which(str_detect(time_vec, '2022'))[1]
            ),
        tickmode = 'array'
    )
)


save_image(fig, 'test-figure.pdf', width=800, height=400)