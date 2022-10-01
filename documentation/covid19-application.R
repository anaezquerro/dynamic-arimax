


data_new <- construct.data(adjust$history, sarscov2, xregs, new_xreg_name = NULL, 
                           optimal_lag = NULL, max_lag=6)
head(data_new)

serie <- data_new[, c(1)]
xregs <- data_new[, -c(1)]

ajuste <- auto.arima(serie, xreg=xregs, d=0, D=0)

ajuste2 <- fit.coefficients(ajuste)

abs(ajuste2$coef) < 1.96*sqrt(diag(ajuste2$var.coef))
