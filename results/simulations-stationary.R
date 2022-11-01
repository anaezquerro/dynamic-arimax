eval(parse("arima-simulation.R", encoding="UTF-8"))
eval(parse("auto-fit.R", encoding="UTF-8"))
eval(parse("auto-select.R", encoding="UTF-8"))

# Librerías de series temporales
library(fpp2)
library(tseries)
library(TSA)
library(seastests)
library(forecast)

# Librerías para los gráficos
library(plotly)
library(forecast)

# Auxiliares
library(prettydoc)
library(stringi)
library(stringr)
library(polynom)
library(parallel)

ics <- c('aicc')
n <- 5000
m <- 20
stationary_methods <- c('auto.arima', 'adf.test')


generate_orders <- function() {
    p <- sample(0:3, size=1, prob=c(0.3, 0.3, 0.2, 0.1))
    q <- sample(0:3, size=1, prob=c(0.3, 0.3, 0.2, 0.1))
    d <- sample(0:2, size=1, prob=c(0.5, 0.4, 0.1))
    return(list(p=p, q=q, d=d))
}

generate_coefficients <- function() {
    sign <- ifelse(runif(4)>0.3, 1, -1)
    betas <- sign*rnorm(4, 2, 0.5) 
    betas[1] <- rnorm(1, 0, 0.1)*ifelse(runif(1)>0.4, 1, -1)
    return(betas)
}

generate_lags <- function() {
    lags <- sample(0:6, size=3, prob=c(0.3, 0.25, 0.25, 0.1, 0.05, 0.025, 0.025))
    return(lags)
}




TP <- 0
TN <- 0
FP <- 0
FN <- 0
lag_right <- 0


ics <- c('aic', 'bic', 'aicc')
stationary_methods <- c('auto.arima', 'adf.test')

results <- list()
for (stmethod in stationary_methods) {
    results[[stmethod]] <- list()
    for (ic in ics) {
        results[[stmethod]][[ic]] <- list(TP=0, TN=0, FP=0, FN=0, lag_right=0)
    }
}
                            
generate_covariate <- function(parinfo) {
    stationary <- parinfo$stationary
    set.seed(parinfo$seed)    
    
    model_orders <- generate_orders()
    if (stationary) {
        model_orders$d <- 0
    }
    xreg <- sim.arima(model_orders, n, F)
    return(xreg)
}

simulate <- function(results) {
    
    
    check_covariate <- function(covariate) {
        if (is.null(ajuste$xreg) || all(is.na(ajuste$history))) {
            return(F)
        }
        
        if (ncol(ajuste$xreg) == 1) {
            return(covariate %in% ajuste$history$var)
        } else {
            return( 
                (covariate %in% ajuste$history$var) &
                    (covariate %in% names(ajuste$coef[ajuste$coef != 0]))
            )
        }
    }
    
    check_lag <- function(covariate, correct_lag) {
        return(as.numeric(
            correct_lag == as.numeric(ajuste$history$lag[ajuste$history$var == covariate])
        )
        )
    }
    
    # generate covariates and residuals
    local_variables <- ls()
    global_variables <- names(.GlobalEnv)
    cl <- makeCluster(detectCores(logical=F))
    clusterExport(cl, local_variables, env=environment())
    clusterExport(cl, global_variables)
    clusterEvalQ(cl, library(polynom))
    clusterEvalQ(cl, library(fpp2))
    clusterEvalQ(cl, library(tseries))
    clusterEvalQ(cl, library(TSA))
    clusterEvalQ(cl, library(seastests))
    clusterEvalQ(cl, library(forecast))
    clusterEvalQ(cl, library(stringr))
    
    covariates <- c('X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'residuals')
    covarinfo <- sapply(covariates, function(x) list(seed=sample(0:10000, 1), stationary=(x=='residuals')), simplify=F)
    
    xregs <- parLapply(cl, covarinfo, generate_covariate)
    stopCluster(cl)
    
    
    # display information about covariates
    for (covar in covariates) {
        cat(paste0(' - ', covar, '~ARIMA(', length(xregs[[covar]]$ar), ',', 
                   xregs[[covar]]$d, ',', length(xregs[[covar]]$ma), ')\n'))
    }
    
    X1 <- xregs[['X1']]$X; X2 <- xregs[['X2']]$X; X3 <- xregs[['X3']]$X
    X4 <- xregs[['X4']]$X; X5 <- xregs[['X5']]$X; X6 <- xregs[['X6']]$X
    residuals <- xregs[['residuals']]$X
    
    
    # generate coefficients and lags
    betas <- generate_coefficients()
    lags <- as.list(setNames(generate_lags(), c('X1', 'X2', 'X3')))
    beta0 <- betas[1]; beta1 <- betas[2]; beta2 <- betas[3]; beta3 <- betas[4]
    lag1 <- lags[['X1']]; lag2 <- lags[['X2']]; lag3 <- lags[['X3']]
    
    # construct the dependent variable
    Y <- beta0 + beta1 * lag(X1, -lag1) + beta2 * lag(X2, -lag2) + beta3 * lag(X3, -lag3) + residuals
    xregs <- cbind(X1, X2, X3, X4, X5, X6)
    
    
    for (ic in ics) {
        for (statmet in stationary_methods) {
            ajuste <- NULL
            while (!is_valid(ajuste)) {
                ajuste <- drm.select(Y, xregs, ic=ic, st_method = statmet, show_info=F)
            }
            cat(paste0('Se ha ajustado un modelo [ic=', ic, ', statmet=', statmet, ']\n'))
            print(ajuste, row.names=F)
            
            
            # update results
            for (covar in c('X1', 'X2', 'X3')) {
                if (check_covariate(covar)) {
                    # update TP values
                    results[[statmet]][[ic]][['TP']] <- results[[statmet]][[ic]][['TP']] + 1
                    
                    # update lag right values
                    results[[statmet]][[ic]][['lag_right']] <- results[[statmet]][[ic]][['lag_right']] + check_lag(covar, lags[[covar]]) 
                    
                } else {
                    results[[statmet]][[ic]][['FN']] <- results[[statmet]][[ic]][['FN']] + 1
                }
            }
            
            for (covar in c('X4', 'X5', 'X6')) {
                if (!check_covariate(covar)) {
                    results[[statmet]][[ic]][['TN']] <- results[[statmet]][[ic]][['TN']] + 1
                } else {
                    results[[statmet]][[ic]][['FP']] <- results[[statmet]][[ic]][['FP']] + 1
                }
            }            
            
        }
    }

    return(results)
    
}



#' Primer caso. Utilizando errores estacionarios
for (i in 1:m) {
    cat(paste0(stri_dup('-', 80), '\n'))
    cat(paste0('Ejecutando la simulación ', i, '/', m, '\n'))
    
    partial_results <- try(simulate(results))
    while (class(partial_results) == 'try-error') {
        
        partial_results <- try(simulate(results), T)
        
    }
    results <- partial_results

}
