PAD <- 86

library(parallel)
library(stringi)
library(stringr)
library(tseries)
library(fpp2)
library(TSA)
library(seastests)
library(forecast)
library(polynom)


#' README:
#' This script contains the necessary functions to automatically select the 
#' ARIMAX covariates from a set of candidates.



#' **Automatic covariates selection in dynamic regression models (DRM)**
#' 
#' This function implements our new approach of covariates selection in 
#' regression models where each covariate is introduced lagged $k\geq 0$ moments. 
#' For more information about this method, check the paper attached with this 
#' file.
#' 
#' @param serie [ts]: Univariate time series which acts as the dependent 
#' variable.
#' @param xregs [mts]: Set of covariate candidates to model the behavior of 
#' `serie`.
#' @param ic [character]: Information criterion used for covariate selection 
#' and model comparisons. The possibilities are: `aic`, `bic` and `aicc`. By 
#' default, th AICc is used.
#' @param alpha [numeric]: Signficance level for hypothesis tests of stationary 
#' and coefficients significance. By default it is set to $0.05$.
#' @param st_method [character]: Method used for checking stationary of a time 
#' series. If it is `auto.arima`, the function `forecast::auto.arima()` is used 
#' and the differentiation order is checked. If it is `adf.test`, the 
#' Dickey-Fuller test is used.
#' @param show_info [boolean]: Displaying or not displaying the historical 
#' of added covariates.
#' @param ndiff [numeric]: Internal argument (**do not use**) to apply 
#' regular differentiations to all data when no ARIMAX model can be fitted with 
#' stationary errors.
#'
#' @returns Fitted ARIMAX model with significant covariates selected. If the 
#' fit does not exists or it cannot be optimized, it returns a ARIMA model 
#' with the dependent variable or `NA` if it cannot be optimized too.
#' 
drm.select <- function(serie, xregs, ic='aicc', alpha=0.05, 
                       st_method='auto.arima', show_info=T, ndiff=0) {

    
    #' --------------------------- AUXILIAR FUNCTIONS ---------------------------
    
    #' Auxiliar function to print the results obtained in the parallel function 
    #' parLapply()
    #' https://www.rdocumentation.org/packages/parallel/versions/3.6.2/topics/makeCluster)
    
    display_results <- function(fitted_models) {
        for (key in names(fitted_models)) {
            if (typeof(fitted_models[[key]]) == 'character') {
                cat(fitted_models[[key]])
            } else {
                cat(paste0('Covariate ', colnames(xregs)[as.numeric(key)], 
                           ' has been tested [ic=', fitted_models[[key]][[ic]], 
                           ', lag=', fitted_models[[key]]$opt_lag, ']\n'))
            }
        }
    }
    
    #' Auxiliar function to obtain the maximum lag (absolute value) where 
    #' maximum significantive correlation occurs in for all covariate candidates.
    #' This function is meant to cut some observations of the data (serie, xregs).
    #' This is needed since for comparing the models via an IC, we need the data 
    #' to have the same number of observations

    get.max.lag <- function() {
        
        # Note that here we use the parallel function parLapply() to distribute 
        # tasks to different processors and we do not restrict the optimal lag 
        # to be less or equal than zero
        opt_lags <- parLapply(
            cl, 
            as.list(1:ncol(xregs)), 
            function(x) get.opt.lag(xregs[, x], response, alpha, NA, st_method, F)
        )

        # Remove those results that are NA
        opt_lags <- unlist(opt_lags[!is.na(opt_lags)])
        
        # It is possible that no lag could be computed for any candidate
        if (is.null(opt_lags)) {
            return(-Inf)
        }
        
        # Return the lag in absolute value
        return(max(abs(opt_lags)))
    }
    
    
    #' Fit of ARIMAX model by adding one covariate to the current set of 
    #' selected covariates.
    #' 
    #' @param j [numeric]: Column index of the `xregs` matrix of the covariate 
    #' to be added to the current `fitted_model`.
    #' @returns Valid fitted ARIMAX model with the new covariate included or a 
    #' string warning that the covariate could not be included because:
    #' a) No significative lag<=0 could be found
    #' b) No model could be optimized.
    #' 
    #' 
    add.regressor <- function(j) {
        xreg <- xregs[, j]                 # covariate
        xreg_name <- colnames(xregs)[j]    # covariate's name
        
        # Obtain the optimal significative correlation lag between the 
        # covariate and the dependent variable. If it is not found, return 
        # warning
        opt_lag <- get.opt.lag(xreg, response, alpha, max_lag, st_method)
        
        if (is.na(opt_lag)) {
            return(
                paste0(
                    'Significative correlation with lag<=0 could not be found ', 
                    'for ', xreg_name, '\n')
            )
        }
        
        # Creation of a new `mts` object with the dependent variable, selected 
        # covariates of `xregs` and the new covariate j
        data_new <- construct.data(
            history, serie, xregs, xreg_name, opt_lag, max_lag)
        
        # FIt ARIMAX model with the new covariate j
        fitted_model <- auto.fit.arima(
            serie=data_new[, c(1)], xregs=data_new[, -c(1)], 
            ic=ic, alpha=alpha, show_info=F
        )
        
        if (!is_valid(fitted_model)) {
            return(paste0('The optimizer could not fit a model for ', 
                          xreg_name, '\n'))
        }
        
        # if (! 
        #     ((xreg_name %in% names(fitted_model$coef)) || 
        #      ('xreg' %in% names(fitted_model$coef)))
        #     ) {
        #     return(paste0('The optimizer could not add ', xreg_name, ' to the model\n'))
        # }
        
        fitted_model$opt_lag <- opt_lag
        return(fitted_model)
    }
    
    
    # ------------------------------- assertions -------------------------------
    if (!any(c('mts', 'ts') %in% class(xregs))) {
        stop('[drm.select] TypeError. Parameter `xregs` must be `ts` or `mts`')
    }
    if (class(serie) != 'ts') {
        stop('[drm.select] TypeError. Parameter `serie` must be `ts`')
    }
    # --------------------------------------------------------------------------
    
    # Global variables initialization
    global_ic <- Inf
    global_fit <- NA
    history <- data.frame(var=NA, lag=NA, ic=NA)
    
    # Note that response is the variable used for computing the optimal lag, 
    # and it will change iteratively by storing the residuals of the current 
    # model
    response <- serie
    
    # Cluster initialization to parallelize the covariates selection
    local_variables <- ls()
    global_variables <- names(.GlobalEnv)
    cl <- makeCluster(detectCores(logical=F))
    clusterExport(cl, local_variables, envir=environment())
    clusterExport(cl, global_variables)
    clusterEvalQ(cl, library(fpp2))
    clusterEvalQ(cl, library(tseries))
    clusterEvalQ(cl, library(TSA))
    clusterEvalQ(cl, library(seastests))
    clusterEvalQ(cl, library(forecast))
    clusterEvalQ(cl, library(stringi))
    clusterEvalQ(cl, library(stringr))
    
    # Get the maximum lag of all covariates
    max_lag <- get.max.lag()
    
    
    # Global loop that adds iteratively regressor variables
    for (i in 1:ncol(xregs)) {
        
        # Consider that no lag could be computed with get.max.lag()
        if (max_lag == -Inf) {
            break
        }
        
        # Construct a list with non-added covariates
        xregs_indexes <- (1:ncol(xregs))[!(colnames(xregs) %in% history$var)]
        xregs_list <- as.list(xregs_indexes)
        names(xregs_list) <- xregs_indexes
        
        # Parallelize the model fitting
        fitted_models <- parLapply(cl, xregs_list, add.regressor)
        
        if (show_info) { 
            display_results(fitted_models) 
        }
        
        # Obtain ICs of each model
        ics <- unlist(lapply(fitted_models, 
                             function(x) ifelse(is.character(x), NA, x[[ic]])))
        
        # If no covariate has been added to the model , the result of 
        # add.regressor() has no "ic" item
        if (all(is.na(ics))) {
            break
        }
        
        # Select the best covariate (minimum IC)
        xreg_selected <- names(ics)[which.min(c(ics))]
        
        # If the new IC does not improce the global IC achieved, then its 
        # corresponding covariate is not selected
        if (global_ic < fitted_models[[xreg_selected]][[ic]]) {
            break
        }
        
        # Otherwise, update the global variables
        global_fit <- fitted_models[[xreg_selected]]
        global_ic <- global_fit[[ic]]
        xreg_selected <- as.numeric(xreg_selected)
        xreg_selected_lag <- global_fit$opt_lag
        
        if (show_info) {
            cat(paste0('Covariate ', colnames(xregs)[xreg_selected], 
                       ' has been added [', ic, '=', global_ic, 
                       ', lag=', xreg_selected_lag, ']\n'))
            print(global_fit, row.names=F)
            cat(paste0(stri_dup('-', PAD), '\n'))
        }
        
        # Add the new added covariate to the historical: (xreg_name, xreg, lag)
        history[i, ] <- c(colnames(xregs)[xreg_selected], 
                          xreg_selected_lag, global_ic)
        
        # Now response is the regression residuals of the global model 
        response <- residuals(global_fit, type='regression')
    }
    
    # Once we have broken the loop, no more covariates will be added
    if (show_info) {
        cat('No more variables will be added\n')
    }
    stopCluster(cl)
    
    
    # If some covariates have been added (at least one) and the residuals of 
    # the global model are not stationary, we need them to be stationary, thus, 
    # fit an ARIMAX model where d=0 and D=0.
    if (!any(is.na(history)) && (sum(global_fit$arma[6:7]) > 0)) {
        
        if (show_info) {
            cat(paste0('The global model does not have stationary errors\n', 
                       'Trying to adjust a model that do have stationary errors\n'))
        }
        
        data_new <- construct.data(history, serie, xregs, new_xreg_name = NULL, 
                                   opt_lag = NULL, max_lag=max_lag)
        global_fit <- auto.fit.arima(data_new[, c(1)], xregs=data_new[, -c(1)],
                                     ic=ic, d=0, D=0, alpha=alpha, show_info=F)    
    }
    
    
    # If the resulting model is valid, return it
    if (is_valid(global_fit)) {
        global_ic <- global_fit[[ic]]
        
        if (show_info) {
            cat(paste0(
                stri_dup('-', PAD), '\n|',
                str_pad(
                    paste0(
                        'Historical of added covariates to the model ', 
                        '(ndiff=', ndiff, ')'), width=PAD-2, side='both',
                        pad=' '), 
                '|\n',
                stri_dup('-', PAD), '\n'
                )
            )
            print(history, row.names=F)
            cat(paste0(stri_dup('-', PAD), '\n'))
            print(global_fit, row.names=F)
        }
        
        
        global_fit$ndiff <- ndiff
        global_fit$history <- history
        return(global_fit)
    }
    
    
    # Here the model is not valid nor could be optiized with stationary errors.
    # We need to:
    # a) Apply regular differentiation to all data and call the function.
    # b) If ndiff=2 it is recomeneded not applying more differentiations, so 
    # fit an ARIMA model with no covariates
    if (ndiff < 3) {
        if (show_info) {
            cat(paste0(
                'No valid model with stationary errors could be optimized\n',
                'Applying regular differentiation (ndiff=', ndiff+1, 
                ') and calling again the function\n'))
            cat(paste0(stri_dup('-', PAD), '\n',
                       stri_dup('-', PAD), '\n'
                       ))
            
        }
        
        # Differentiate all data and call the function again
        serie <- diff(serie)
        xregs <- diff(xregs)
        return(
            drm.select(serie, xregs, ic, alpha, st_method, show_info, ndiff+1)
        )
    }
    
    # Othewise fit an ARIMA model with no covariates
    global_fit <- auto.fit.arima(serie, NULL, ic, NA, NA, alpha, F, F)
    
    if (!is_valid(global_fit)) {
        warning('No valid model could be optimized')
        return(NA)
    } 
    
    global_ic <- global_fit[[ic]]
    
    if (show_info) {
        cat(paste0(stri_dup('-', PAD), '\n'))
        cat(paste0('Model without covariates  (ndiff=', ndiff, ') [', 
                   ic, '=', global_ic, ']\n'))
        cat(paste0(stri_dup('-', PAD), '\n'))
        print(global_fit, row.names=F)
        cat(paste0(stri_dup('-', PAD), '\n'))
    }
    
    global_fit$ndiff <- ndiff
    
    return(global_fit)
}





#' Selection of the maximum significative correlation lag for a time series and 
#' its regressor variable. Optionally, a maximum magnitude can be fixed and it 
#' is possible to constraint the lagging to be negative or null.
#' 
#' 
#' @param xreg [ts]: Univariate time series (covariate).
#' @param serie [ts]: Univariate time series (dependent).
#' @param alpha [numeric]: Significance level for hypothesis tests of trend and 
#' seasonality.
#' @param max_lag [numeric]: Maximum magnitud of the selected lag.
#' @param method [character]: Method used to check stationary of the series 
#' (note that we are asumming trend approximates stationary).
#' 
#' @returns Optimal correlation lag between the covariate and the dependent 
#' variable based of the method proposed by @cryer2008time. If it is not 
#' possible to obtain a lag, it returns `NA`.
get.opt.lag <- function(xreg, serie, alpha=0.05, max_lag=NA, 
                        method='auto.arima', less0=T) {
    if (!method %in% c('auto.arima', 'adf.test')) {
        stop('[get.opt.lag] ArgumentError: `method`  must be "auto.arima" or "adf.test"')
    }
    
    # Check if some variable has trend component
    y_trend <- has_trend(serie, method=method)
    x_trend <- has_trend(xreg, method=method)
    
    # Regular differentiation to remove trend component
    while (y_trend || x_trend) {
        serie <- diff(serie)
        xreg <- diff(xreg)
        y_trend <- has_trend(serie, method=method, alpha=alpha)
        x_trend <- has_trend(xreg, method=method, alpha=alpha)
    }
    
    # Check if some variable has seasonal component
    y_seasonal <- frequency(serie) > 1 && isSeasonal(serie)
    x_seasonal <- frequency(xreg) > 1 && isSeasonal(xreg)
    
    # Seasonal differentiation to remove the seasonal component
    while (y_seasonal || x_seasonal) {
        serie <- diff(serie, lag=frequency(serie))
        xreg <- diff(xreg, lag=frequency(xreg))
        
        y_seasonal <- isSeasonal(serie)
        x_seasonal <- isSeasonal(xreg)
    }
    
    # Prewhitening
    series_prewhiten <- TSA::prewhiten(xreg, serie, plot=F)
    
    # Obtain lags and criss-correlations values
    series_lags <- series_prewhiten$ccf$lag
    series_acfs <- series_prewhiten$ccf$acf
    if (less0) {
        series_acfs <- series_acfs[series_lags <= 0]
        series_lags <- series_lags[series_lags <= 0]
    }
    
    # Filter by the max_lag
    if (!is.na(max_lag)) {
        series_acfs <- series_acfs[abs(series_lags) <= max_lag]
        series_lags <- series_lags[abs(series_lags) <= max_lag]
    }
    
    # Select significative values
    n_used <- series_prewhiten$ccf$n.used
    significative_lags <- abs(series_acfs) >= qnorm(1-alpha/2)/sqrt(n_used)
    series_lags <- series_lags[significative_lags]
    series_acfs <- series_acfs[significative_lags]
    
    # Select maximum correlation value (absolute) and return it
    opt_lag <- series_lags[which.max(abs(series_acfs))]
    
    if (length(opt_lag) == 0) {
        return(NA)
    }
    
    return(opt_lag*frequency(serie))
}






#' ----------------------------- AUXILIAR FUNCTIONS-----------------------------
#' In this section of the script we provide some auxiliar functions that are not 
#' related to covariates selection, but were needed to implement the selection 
#' method
#' -----------------------------------------------------------------------------


#' Auxiliar function to construct a `mts` matrix given the historical of the 
#' selection method and a new covariate to add to the model. 
#' 
#' @param history [data.frame]: Historical of added covariates (name, lag, IC 
#' obtained in the fitting).
#' @param serie [ts]: Univariate time series that act as the dependent variable.
#' @param xregs [mts]: Set of covariates.
#' @param new_xreg_name [character]: Name of the new covariate to be added to 
#' the model. If it is `NULL` we assume that no covariate is being added and 
#' the covariates matrix will be constructed with those that are mentioned in 
#' the historical.
#' @param opt_lag [numeric]: Only if `new_xreg_name` is not `NULL`. Lag where 
#' the significative correlation occurs in with the dependen variable.
#' @param max_lag [numeric]: Maximum lag computed at the initialization of the 
#' selection method. This function must cut all data considering this value.
#'
#' @returns Matrix of lagged and cut covariates to introduce them in the 
#' `auto.fit.arima()` function.
#'
construct.data <- function(history, serie, xregs, new_xreg_name, opt_lag, max_lag) {
    
    # Cut the first max_lag values of the dependent varialbe
    serie <- window(serie, start=add_t(start(serie), max_lag, frequency(serie)))
    data <- matrix(serie)
    
    # If the historical is not empy, add covariates to data
    if (!all(is.na(history))) {
        
        for (j in 1:nrow(history)) {
            xreg <- xregs[, c(history$var[j])]                        # covariate
            xreg <- stats::lag(xreg, as.integer(history$lag[j]))      # lagged covariate
            xreg <- window(xreg, start=start(serie), end=end(serie))  # cut covariate
            data <- cbind(data, xreg)
            
        }
        colnames(data) <- c('serie', history$var)
    } 
    
    # Add the new covariate
    if (!is.null(new_xreg_name)) {
        
        new_xreg <- stats::lag(xregs[, c(new_xreg_name)], opt_lag)        # lagged new covariate
        new_xreg <- window(new_xreg, start=start(serie), end=end(serie))  # cut new covariate
        data <- cbind(data, new_xreg)
        
        if (!all(is.na(history))) {
            colnames(data) <- c('serie', history$var, new_xreg_name)
        } else { colnames(data) <- c('serie', new_xreg_name) }
    }
    
    return(data)
}


#' Auxiliar function to add $t$ moments to a given moment of the `ts` class 
#' (handling seasonal data).
#' 
#' @param moment [vector]: Vector of two values that give information about 
#' the current moment. Note that the meaning of each value is the same of `ts` 
#' indexes.
#' @param t [numeric]: Number of moments that we want to add.
#' @param freq [numeric]: Frequency of the time series. It is 1 if there is 
#' not a seasonal component. 
#' @returns Moment resulting from adding `t` moments to the current moment.
#' @examples 
#' --------------------------------------
#' > add_t(c(6, 3), 6, 7)
#' c(7,2)
#' --------------------------------------
#' > add_t(c(6, 1), 6, 1)
#' c(12, 1)
#' --------------------------------------

add_t <- function(moment, t, freq=1) {
    if (freq == 1) {
        return(c(moment[1] + t, moment[2]))
    } else {
        raw_moment <- moment[2] + t
        return(c(moment[1] + floor(raw_moment/freq), raw_moment %% freq))
    }
}


#' Check the trend of a time series
#' 
#' @param x [ts]: Univariate time series.
#' @param method [character]: Method used to check if the serie has trend. If 
#' it is `auto.arima`, uses `auto.arima()` function to fit an ARIMA model and 
#' check the $d$ order. If it is `adf.test`, the Dickey-Fuller test is used.
#' @param alpha [numeric]: Signficance level for the Dickey-Fuller test. By 
#' default, $0.05$. 

#' @returns `True` if `x` has trend, `False` otherwise
#'
has_trend <- function(x, method='auto.arima', alpha=0.05) {
    if (!method %in% c('auto.arima', 'adf.test')) {
        stop('[has_trend] ArgumentError: `method` must be "auto.arima" or "adf.test"')
    }
    
    if (method == 'auto.arima') {
        fitted_model <- suppressWarnings(auto.arima(x, max.d=3))
        d <- fitted_model$arma[6]
        return(d>0)
    }
    if (method == 'adf.test') {
        pvalue <- suppressWarnings(adf.test(x[!is.na(x)])$p.value)
        return(alpha <= pvalue)
    }
}
