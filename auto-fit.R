
#' README:
#' This script contains the necessary functions to automatically fit ARIMA and 
#' ARIMAX models with statistically significative coefficients with valid 
#' residuals (i.e. with zero mean and independent)



PAD <- 86  # display-info parameter 



#' **Automatic fitting of ARIMA or ARIMAX model**
#' 
#' Implementaton of the ARIMA or ARIMAX model selection that optimizes selected 
#' information criterion and satisifies that:
#' 
#' 1) All coefficients are statistically significative.
#' 2) Model residuals have zero mean and are independent.
#' 
#' 
#' @param serie [ts]: Univariate time series to fit an ARIMA model with. If 
#' `xregs` is not `NULL` it will act as the dependent variable.
#' @param xregs [mts]: Matrix of time series that act as covariates in an 
#' ARIMAX model. By default, `NULL`, meaning that an ARIMA will be fitted with 
#' `serie`.
#' @param ic [character]: Information criterion to be used in ARIMA orders' 
#' selection. Options available are the same as in `forecast::auto.arima()` 
#' function: `ic`, `bic` or `aicc`.
#' @param d [numeric]: Value of the diferentiation order. By default `NA`, thus 
#' there is no restriction about the $d$ order value.
#' @param D [numeric]: Value fo the seasonal diferentiation order. By default 
#' `NA`, thus there is no restriction about the $D$ order value.
#' @param alpha [numeric]: Significance level of hypothesis tests used for 
#' checking independence, zero mean and normality of residuals, and significance 
#' of estimated coefficients.
#' @param show_info [boolean]: Displaying or not displaying the historical of 
#' fitted models.
#' @param plot_result [boolean]: Returning or not with the fitted model a 
#' sequential plot of time series and residuals.
#' 
#' @returns [Arima] Fitted ARIMA/ARIMAX model where all estimated coefficients are 
#' significative and residuals are independent and with zero mean. In case it 
#' is no possible to optimize a model, the function returns `NA`.
#'
auto.fit.arima <- function(serie, xregs=NULL, ic='aicc', d=NA, D=NA, alpha=0.05, 
                           show_info=T, plot_result=F) {
    
    # ------------------------------- assertions -------------------------------
    if (class(serie) != 'ts') {
        stop('[auto.fit.arima] TypeError. Parameter `serie` must be `ts`')
    }
    if (!is.null(xregs) && !any(c('mts', 'ts') %in% class(xregs))) {
        stop(paste0('[auto.fit.arima] ', 
                    'TypeError. Parameter `xregs` must be `ts` or `mts`'))
    }
    # --------------------------------------------------------------------------
    
    # Save the output of auto.arima() function
    trace <- capture.output({
        aux <- suppressWarnings(
            auto.arima(serie, xreg=xregs, d=d, D=D, seasonal=frequency(serie)>1, ic=ic, max.d=4, 
                       max.D=3, stepwise=F, approximation=F, trace=T) 
        )
    })
    con    <- textConnection(trace)
    models <- read.table(con, sep=':')
    models <- models[1:nrow(models)-1,]
    names(models) <- c('Model', ic)
    models <- sapply(models, trimws)
    models <- as.data.frame(models)
    
    # Global loop. It will be fitting models following the order (ascendent) of the IC
    while (nrow(models) > 0) {
        
        # Get the model with the minimum IC 
        best_model <- models[which.min(models[[ic]]), 1]
        
        # Get its orders from the string with parse.orders()
        orders <- parse.orders(best_model, seasonal=frequency(serie)>1, 
                               regressor=!is.null(xregs))
        
        # Try to fit ARIMA(X) model with fit.model()
        fitted_model <- fit.model(serie, xreg=xregs, orders)
        
        # If the result is not valid (with is_valid() function), remove the
        # corresponding description of `models` and continue the selection
        if (!is_valid(fitted_model)) {
            models <- models[-c(which.min(models[[ic]])), ]
            next
        }
        if (show_info) {
            cat(paste0(stri_dup('-', PAD), '\n'))
            print(fitted_model, row.names=F)
            cat(paste0(stri_dup('-', PAD), '\n'))
        }
        
        # Otherwise, continue by deleting non-significative coefficients with 
        # fit.coefficients() function
        fitted_model <- fit.coefficients(fitted_model, alpha, show_info)
        
        # Again, the model might not be "valid"
        if (!is_valid(fitted_model)) {
            models <- models[-c(which.min(models[[ic]])),]
            next
        }
        
        # ----------------------------------------------------------------------
        # Here, if no optimization problem has been found, the fitted ARIMA(X) 
        # is valid, so we must check the residuals' properties
        # 1) Independence (LjungBox)
        H <- ljungbox_lag(fitted_model$residuals)
        ljungbox <- Box.test(fitted_model$residuals, lag=H, type='Ljung-Box', 
                             fitdf=sum(fitted_model$coef!=0)
                    )$p.value
        
        # 2) Zero mean (t-test)
        ttest    <- t.test(fitted_model$residuals, mu=0)$p.value
        
        # 3) Normality (JarqueBera and Shapiro-Wilks)
        testJB   <- jarque.bera.test(
            fitted_model$residuals[!is.na(fitted_model$residuals)])$p.value
        testSW   <- shapiro.test(fitted_model$residuals)$p.value
        
        if (ljungbox < alpha) {
            models <- models[-c(which.min(models[[ic]])),] 
            if (show_info) {
                cat(paste0(
                    'Independence hypothesis of residuals is rejeceted\n',
                    'Model is not valid. Trying with the next model following ', 
                    ic, '\n', stri_dup('-', PAD), '\n'
                    )
                )
            }
            next
        }
        if (ttest < alpha) {
            models <- models[-c(which.min(models[[ic]])),]    # eliminamos el modelo del historial
            if (show_info) {
                cat(paste0(
                    'Zero mean hypothesis of residuals is rejected\n',
                    'Model is not valid. Trying with the next model following ',
                    ic, '\n', stri_dup('-', PAD), '\n'
                ))
            }
            next
        }
        
        if (any(c(testJB, testSW) < alpha)) {
            if (show_info) {
                cat(paste0(
                    'Normality hypothesis is rejected\nModel is valid ',
                    'but forecasting asuming normality is not available\n', 
                    stri_dup('-', PAD), '\n'
                ))
            }
        }
        
        # Here, the fitted model and its residuals are valid, so exit
        break
    }
    
    # Check that the model is valid. These lines are needed since it is possible 
    # that no fitted valid model can be optimized
    if (!is_valid(fitted_model))  {
        if (show_info) {
            warning(
                'No valid model could be fitted for this time serie'
            )
        }
        return(NA)   # in this case, return NA
    }
    
    
    # --------------------------------------------------------------------------
    # Here the function states that the fitted model is valid
    if (show_info) {
        cat(paste0(stri_dup('-', PAD), '\n'))
        cat(paste0('|', str_pad('FINAL MODEL', width=PAD-2, side='both', pad=' '), '|\n'))
        cat(paste0(stri_dup('-', PAD), '\n'))
        print(fitted_model, row.names=F)
    }
    
    # It is possible to add to the result some plots (time series and residuals)
    if (plot_result) {
        fitted_model$fig_serie <- suppressWarnings(plot_serie(serie, alpha=alpha))
        fitted_model$fig_residuals <- suppressWarnings(
            plot_residuals(fitted_model, alpha=alpha)
            )
    }
    return(fitted_model)
}




#' **Fitting of ARIMA or ARIMAX model trying with multiple optimizators**
#'
#' @param serie [ts]: Univariate time series to fit an ARIMA model with. If 
#' `xregs` is not `NULL` it will act as the dependent variable.
#' @param orders [list]: Information about the orders of the ARIMA model. It 
#' corresponds to the output of `parse.orders()`. It must have the following 
#' items:
#' * `regular`: vector of regular orders (p, d, q).
#' * `seasonal`: vector of seasonal orders (P, D, Q).
#' * `include_mean`: boolean value indicating if a constant value is included 
#' in the definition of the model.
#' @param xregs [mts]: Matrix of time series that act as covariates in an 
#' ARIMAX model. By default, `NULL`, meaning that an ARIMA will be fitted with 
#' `serie`.
#' @param fixed [vector]: Boolean vector which masks the coefficients set to 0.
#' @param show_info [boolean]: Indicates if information about the fitted models 
#' is displayed in the console. 
#'
#' @returns Fitted model or `NA` in case no optimizer could fit coefficients of 
#' the model.
#'
fit.model <- function(serie, orders, xregs=NULL, fixed=NULL, show_info=F) {
    
    # ------------------------------- ASSERTIONS -------------------------------
    total_params <- get.total.params(orders, xregs)
    
    if (!is.null(fixed) && (length(fixed) != total_params)) { 
        stop('`fixed` size is not the same as the number of parameters')
    }
    # --------------------------------------------------------------------------
    
    optimizers <- c('BFGS', 'Nelder-Mead', 'CG', 'L-BFGS-B', 'SANN', 'Brent')
    
    for (opt in optimizers) {
        fitted_model <- try(
            suppressWarnings(
                Arima(serie, xreg=xregs, order=orders$regular, 
                      seasonal=orders$seasonal, fixed=fixed, 
                      include.mean=orders$include_mean)),
            silent=T)
        
        if (is_valid(fitted_model)) {
            fitted_model$fixed <- fixed
            return(fitted_model)
        }
    }
    
    # Here no optimizer could estimate the parameters of the ARIMA(X) model
    if (show_info) {
        warning('No model could be optimized\n') 
    }
    return(NA)
}


#' **Coefficients estimation in ARIMA or ARIMAX model**
#' 
#' Checks the significance of all ARIMA(X) coefficients and iteratively sets 
#' those non significative to zero and readjusts the model.
#' 
#' @param fitted_model [Arima]: Initial fitted ARIMA model. Its coefficients 
#' might be non significative.
#' @param orders [list]: Model orders. Returned value of `parse.orders()`.
#' @param alpha [numeric]: Significance level used for hypothesis tests of 
#' coefficients.
#' @param show_info [boolean]: Indicates showing or not showing the parameters 
#' that are set to zero and the resulting fit.
#'
#' @return Fitted ARIMA(X) model (if it exists) with significative coefficients.
#'
fit.coefficients <- function(fitted_model, alpha=0.05, show_info=T) {
    
    stat <- qnorm(1-alpha/2)                  # pivot
    fixed <- rep(NA, length(fitted_model$coef))     # initial fixed vector
    
    
    # Start an endless loop to search non-significative coefficients
    while (TRUE) {
        
        # Obtain the "first" non significative coefficient and set its 
        # corresponding value in `fixed` to zero
        remov <- get.nonsignificative(fitted_model, alpha)
        
        if (is.na(remov)) {         # all coefficients are significative
            return(fitted_model)
        }
        
        fixed[names(fitted_model$coef) == remov] <- 0 
        
        if (show_info) {
            cat(paste0('Removing parameter ', remov, 
                       ' since it is not significative\n'))
        }
        
        # Update the orders of the model add the fixed vector. This lines are 
        # needed since, if some coefficients are set to zero, some orders of 
        # the ARIMA(X) might change
        orders_fixed_update <- update.orders(fitted_model, fixed)
        orders <- orders_fixed_update$orders
        fixed <- orders_fixed_update$fixed
        
        # Fit a new ARIMA(X) model
        fitted_model <- fit.model(fitted_model$x, xregs=fitted_model$xreg, 
                            orders=orders, fixed=fixed)
        
        # Check if the new fitted model is valid
        if (!is_valid(fitted_model)) { 
            if (show_info) {
                cat(paste0('No ARIMA model could be optimized\n'))
            }
            return(NA) 
        }
        
        if (show_info) {
            cat(paste0(stri_dup('-', PAD), '\n'))
            print(fitted_model, row.names=F)
            cat(paste0(stri_dup('-', PAD), '\n'))
        }
    }
    
    # Once no more coefficients are set to zero, return the resulting model
    return(fitted_model)
}



#' ----------------------------- AUXILIAR FUNCTIONS-----------------------------
#' In this section of the script we provide some auxiliar functions that are not 
#' related to ARIMA(X) model fitting, but were needed to implement the selection 
#' method
#' -----------------------------------------------------------------------------

#' Obtains ARIMA orders given the description of the trace of 
#' `forecast::auto.arima()`. Test with
#' 
#' > parse.orders('ARIMA(2,1,0)(0,1,3)[12] with non-zero mean', seasonal=T)
#' --------------------------
#' $regular
#' [1] 2 1 0
#' 
#' $seasonal
#' [1] 0 1 3
#' 
#' $include_mean
#' [1] TRUE
#' --------------------------
#' 
#' @param description [character]: Description of the ARIMA model (provided by 
#' `forecast::auto.arima()`).
#' @param seasonal [boolean]: Indicates if the model has a seasonal component.
#' @param regressor [boolean]: Indicates if the model has covariates (i.e. if 
#' it is an ARIMAX).
#' @param returns A List object with the vector of regular orders, seasonal 
#' orders and if it includes a constant value (only if it has not seasonal 
#' component and no regressor variables).
parse.orders <- function(description, seasonal=F, regressor=F) {
    info <- gsub('[\\(\\)]', '', 
                 regmatches(description, gregexpr('\\(.*?\\)', description))[[1]])
    regular_order <- unlist(lapply(strsplit(info[1], ','), as.double))
    if (length(info) > 1) {
        seasonal_order <- unlist(lapply(strsplit(info[2], ','), as.double))
    } else {seasonal_order <- c(0, 0, 0)}
    
    
    if (grepl('non-zero mean', description, fixed=T) || 
        grepl('drift', description, fixed=T) || regressor) {
        include_mean <- TRUE
    } else {
        include_mean <- FALSE
    }
    
    if (regressor && seasonal) {
        include_mean <- FALSE
    }
    
    orders <- list(regular=regular_order, seasonal=seasonal_order, 
                   include_mean=include_mean)
    return(orders)
}


#' Calculates the number of parameters/coefficients needed given the orders of 
#' the ARIMA(X) model and the matrix of regressor variables.
#' @param orders [list]: List of ARIMA(X) model (object returned by 
#' `parse.orders()`).
#' @param xregs [mts]: Matrix of covariates of the ARIMAX. By default `NULL`, 
#' meaning that there are no covariates.
#' @return The number of coefficients/parameters of the ARIMA or ARIMAX model.
get.total.params <- function(orders, xregs=NULL) {
    total_params <- sum(
        orders$regular[1], orders$regular[3],
        orders$seasonal[1], orders$seasonal[3],
        ifelse(orders$include_mean, 1, 0),
        ifelse(is.null(xregs), 0, ifelse(is.null(ncol(xregs)), 1, ncol(xregs)))
    )
    return(total_params)
}



#' Updates the orders of a fitted ARIMA/ARIMAX model given the vector `fixed` 
#' that defines which coefficients are set to zero.
#' 
#' @param fitted_model [Arima]: Fitted ARIMA model.
#' @param fixed [vector]: Vector mask that indicates which parameters of the 
#' ARIMA(X) model are set to zero.
#' @returns List object with new orders (item `$orders`) and an updated vector 
#' `fixed` (item `$fixed`).
#' 
update.orders <- function(fitted_model, fixed) {
    
    # Consider no coefficient has been set to zero
    if (is.null(fixed) || all(is.na(fixed))) {
        return(list(orders=get.orders(fitted_model), fixed=fixed))
    }
    
    # Obtain coefficients names that are set to zero
    remov_coefs <- names(fitted_model$coef[!is.na(fixed)])
    arma_orders <- fitted_model$arma[1:4]
    patterns <- c('^ar\\d+', '^ma\\d+', '^sar\\d+', '^sma\\d+')
    
    # Variables where new values will be stored
    new_fixed <- c()
    new_orders <- arma_orders
    
    # For each order (AR, MA, SAR, SMA), update its value and elements in the 
    # vector fixed
    for (i in 1:4) {
        new_order <- update.order(arma_orders[i], str_subset(remov_coefs, patterns[i]))
        new_orders[i] <- new_order
        if (new_order > 0) {
            start_index <- ifelse(i>1, sum(arma_orders[1:i-1]), 0)
            new_fixed <- c(new_fixed, fixed[(start_index+1):(start_index+new_order)])   
        }
    }
    
    # Detect if the ARIMA(X) constant value has been set to zero, update the 
    # order and the vector fixed
    if (any(str_detect(remov_coefs, 'mean|drift|intercept')) ) {
        include_mean <- FALSE
        if (sum(arma_orders, 2) <= length(fixed)) {
            new_fixed <- c(new_fixed, fixed[sum(arma_orders, 2):length(fixed)])
        }
    } else  {
        include_mean <- any(str_detect(names(fitted_model$coef), 'mean|drift|intercept'))
        if (sum(arma_orders) < length(fixed)) {
            new_fixed <- c(new_fixed, fixed[sum(arma_orders,1):length(fixed)])
        }
    }
    
    new_orders <- list(
        regular=c(new_orders[1], fitted_model$arma[6], new_orders[2]),
        seasonal=c(new_orders[3], fitted_model$arma[7], new_orders[4]),
        include_mean=include_mean)
    
    return(list(orders=new_orders, fixed=new_fixed))
    
}





#' 
#' Obtain the orders of the ARIMA or ARIMAX model from the Arima object.
#'
#' @param fitted_model [Arima]: Fitted ARIMA(X) model.
#' @returns Object List with the regular (`$regular`) and seasonal (`$seasonal`) 
#' orders and a boolean variable indicating if a constant parameter is included 
#' (`$include_mean`).
#'
get.orders <- function(fitted_model) {
  if (!('Arima' %in% class(fitted_model))) {
    stop('El argumento `fitted_model` debe ser un objeto Arima')
  }
  
  p <- fitted_model$arma[1]
  q <- fitted_model$arma[2]
  P <- fitted_model$arma[3]
  Q <- fitted_model$arma[4]
  d <- fitted_model$arma[6]
  D <- fitted_model$arma[7]
  
  include_mean <- any(c('mean', 'drift', 'intercept') %in% names(fitted_model$coef))
  orders <- list(regular=c(p, d, q), seasonal=c(P, D, Q), include_mean=include_mean)
  return(orders)
}



#' Updates an order given the names of the removed coefficients of that order.
#' 
#' @param order [numeric]: Order value of the ARIMA(X) model.
#' @param remov_coefs [vector]: Vector of the removed coefficients of the given 
#' order.
#' @returns The updated order.
#' 
#' > update.order(4, 'ma4')
#' --------------------------
#' 3
#' --------------------------
#'

update.order <- function(order, remov_coefs) {
    if (length(remov_coefs) == 0) {return(order)}
    while (as.character(order) %in% str_extract(remov_coefs, '\\d+')) {
        order <- order - 1
    }
    return(order)
}



#' Obtains a non-significative coefficient (considering the pvalue of the 
#' significance test)
#' 
#' @param fitted_model [Arima]: Fitted ARIMA or ARIMAX model.
#' @param alpha [numeric]: Significance level for hypothesis tests.
#'
get.nonsignificative <- function(fitted_model, alpha) {
    stat <- qnorm(1-alpha/2)
    
    # Obtain those those coefficients that are not set to zero
    coefs_mask <- fitted_model$coef != 0
    coefs <- fitted_model$coef[coefs_mask]
    
    # If there is no coefficient different than zero, return NA
    if (length(names(coefs)) == 0) {
        return(NA)
    }
    # Obtain a mask of those non-significative coefficients
    coefs_sd <- suppressWarnings(sqrt(diag(fitted_model$var.coef)))
    nonsig <- abs(coefs) < stat*coefs_sd        # TRUE -> non-significative
    
    # Consider all coefficients are significative
    if (!any(nonsig)) {
        return(NA)
    }
    
    # Obtain the coefficient that must be set to zero (smaller pvalue)
    remov <- names(which.min(abs(coefs)/(stat*coefs_sd)))
    return(remov)
}




#' Check if the fitted ARIMA/ARIMAX model has been correctly estimated.
#' 
#' @param fitted_model [Arima]: Fitted ARIMA or ARIMAX model
#' @returns `TRUE` if the model is valid, `FALSE` otherwise.
is_valid <- function(fitted_model) {
    
    if (length(fitted_model) == 1 && class(fitted_model) == 'try-error') {
        return(FALSE)
    }
    else if (all(is.na(fitted_model))) {
        return(FALSE)
    }
    else if (suppressWarnings(any(is.na(sqrt(diag(fitted_model$var.coef)))))) {
        return(FALSE)
    } else if (suppressWarnings(any(is.na(fitted_model$var.coef)))) {
        return(FALSE)
    }
    return(TRUE)
}


#' Obtains the optimal lag for the LjungBox tests
#'
#' @param serie [ts]: Time series to check LjungBox test with.
#'
ljungbox_lag <- function(serie) {
  m <- frequency(serie)
  n <- length(serie)
  if (m > 1) {
    return(floor(min(2*m, n/5)))
  } else {
    return(floor(min(10, n/5)))
  }
}




