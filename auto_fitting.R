
#' README:
#' This script contains the necessary functions to automatically fit ARIMA and 
#' ARIMAX models with statistically significative coefficients with valid 
#' residuals (i.e. with zero mean and independent)



PAD <- 86  # display-info parameter 



#' **Automatic fitting of ARIMA or ARIMAX model**
#' 
#' Implementaton of the ARIMA or ARIMAX model selection that optimizes an 
#' information criterion selected and satisifies that:
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
#' @param show_info [boolean]: Boolean value for displaying or not displaying 
#' the historical of fitted models.
#' 
#' @returns Fitted ARIMA/ARIMAX model where all estimated coefficients are 
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
                    ic, '\n', stri_dup('-', PAD, '\n')
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
                    ic, '\n', stri_dup('-', PAD, '\n')
                ))
            }
            next
        }
        
        if (any(c(testJB, testSW) < alpha)) {
            if (show_info) {
                cat(paste0(
                    'Normality hypothesis is rejected\nModel is valid ',
                    'but forecasting asuming normality is not available'
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
fit.coefficients <- function(ajuste, alpha=0.05, show_info=T) {
    
    stat <- qnorm(1-alpha/2)                  # estadístico de contraste
    fixed <- rep(NA, length(ajuste$coef))     # valor inicial de los coeficientes
    
    
    while (TRUE) {
        
        # Obtención de los nombres de los coeficientes ARMA
        # arma_coefs_names <- get_arma_coefs_names(ajuste)
        
        # En caso de que no haya coeficientes ARIMA, se devuelve el ajuste 
        # if (is.null(arma_coefs_names)) { return(ajuste) }
        
        # Obtenemos los valores de dichos coeficientes y sus desviaciones típicas
        # arma_coefs <- ajuste$coef[arma_coefs_names]
        # arma_coefs_sd <- suppressWarnings(sqrt(diag(ajuste$var.coef[arma_coefs_names, arma_coefs_names])))
        # res <- abs(arma_coefs) < stat*arma_coefs_sd
        
        # if (!any(res)) { # en caso de que todos sean significativos, se para el bucle
        #   break
        # }
        # 
        # # Configuración del vector fixed
        # remov <- names(which.min(abs(arma_coefs)/(stat*arma_coefs_sd)))
        # fixed[names(ajuste$coef) == remov] <- 
        
        
        remov <- get_nonsig(ajuste, alpha)
        fixed[names(ajuste$coef) == remov] <- 0 
        if (is.na(remov)) { # no hay ningún coeficiente que retirar
            return(ajuste)
        }
        
        if (show_info) {
            cat(paste0('Es necesario retirar del modelo el parámetro: ', remov, '\n'))
        }
        
        # Actualización de los órdenes
        orders_fixed_update <- update_orders(ajuste, fixed)
        orders <- orders_fixed_update$orders
        fixed <- orders_fixed_update$fixed
        
        # Actualización de los coeficientes de regresión
        # xregs_update <- update_xregs(ajuste, fixed)
        # 
        # if (!is.null(xregs_update)) {
        #     fixed <- fixed[!(names(ajuste$coef) == remov)]
        #     xregs <- xregs_update
        # }
        
        
        # Ajuste del nuevo modelo
        ajuste <- fit.model(ajuste$x, xregs=ajuste$xreg, orders=orders, fixed=fixed)
        
        # Si no se puede optimizar el modelo se devuelve NA
        if (!is_valid(ajuste)) { 
            if (show_info) {
                cat(paste0('No se ha podido optimizar el modelo ARIMA(', string_concat(orders$regular), ')[',
                           string_concat(orders$seasonal), '] con esta configuración\n'))
            }
            return(NA) 
        }
        
        if (show_info) {
            cat(paste0(stri_dup('-', PAD), '\n'))
            print(ajuste, row.names=F)
            cat(paste0(stri_dup('-', PAD), '\n'))
        }
    }
    return(ajuste)
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
    
    
    if (grepl('non-zero mean', description, fixed=T) || grepl('drift', description, fixed=T) || regressor) {
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



#' Obtener los órdenes del modelo ARIMA a partir de su ajuste
#'
#' @param ajuste Ajuste del modelo ARIMA.
#'
#' @return Lista con los órdenes regulares (`regular`), estacionales (`seasonal`) y una 
#' variable booleana para indicar si se incluye media/constante (`include_mean`).
#' @export
#'
get_orders <- function(ajuste) {
  if (!('Arima' %in% class(ajuste))) {
    stop('El argumento `ajuste` debe ser un objeto Arima')
  }
  
  p <- ajuste$arma[1]
  q <- ajuste$arma[2]
  P <- ajuste$arma[3]
  Q <- ajuste$arma[4]
  d <- ajuste$arma[6]
  D <- ajuste$arma[7]
  
  include_mean <- any(c('mean', 'drift', 'intercept') %in% names(ajuste$coef))
  orders <- list(regular=c(p, d, q), seasonal=c(P, D, Q), include_mean=include_mean)
  return(orders)
}


update_orders <- function(ajuste, fixed) {
    
    if (is.null(fixed)) {return(list(orders=get_orders(ajuste), fixed=NULL))}
    if (all(is.na(fixed))) {return(list(orders=get_orders(ajuste), fixed=fixed))}
    
    remov_coefs <- names(ajuste$coef[!is.na(fixed)])
    arma_orders <- ajuste$arma[1:4]
    patterns <- c('^ar\\d+', '^ma\\d+', '^sar\\d+', '^sma\\d+')
    new_fixed <- c()
    new_orders <- arma_orders
    
    for (i in 1:4) {
        new_order <- update_order(arma_orders[i], str_subset(remov_coefs, patterns[i]))
        new_orders[i] <- new_order
        if (new_order > 0) {
            start_index <- ifelse(i>1, sum(arma_orders[1:i-1]), 0)
            new_fixed <- c(new_fixed, fixed[(start_index+1):(start_index+new_order)])   
        }
    }
    
    if (any(str_detect(remov_coefs, 'mean|drift|intercept')) ) {
        include_mean <- FALSE
        if (sum(arma_orders, 2) <= length(fixed)) {
            new_fixed <- c(new_fixed, fixed[sum(arma_orders,2):length(fixed)])
        }
    } else  {
        include_mean <- any(str_detect(names(ajuste$coef), 'mean|drift|intercept'))
        if (sum(arma_orders) < length(fixed)) {
            new_fixed <- c(new_fixed, fixed[sum(arma_orders,1):length(fixed)])
        }
    }
    
    new_orders <- list(regular=c(new_orders[1], ajuste$arma[6], new_orders[2]),
                       seasonal=c(new_orders[3], ajuste$arma[7], new_orders[4]),
                       include_mean=include_mean)
    
    return(list(orders=new_orders, fixed=new_fixed))
    
}

update_order <- function(order, remov_coefs) {
    if (length(remov_coefs) == 0) {return(order)}
    while (as.character(order) %in% str_extract(remov_coefs, '\\d+')) {
        order <- order - 1
    }
    return(order)
}






#' Obtención de los coeficientes ARMA a chequear de un ajuste
#'
#' @param ajuste Ajuste ARIMA sobre el que se quieren obtener los nombres de los 
#' coeficientes correspondientes a la parte ARMA (se excluyen coeficientes que se 
#' correspondan con las variables regresoras en un mdoelo ARIMAX). Esto es, `ar`, `ma`, 
#' `sar`, `sma` y, en caso de que d=0, se escoge `mean`, y en caso de que d=1, se escoge 
#' `drift`.
#'
#' @return Vector con los nombres de los coeficientes que se deben chequear. En caso de 
#' que no haya ninguno, se devuelve `NULL`.
#' @export
#'
get_arma_coefs_names <- function(ajuste) {
  if (!('Arima' %in% class(ajuste))) {
    stop('El argumento `ajuste` debe ser un objeto Arima')
  }
  
  # Obtenemos la máscara de los coeficientes ARMA
  arma_coefs_mask <- sapply(names(ajuste$coef), str_detect, pattern=c('s?ma\\d+', 's?ar\\d+', 'mean', 'drift', 'intercept'))
  
  # En caso de que la máscara esté vacía esto indica que el modelo es un ARIMA(0, 0, 0)
  if (length(arma_coefs_mask) == 0) { return(NULL) }
  
  # La máscara se calcula sobre cualquier nombre del coeficiente que cumpla algún pattern
  arma_coefs_mask <- apply(arma_coefs_mask, 2, sum) == 1
  
  # Se retiran aquellos que ya sean iguales a 0
  arma_coefs_mask <- arma_coefs_mask & (ajuste$coef != 0)
  
  # En caso de que ya no haya ningún coeficiente (todos están fijados a 0), se devuelve NULL
  if (sum(arma_coefs_mask) == 0) { return(NULL) }
  
  # Se obtienen sus nombres
  arma_coefs_names <- names(ajuste$coef[arma_coefs_mask])
  
  return(arma_coefs_names)
}

#' Obtención del coeficiente no significativo (considerando el p-valor)
#' 
#'
get_nonsig <- function(ajuste, alpha) {
    stat <- qnorm(1-alpha/2)
    
    # obtenemos los coeficientes que no están fijados a 0
    coefs_mask <- ajuste$coef != 0
    coefs <- ajuste$coef[coefs_mask]
    
    # si no hay ningún coeficiente con el primer requisito, devolver NA
    if (length(names(coefs)) == 0) {
        return(NA)
    }
    # obtenemos una máscara de los coeficientes NO significativos
    coefs_sd <- suppressWarnings(sqrt(diag(ajuste$var.coef)))
    nonsig <- abs(coefs) < stat*coefs_sd        # TRUE -> no significativo
    
    if (!any(nonsig)) {
        return(NA)
    }
    
    
    # obtenemos el coeficiente a fijar a 0
    remov <- names(which.min(abs(coefs)/(stat*coefs_sd)))
    
    return(remov)
}

update_xregs <- function(ajuste, fixed) {
    
    xregs_coefs_mask <- sapply(names(ajuste$coef), str_detect, pattern=c('s?ma\\d+', 's?ar\\d+', 'mean', 'drift', 'intercept'))
    xregs_coefs_mask <- !(apply(xregs_coefs_mask, 2, sum) == 1)
    remov_mask <- (!is.na(fixed)) & xregs_coefs_mask
    
    if (!any(remov_mask)) {
        return(NULL)
    }
    
    xregs <- ajuste$xreg[, names(ajuste$coef)[xregs_coefs_mask & (!remov_mask)]]
    return(xregs)
}




string_concat <- function(arr, delimiter=',') {
    result <- ''
    for (item in arr[1:(length(arr)-1)]) {
        result <- paste0(result, item, delimiter)
    }
    result <- paste0(result, arr[length(arr)])
    return(result)
}



#' Actualización correcta de los órdenes del modelo cuando alguno se fija a 0
#' Por ejemplo, cuando el parámetro ma2 de un ARMA(1, 2) se fija a cero, éste 
#' pasa a ser un ARMA(1, 1), por tanto el objeto `orders` y `fixed` deben ser actualizados.
#'
#' @param orders Antiguos órdenes que deben ser actualizados.
#' @param coefs_names Nombres de los coeficientes del modelo (es necesario saber esto 
#' para poder interpretar bien el vector fixed y actualizarlo).
#' @param fixed Máscara para indicar qué coeficientes (`coefs_names`) se van a fijar a 0.
#' 
#' @returns Lista con los órdenes (`orders`) y vector máscara (`fixed`) actualizados.
update.orders <- function(orders, coefs_names, fixed=NULL) {
  
  if (is.null(fixed)) {return(list(orders=orders, fixed=NULL))}
  if (all(is.na(fixed))) {return(list(orders=orders, fixed=fixed))}
  
  
  ar <- str_detect(coefs_names, '^ar\\d+')
  ar_updated <- update.order(orders$regular[1], fixed, ar)
  orders$regular[1] <- ar_updated$order
  ar_fixed <- ar_updated$fixed
  
  
  ma <- str_detect(coefs_names, '^ma\\d+')
  ma_updated <- update.order(orders$regular[3], fixed, ma)
  orders$regular[3] <- ma_updated$order
  ma_fixed <- ma_updated$fixed
  
  
  sar <- str_detect(coefs_names, '^sar\\d+')
  sar_updated <- update.order(orders$seasonal[1], fixed, sar)
  orders$seasonal[1] <- sar_updated$order
  sar_fixed <- sar_updated$fixed
  
  sma <- str_detect(coefs_names, '^sma\\d+')
  sma_updated <- update.order(orders$seasonal[3], fixed, sma)
  orders$seasonal[3] <- sma_updated$order
  sma_fixed <- sma_updated$fixed
  
  
  if (orders$include_mean) {
    mask <- sapply(coefs_names, str_detect, c('mean', 'drift', 'intercept'))
    mask <- apply(mask, 2, sum) == 1
    const_value <- fixed[mask]
    if (length(const_value) == 0) {
        orders$include_mean <- FALSE
    } else if (!is.na(const_value)) {
      orders$include_mean <- FALSE
    }
  }
  
  cnt <- sum(ar, ma, sar, sma)
  
  fixed <- if (cnt>=length(fixed)) c(ar_fixed, ma_fixed, sar_fixed, sma_fixed) else 
    c(ar_fixed, ma_fixed, sar_fixed, sma_fixed, fixed[(1+cnt):length(fixed)])
  return(list(orders=orders, fixed=fixed))
}


update.order <- function(order, fixed, mask) {
  if (sum(mask) == 0 || order==0) {
    return(list(order=order, fixed=c()))
  } 
  
  fixed <- fixed[mask]
  
  if (all(is.na(fixed))) {
    return(list(order=order, fixed=fixed))
  }
  
  for (val in rev(fixed)) {
    if (!is.na(val) && val==0) {
      order <- order - 1
    } else {
      break
    }
  }
  fixed <- if (order==0) c() else fixed[1:order]
  return(list(order=order, fixed=fixed))
}



#' Comprobación de que un ajuste ha sido correctamente optimizado
#' 
#' @param ajuste Ajuste a comprobar.
#' @returns `TRUE` si el modelo es válido, `FALSE` en otro caso.
is_valid <- function(ajuste) {
    
    if (length(ajuste) == 1 && class(ajuste) == 'try-error') {
        return(FALSE)
    }
    else if (all(is.na(ajuste))) {
        return(FALSE)
    }
    else if (suppressWarnings(any(is.na(sqrt(diag(ajuste$var.coef)))))) {
        return(FALSE)
    } else if (suppressWarnings(any(is.na(ajuste$var.coef)))) {
        return(FALSE)
    }
    return(TRUE)
}


#' Obtención del lag óptimo de Ljung-Box
#'
#' @param serie 
#'
#' @return
#' @export
#'
#' @examples
ljungbox_lag <- function(serie) {
  m <- frequency(serie)
  n <- length(serie)
  if (m > 1) {
    return(floor(min(2*m, n/5)))
  } else {
    return(floor(min(10, n/5)))
  }
}




