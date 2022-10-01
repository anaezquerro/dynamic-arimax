

#' En este archivo se encuentran las funciones iniciales que se implementaron 
#' en una primera aproximación secuencial (ya no se usan en los notebooks). 
#' Se han mantenido para realizar comparaciones en tiempos de ejecución con la 
#' versión final paralelizada.


auto.fit.arima.regression <- function(serie, xregs, ic='aicc', alpha=0.05, 
                                      stationary_method='auto.arima', show_info=T, 
                                      ndiff=0) {
    
    # Comprobación de argumentos
    if (!is.null(xregs) && !any(c('mts', 'ts') %in% class(xregs))) {
        stop('El argumento `xregs` debe ser de tipo mts o ts')
    }
    if (class(serie) != 'ts') {
        stop('El arguemnto `serie` debe ser de tipo ts')
    }
    
    # Inicialización de variables que se utilizarán a lo largo del algoritmo
    global_ic <- Inf
    model_history <- data.frame(var=NA, lag=NA, ic=NA)
    
    max_lag <- get.maximum.lag(serie, xregs, alpha=alpha, method=stationary_method)
    response <- serie

        # Inicio del bucle para añadir variables regresoras
    for (i in 1:ncol(xregs)) {
        if (max_lag == -Inf) {break}
        
        # Inicialzación de parámetros
        best_xreg <- NA            # nombre de la variable regresora que se añadirá
        best_xreg_lag <- NA        # retardo de la variable regresora que se añadirá
        added <- FALSE             # indicación de si se ha añadido una regresora más
        
        # Inicialización del bucle anidado para seleccionar la "mejor" regresora
        
        for (j in 1:ncol(xregs)) {
            
            xreg <- xregs[, j]              # variable que se está evaluando
            xreg_name <- colnames(xregs)[j]    # nombre de la variable que se está evaluando
            
            # Si la variable ya está registrada en el historial del modelo, la ignoramos
            if (xreg_name %in% model_history$var) {
                if (show_info) {
                    cat(paste0('Saltamos ', xreg_name, '\n'))
                }
                next
            }
            
            # Cálculo del retardo óptimo entre la variable regresora y la respuesta
            optimal_lag <- select.optimal.lag(xreg, response, alpha=0.05, max_lag=max_lag, method=stationary_method)
            
            
            # Si no hay retardo significativo, se pasa a evaluar la siguiente variable
            if (is.na(optimal_lag)) {
                if (show_info) {
                    cat(paste0('No se ha podido encontrar un retardo significativo para ', xreg_name, '\n'))}
                next
            }
            
            data_new <- construct.data(model_history, serie, xregs, xreg_name, optimal_lag, max_lag)
            ajuste <- auto.fit.arima(data_new[, c(1)], xregs=data_new[, -c(1)], 
                                     ic=ic, d=NA, D=NA, alpha=alpha, show_info=F)
            
            # Si no se consigue un ajuste válido, se prueba con la siguiente variable
            if (all(is.na(ajuste))) {   
                if (show_info) {
                    cat(paste0('No se ha podido ajustar un modelo para ', xreg_name, '\n'))
                }
                next
            }
            
            if (show_info) {
                cat(paste0('Se ha probado con la variable ', xreg_name, ' [ic=', ajuste[[ic]], ', lag=', optimal_lag, ']\n'))
            }
            
            # Si se consigue un modelo válido, se comprueba si se ha mejorado el ic con el nuevo ajuste
            if (ajuste[[ic]] < global_ic) {
                global_fit <- ajuste
                global_ic <- ajuste[[ic]]
                best_xreg <- j
                best_xreg_lag <- optimal_lag
                added <- TRUE
            }
        }
        
        # Comprobación: Se tienen que haber añadido alguna variable al modelo
        if (!added) {
            break
        }
        
        # Se anota la variable regresora que se ha añadido al modelo
        if (show_info) {
            cat(paste0('Se ha añadido la variable regresora ', colnames(xregs)[best_xreg], ' [', ic, '=', global_ic, ']\n'))
            print(global_fit, row.names=F)
            cat(paste0(stri_dup('-', PAD), '\n'))
        }
        
        # Se añade al historial del modelo la nueva variable añadida
        model_history[i, ] <- c(colnames(xregs)[best_xreg], best_xreg_lag, global_ic)
        
        # La "respuesta" pasa a ser los residuos del ajuste
        response <- residuals(global_fit, type='regression')
    }
    
    if (show_info) { cat('No se añaden más variables\n') }
    
    
    # Si ya no entran más variables (y ha entrado al menos una) se ajusta el 
    # modelo completo donde las innovaciones sigan un proceso ARMA
    if (!any(is.na(model_history)) && (sum(global_fit$arma[6:7]) > 0)) {
        
        if (show_info) {
            cat(paste0('El modelo global no tiene errores estacionarios\n', 
                       'Se intenta ajustar uno que sí los tenga\n'))
        }
        data_new <- construct.data(model_history, serie, xregs, new_xreg_name = NULL, 
                                   optimal_lag = NULL, max_lag=max_lag)
        global_fit <- auto.fit.arima(data_new[, c(1)], xregs=data_new[, -c(1)],
                                     ic=ic, d=0, D=0, alpha=alpha, show_info=F)    
    }
    
    # Si se consigue un modelo válido, se devuelve
    if (is_valid(global_fit)) {
        global_ic <- global_fit[[ic]]
        
        if (ndiff > 1) {
            warning(paste0('Se han aplicado ', ndiff, ' diferencias regulares'))
        }
        
        cat(paste0(stri_dup('-', PAD), '\n'))
        cat(paste0('|', 
                   str_pad(paste0('Histórico de variables añadidas al modelo (ndiff=', ndiff, ')'), width=PAD-2, 
                           side='both', pad=' '), '|\n'))
        cat(paste0(stri_dup('-', PAD), '\n'))
        print(model_history, row.names=F)
        cat(paste0(stri_dup('-', PAD), '\n'))
        print(global_fit, row.names=F)
        global_fit$ndiff <- ndiff
        global_fit$history <- model_history
        return(global_fit)
    }
    
    
    if (ndiff < 3) {
        if (show_info) {
            cat(paste0('No se ha podido encontrar un modelo válido con errores estacionarios\n',
                       'Se aplica una diferenciación regular (ndiff=', ndiff+1, 
                       ') y se vuelve a llamar a la función\n'))            
        }
        
        serie <- diff(serie)
        xregs <- diff(xregs)
        return(
            auto.fit.arima.regression(serie, xregs, ic, alpha, stationary_method, 
                                      show_info, ndiff+1)
        )
    }
    
    # En otro caso se intenta ajustar un ARIMA sin variables regresoras
    global_fit <- auto.fit.arima(serie, ic=ic, d=NA, D=NA, alpha=alpha, show_info=F)
    
    # En el caso de que tampoco haya un ARIMA válido sin regresoras
    if (!is_valid(global_fit)) {
        warning('Ningún modelo es válido')
        return(NA)
    } 
    
    # En caso de que sí lo haya, se devuelve éste
    global_ic <- global_fit[[ic]]
    if (show_info) {
        cat(paste0(stri_dup('-', PAD), '\n'))
        cat(paste0('Modelo (sin variables regresoras cond ndiff=', ndiff, ') [', ic, '=', global_ic, ']\n'))
        cat(paste0(stri_dup('-', PAD), '\n'))
        print(global_fit, row.names=F)
        cat(paste0(stri_dup('-', PAD), '\n'))
    }
    global_fit$ndiff <- ndiff
    return(global_fit)
}


get.maximum.lag <- function(serie, xregs, alpha=0.05, method='auto.arima') {
    
    max_lag <- -Inf
    
    for (i in 1:ncol(xregs)) {
        xreg <- xregs[, i]
        optimal_lag <- select.optimal.lag(xreg, serie, alpha=alpha, method=method, less0=F)
        if (is.na(optimal_lag)) {
            next
        }
        if (abs(optimal_lag) > max_lag) {
            max_lag <- abs(optimal_lag)
        }
    }
    return(max_lag)
}





