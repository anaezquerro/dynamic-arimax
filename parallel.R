

auto.fit.arima.regression <- function(serie, xregs, ic='aicc', alpha=0.05, 
                                      stationary_method='auto.arima', show_info=T, 
                                      ndiff=0) {
    
    display_results <- function(ajustes) {
        for (key in names(ajustes)) {
            if (typeof(ajustes[[key]]) == 'character') {
                cat(ajustes[[key]])
            } else {
                cat(paste0('Se ha probado con la variable ', names(xregs)[as.numeric(key)], 
                           ' [ic=', ajustes[[key]][[ic]], ', lag=', ajustes[[key]]$optimal_lag, ']\n'))
            }
        }
    }
    
    
    fit.simple.regression <- function(j) {
        xreg <- xregs[, j]              # variable respuesta
        xreg_name <- names(xregs)[j]    # nombre de la variable respuesta
        
        # Obtención del retardo óptimo entre la variable regresora y la respuesta
        optimal_lag <- select.optimal.lag(response, xreg, alpha=alpha, max_lag=max_lag, method=stationary_method)
        
        if (is.na(optimal_lag)) {
            return(paste0('No se ha podido encontrar un retardo significativo para ', xreg_name, '\n'))
        }
        
        # Creación de la matriz de tipo `ts` con la variables respuesta, variables 
        # regreosras anteriormente añadidas y nueva variable regresora
        data_new <- construct.data(model_history, response, xregs, xreg_name, optimal_lag, max_lag)
        
        # Ajuste del modelo de regresión dinámica con la nueva variable regresora
        ajuste <- auto.fit.arima(data_new[, c(1)], xregs=data_new[, -c(1)], 
                                 d=NA, D=NA, alpha=alpha, show_info=F)
        
        if (!is_valid(ajuste)) {
            return(paste0('No se ha podido ajustar un modelo para ', xreg_name, '\n'))
        }
        
        ajuste$optimal_lag <- optimal_lag
        return(ajuste)
    }
    
    
    # Comprobación de argumentos
    if (class(xregs) != 'data.frame') {
        stop('El argumento `xregs` debe ser un data.frame')
    }
    if (!all(unlist(lapply(xregs, class)) == "ts")) {
        stop('Las variables regresoras del data.frame `xregs` deben ser de tipo ts')
    }
    if (class(serie) != 'ts') {
        stop('El arguemnto `serie` debe ser de tipo ts')
    }
    if (frequency(serie) == 1) { 
        seasonal <- F
    }
    
    # Inicialización de variables que se utilizarán a lo largo del algoritmo
    best_ic <- Inf
    model_history <- data.frame(var=NA, lag=NA, ic=NA)
    
    # Cálculo del retardo óptimo máximo de cada variable regresora sobre la variable respuesta. 
    # Se utilizará posteriormente para recortar las series y poder comparar modelos construidos 
    # con series del mismo tamaño
    max_lag <- get.maximum.lag(serie, xregs, alpha=alpha, method=stationary_method)
    
    # Recortamos todas las series
    response <- serie       # inicialmente, response = serie
    
    # Inicio del bucle para añadir variables regresoras
    for (i in 1:ncol(xregs)) {
        
        # Inicialzación de parámetros
        best_xreg <- NA            # nombre de la variable regresora que se añadirá
        best_xreg_lag <- NA        # retardo de la variable regresora que se añadirá
        
        # Selección de la mejor variable regresora
        xregs_indexes <- (1:ncol(xregs))[!(names(xregs) %in% model_history$var)]
        xregs_list <- as.list(xregs_indexes)
        names(xregs_list) <- xregs_indexes
            
        ajustes <- mclapply(xregs_list, fit.simple.regression, mc.cores=detectCores(logical=F))
        if (show_info) { display_results(ajustes) }
        

        results <- unlist(lapply(ajustes, function(x) ifelse(is.character(x), NA, x[[ic]])))
        if (all(is.na(results))) {    # no se ha podido añadir ninguna variable más
            break
        }
        
        best_xreg <- names(results)[which.min(c(results))]
        
        # Se comprueba que se ha mejorado el IC respecto al mejor obtenido
        if (best_ic < ajustes[[best_xreg]][[ic]]) {
            break
        }
        
        global_fit <- ajustes[[best_xreg]]
        best_ic <- global_fit[[ic]]
        best_xreg <- as.numeric(best_xreg)
        best_xreg_lag <- global_fit$optimal_lag
            
        if (show_info) {
            cat(paste0('Se ha añadido la variable regresora ', names(xregs)[best_xreg], ' [', ic, '=', best_ic, ']\n'))
            print(global_fit, row.names=F)
            cat(paste0(stri_dup('-', PAD), '\n'))
        }
        
        # Se añade al historial del modelo la nueva variable añadida
        model_history[i, ] <- c(names(xregs)[best_xreg], best_xreg_lag, best_ic)
        
        # La "respuesta" pasa a ser los residuos del ajuste
        response <- residuals(global_fit, type='regression')
    }
    
    
    # Si ya no entran más variables (y ha entrado al menos una) se ajusta el 
    # modelo completo donde las innovaciones sigan un proceso ARMA
    if (!any(is.na(model_history)) && (sum(global_fit$arma[6:7]) > 0)) {
        data_new <- construct.data(model_history, serie, xregs, new_xreg_name = NULL, 
                                   optimal_lag = NULL)
        global_fit <- auto.fit.arima(data_new[, c(1)], xregs=data_new[, -c(1)],
                                     ic=ic, d=0, D=0, alpha=alpha, show_info=F)    
    }
        
    # Si se consigue un modelo válido, se devuelve
    if (is_valid(global_fit)) {
        best_ic <- global_fit[[ic]]
        
        if (show_info) {
            cat('No se añaden más variables\n\n')
        }
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
    
    
    # Si el model_history está vacío  (es decir, no hay variables regresoras que 
    # aporten algo a la variable respuesta), o no se ha podido ajustar uno con 
    # errores ARMA, se trata de de:
    # a) Diferenciar todas las variables y volver a llamar a la función auto.fit.arima.regression
    # b) Si ndiff=2 no se recomienda diferenciar más y se ajusta un ARIMA sin regresoras
    if (ndiff < 3) {
        cat(paste0('No se ha podido encontrar un modelo válido con errores estacionarios\n',
                   'Se aplica una diferenciación regular (ndiff=', ndiff+1, 
                   ') y se vuelve a llamar a la función\n'))
        cat(paste0(stri_dup('-', PAD), '\n'))
        
        serie <- diff(serie)
        xregs <- as.data.frame(lapply(xregs, diff))
        return(
            auto.fit.arima.regression(serie, xregs, ic, alpha, stationary_method, 
                                      show_info, ndiff+1)
        )
    }
    
    # En otro caso se intenta ajustar un ARIMA sin variables regresoras
    global_fit <- auto.fit.arima(serie, seasonal=seasonal, ic=ic, d=NA, D=NA, alpha=alpha, show_info=F)
    
    # En el caso de que tampoco haya un ARIMA válido sin regresoras
    if (!is_valid(global_fit)) {
        warning('Ningún modelo es válido')
        return(NA)
    } 
    
    # En caso de que sí lo haya, se devuelve éste
    best_ic <- global_fit[[ic]]
    
    if (show_info) {
        cat(paste0(stri_dup('-', PAD), '\n'))
        cat(paste0('Modelo (sin variables regresoras cond ndiff=', ndiff, ') [', ic, '=', best_ic, ']\n'))
        cat(paste0(stri_dup('-', PAD), '\n'))
        print(global_fit, row.names=F)
        cat(paste0(stri_dup('-', PAD), '\n'))
    }
    global_fit$ndiff <- ndiff
    return(global_fit)
}




