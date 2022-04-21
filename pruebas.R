eval(parse("automatic_selection.R", encoding="UTF-8"))
eval(parse("arima_simulation.R", encoding="UTF-8"))

auto.fit.arima.regression <- function(serie, xregs, seasonal=T, ic='aicc', alpha=0.05, 
                                      stationary_method='auto.arima', show_info=T, 
                                      ndiff=0) {
    
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
    
    # response: serie sobre la que se calcula el retardo óptimo en cada iteración (residuos)
    response <- serie       # inicialmente, response = serie
    
    # Cálculo del retardo óptimo máximo de cada variable regresora sobre la variable respuesta. 
    # Se utilizará posteriormente para recortar las series y poder comparar modelos construidos 
    # con series del mismo tamaño
    max_lag <- get.maximum.lag(serie, xregs, alpha=alpha, method=stationary_method)
    
    # Inicio del bucle para añadir variables regresoras
    for (i in 1:ncol(xregs)) {
        if (show_info) {
            cat('----------------------------------------------------------------------\n')
            cat(paste0('Iteración ', i, ' del algoritmo\n'))
            cat('----------------------------------------------------------------------\n')
        }
        
        # Inicialzación de parámetros
        best_xreg <- NA            # nombre de la variable regresora que se añadirá
        best_xreg_lag <- NA        # retardo de la variable regresora que se añadirá
        added <- FALSE             # indicación de si se ha añadido una regresora más
        
        # Inicialización del bucle anidado para seleccionar la "mejor" regresora
        for (j in 1:ncol(xregs)) {
            
            xreg <- xregs[, j]              # variable que se está evaluando
            xreg_name <- names(xregs)[j]    # nombre de la variable que se está evaluando
            
            # Si la variable ya está registrada en el historial del modelo, la ignoramos
            if (xreg_name %in% model_history$var) {
                if (show_info) {
                    cat(paste0('Saltamos ', xreg_name, '\n'))
                }
                next
            }
            
            # Cálculo del retardo óptimo entre la variable regresora y la respuesta
            optimal_lag <- select.optimal.lag(response, xreg, max_lag=max_lag, method=stationary_method)
            
            # Si no hay retardo significativo, se pasa a evaluar la siguiente variable
            if (is.na(optimal_lag)) {
                if (show_info) {
                    cat(paste0('No se ha podido encontrar un retardo significativo para ', xreg_name, '\n'))}
                next
            }
            
            # Creación de la matriz de tipo `ts` con la variable respuesta, variables 
            # regresoras anteriormente añadidas y nueva variable regresora
            data_new <- construct.data(model_history, serie, xregs, xreg_name, optimal_lag, max_lag)
            
            # Ajuste del modelo de regresión dinámica con la nueva variable regresora
            ajuste <- auto.fit.arima(data_new[, c(1)], xregs=data_new[, -c(1)], 
                                     seasonal=seasonal, ic=ic, d=NA, D=NA, alpha=alpha, show_info=F)
            
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
            if (ajuste[[ic]] < best_ic) {
                global_fit <- ajuste
                best_ic <- ajuste[[ic]]
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
            cat('----------------------------------------------------------------------\n')
            cat(paste0('Se ha añadido la variable regresora ', names(xregs)[best_xreg], ' [', ic, '=', best_ic, ']\n'))
            cat('----------------------------------------------------------------------\n')
            print(global_fit, row.names=F)
        }
        
        # Se añade al historial del modelo la nueva variable añadida
        model_history[i, ] <- c(names(xregs)[best_xreg], best_xreg_lag, best_ic)
        
        # La "respuesta" pasa a ser los residuos del ajuste
        response <- residuals(global_fit, type='regression')
    }
    
    
    # Si ya no entran más variables (y ha entrado al menos una) se ajusta el 
    # modelo completo donde las innovaciones sigan un proceso ARMA
    if (!any(is.na(model_history))) {
        data_new <- construct.data(model_history, serie, xregs, new_xreg_name = NULL, 
                                   optimal_lag = NULL, max_lag = max_lag)
        global_fit <- auto.fit.arima(data_new[, c(1)], xregs=data_new[, -c(1)],
                                     seasonal=seasonal, ic=ic, d=0, D=0, alpha=alpha, show_info=F)    
        
        # Si se consigue un modelo válido, se devuelve
        if (is_valid(global_fit)) {
            best_ic <- global_fit[[ic]]
            
            if (show_info) {
                cat('No se añaden más variables\n\n')
            }
            
            cat('----------------------------------------------------------------------\n')
            cat('|            Histórico de variables añadidas al modelo                |\n')
            cat('----------------------------------------------------------------------\n')
            print(model_history, row.names=F)
            return(global_fit)
        }
    }
   
    
    # Si el model_history está vacío  (es decir, no hay variables regresoras que 
    # aporten algo a la variable respuesta), o no se ha podido ajustar uno con 
    # errores ARMA, se trata de de:
    # a) Diferenciar todas las variables y volver a llamar a la función auto.fit.arima.regression
    # b) Si ndiff=2 no se recomienda diferenciar más y se ajusta un ARIMA sin regresoras
    if (ndiff < 3) {
        warning(paste0('No se ha podido encontrar un modelo válido\n',
                       'Se aplica una diferenciación regular (ndiff=', ndiff+1, 
                       ') y se vuelve a llamar a la función\n'))
        
        serie <- diff(serie)
        xregs <- as.data.frame(lapply(xregs, diff))
        return(
            auto.fit.arima.regression(serie, xregs, seasonal, ic, alpha, stationary_method, 
                                      show_info, ndiff+1)
        )
    }
    
    # En otro caso se intenta ajustar un ARIMA sin variables regresoras
    global_fit <- auto.fit.arima(serie, seasonal=seasonal, ic=ic, d=NA, D=NA, alpha=alpha, show_info=F)

    # En el caso de que tampoco haya un ARIMA válido sin regresoras
    if (!is_valid(global_fit)) {
        stop('Ningún modelo es válido')
    } 
    
    # En caso de que sí lo haya, se devuelve éste
    best_ic <- global_fit[[ic]]
    
    if (show_info) {
        cat('----------------------------------------------------------------------\n')
        cat(paste0('Modelo (sin variables regresoras cond d=', ndiff, ') [', ic, '=', best_ic, ']\n'))
        cat('----------------------------------------------------------------------\n')
        print(global_fit, row.names=F)
        cat('----------------------------------------------------------------------\n')
    }
    return(global_fit)
}