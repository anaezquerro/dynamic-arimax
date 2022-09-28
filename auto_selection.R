PAD <- 86



#' Selección automática de variables regresoras en un modelo de regresión dinámica
#' 
#' A partir del proceso de selección de retardos óptimos para cada variable 
#' regresora y el ajuste automático de modelos ARIMAX, se añaden de forma 
#' automática aquellas variables de un conjunto dado que inciden en la variable 
#' respuesta siguiendo un procedimiento iterativo en base a @cryer2008time. Como 
#' los errores del modelo deben ser estacionarios, en caso de que no sea posible 
#' optimizar un modelo con este tipo de errores, se procede a diferenciar todas 
#' los datos y repetir el proceso de selección de variables regresoras.
#'
#' @param serie [`ts`] Serie de tiempo que funciona como variable respuesta.
#' @param xregs [`mts`] Conjunto de variables candidatas a ser variables 
#' regresoras de `serie`.
#' @param ic [`character`] Criterio de información sobre el que se seleccionarán 
#' variables regresoras  y se comparará la efectividad de un ajuste. Puede ser 
#' el AIC (`"aic"`), BIC (`"bic"`) o AICc (`"aicc"`). Por defecto se usa el AICc.
#' @param alpha [`numeric`] Nivel de significación para los contrastes de 
#' hipótesis de estacionariedad y significancia. Por defecto es igual a $0.05$.
#' @param stationary_method [`character`] Método utilizado para chequear la 
#' estacionariedad de una serie de tiempo. Si es `"auto.arima"`, se utiliza la 
#' función `forecast::auto.arima()` para ajustar un modelo ARIMA(p,d,q) y 
#' chequear si $d > 0$ (si se cumple esta condición se asume que la serie no es 
#' estacionaria). Si es `"adf.test"` se  usa el test Dickey-Fuller 
#' (`tseries:.adf.test()`).
#' @param show_info Valor booelano que indica si mostrar o no por pantalla los 
#' modelos que se van ajustando, así como información relacionada con la 
#' selección de variables.
#' @param ndiff Parámetro interno del programa (*no utilizar*) para diferenciar 
#' todas las variables cuando no se pueda ajustar un modelo válido y mantener un 
#' registro del número de diferenciaciones que se están realizando.
#'
#' @returns Ajuste del modelo ARIMAX con las variables regresoras seleccionadas. 
#' Si dicho ajuste no existe (o no es posible de optimizar), devuelve el ajuste 
#' válido de un modelo ARIMA sobre la variable respuesta. Si éste último ajuste 
#' tampoco existe (no es posible de optimizar), devuelve `NA`. 
#' 
auto.fit.arima.regression <- function(serie, xregs, ic='aicc', alpha=0.05, 
                                      stationary_method='auto.arima', 
                                      show_info=T, ndiff=0) {
    
    
    # Función auxiliar para imprimir los resultados obtenidos en la función 
    # parLapply() ejecutada en el cluster de procesadores
    display_results <- function(ajustes) {
        for (key in names(ajustes)) {
            if (typeof(ajustes[[key]]) == 'character') {
                cat(ajustes[[key]])
            } else {
                cat(paste0('Covariate ', colnames(xregs)[as.numeric(key)], 
                           ' has been tested [ic=', ajustes[[key]][[ic]], 
                           ', lag=', ajustes[[key]]$optimal_lag, ']\n'))
            }
        }
    }
    
    
    # Obtención del retardo máximo (en valor absoluto) sobre el conjunto de todas 
    # las variables regresoras de forma paralelizada
    get.maximum.lag <- function() {
        
        optimal_lags <- parLapply(cl, as.list(1:ncol(xregs)), 
            function(x) select.optimal.lag(xregs[, x], response, alpha, 
                             method=stationary_method, less0=F))
        optimal_lags <- unlist(optimal_lags[!is.na(optimal_lags)])
        
        if (is.null(optimal_lags)) {
            return(-Inf)
        }
        
        return(max(abs(optimal_lags)))
    }
    
    
    # Ajuste de un modelo de regresión añadiendo la variable j al modelo global
    fit.simple.regression <- function(j) {
        xreg <- xregs[, j]                 # variable respuesta
        xreg_name <- colnames(xregs)[j]    # nombre de la variable respuesta
        
        # Obtención del retardo óptimo entre la variable regresora y la respuesta
        optimal_lag <- select.optimal.lag(xreg, response, 
            alpha=alpha, max_lag=max_lag, method=stationary_method)
        
        if (is.na(optimal_lag)) {
            return(paste0('Significative correlation with lag<=0 could not be found for ', xreg_name, '\n'))
        }
        
        # Creación de la matriz de tipo `ts` con la variables respuesta, variables 
        # regresoras anteriormente añadidas y nueva variable regresora
        data_new <- construct.data(model_history, serie, xregs, xreg_name, optimal_lag, max_lag)
        
        # Ajuste del modelo de regresión dinámica con la nueva variable regresora
        ajuste <- auto.fit.arima(data_new[, c(1)], xregs=data_new[, -c(1)], 
                                 d=NA, D=NA, alpha=alpha, show_info=F)
        
        if (!is_valid(ajuste)) {
            return(paste0('The optimizer could not fit a model for ', xreg_name, '\n'))
        }
        
        if (! 
            ((xreg_name %in% names(ajuste$coef)) || 
             ('xreg' %in% names(ajuste$coef)))
            ) {
            return(paste0('The optimizer could not add ', xreg_name, ' to the model\n'))
        }
        
        ajuste$optimal_lag <- optimal_lag
        return(ajuste)
    }
    
    
    # Comprobación de argumentos
    if (!any(c('mts', 'ts') %in% class(xregs))) {
        stop('El argumento `xregs` debe ser de tipo mts o ts')
    }
    if (class(serie) != 'ts') {
        stop('El arguemnto `serie` debe ser de tipo ts')
    }
    
    # Inicialización de variables que se utilizarán a lo largo del algoritmo
    global_ic <- Inf
    global_fit <- NA
    model_history <- data.frame(var=NA, lag=NA, ic=NA)
    
    # Creación del cluster de procesos para paralelizar
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
    
    # Cálculo del retardo óptimo máximo de cada variable regresora sobre la 
    # variable respuesta. Se utilizará posteriormente para recortar las series y 
    # poder comparar modelos construidos con series del mismo tamaño
    response <- serie               # inicialmente, response = serie
    max_lag <- get.maximum.lag()
    
    # Inicio del bucle para añadir variables regresoras
    for (i in 1:ncol(xregs)) {
        
        # Consideramos el caso en el que no se ha conseguido obtener ningún 
        # retardo en el que se produzca una correlación significativa para 
        # ninguna de las variables
        if (max_lag == -Inf) {
            break
        }
        
        # Selección de la "mejor" variable regresora
        xregs_indexes <- (1:ncol(xregs))[!(colnames(xregs) %in% model_history$var)]
        xregs_list <- as.list(xregs_indexes)
        names(xregs_list) <- xregs_indexes
        ajustes <- parLapply(cl, xregs_list, fit.simple.regression)
        
        if (show_info) { display_results(ajustes) }
        
        # Obtenemos los ICs de cada ajuste
        results <- unlist(lapply(ajustes, function(x) ifelse(is.character(x), NA, x[[ic]])))
        
        if (all(is.na(results))) {# no se ha podido añadir ninguna variable más
            break
        }
        
        # La variable que se introduce al modelo es aquella que minimiza el IC
        best_xreg <- names(results)[which.min(c(results))]
        
        # Si no se ha mejorado el IC con respecto al IC del modelo global, se 
        # detiene el algoritmo
        if (global_ic < ajustes[[best_xreg]][[ic]]) {
            break
        }
        
        # Actualización de las variables globales
        global_fit <- ajustes[[best_xreg]]
        global_ic <- global_fit[[ic]]
        best_xreg <- as.numeric(best_xreg)
        best_xreg_lag <- global_fit$optimal_lag
        
        if (show_info) {
            cat(paste0('Covariate ', colnames(xregs)[best_xreg], 
                       ' has been added [', ic, '=', global_ic, ', lag=', best_xreg_lag, ']\n'))
            print(global_fit, row.names=F)
            cat(paste0(stri_dup('-', PAD), '\n'))
        }
        
        # Se añade al historial del modelo la nueva variable añadida
        model_history[i, ] <- c(colnames(xregs)[best_xreg], best_xreg_lag, global_ic)
        
        # La "respuesta" pasa a ser los residuos del ajuste
        response <- residuals(global_fit, type='regression')
    }
    
    if (show_info) { cat('No more variables will be added\n') }
    
    
    # Una vez finalizado el algoritmo de selección se puede detener el cluster 
    # de procesos
    stopCluster(cl)
    
    
    # Si ya no entran más variables (y ha entrado al menos una) se ajusta el 
    # modelo completo donde los errores sigan un proceso ARMA
    if (!any(is.na(model_history)) && (sum(global_fit$arma[6:7]) > 0)) {
        
        if (show_info) {
            cat(paste0('The global model does not have stationary errors\n', 
                       'Trying to adjust a model that do have stationary errors\n'))
        }
        
        data_new <- construct.data(model_history, serie, xregs, new_xreg_name = NULL, 
                                   optimal_lag = NULL, max_lag=max_lag)
        global_fit <- auto.fit.arima(data_new[, c(1)], xregs=data_new[, -c(1)],
                                     ic=ic, d=0, D=0, alpha=alpha, show_info=F)    
    }
    
    # Si se consigue un modelo válido, se devuelve
    if (is_valid(global_fit)) {
        global_ic <- global_fit[[ic]]
        
        # Se imprime por pantalla los resultados
        if (show_info) {
            cat(paste0(stri_dup('-', PAD), '\n'))
            cat(paste0('|', 
                       str_pad(paste0('Historical of added covariates to the model (ndiff=', ndiff, ')'), width=PAD-2, 
                               side='both', pad=' '), '|\n'))
            cat(paste0(stri_dup('-', PAD), '\n'))
            print(model_history, row.names=F)
            cat(paste0(stri_dup('-', PAD), '\n'))
            print(global_fit, row.names=F)
        }
        
        # Se devuelve junto con el modelo el número de diferenciaciones 
        # aplicadas a los datos y el historial de las variables regresoras 
        # añadidads
        global_fit$ndiff <- ndiff
        global_fit$history <- model_history
        
        return(global_fit)
    }
    
    
    # Si el model_history está vacío  (es decir, no hay variables regresoras que 
    # aporten algo a la variable respuesta), o no se ha podido ajustar uno con 
    # errores ARMA, se trata de de:
    # a) Diferenciar todas las variables y volver a llamar a la función.
    # b) Si ndiff=2 no se recomienda diferenciar más y se ajusta un ARIMA sin regresoras
    if (ndiff < 3) {
        
        if (show_info) {
            cat(paste0('No valid model with stationary errors could be optimized\n',
                       'Applying regular differentiation (ndiff=', ndiff+1, 
                       ') and calling again the function\n'))
            cat(paste0(stri_dup('-', PAD), '\n'))
            cat(paste0(stri_dup('-', PAD), '\n'))
            
        }
        
        # Se diferencian todos los datos
        serie <- diff(serie)
        xregs <- diff(xregs)
        return(
            auto.fit.arima.regression(serie, xregs, ic, alpha, stationary_method, 
                                      show_info, ndiff+1)
        )
    }
    
    # En otro caso se intenta ajustar un ARIMA sin variables regresoras
    global_fit <- auto.fit.arima(serie, seasonal=seasonal, ic=ic, d=NA, D=NA, alpha=alpha, show_info=F)
    
    # En el caso de que tampoco haya un ARIMA válido sin regresoras
    if (!is_valid(global_fit)) {
        warning('No valid model could be optimized')
        return(NA)
    } 
    
    # En caso de que sí lo haya, se devuelve éste
    global_ic <- global_fit[[ic]]
    
    if (show_info) {
        cat(paste0(stri_dup('-', PAD), '\n'))
        cat(paste0('Model without covariates  (ndiff=', ndiff, ') [', ic, '=', global_ic, ']\n'))
        cat(paste0(stri_dup('-', PAD), '\n'))
        print(global_fit, row.names=F)
        cat(paste0(stri_dup('-', PAD), '\n'))
    }
    
    global_fit$ndiff <- ndiff
    
    return(global_fit)
}



#' Función auxiliar para la construcción de una matriz de series temporales en 
#' base al histórico del modelo y una nueva variable a añadir con su respectivo 
#' retardo óptimo y retardo máximo del conjunto total. La matriz resultante 
#' debe satisfacer que:
#' 
#' - Su primera columna sea la variable respuesta.
#' - Las otras columnas se correspondan con las variables regresoras *en orden* 
#' del modelo ya retardadas por su respectivo retardo.
#' - A todas las variables se les haya recortado los primeros `max_lag` valores 
#' (para tener modelos construidos con el mismo número de muestras).
#' 
#' @param model_history [`data.frame`] Histórico de las variables que ya se han 
#' añadido al modelo (nombre de la variable, su retardo y su valor en el 
#' criterio de información seleccionado).
#' @param serie [`ts`] Serie temporal que funciona como variable respuesta.
#' @param xregs [`mts`] Conjunto de variables candidatas a ser variables 
#' regresoras de `serie`.
#' @param new_xreg_name [`character`] Nombre de la nueva variable regresora que 
#' se quiere añadir al modelo. Si es `NULL` se asume que no se quiere añadir 
#' ninguna variable nueva y la matriz se contruye en base a las que ya están en 
#' el historial.
#' @param optimal_lag [`numeric` ] Retardo óptimo de la nueva variable a incluir. 
#' Puede ser `NULL` en caso de que `new_xreg_name` también lo sea.
#' @param max_lag [`numeric`] Retardo máximo de todos los retardos óptimos (en 
#' valor absoluto) del conjunto de variables regresoras.
#'
#' @returns Matriz con las variables retardadas y recortadas de forma adecuada 
#' para introducir a la función de ajuste automático de ARIMAXs.
#'
construct.data <- function(model_history, serie, xregs, new_xreg_name, 
                           optimal_lag, max_lag) {
    
    # Primero recortamos los primeros max_lag valores de la variable respuesta 
    serie <- window(serie, start=add_t(start(serie), max_lag, frequency(serie)))
    data <- matrix(serie)
    
    # Si hay más variables en el historial del modelo, las añadimos
    if (!all(is.na(model_history))) {
        
        for (j in 1:nrow(model_history)) {
            xreg <- xregs[, c(model_history$var[j])]                    # variable regresora
            xreg <- stats::lag(xreg, as.integer(model_history$lag[j]))  # variable regresora retardada
            xreg <- window(xreg, start=start(serie), end=end(serie))    # variable recortada
            
            data <- cbind(data, xreg)
            
        }
        colnames(data) <- c('serie', model_history$var)
    } 
    
    # Se añade la nueva variable regresora en caso de que sea necesario
    if (!is.null(new_xreg_name)) {
        
        new_xreg <- stats::lag(xregs[, c(new_xreg_name)], optimal_lag)           # nueva variable retardada
        new_xreg <- window(new_xreg, start=start(serie), end=end(serie))  # nueva variable recortada
        data <- cbind(data, new_xreg)
        
        # Actualización de el encabezado de la matriz
        if (!all(is.na(model_history))) {
            colnames(data) <- c('serie', model_history$var, new_xreg_name)
        } else { colnames(data) <- c('serie', new_xreg_name) }
    }
    
    return(data)
}


#' Función auxiliar para añadir $t$ tiempos a un momento dado de la clase ts.
#' 
#' @param moment [`vector`] Vector de dos valores que dan información acerca del 
#' momento actual al que se le quiere añadir $t$ momentos. El significado de 
#' cada valor es el mismo que el que tienen los índices de la clase ts.
#' @param t [`numeric`] Número de momentos que se le quieren añadir al momento.
#' actual.
#' @param freq [`freq`] Frecuencia estacional del tiempo.
#' 
#' @examples 
#' 
#' - add_t(c(6, 3), 6, 7): c(7,2)
#' - add_t(c(6, 1), 6, 1): c(12, 1)
#' 
add_t <- function(moment, t, freq=1) {
    if (freq == 1) {
        return(c(moment[1] + t, moment[2]))
    } else {
        raw_moment <- moment[2] + t
        return(c(moment[1] + floor(raw_moment/freq), raw_moment %% freq))
    }
}


#' Función auxiliar para restar $t$ tiempos a un momento dado en la clase ts.
#' 
#' @param moment [`vector`] Vector de dos valores que dan información acerca del 
#' momento actual al que se le quiere restar $t$ momentos. El significado de 
#' cada valor es el mismo que el que tienen los índices de la clase ts.
#' @param t [`numeric`] Número de momentos que se le quieren restar al momento.
#' actual.
#' @param freq [`freq`] Frecuencia estacional del tiempo.
#' 
#' @examples 
#' 
#' - substract_t(c(6, 3), 6, 7): c(5, 4)
#' - substract_t(c(6, 1), 6, 1): c(0, 1)
#' 
substract_t <- function(moment, t, freq=1) {
    if (length(moment) != 2) {stop('El momento introducido no es correcto')}
    
    if (freq == 1) {
        return(c(moment[1]-t, moment[2]))
    } else {
        raw_moment <- (moment[1]*freq + moment[2]) - t
        return(c(floor(raw_moment/freq), raw_moment %% freq))
    }
}


#' Comprobación de la tendencia de una serie temporal
#' 
#' @param x [`ts`] Serie temporal de la que se quiere chequear la tendencia.
#' @param test [`character`] Método usado para chequear la tendencia Si es 
#' `"auto.arima"` se intenta ajustar a la serie `x` un modelo ARIMA. Si el orden 
#' de diferenciación $d$ es mayor a 0, se considera que la serie tiene tendencia. 
#' Si es `"adf.test"` se utiliza el test Dickey-Fuller.
#' @param alpha [`numeric`] Nivel de significación para el contraste de 
#' hipótesis de Dickey-Fuller. Por defecto $0.05$..
#'
#' @returns `True` si la serie tiene tendencia, `False` en otro caso.
#'
has_trend <- function(x, test='auto.arima', alpha=0.05) {
    if (!test %in% c('auto.arima', 'adf.test')) {
        stop('El parámetro `test` debe ser "auto.arima" o "adf.test"')
    }
    
    if (test == 'auto.arima') {
        ajuste <- suppressWarnings(auto.arima(x, max.d=3))
        d <- ajuste$arma[6]
        return(d>0)
    }
    if (test == 'adf.test') {
        pvalue <- suppressWarnings(adf.test(x[!is.na(x)])$p.value)
        return(alpha <= pvalue)
    }
}



#' Selección del retardo óptimo (siempre menor o igual a 0) para una serie 
#' temporal y su variable regresora. De forma opcional se puede fijar una 
#' magnitud máxima (`max_lag`) del retardo a escoger.
#' @param xreg [`ts`] Variable regresora.
#' @param serie [`ts`] Variable respuesta.
#' @param alpha [`numeric`] Nivel de significación para los contrastes de 
#' hipótesis de tendencia, estacionalidad y significación.
#' @param max_lag [`numeric`] Magnitud máxima del retardo a escoger (filtrado de 
#' retardos).
#' @param method [`character`] Método para chequear estacionariedad de las 
#' series.
#' 
#' @returns Retardo óptimo para la variable respuesta y su variable regresora 
#' basada en el método propuesto por @cryer2008time. Si no es posible obtener 
#' este retardo, se devuelve `NA`.
#' @export
#'
select.optimal.lag <- function(xreg, serie, alpha=0.05, max_lag=NA, method='auto.arima', less0=T) {
    if (!method %in% c('auto.arima', 'adf.test')) {
        stop('El parámetro `method` debe ser "auto.arima" o "adf.test"')
    }
    
    # Chequeamos si alguna de las variables tiene tendencia
    y_trend <- has_trend(serie, test=method)
    x_trend <- has_trend(xreg, test=method)
    
    # Si alguna de ellas la tiene, aplicamos a ambas una diferenciación regular
    while (y_trend || x_trend) {
        serie <- diff(serie)
        xreg <- diff(xreg)
        y_trend <- has_trend(serie, test=method, alpha=alpha)
        x_trend <- has_trend(xreg, test=method, alpha=alpha)
    }
    
    # Chequeamos si las series tienen componente estacional (nota, frequency > 1)
    y_seasonal <- frequency(serie) > 1 && isSeasonal(serie)
    x_seasonal <- frequency(xreg) > 1 && isSeasonal(xreg)
    
    # Si alguna de ellas la tiene, aplicamos a ambas una diferenciación estacional
    while (y_seasonal || x_seasonal) {
        serie <- diff(serie, lag=frequency(serie))
        xreg <- diff(xreg, lag=frequency(xreg))
        
        y_seasonal <- isSeasonal(serie)
        x_seasonal <- isSeasonal(xreg)
    }
    
    # Realizamos el preblanqueo de las serie con su variable regresora
    if (any(is.na(xreg))) {
        print(head(xreg))
    }
    series_prewhiten <- TSA::prewhiten(xreg, serie, plot=F)
    
    # Obtenemos del objeto resultante las correlaciones y los retardos menores o iguales que 0
    series_lags <- series_prewhiten$ccf$lag
    series_acfs <- series_prewhiten$ccf$acf
    if (less0) {
        series_acfs <- series_acfs[series_lags <= 0]
        series_lags <- series_lags[series_lags <= 0]
    }
    
    # Si tenemos un límite de magnitud para el retardo, filtramos aquellos que tengan 
    # magnitud mayor
    if (!is.na(max_lag)) {
        series_acfs <- series_acfs[abs(series_lags) <= max_lag]
        series_lags <- series_lags[abs(series_lags) <= max_lag]
    }
    
    # Seleccionamos los retardos significativos
    n_used <- series_prewhiten$ccf$n.used
    significative_lags <- abs(series_acfs) >= qnorm(1-alpha/2)/sqrt(n_used)
    series_lags <- series_lags[significative_lags]
    series_acfs <- series_acfs[significative_lags]
    
    # Seleccionamos el retardo "más significativo" 
    optimal_lag <- series_lags[which.max(abs(series_acfs))]
    
    # Si no se ha conseguido ningún retardo con estas condiciones, se devuelve NA
    if (length(optimal_lag) == 0) {
        return(NA)
    }
    
    # En otro caso, se devuelve el retardo obtenido
    return(optimal_lag*frequency(serie))
}


