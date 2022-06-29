PAD <- 90

#' Selección automática de variables regresoras en un modelo de regresión dinámica
#' 
#' \deqn{ Y_t = \beta_0 + \beta_1 X_{t-r_1}^{(1)} + \beta_2 X_{t-r_2}^{(2)} + \cdots 
#' + X_{t-r_p}^{(p)} + \eta_t}
#' 
#' donde \eqn{Y_t} es la variable respuesta; \eqn{\{ X_{t-r_i}^{(i)}, \ i=1,...p \} } es 
#' el conjunto de variables explicativas con sus respectivos retardos \eqn{i=1,...,p}; 
#' \eqn{\beta_0,\beta_1,...,\beta_p} son los coeficientes de regresión lineal; y 
#' \eqn{\eta_t} \sim ARMA(p,q) son las innovaciones del modelo.
#' 
#' A partir del proceso de selección de retardos óptimos para cada variable regresora 
#' y el ajuste automático de modelos ARIMAX, se añaden de forma automática aquellas 
#' variables de un conjunto dado que inciden en la variable respuesta siguiendo un 
#' procedimiento iterativo en base a @cryer2008time.
#'
#' @param serie Serie temporal que funciona como variable respuesta.
#' @param xregs Conjunto de variables candidatas a ser variables regresoras de `serie`.
#' @param ic Criterio de información sobre el que se seleccionarán variables regresoras 
#' y se comparará la efectividad de un ajuste.
#' @param alpha Nivel de significación para los contrastes de hipótesis de estacionariedad 
#' y significancia.
#' @param stationary_method Método utilizado para chequear la estacionariedad de una 
#' serie temporal. Si `stationary_method = 'auto.arima'`, se utiliza la función 
#' \code{\link[forecast]{auto.arima}} para ajustar un modelo ARIMA(p,d,q) y chequear si 
#' \eqn{d > 0} (si se cumple esta condición se asume que la serie no es estacionaria). 
#' Si `stationary_method = 'adf.test'` se usa el \link[tseries:adf.test]{test Dickey-Fuller} 
#' para chequear la estacionariedad de una serie temporal.
#' @param show_info Valor booelano que indica si mostrar o no por pantalla los modelos 
#' que se van ajustando, así como información relacionada con la selección de variables.
#' @param ndiff Parámetro interno del programa (no utilizar!) para diferenciar todas las variables 
#' cuando no se pueda ajustar un modelo válido y mantener un registro del número de diferenciaciones 
#' que se están realizando.
#'
#' @returns Ajuste del modelo ARIMAX con las variables regresoras seleccionadas. Si dicho 
#' ajuste no existe (o no es posible de optimizar), devuelve el ajuste válido de un 
#' modelo ARIMA sobre la variable respuesta. Si éste último ajuste tampoco existe (no es 
#' posible de optimizar), devuelve `NA`.
#' @export
#'

auto.fit.arima.regression <- function(serie, xregs, ic='aicc', alpha=0.05, 
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
    
    # Inicialización de variables que se utilizarán a lo largo del algoritmo
    best_ic <- Inf
    model_history <- data.frame(var=NA, lag=NA, ic=NA)
    
    # Cálculo del retardo óptimo máximo de cada variable regresora sobre la variable respuesta. 
    # Se utilizará posteriormente para recortar las series y poder comparar modelos construidos 
    # con series del mismo tamaño
    max_lag <- get.maximum.lag(serie, xregs, alpha=alpha, method=stationary_method)
    
    response <- serie       # inicialmente, response = serie
    
    # Inicio del bucle para añadir variables regresoras
    for (i in 1:ncol(xregs)) {
        
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
            optimal_lag <- select.optimal.lag(response, xreg, alpha=0.05, max_lag=max_lag, method=stationary_method)
                
            
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
    if (!any(is.na(model_history))) {
        data_new <- construct.data(model_history, serie, xregs, new_xreg_name = NULL, 
                                   optimal_lag = NULL, max_lag=max_lag)
        global_fit <- auto.fit.arima(data_new[, c(1)], xregs=data_new[, -c(1)],
                                     ic=ic, d=0, D=0, alpha=alpha, show_info=F)    
        
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
    }
    
    
    # Si el model_history está vacío (es decir, no hay variables regresoras que 
    # aporten algo a la variable respuesta), o no se ha podido ajustar uno con 
    # errores ARMA, se trata de de:
    # a) Diferenciar todas las variables y volver a llamar a la función auto.fit.arima.regression
    # b) Si ndiff=2 no se recomienda diferenciar más y se ajusta un ARIMA sin regresoras
    if (ndiff < 3) {
        cat(paste0('No se ha podido encontrar un modelo válido con errores estacionarios\n',
            'Se aplica una diferenciación regular (ndiff=', ndiff+1, 
                       ') y se vuelve a llamar a la función\n'))
        
        serie <- diff(serie)
        xregs <- as.data.frame(lapply(xregs, diff))
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




#' Función auxiliar para la construcción de una matriz de series temporales en base al 
#' histórico del modelo y una nueva variable a añadir con su respectivo retardo óptimo y 
#' retardo máximo del conjunto total. La matriz resultante debe satisfacer que:
#' 
#' - Su primera columna sea la variable respuesta.
#' - Las otras columnas se correspondan con las variables regresoras del modelo ya 
#' retardadas por su retardo óptimo.
#' - A todas las variables se les haya recortado los primeros `max_lag` valores (para 
#' tener modelos construidos con el mismo número de muestras).
#' 
#' @param model_history DataFrame con el histórico de las variables que ya se han 
#' añadido al modelo (nombre de la variable, su retardo y su valor en el criterio de 
#' información seleccionado).
#' @param serie Serie temporal que funciona como variable respuesta.
#' @param xregs Conjunto de variables candidatas a ser variables regresoras de `serie`.
#' @param new_xreg_name Nombre de la nueva variable regresora que se quiere añadir al 
#' modelo. Si es `NULL` se asume que no se quiere añadir ninguna variable nueva y la 
#' matriz se contruye en base a las que ya están en el historial.
#' @param optimal_lag Retardo óptimo de la nueva variable a incluir. Puede ser `NULL` en 
#' caso de que `new_xreg_name` también lo sea.
#' @param max_lag Retardo máximo de todos los retardos óptimos (en valor absoluto) del 
#' conjunto de variables regresoras.
#'
#' @returns Matriz con las variables retardadas y recortadas de forma adecuada para 
#' introducir a la función de ajuste automático de ARIMAXs.
#'
construct.data <- function(model_history, serie, xregs, new_xreg_name, optimal_lag, max_lag) {
    # Primero recortamos los primeros max_lag valores de la variable respuesta 
    serie <- window(serie, start=add_t_start(serie, max_lag))
    data <- cbind(serie)
    
    # Si hay más variables en el historialdel modelo, las añadimos
    if (!all(is.na(model_history))) {
        for (j in 1:nrow(model_history)) {
            data <- cbind(data, 
                    window(lag(xregs[[model_history$var[j]]], as.integer(model_history$lag[j])),
                           start=start(serie), end=end(serie)))
        }
        colnames(data) <- c('serie', model_history$var)
    } 
    
    if (!is.null(new_xreg_name)) {
        new_xreg <- lag(xregs[[new_xreg_name]], optimal_lag)
        data <- cbind(data, 
                      window(new_xreg, start=start(serie), end=end(serie)))
        
        if (!all(is.na(model_history))) {
            colnames(data) <- c('serie', model_history$var, new_xreg_name)
        } else {  colnames(data) <- c('serie', new_xreg_name) }
    }
    return(data)
}

add_t_start <- function(serie, t) {
    if (frequency(serie) == 1) {
        return( start(serie)[1] + t )
    } else {
        freq <- frequency(serie)
        raw_start <- start(serie)[2] + t
        return( 
            c(start(serie)[1] + floor(raw_start/freq), t %% freq ) 
        )
    }
}

add_t_end <- function(serie, t) {
    if (frequency(serie)==1) {
        return(c(end(serie)[1] + t, 1))
    } else {
        freq <- frequency(serie)
        raw_end <- end(serie)[2] + t
        return(
            c( end(serie)[1] + floor(raw_end/freq), raw_end%%freq)
        )
    }
}

#' Comprobación de la tendencia de una serie temporal
#' 
#' @param x Serie temporal de la que se quiere chequear la tendencia.
#' @param test Método usado para chequear la tendencia Si `test = 'auto.arima'` se 
#' intenta ajustar a la serie `x` un modelo ARIMA. Si el orden de diferenciación $d$ es 
#' mayor a 0, se considera que la serie tiene tendencia. Si `test = 'adf.test'` se 
#' utiliza el test Dickey-Fuller.
#' @param alpha Nivel de significación para el contraste de hipótesis de Dickey-Fuller. 
#' Por defecto 0.05.
#'
#' @returns `True` si la serie tiene tendencia, `False` en otro caso.
#' @export
#'
has_trend <- function(x, test='auto.arima', alpha=0.05) {
    if (!test %in% c('auto.arima', 'adf.test')) {
        stop('El parámetro `test` debe ser "auto.arima" o "adf.test"')
    }
    
    if (test == 'auto.arima') {
        ajuste <- suppressWarnings(auto.arima(x, stepwise=F, approximation=F, max.d=3))
        d <- ajuste$arma[6]
        return(d>0)
    }
    if (test == 'adf.test') {
        pvalue <- suppressWarnings(adf.test(x[!is.na(x)])$p.value)
        return(alpha <= pvalue)
    }
}



#' Selección del retardo óptimo (siempre menor o igual a 0) para una serie temporal y su 
#' variable regresora. De forma opcional se puede fijar una magnitud máxima (`max_lag`) 
#' del retardo a escoger.
#' 
#' @param serie Variable respuesta.
#' @param xreg Variable regresora.
#' @param alpha Nivel de significación para los contrastes de hipótesis de tendencia, 
#' estacionalidad y significación.
#' @param max_lag Magnitud máxima del retardo a escoger (filtrado de retardos).
#' @param method Método para chequear estacionariedad de las series.
#' 
#' @returns Retardo óptimo para la variable respuesta y su variable regresora basada 
#' en el método propuesto por @cryer2008time. Si no es posible obtener este retardo, 
#' se devuelve `NA`.
#' @export
#'
select.optimal.lag <- function(serie, xreg, alpha=0.05, max_lag=NA, method='auto.arima', less0=T) {
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
    series_prewhiten <- TSA::prewhiten(ts(xreg), ts(serie), plot=F)
    
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
    
    # Si no se ha conseguido ningún retardo con estas condiciones, se devuleve NA
    if (length(optimal_lag) == 0) {
        return(NA)
    }
    
    # En otro caso, se devuelve el retardo obtenido
    return(optimal_lag)
}



#' Obtención del retardo máximo (en valor absoluto) sobre un conjunto de variables 
#' regresoras y su variable respuesta.
#' 
#' @param serie Variable respuesta.
#' @param xregs Conjunto de variables regresoras.
#' @param alpha Nivel de significación para los contrastes de hipótesis de tendencia, 
#' estacionalidad y significación.
#' @param method Método para chequear estacionariedad de las series.
#' 
#' @returns Máximo retardo (en valor absoluto) de todos los retardos que se obtengan sobre 
#' el conjunto de varaibles regresoras.
#' @export
#'
get.maximum.lag <- function(serie, xregs, alpha=0.05, method='auto.arima') {
    if (class(xregs) != 'data.frame') {
        stop('El argumento `xregs` debe ser un data.frame')
    }
    if (class(serie) != 'ts') {
        stop('El arguemnto `serie` debe ser de tipo ts')
    }
    
    max_lag <- -Inf
    
    for (i in 1:ncol(xregs)) {
        xreg <- xregs[, i]
        optimal_lag <- select.optimal.lag(serie, xreg, alpha=alpha, method=method, less0=F)
        if (is.na(optimal_lag)) {
            next
        }
        if (abs(optimal_lag) > max_lag) {
            max_lag <- abs(optimal_lag)
        }
    }
    return(max_lag)
}
