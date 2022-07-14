


#' Predicciones puntuales a horizonte $h$ e intervalos de predicción en base a 
#' un modelo de regresión dinámica.
#' 
#' @param serie [`ts`] Serie de tiempo que funciona como variable respuesta. 
#' Debe ser la serie original que se utilizó para la selección de variables 
#' regresoras.
#' @param xregs [`mts`] Conjunto original de todas las variables regresoras. 
#' @param ajuste [`Arima`] Ajuste del modelo de regresión dinámica obtenido y 
#' sobre el que se hacen las predicciones puntuales e intervalos de predicción.
#' @param h [`numeric`] Valor horizonte de las predicciones.
#' @param mode [`character`] Modo de realizar los intervalos de predicción: 
#' basados en normalidad sobre los residuos (`norm`) o a través de bootstrap 
#' (`bootstrap`). Por defecto se realizan a través de bootstrap.
#' @param levels [`vector`] Vector numérico de los niveles a los que se quieren 
#' hacer los intervalos de predicción.
#' 
#' @returns Predicciones puntuales e intervalos de confianza en unidades 
#' originales.
#' 
forecast_model <- function(serie, xregs, ajuste, h, mode='bootstrap', levels=c(80, 90)) {
    
    # Comprobación de los parámetros     
    if (!any(c('mts', 'ts') %in% class(xregs))) {
        stop('El argumento `xregs` debe ser de tipo mts o ts')
    }
    if (class(serie) != 'ts') {
        stop('El arguemnto `serie` debe ser de tipo ts')
    }
    if (class(mode) != 'character') {stop('Argumento `mode` debe ser una cadena de caracteres')}
    if (!(mode %in% c('bootstrap', 'norm'))) {stop('Modo no válido')}
    
    
    # Si el ajuste se realizado con los datos diferenciados, es necesario 
    # diferenciarlos los originales también (aunque se guardan sin diferenciar)
    
    serie_original <- serie; xregs_original <- xregs  # datos sin diferenciar
    
    if (ajuste$ndiff > 0) {
        serie <- mul_diff(serie, ajuste$ndiff)
        xregs <- mul_diff(xregs, ajuste$ndiff)
    }
    
    # Matriz donde se guardarán las predicciones a horizonte h de las variables 
    # regresoras
    xregs_pred <- matrix(NA, nrow=h, ncol=ncol(ajuste$xreg))
    colnames(xregs_pred) <- ajuste$history$var
    
    
    # Para cada variable regresora obtenemos predicciones puntuales (o bien las 
    # podemos sacar de los valores de la serie original)
    for (xreg_name in ajuste$history$var) {
        if (mode == 'bootstrap') {bootst=T} else if (mode == 'norm') {bootst <- F}
        
        xreg <- xregs[, c(xreg_name)]   # cogemos la variable completa
        xreg_ajuste <- auto.fit.arima(xreg, show_info = F)
        xreg_lag <- as.numeric(ajuste$history$lag[ajuste$history$var == xreg_name])
        
        # cabe la posibilidad de que ya tengamos valores en la serie original 
        # (pero que han sido recortados en el objeto ajuste)
        if (abs(xreg_lag) >= h) {
            xregs_pred[, c(xreg_name)] <- xreg[(length(xreg) - h + 1):length(xreg)]
            next
        } 
          
         
        if (!bootst && !is_norm(xreg_ajuste)) { bootst <- TRUE }
        
        # realizamos las suficientes predicciones a horizonte h

        xreg_pred <- forecast(xreg_ajuste, bootstrap=bootst, h=(h + xreg_lag))$mean
        
        if (xreg_lag == 0) {
            xregs_pred[, c(xreg_name)] <- xreg_pred
        } else {
            xregs_pred[, c(xreg_name)] <- c(xreg[(length(xreg) + xreg_lag + 1):length(xreg)], xreg_pred)
        }
    }
    
    # Si sólo hay una variable regresora es necesario renombrarla para utilizarla 
    # con la función `forecast()`
    if (ncol(ajuste$xreg) == 1) {
        colnames(xregs_pred) <- c('xreg')
    }
    
    # Realizamos las predicciones
    preds <- forecast(ajuste, h=h, bootstrap=(!is_norm(ajuste) || (mode=="bootstrap")), 
                      xreg=xregs_pred, level=levels)
    
    
    # Si se han aplicado diferenciaciones, pasamos las predicciones y los 
    # intervalos de predicción a unidades originales (con los datos originales 
    # sin diferenciar que se habían guardado)
    if (ajuste$ndiff > 0) {
        cat('Se devuelven las predicciones en unidades originales...\n')
        apply_undiff <- function(x) {
            return(
                mul_undiff(ts(c(serie, x), start=start(serie)), 
                           previous=window(serie_original, 
                                           start=substract_t(start(serie), ajuste$ndiff, frequency(serie)),
                                           end=substract_t(start(serie), 1, frequency(serie))),
                           ajuste$ndiff)
        )}
        
        
        serie_mean <- apply_undiff(preds$mean)
        serie_lower <- ts(as.matrix(apply(preds$lower, 2, apply_undiff)), start=start(serie_original))
        serie_upper <- ts(as.matrix(apply(preds$upper, 2, apply_undiff)), start=start(serie_original))
        
        preds$mean <- window(serie_mean, start=start(preds$mean), end=end(preds$mean))
        preds$lower <- window(serie_lower, start=start(preds$lower))
        preds$upper <- window(serie_upper, start=start(preds$upper))
        preds$x <- window(serie_mean, start=start(preds$x), end=end(preds$x))
    }
    
    return(preds)
}


# Función auxiliar para comprobar si los residuos de un ajuste son normales
is_norm <- function(ajuste, alpha=0.05) {
    testJB <- jarque.bera.test(ajuste$residuals[!is.na(ajuste$residuals)])$p.value
    testSW <- shapiro.test(ajuste$residuals)$p.value
    
    if (any(c(testJB, testSW) < alpha)) {
        return(FALSE)
    }
    return(TRUE)
} 

# Función auxiliar para realizar el proceso "inverso" de diferenciar (partiendo 
# de que se cuenta con el primer valor de la serie de tiempo que se "perdió" 
# en la diferenciación).
undiff <- function(x, previous) {
    
    if (class(x) != 'ts') {stop('El objeto x debe ser de tipo ts')}
    
    x_undiff <- rep(NA, length(x) + 1)
    x_undiff[1] <- previous
    
    for (i in 2:length(x_undiff)) {
        x_undiff[i] <- x_undiff[i-1] + x[i-1]
    }
    
    return(ts(x_undiff, end=end(x)))
}

# Función auxiliar que aplicar el "undiff" de múltiples diferenciaciones 
mul_undiff <- function(x, previous, ndiff) {
    if (class(x) != 'ts') {stop('El objeto x debe ser de tipo ts')}
    if (length(previous) != ndiff) {stop('No hay datos suficientes')}
    
    
    # si ndiff==1 se puede utilizar directamente la función simple
    if (ndiff == 1) {
        return(undiff(x, previous))
    }
    
    return(mul_undiff(undiff(x, mul_diff(previous, ndiff-1)), previous[1:(ndiff-1)], ndiff-1))
    
}


# Función auxiliar para aplicar múltiples diferenciaciones regulares
mul_diff <- function(x, ndiff) {
    x_diff <- x
    for (i in 1:ndiff) {
        x_diff <- diff(x_diff)
    }
    return(x_diff)
}


plot_forecast <- function(preds, rang=NULL) {
    colors <- c('rgba(253,205,172,0.6)', 'rgba(203,213,232,0.6)', '
                rgba(244,202,228,0.6)', 'rgba(230,245,201,0.6)', 'rgba(255,242,174,0.6)', 
                'rgba(241,226,204,0.6)', 'rgba(204,204,204,0.6)')
    
    if (!is.null(rang)) {
        preds$x <- suppressWarnings(window(preds$x, start=rang[1], end=rang[2]))
        preds$upper <- suppressWarnings(window(preds$upper, start=rang[1], end=rang[2]))
        preds$lower <- suppressWarnings(window(preds$lower, start=rang[1], end=rang[2]))
        preds$mean <- suppressWarnings(window(preds$mean, start=rang[1], end=rang[2]))
    }
    
    fig <- plot_ly(type='scatter', mode='lines') %>%
        add_trace(x=time(preds$x), y=preds$x, line=list(color='grey'),
                  name='serie temporal', showlegend=T)
    
    for (i in rev(1:length(preds$level))) {
        fig <- fig %>% 
            add_trace(x=time(preds$upper[, i]), y=preds$upper[, i],
                      line=list(color=colors[i]), name=paste0('IC ', preds$level[i], ' % conf.')) %>%
            add_trace(x=time(preds$lower[, i]), y=preds$lower[, i], 
                      line=list(color=colors[i]), name=paste0('IC ', preds$level[i], ' % conf.'),
                      fill='tonexty', fillcolor=colors[i], showlegend=F)
        
    }
    fig <- fig %>% add_trace(x=time(preds$mean), y=preds$mean, line=list(color='#2F63D4'),
                     name='predicciones puntuales', showlegend=T)
    
    return(fig)
}
    