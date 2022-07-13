


forecast_model <- function(serie, xregs, ajuste, h, mode='bootstrap', levels=c(80, 90)) {
    if (class(xregs) != 'data.frame') { stop('El argumento `xregs` debe ser un data.frame') }
    if (!all(unlist(lapply(xregs, class)) == "ts")) {
        stop('Las variables regresoras del data.frame `xregs` deben ser de tipo ts')
    }
    if (class(serie) != 'ts') { stop('El argumento `serie` debe ser de tipo ts') }
    if (class(mode) != 'character') {stop('Argumento `mode` debe ser una cadena de caracteres')}
    if (!(mode %in% c('bootstrap', 'norm'))) {stop('Modo no válido')}
    
    serie_original <- serie; xregs_original <- xregs
    if (ajuste$ndiff > 0) {
        for (i in 1:ajuste$ndiff) {
            serie <- diff(serie)
            xregs <- as.data.frame(lapply(xregs, diff))
        }
    }
    
    h_moment <- add_t_end(ajuste$x, h)
    xregs_pred <- matrix(NA, nrow=h, ncol=ncol(ajuste$xreg))
    colnames(xregs_pred) <- ajuste$history$var
    
    
    # Para cada variable regresora obtenemos predicciones puntuales (o bien las podemos 
    # sacar de los valores de la serie original)
    for (xreg_name in ajuste$history$var) {
        if (mode == 'bootstrap') {bootst=T} else if (mode == 'norm') {bootst <- F}
        
        xreg <- xregs[[xreg_name]]   # cogemos la variable completa
        xreg_ajuste <- auto.fit.arima(xreg, show_info = F)
        xreg_lag <- as.numeric(ajuste$history$lag[ajuste$history$var == xreg_name])
        
        # cabe la posibilidad de que ya tengamos valores en la serie original 
        # (pero que han sido recortados en el objeto ajuste)
        if (abs(xreg_lag) >= h) {
            xregs_pred[, c(xreg_name)] <- xreg[(length(xreg) - h + 1):length(xreg)]
            next
        } 
          
         
        if (!bootst && !is_norm(xreg_ajuste)) { bootst <- TRUE }
        
        # realizamos las suficientes predicciones a horizonte (h_moment- end(xreg)
        
        xreg_pred <- forecast(xreg_ajuste, bootstrap=bootst, h=(h + xreg_lag))$mean
        
        if (xreg_lag == 0) {
            xregs_pred[, c(xreg_name)] <- xreg_pred
        } else {
            xregs_pred[, c(xreg_name)] <- c(xreg[(length(xreg) + xreg_lag + 1):length(xreg)], xreg_pred)
        }
    }
    
    
    # Calculamos las predicciones
    if (ncol(ajuste$xreg) == 1) {
        colnames(xregs_pred) <- c('xreg')
    }
    
    preds <- forecast(ajuste, h=h, bootstrap=(!is_norm(ajuste) || (mode=="bootstrap")), 
                      xreg=xregs_pred, level=levels)
    return(preds)
    
}

moments_diff <- function(moment1, moment2, freq) {
    if (freq==1) {
        return(moment1[1]-moment2[1])
    } else {
        return(
            (moment1[1]*freq + moment1[2]) - (moment2[1]*freq + moment2[2])
        )
    }
}


is_norm <- function(ajuste, alpha=0.05) {
    testJB <- jarque.bera.test(ajuste$residuals[!is.na(ajuste$residuals)])$p.value
    testSW <- shapiro.test(ajuste$residuals)$p.value
    
    if (any(c(testJB, testSW) < alpha)) {
        return(FALSE)
    }
    return(TRUE)
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
    fig <- fig %>% add_trace(x=time(preds$mean), y=preds$mean, line=list(color='grey'),
                     name='predicciones puntuales', showlegend=T)
    
    return(fig)
}
    