
forecast_model <- function(ajuste, h, mode=c('bootstrap', 'norm'), levels=c(80, 90)) {
    if (class(mode) != 'character') {stop('Argumento `mode` debe ser una cadena de caracteres')}
    if (!(mode %in% c('bootstrap', 'norm'))) {stop('Modo no válido')}
    
    finalT <- end(na.remove(ajuste$x))[1]
    
    xregs_pred <- matrix(NA, nrow=h, ncol=ncol(ajuste$xreg))
    for (j in 1:ncol(ajuste$xreg)) {
        if (mode == 'bootstrap') {bootst=T} else if (mode == 'norm') {bootst <- F}
        
         # Predicciones puntuales para cada variable
         xreg <- na.remove(ajuste$xreg[, j])
         xreg_ajuste <- auto.fit.arima(xreg, show_info = F)
         
         if (!bootst && !is_norm(xreg_ajuste)) {
             bootst <- TRUE
         }
         if (end(xreg)[1] == finalT) {
             xreg_pred <- forecast(xreg_ajuste, bootstrap=bootst, h=h)$mean
         } else if (end(xreg)[1] > finalT) {   # se pueden aprovechar valores futuros
             xreg_pred <- window(xreg, start=finalT+1, end=min(finalT + h, end(xreg)[1]))
             if (length(xreg_pred) < h) {
                 xreg_pred <- ts(c(xreg_pred, forecast(xreg_ajuste, bootstrap=bootst, h=h-length(xreg_pred))$mean))
             }
         }
         xregs_pred[, j] <- xreg_pred
    }
    colnames(xregs_pred) <- colnames(ajuste$xreg)
    
    # Calculamos las predicciones
    ajuste$x <- na.remove(ajuste$x)
    ajuste$fitted <- suppressWarnings(window(ajuste$fitted, start=start(ajuste$x)[1], end=end(ajuste$x)[1]))
    ajuste$residuals <- suppressWarnings(window(ajuste$residuals, start=start(ajuste$x)[1], end=end(ajuste$x)[1]))
    ajuste$xreg <- suppressWarnings(window(ajuste$xreg, start=start(ajuste$x)[1], end=end(ajuste$x)[1]))
    preds <- forecast(ajuste, h=h, bootstrap=(!is_norm(ajuste) || (mode=="bootstrap")), 
                      xreg=xregs_pred, level=levels)
    return(preds)
    
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
    