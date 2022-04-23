
forecast_model <- function(ajuste, h, mode=c('bootstrap', 'norm'), levels=c(80, 90)) {
    if (class(mode) != 'character') {stop('Argumento `mode` debe ser una cadena de caracteres')}
    if (!(mode %in% c('bootstrap', 'norm'))) {stoop('Modo no válido')}
    
    # Quitamos los valores nulos
    ajuste <- na.remove(ajuste)
    
    xregs_pred <- matrix(NA, nrow=h, ncol=ncol(ajuste$xreg))
    for (j in 1:ncol(ajuste$xreg)) {
        if (mode == 'bootstrap') {bootst=F} else if (mode == 'norm') {bootst <- F}
        
         # Predicciones puntuales para cada variable
         xreg_ajuste <- auto.fit.arima(ajuste$xreg[, j], show_info = F)
         if (!bootst && !is_norm(xreg_ajuste)) {
             bootst <- TRUE
         }
         xreg_pred <- forecast::forecast(xreg_ajuste, bootstrap=bootst, h=h, level=levels)
         xregs_pred[, j] <- xreg_pred$mean
    }
    colnames(xregs_pred) <- colnames(ajuste$xreg)
    
    # Calculamos las predicciones
    preds <- forecast(ajuste, h=h, bootstrap=(!is_valid(ajuste) || (mode=="bootstrap")), xreg=xregs_pred)
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


plot_forecast <- function(preds) {
    colors <- c('rgba(253,205,172,0.6)', 'rgba(203,213,232,0.6)', '
                rgba(244,202,228,0.6)', 'rgba(230,245,201,0.6)', 'rgba(255,242,174,0.6)', 
                'rgba(241,226,204,0.6)', 'rgba(204,204,204,0.6)')
    
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
    