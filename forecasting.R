
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
         xreg_pred <- forecast(xreg_ajuste, bootstrap=bootst, h=h, level=levels, xreg=NULL)
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

forecast_model(ajuste, h=10, mode='bootstrap')

