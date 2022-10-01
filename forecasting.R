

#' README:
#' This script contains the necessary functions to automatically apply 
#' forecasting to dynamic regression models.


#' 
#' Given an horizon $h$, compute point predictions and prediction intervals 
#' based in a DRM computed by the function `drm.select()`.
#' 
#' @param serie [ts]: Univariate time series which acts as the dependent 
#' variable. It must be the same variable used in `drm.select()`.
#' @param xregs [mts]: Matrix time series with the selected covariates. 
#' It must be the same matrix used in `drm.select()`.
#' @param fitted_model [Arima]: Fitted ARIMAX model obtained.
#' @param h [numeric]: Horizon value for forecasting.
#' @param mode [character]: Type of prediction intervals.: based on 
#' residuals normality (`norm`) or via bootstrap (`bootstrap`). By 
#' default `bootstrap` .
#' @param levels [vector]: Numerical vector of the confidence intervals.
#' 
#' @returns Point predictions and prediction intervals in original scale.
forecast_model <- function(serie, xregs, fitted_model, h, mode='bootstrap', levels=c(80, 90)) {
    
    # ------------------------------- assertions -------------------------------
    if (!any(c('mts', 'ts') %in% class(xregs))) {
        stop('[forecast_model] TypeError: Parameter `xregs` must be `ts` or `mts`')
    }
    if (class(serie) != 'ts') {
        stop('[forecast_model] TypeError: Parameter `serie` must be `ts`')
    }
    if (class(mode) != 'character') {
        stop('[forecast_model] TypeError: Parameter `mode` must be `character`')}
    if (!(mode %in% c('bootstrap', 'norm'))) {
        stop('[forecast_model] ArgumentError: Parameter `mode` must be `bootstrap` or `norm`')
    }
    # --------------------------------------------------------------------------

    # If the fitted model has been computed with differentiated data, it is 
    # necessary to differentiate the original data too
    serie_original <- serie; xregs_original <- xregs # store original data
    
    if (fitted_model$ndiff > 0) {
        serie <- mul_diff(serie, fitted_model$ndiff)
        xregs <- mul_diff(xregs, fitted_model$ndiff)
    }
    
    # Prediction matrix with horizon h of covariates
    xregs_pred <- matrix(NA, nrow=h, ncol=ncol(fitted_model$xreg))
    colnames(xregs_pred) <- fitted_model$history$var
    
    
    # For each covariate, obtain its point predictions or obtain future 
    # values from the original covariate data
    for (xreg_name in fitted_model$history$var) {
        if (mode == 'bootstrap') {bootst=T} else if (mode == 'norm') {bootst <- F}
        
        xreg <- xregs[, c(xreg_name)]           # this has complete data
        xreg_fitted_model <- auto.fit.arima(xreg, show_info = F)
        xreg_fitted_model$call <- NULL
        xreg_lag <- as.numeric(fitted_model$history$lag[fitted_model$history$var == xreg_name])
        
        # We might have future values from original data
        if (abs(xreg_lag) >= h) {
            xregs_pred[, c(xreg_name)] <- xreg[(length(xreg) - h + 1):length(xreg)]
            next
        } 
          
        # We need to check if the residuals of the covariate model satisfy 
        # normality
        if (!bootst && !is_norm(xreg_fitted_model)) {
            bootst <- TRUE 
        }
        
        # Compute the necessary predictions to horizon h
        xreg_pred <- forecast::forecast(xreg_fitted_model, bootstrap=bootst, h=(h + xreg_lag))$mean
        
        if (xreg_lag == 0) {
            xregs_pred[, c(xreg_name)] <- xreg_pred
        } else {
            xregs_pred[, c(xreg_name)] <- c(xreg[(length(xreg) + xreg_lag + 1):length(xreg)], xreg_pred)
        }
    }

    # If only one covariate is in the model, rename it to use it with the function 
    # forecast()
    if (ncol(fitted_model$xreg) == 1) {
        colnames(xregs_pred) <- c('xreg')
    }
    
    preds <- forecast(fitted_model, h=h, 
                      bootstrap=(!is_norm(fitted_model) || (mode=="bootstrap")), 
                      xreg=xregs_pred, level=levels)
    
    
    # If differentiation has been applied, convert predictions to original scale
    if (fitted_model$ndiff > 0) {
        cat('Returning predictions in original scale...\n')
        apply_undiff <- function(x) {
            return(
                mul_undiff(ts(c(serie, x), start=start(serie)), 
                           previous=window(serie_original, 
                                           start=substract_t(start(serie), fitted_model$ndiff, frequency(serie)),
                                           end=substract_t(start(serie), 1, frequency(serie))),
                           fitted_model$ndiff)
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


#' Auxiliar function to check if the residuals of a fitted model statisfy normality
#' 
#' @param fitted_model [Arima]: Fitted ARIMA(X) model.
#' @param alpha [numeric]: Significance level.
#' @returns `TRUE` if the residuals of the model satisfy normality. `FALSE`, 
#' otherwise.
is_norm <- function(fitted_model, alpha=0.05) {
    testJB <- jarque.bera.test(fitted_model$residuals[!is.na(fitted_model$residuals)])$p.value
    testSW <- shapiro.test(fitted_model$residuals)$p.value
    
    if (any(c(testJB, testSW) < alpha)) {
        return(FALSE)
    }
    return(TRUE)
} 


#' Auxiliar function which makes the inverse of the diff() function (given the 
#' first original value which we lose after the differentiation).
#' 
#' @param x [ts]: Differentiated time series.
#' @param previous [numeric]: First value of the original serie.
#' 
#' @return Original data.
undiff <- function(x, previous) {
    
    if (class(x) != 'ts') {
        stop('[undiff] TypeError: Parameter `x` must be `ts`')
    }
    
    x_undiff <- rep(NA, length(x) + 1)
    x_undiff[1] <- previous
    
    for (i in 2:length(x_undiff)) {
        x_undiff[i] <- x_undiff[i-1] + x[i-1]
    }
    
    return(ts(x_undiff, end=end(x)))
}

#' Auxiliar function that applies multiple `undiffs`.
#' @param x [ts]: Differentiated time series.
#' @param previous [numeric]: First value of the original serie.
#' @param ndiff [numeric]: Number of undiffs that are applied.
#' 
#' @returns Origial data.
mul_undiff <- function(x, previous, ndiff) {
    if (class(x) != 'ts') {stop('El objeto x debe ser de tipo ts')}
    if (length(previous) != ndiff) {stop('No hay datos suficientes')}
    
    if (ndiff == 1) {
        return(undiff(x, previous))
    }
    
    return(mul_undiff(undiff(x, mul_diff(previous, ndiff-1)), previous[1:(ndiff-1)], ndiff-1))
    
}


#' Auxiliar function that applies multiple regular differentiations.
#' 
#' @param x [ts]: Univariate time series.
#' @param ndiff [numeric]: Number of regular differentiations that are applied.
#' @returns: Differentiated data.
#' 
mul_diff <- function(x, ndiff) {
    x_diff <- x
    for (i in 1:ndiff) {
        x_diff <- diff(x_diff)
    }
    return(x_diff)
}


#' Auxiliar function to substract $t$ moments to a given moment of the `ts` 
#' class (handling seasonal data).
#' 
#' @param moment [vector]: Vector of two values that give information about 
#' the current moment. Note that the meaning of each value is the same of `ts` 
#' indexes.
#' @param t [numeric]: Number of moments that we want to substract.
#' @param freq [numeric]: Frequency of the time series. It is 1 if there is 
#' not a seasonal component. 
#' @returns Moment resulting from substracting `t` moments to the current moment.
#' 
#' @examples
#' --------------------------------------
#' > substract_t(c(6, 3), 6, 7)
#' c(5, 4)
#' --------------------------------------
#' > substract_t(c(6, 1), 6, 1)
#' c(0, 1)
#' --------------------------------------
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



#' Auxiliar plotting function to display the predictions of the 
#' `forecast_model()` function
#' @param preds: Output of the `forecast_model()`  function.
#' @param rang [vector]: Two-size vector that indicates the range of the x-axis.
#' @returns Plotly figure.
#' 
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
                  name='time serie', showlegend=T)
    
    for (i in rev(1:length(preds$level))) {
        fig <- fig %>% 
            add_trace(x=time(preds$upper[, i]), y=preds$upper[, i],
                      line=list(color=colors[i]), name=paste0('IC ', preds$level[i], ' % conf.')) %>%
            add_trace(x=time(preds$lower[, i]), y=preds$lower[, i], 
                      line=list(color=colors[i]), name=paste0('IC ', preds$level[i], ' % conf.'),
                      fill='tonexty', fillcolor=colors[i], showlegend=F)
        
    }
    fig <- fig %>% add_trace(x=time(preds$mean), y=preds$mean, line=list(color='#2F63D4'),
                     name='forecasting', showlegend=T)
    
    
    return(fig)
}
    