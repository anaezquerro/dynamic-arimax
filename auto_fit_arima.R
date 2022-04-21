
#' Ajuste automático de un modelo ARIMA o un modelo ARIMAX
#' 
#'
#' @param serie : Serie temporal sobre la que se quiere ajustar un modelo ARIMA.
#' @param xregs : Conjunto de variables regresoras que inciden en `serie`. 
#' @param ic : Criterio de información para evaluar los ajustes obtenidos sobre la 
#' serie de tiempo.
#' @param d : Valor del orden de diferenciación regular. Por defecto `NA` (se escoge de 
#' forma automática siguiendo el criterio de información).
#' @param D : Valor del orden de diferenciación estacional. Por defecto `NA` (se escoge 
#' de forma automática siguiendo el criterio de información).
#' @param alpha : Nivel de significación para los tests estadísticos (chequear 
#' independencia, media nula y normalidad de los residuos y significancia de los 
#' coeficientes ajustados).
#' @param show_info : Valor booleano que indica si se muestra o no por pantalla los 
#' modelos que se van ajustando.
#'
#' @return Ajuste del modelo ARIMA cuyos coeficientes son significativos y cuyos residuos 
#' cumplen las hipótesis de independencia y media nula. En caso de que no sea posible  
#' ajustar este modelo, se devuelve `NA`.
#' @export
#' 
auto.fit.arima <- function(serie, xregs=NULL, ic='aicc', d=NA, D=NA, alpha=0.05, show_info=T) {
    
    if (class(serie) != 'ts') {
        stop('El argumento `serie` debe ser de tipo ts')
    }
    if (!is.null(xregs) && !all(c(unlist(lapply(xregs, class))) == 'ts')) {
        stop('El argumento `xregs` debe ser de tipo ts')
    }

    # Selección del mejor modelo con la función auto.arima guardando el historial 
    trace <- capture.output({
        aux <- suppressWarnings(
            auto.arima(serie, xreg=xregs, d=d, D=D, seasonal=frequency(serie) > 1, 
                       ic=ic, max.d=3, stepwise=FALSE, approximation=FALSE, trace = TRUE) 
        )
    })
    con    <- textConnection(trace)
    models <- read.table(con, sep=':')
    models <- models[1:nrow(models)-1,]
    names(models) <- c('Modelo', ic)
    models <- sapply(models, trimws)
    models <- as.data.frame(models)
    
    # Bucle global que irá ajustando los modelos siguiendo el orden del criterio escogido
    while (nrow(models) > 0) {
        
        # Tomamos el modelo con el mejor valor para el criterio escogido (ic)
        best_model <- models[which.min(models[[ic]]), 1]
        
        # Obtenemos sus órdenes de la cadena de caracteres
        orders <- parse.orders(best_model, reg=!is.null(xregs), seas=frequency(serie)>1)
        
        # Intentamos estimar el modelo. Si no es posible optimizarlo (ajuste == NA), entonces 
        # lo eliminamos del trace recogido y probamos con el siguiente.
        ajuste <- fit.model(serie, xreg=xregs, orders)
        
        if (all(is.na(ajuste))) {
            models <- models[-c(which.min(models[[ic]])), ]
            next
        }
        
        if (show_info) {
            cat('----------------------------------------------------------------------\n')
            print(ajuste, row.names=F)
            cat('----------------------------------------------------------------------\n')
        }
        
        # Chequeamos si en el modelo hay parámetros que se pueden sacar
        ajuste <- fit.coefficients(ajuste, orders, alpha=alpha, show_info=show_info)
        
        # 2.5. Chequeamos si no ha habido problemas en el ajuste del modelo.
        if (all(is.na(ajuste)) || any(is.na(ajuste$var.coef))){
            ajuste <- NA
            models <- models[-c(which.min(models[[ic]])),]
            next
        }
        
        # 2.6 Si no ha habido problemas de optimización, realiazmos el 
        # análisis de residuos
        H <- ljungbox_lag(ajuste$residuals)
        ljungbox <- Box.test(ajuste$residuals, lag=H, type='Ljung-Box', fitdf=sum(ajuste$coef!=0))$p.value
        ttest <- t.test(ajuste$residuals, mu=0)$p.value
        testJB <- jarque.bera.test(ajuste$residuals[!is.na(ajuste$residuals)])$p.value
        testSW <- shapiro.test(ajuste$residuals)$p.value
        
        if (ljungbox < alpha) {
            models <- models[-c(which.min(models[[ic]])),]  # eliminamos el modelo del historial
            if (show_info) {
                cat(paste0('Falla la hipótesis de independencia de los residuos: p-valor=', ljungbox, 
                           '\nModelo no válido. Probamos con el siguiente modelo vía criterio ', ic, '\n'))
                cat('**********************************************************************\n') 
            }
            next
        }
        if (ttest < alpha) {
            models <- models[-c(which.min(models[[ic]])),]    # eliminamos el modelo del historial
            if (show_info) {
                cat(paste0('Falla la hipótesis de media nula de los residuos: p-valor=', ttest,
                           '\nModelo no válido. Probamos con el siguiente modelo vía criterio ', ic, '\n'))
                cat('**********************************************************************\n') 
            }
            next
        }
        
        if (any(c(testJB, testSW) < alpha)) {
            if (show_info) {
                cat(paste0('Falla la hipótesis de normalidad sobre los residuos.\n', 
                           'El modelo es válido pero los intervalos de predicción basados en la\n',
                           'dist. asintótica no son válidos\n')) 
            }
        }
        
        # Si se llega al final del bucle, el modelo es válido
        break
    }
    
    # 3. Comprobamos que el ajuste que hemos obtenido es válido
    if (all(is.na(ajuste)))  {
        if (show_info) {
            warning('No se ha podido encontrar ningún modelo para la serie temporal')
        }
        return(NA)
    }
    
    # 4. En caso de que si se haya obtenido un modelo válido, se muestra en pantalla 
    # y se devuelve
    if (show_info) {
        cat('\n----------------------------------------------------------------------\n')
        cat('|                             MODELO FINAL                           |\n')
        cat('----------------------------------------------------------------------\n')
        print(ajuste, row.names=F)
    }
    return(ajuste)
}



#' Obtener los órdenes del modelo ARIMA a partir de su ajuste
#'
#' @param ajuste Ajuste del modelo ARIMA.
#'
#' @return Lista con los órdenes regulares (`regular`), estacionales (`seasonal`) y una 
#' variable booleana para indicar si se incluye media/constante (`include_constant`).
#' @export
#'
get_orders <- function(ajuste) {
  if (!('Arima' %in% class(ajuste))) {
    stop('El argumento `ajuste` debe ser un objeto Arima')
  }
  
  p <- ajuste$arma[1]
  q <- ajuste$arma[2]
  P <- ajuste$arma[3]
  Q <- ajuste$arma[4]
  d <- ajuste$arma[6]
  D <- ajuste$arma[7]
  
  include_constant <- any(c('mean', 'drift') %in% names(ajuste$coef))
  orders <- list(regular=c(p, d, q), seasonal=c(P, D, Q), include_constant=include_constant)
  return(orders)
}



#' Obtención de los coeficientes ARMA a chequear de un ajuste
#'
#' @param ajuste Ajuste ARIMA sobre el que se quieren obtener los nombres de los 
#' coeficientes correspondientes a la parte ARMA (se excluyen coeficientes que se 
#' correspondan con las variables regresoras en un mdoelo ARIMAX). Esto es, `ar`, `ma`, 
#' `sar`, `sma` y, en caso de que d=0, se escoge `mean`, y en caso de que d=1, se escoge 
#' `drift`.
#'
#' @return Vector con los nombres de los coeficientes que se deben chequear. En caso de 
#' que no haya ninguno, se devuelve `NULL`.
#' @export
#'
get_arma_coefs_names <- function(ajuste) {
  if (!('Arima' %in% class(ajuste))) {
    stop('El argumento `ajuste` debe ser un objeto Arima')
  }
  
  # Obtenemos la máscara de los coeficientes ARMA
  arma_coefs_mask <- sapply(names(ajuste$coef), str_detect, pattern=c('s?ma\\d+', 's?ar\\d+', 'mean', 'drift'))
  
  # En caso de que la máscara esté vacía esto indica que el modelo es un ARIMA(0, 0, 0)
  if (length(arma_coefs_mask) == 0) { return(NULL) }
  
  # La máscara se calcula sobre cualquier nombre del coeficiente que cumpla algún pattern
  arma_coefs_mask <- apply(arma_coefs_mask, 2, sum) == 1
  
  # Se retiran aquellos que ya sean iguales a 0
  arma_coefs_mask <- arma_coefs_mask & (ajuste$coef != 0)
  
  # En caso de que ya no haya ningún coeficiente (todos están fijados a 0), se devuelve NULL
  if (sum(arma_coefs_mask) == 0) { return(NULL) }
  
  # Se obtienen sus nombres
  arma_coefs_names <- names(ajuste$coef[arma_coefs_mask])
  
  return(arma_coefs_names)
}

#' Ajuste de los coeficientes de un modelo ARIMA
#'
#' @param ajuste Ajuste inicial de un modelo ARIMA. Sus coeficientes pueden no ser 
#' significativos.
#' @param orders Órdenes del modelo.
#' @param alpha Nivel de significación estadística para los tests de significancia de los 
#' coeficientes.
#' @param show_info Valor booleano que indica si mostrar o no los parámetros que se vayan 
#' fijando a cero y los ajustes resultantes.
#'
#' @return Modelo ajustado (si existe) con coeficientes significativos.
#' @export
#'
fit.coefficients <- function(ajuste, orders, alpha=0.05, show_info=T) {
  
  stat <- qnorm(1-alpha/2)                  # estadístico de contraste
  fixed <- rep(NA, length(ajuste$coef))     # valor inicial de los coeficientes
  
  
  while (TRUE) {
    
    # Obtención de los nombres de los coeficientes ARMA
    arma_coefs_names <- get_arma_coefs_names(ajuste)
    
    # En caso de que no haya coeficientes ARIMA, se devuelve el ajuste 
    if (is.null(arma_coefs_names)) { return(ajuste) }
    
    # Obtenemos los valores de dichos coeficientes y sus desviaciones típicas
    arma_coefs <- ajuste$coef[arma_coefs_names]
    arma_coefs_sd <- suppressWarnings(sqrt(diag(ajuste$var.coef[arma_coefs_names, arma_coefs_names])))
    res <- abs(arma_coefs) < stat*arma_coefs_sd
    
    if (!any(res)) { # en caso de que todos sean significativos, se para el bucle
      break
    }
    
    # Configuración del vector fixed
    remov <- names(which.min(abs(arma_coefs)/(stat* arma_coefs_sd)))
    fixed[names(ajuste$coef) == remov] <- 0
    if (show_info) {
      cat(paste0('Es necesario retirar del modelo el parámetro: ', remov, '\n'))
    }
    
    # Actualización de los órdenes
    orders_fixed_update <- update.orders(orders, coefs_names=names(ajuste$coef), fixed=fixed)
    orders <- orders_fixed_update$orders
    fixed <- orders_fixed_update$fixed
    
    # Ajuste del nuevo modelo
    ajuste <- fit.model(ajuste$x, xreg=ajuste$xreg, orders=orders, fixed=fixed)
    
    # Si no se puede optimizar el modelo se devuelve NA
    if (!is_valid(ajuste)) { 
      if (show_info) {warning('No se ha podido optimizar el modelo')}
      return(NA) 
    }
    
    if (show_info) {
      cat('----------------------------------------------------------------------\n')
      print(ajuste, row.names=F)
      cat('\n**********************************************************************\n')
    }
  }
  return(ajuste)
}


#' Actualización correcta de los órdenes del modelo cuando alguno se fija a 0
#' Por ejemplo, cuando el parámetro ma2 de un ARMA(1, 2) se fija a cero, éste 
#' pasa a ser un ARMA(1, 1), por tanto el objeto `orders` y `fixed` deben ser actualizados.
#'
#' @param orders Antiguos órdenes que deben ser actualizados.
#' @param coefs_names Nombres de los coeficientes del modelo (es necesario saber esto 
#' para poder interpretar bien el vector fixed y actualizarlo).
#' @param fixed Máscara para indicar qué coeficientes (`coefs_names`) se van a fijar a 0.
#' 
#' @returns Lista con los órdenes (`orders`) y vector máscara (`fixed`) actualizados.
update.orders <- function(orders, coefs_names, fixed=NULL) {
  
  if (is.null(fixed)) {return(list(orders=orders, fixed=NULL))}
  if (all(is.na(fixed))) {return(list(orders=orders, fixed=fixed))}
  
  
  ar <- str_detect(coefs_names, '^ar\\d+')
  ar_updated <- update.order(orders$regular[1], fixed, ar)
  orders$regular[1] <- ar_updated$order
  ar_fixed <- ar_updated$fixed
  
  
  ma <- str_detect(coefs_names, '^ma\\d+')
  ma_updated <- update.order(orders$regular[3], fixed, ma)
  orders$regular[3] <- ma_updated$order
  ma_fixed <- ma_updated$fixed
  
  
  sar <- str_detect(coefs_names, '^sar\\d+')
  sar_updated <- update.order(orders$seasonal[1], fixed, sar)
  orders$seasonal[1] <- sar_updated$order
  sar_fixed <- sar_updated$fixed
  
  sma <- str_detect(coefs_names, '^sma\\d+')
  sma_updated <- update.order(orders$seasonal[3], fixed, sma)
  orders$seasonal[3] <- sma_updated$order
  sma_fixed <- sma_updated$fixed
  
  
  if (orders$include_constant) {
    mask <- sapply(coefs_names, str_detect, c('mean', 'drift', 'intercept'))
    mask <- apply(mask, 2, sum) == 1
    const_value <- fixed[mask]
    if (length(const_value) == 0) {
        orders$include_constant <- FALSE
    } else if (!is.na(const_value)) {
      orders$include_constant <- FALSE
    }
  }
  
  cnt <- sum(ar, ma, sar, sma)
  
  fixed <- if (cnt>=length(fixed)) c(ar_fixed, ma_fixed, sar_fixed, sma_fixed) else 
    c(ar_fixed, ma_fixed, sar_fixed, sma_fixed, fixed[(1+cnt):length(fixed)])
  return(list(orders=orders, fixed=fixed))
}


update.order <- function(order, fixed, mask) {
  if (sum(mask) == 0 || order==0) {
    return(list(order=order, fixed=c()))
  } 
  
  fixed <- fixed[mask]
  
  if (all(is.na(fixed))) {
    return(list(order=order, fixed=fixed))
  }
  
  for (val in rev(fixed)) {
    if (!is.na(val) && val==0) {
      order <- order - 1
    } else {
      break
    }
  }
  fixed <- if (order==0) c() else fixed[1:order]
  return(list(order=order, fixed=fixed))
}



#' Ajuste de un modelo ARIMA probando de forma iterativa con múltiples optimizadores 
#' para encontrar un ajuste
#'
#' @param serie Serie temporal sobre la que se ajusta el ARIMA.
#' @param orders Órdenes del modelo ARIMA.
#' @param xregs Variables regresoras del modelo.
#' @param fixed Vector máscara con los coeficientes que se deben fijar a 0.
#'
#' @returns Ajuste o `NA` en caso de que no sea posible optimizar un modelo.
#' @export
#'
fit.model <- function(serie, orders, xregs=NULL, fixed=NULL, show_info=F) {
  
  optimizers <- c('BFGS', 'Nelder-Mead', 'CG', 'L-BFGS-B', 'SANN', 'Brent')
  
  for (opt in optimizers) {
    ajuste <- try(
      suppressWarnings(
        Arima(serie, xreg=xregs, order=orders$regular, seasonal=orders$seasonal, 
            fixed=fixed, include.constant = orders$include_constant)),
      silent=T)
    if (is_valid(ajuste)) {
      return(ajuste)
    }
  }
  if (show_info) { warning('No se ha podido optimizar un modelo\n') }
  return(NA)
}


is_valid <- function(ajuste) {
    
    if (length(ajuste) == 1 && class(ajuste) == 'try-error') {
        return(FALSE)
    }
    else if (all(is.na(ajuste))) {
        return(FALSE)
    }
    else if (suppressWarnings(any(is.na(sqrt(diag(ajuste$var.coef)))))) {
        return(FALSE)
    }
    return(TRUE)
}


#' Obtención del lag óptimo de Ljung-Box
#'
#' @param serie 
#'
#' @return
#' @export
#'
#' @examples
ljungbox_lag <- function(serie) {
  m <- frequency(serie)
  n <- length(serie)
  if (m > 1) {
    return(floor(min(2*m, n/5)))
  } else {
    return(floor(min(10, n/5)))
  }
}


#' Obtener los órdenes de un modelo ARIMA a través de una cadena de caracteres
#'
#' @param cadena Cadena de caracteres describiendo el modelo ARIMA (encabezado del 
#' output del objeto Arima). Por ejemplo: ARIMA(1, 0, 2) with non-zero mean.
#' @param reg Valor booleano indicando si se trata de un modelo con o sin variables 
#' regresoras.
#' @param seas Valor booleano indicando si se trata de un modelo con o sin componente 
#' estacional.
parse.orders <- function(cadena, reg=F, seas=F) {
    info <- gsub('[\\(\\)]', '', regmatches(cadena, gregexpr('\\(.*?\\)', cadena))[[1]])
    regular_order <- unlist(lapply(strsplit(info[1], ','), as.double))
    if (length(info) > 1) {
        seasonal_order <- unlist(lapply(strsplit(info[2], ','), as.double))
    } else {seasonal_order <- c(0, 0, 0)}
    
    
    if (grepl('non-zero mean', cadena, fixed=T) || grepl('drift', cadena, fixed=T) || reg) {
        include_constant <- TRUE
    } else {
        include_constant <- FALSE
    }
    
    if (reg && seas) {
        include_constant <- FALSE
    }
    
    orders <- list(regular=regular_order, seasonal=seasonal_order, 
                   include_constant=include_constant)
    return(orders)
}

