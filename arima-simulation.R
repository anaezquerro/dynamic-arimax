library(polynom)
eval(parse("auto-fit.R", encoding="UTF-8"))

generate.ar.coef <- function(order) {
    if (order <= 0) {stop('Order coefficients must be greater than 0')}
    
    valid <- FALSE 
    
    # Se inicializa un bucle para encontrar de forma aleatoria coeficientes que 
    # cumplan las condiciones de estacionariedad y causalidad
    while (!valid) {
        # Primero se generan los valores (media cero porque pueden ser no significativos)
        if (order > 1) {
            sign <- ifelse(runif(order-1) < 0.2, -1, 1)
            first_values <- sign*rnorm(order-1, 0.3, 0.3)
        } else {
            first_values <- c()
        }
        
        # Después se genera el valor del último coeficiente (que debe ser significativo)
        sign <- ifelse(runif(1) < 0.1, -1, 1)
        last_value <- sign*rnorm(1, 0.8, 0.05) 
        
        # Se concatenan los dos vectores
        ar <- append(first_values, last_value)
        
        # Se comprueban las condiciones de estacionariedad y causalidad
        z <- solve(polynomial(coef=append(1, -ar)), 0)
        
        if (!any(abs(z) <= 1)) {
            valid <- TRUE
        }
    }
    return(ar)  
}

generate.ma.coef <- function(order) {
    if (order <= 0) {stop('MA order must be grater than 0')}
    valid <- FALSE
    
    while (!valid) {
        if (order > 1) {
            sign <- ifelse(runif(order-1) < 0.8, -1, 1)
            first_values <- sign*rnorm(order-1, 0.3, 0.07)
        } else {
            first_values <- c()
        }
        
        sign <- ifelse(runif(1) < 0.1, -1, 1)
        last_value <- sign*rnorm(1, 0.6, 0.05)
        
        ma <- append(first_values, last_value)
        
        z <- solve(polynomial(coef=append(1, ma)), 0)
        
        if (!any(abs(z) <= 1)) {
            valid <- TRUE
        }
    }
    return(ma)
}


comb <- function(n, x) {
    return(factorial(n) / factorial(n-x) / factorial(x))
}

B.binom <- function(X, d, t) {
    binom_coeff <- comb(d, 1:d)
    Xwindow <- rev(window(X, start=t-d, end=t-1))
    return(sum(binom_coeff*((-Xwindow)**(1:d))))
}


sim.arima <- function(model=list(p=0, d=0, q=0), n=1000, with.constant=FALSE, 
                      x_ini_mean=0) {
    # Corrección de los valores de model
    if ('ar' %in% names(model)) {
        model$p <- length(model[['ar']])
    }
    if ('ma' %in% names(model)) {
        model$q <- length(model[['ma']])
    }
    if (!'d' %in% names(model)) {
        model$d <- 0
    }
    
    if (with.constant && (model$d > 0)) {stop('Si d > 0 no puede haber constante')}
    
    # Generación del proceso de ruido blanco
    a <- ts(rnorm(n+model$q, 0, 0.05), start=-model$q, end=n-1)
    
    # Generación de una constante
    c <- ifelse(with.constant, rnorm(1, x_ini_mean, 0.1), 0)
    
    # En caso de que no haya órdenes, devolver ruido gaussiano
    if ((model$p==0) && (model$q==0) && (model$d==0)) {
        x <- c + a
        return(list(X=x, ar=NULL, ma=NULL, c=c, d=0))
    } 
    
    # Generación de las muestras iniciales de X (aleatorias)
    if (model$p > 0) {
        X <- ts(rnorm(max(c(model$p, model$d)), 0, 0.5), end=-1)
        # X <- ts(rep(0, max(c(model$p, model$d))), end=-1)
    } else {
        X <- c()
    }
    
    # Si order no tiene el parámetro ar (coeficientes), se generan
    if (!'ar' %in% names(model) && (model$p>0)) {
        model$ar <- generate.ar.coef(model$p)
    }
    
    if (model$p > 0) {
        apply.ar <- function(X, t) {
            Xwindow <- rev(window(X, start=t-model$p, end=t-1))
            return(sum(model$ar*Xwindow))
        }
    } else { apply.ar <- function(X, t) {return(0)} }
    
    # Si order no tiene el parámetro ma (coeficientes), se generan
    if (!'ma' %in% names(model) && (model$q > 0)) {
        model$ma <- generate.ma.coef(model$q)
    }
    
    
    if (model$q > 0) {
        apply.ma <- function(a, t) {
            awindow <- rev(window(a, start=t-model$q, end=t-1))
            return(sum(model$ma*awindow))
        }
    } else { apply.ma <- function(x, t) {return(0)} }
    
    
    # Se van actualizando los valores
    for (t in 0:(n-1)) {
        ar_part <- apply.ar(X, t)
        ma_part <- apply.ma(a, t)
        
        Xt <- c +  ar_part + ma_part + window(a, start=t, end=t)
        X <- ts(append(X, Xt), end=t)
    }
    
    X <- window(X, start=0, end=n-1)
    
    # Como X es un proceso ARMA(p, q), tenemos que obtener Z tal que X_t = Z_t - Z_{t-1}
    if (model$d > 0) {
        for (i in 1:model$d) {
            Z <- ts(rnorm(1, 0, 0.01), start=-1)
            for (t in 0:(n-1)) {
                Zt <-  window(X, start=t, end=t) - B.binom(Z, 1, t)
                Z <- ts(append(Z, Zt), start=-1, end=t)
            }
            X <- window(Z, start=0, end=n-1) 
        }
    }
    
    result <- list(ar=model$ar, ma=model$ma, c=c, X=X, d=model$d)
    return(result)
    
    
    # Comprobación
    # ajuste <- auto.fit.arima(result$X, show_info=F)
    # 
    # if (is_valid(ajuste)) {
    #     valid <- (model$p == ajuste$arma[1]) & (model$q == ajuste$arma[2]) & (model$d == ajuste$arma[6]) & 
    #         (('intercept' %in% names(ajuste$coef)) == with.constant)
    # } else {
    #     valid <- F
    # }
    # 
    # 
    # if (valid) {
    #     return(result)
    # } else {
    #     return(sim.arima(model, n, with.constant, x_ini_mean))
    # }
    
}