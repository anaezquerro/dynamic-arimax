update.xregs <- function(fitted_model) {
    
    xregs_names <- colnames(fitted_model$xreg)
    xregs_in <- xregs_names[fitted_model$coef[xregs_names] != 0]
    xregs_out <- xregs_names[fitted_model$coef[xregs_names] == 0]
    
    
    fixed <- fitted_model$fixed[
        !str_detect(names(fitted_model$coef), 
                    str_c(xregs_out, collapse='|'))]
    return(
        fit.model(
            fitted_model$x, orders=get.orders(fitted_model), 
            xregs=fitted_model$xreg[, xregs_in], 
            fixed=fixed, show_info=F
        )
    )
    
}