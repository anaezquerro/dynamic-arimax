library(plotly)




#' Visualización de la información básica sobre series temporales.
#'
#' @param serie : Serie temporal de entrada.
#' @param title : Título del gráfica secuencial (por defecto, "Gráfico secuencial")
#' @param alpha : Significance level of ACF and PACF barplot.
#'
#' @return `plotly` figure with sequential plot and ACF and PACF barplots.
#' @export
plot_serie <- function(serie, title='Gráfico secuencial', alpha=0.05) {
  if (class(serie) != 'ts') { stop('El parámetro serie debe ser un objeto de tipo ts.') }
  
  # Gráfico secuencial
  seq_plot <- plot_ly(type='scatter', mode='lines') %>%
    add_trace(x=time(serie), y=serie, 
              line=list(color='grey'), name=title, showlegend=F) %>%
    layout(xaxis=list(zeroline=F), yaxis=list(zeroline=F), 
           plot_bgcolor='#e5ecf6', separators='.',
           annotations=list(
             list(x=0.5, y=1.1, text=title, showarrow=F, margin=list(t=50),
                  xref='paper', yref='paper', font=list(size=16))
           )
    )
  
  # Gráfico de correlaciones simples
  acfs <- acf(serie, plot=F)
  stat <- qnorm(1-alpha/2)/sqrt(acfs$n.used)
  sig_lines <- data.frame(x=c(0, 0), xend=rep(max(acfs$lag), 2), 
                          y=c(stat, -stat), yend=c(stat, -stat))
  
  acf_plot <- plot_ly() %>% 
    add_bars(x=c(acfs$lag), y=c(acfs$acf), type='bar', 
             width=(max(acfs$lag) - min(acfs$lag))/100, 
             marker=list(color='grey')) %>%
    add_segments(data=sig_lines, x=~x, xend=~xend, y=~y, yend=~yend,
                 line=list(color='blue', dash='dot'), showlegend=F,
                 name='acf') %>% 
    layout(plot_bgcolor='#e5ecf6',
           annotations=list(
             list(x=max(acfs$lag)/2, y=1.1, text='Autocorrelaciones', 
                  showarrow=F, yref='paper', font=list(size=16))
           )
    )
  
  # Gráfico de correlaciones parciales
  pacfs <- pacf(serie, plot=F)
  pstat <- qnorm(1-alpha/2)/sqrt(pacfs$n.used)
  psig_lines <- data.frame(x=c(0, 0), xend=rep(max(pacfs$lag), 2),
                           y=c(pstat, -pstat), yend=c(pstat, -pstat))
  
  pacf_plot <- plot_ly() %>% 
    add_bars(x=c(pacfs$lag), y=c(pacfs$acf), type='bar', 
             width=(max(pacfs$lag) - min(pacfs$lag))/100, 
             marker=list(color='grey')) %>%
    add_segments(data=psig_lines, x=~x, xend=~xend, y=~y, yend=~yend,
                 line=list(color='blue', dash='dot'), showlegend=F,
                 name='pacf')  %>% 
    layout(plot_bgcolor='#e5ecf6',
           annotations=list(
             list(x=max(pacfs$lag)/2, y=1.1, text='Autocorrrelaciones parciales', 
                  showarrow=F, yref='paper', font=list(size=16))
           )
    )
  
  subfig <- subplot(acf_plot, pacf_plot, margin=0.07)
  
  fig <- subplot(seq_plot, subfig, nrows=2, margin=0.07) %>%
    layout(showlegend=FALSE, showlegend2=FALSE)
  
  return(fig)
}




plot_residuals <- function(ajuste, title='Gráfico secuencia de los residuos', alpha=0.05) {
  
  # Gráfico distribución de los residuos
  bounds <- sqrt(ajuste$sigma2)*3
  x_grid <- seq(from=-bounds, to=bounds, by=2*bounds/1000)
  
  hist_residuals <- ggplot(data=NULL, aes(x=c(ajuste$residuals))) + 
    geom_histogram(aes(y=..density..), fill='lightblue', color='grey', binwidth=0.1) +
    geom_rug(color='grey') + 
    geom_line(aes(x=x_grid, y=dnorm(x_grid, sd=sqrt(ajuste$sigma2))), 
              color='orange')
  
  hist_residuals <- ggplotly(hist_residuals) %>%
    layout(xaxis=list(title='residuos', zeroline=F), 
           yaxis=list(title='densidad', zeroline=F),  plot_bgcolor='#e5ecf6',
           annotations=list(
             list(x=0, y=1.1, text='Test de normalidad', showarrow=F, yref='paper')
           )
    )
  
  # Gráfico secuencial de los residuos
  seq_plot <-  plot_ly(type='scatter', mode='lines') %>%
    add_trace(x=time(ajuste$residuals), y=ajuste$residuals, 
              line=list(color='grey'), name=title, showlegend=F) %>%
    layout(xaxis=list(zeroline=F), yaxis=list(zeroline=F), plot_bgcolor='#e5ecf6',
           annotations=list(
             list(x=0.5, y=1.1, text=title, showarrow=F, xref='paper', yref='paper',
                  margin=list(t=50))
           )
    )
  
  # Gráfico de las autocorrelaciones simples
  acfs <- acf(ajuste$residuals, plot=F)
  stat <- qnorm(1-alpha/2)/sqrt(acfs$n.used)
  sig_lines <- data.frame(x=c(0, 0), xend=rep(max(acfs$lag), 2),
                          y=c(stat, -stat), yend=c(stat, -stat))
  acf_plot <- plot_ly() %>% 
    add_bars(x=c(acfs$lag), y=c(acfs$acf), type='bar', 
             width=(max(acfs$lag) - min(acfs$lag))/200, 
             marker=list(color='grey')) %>%
    add_segments(data=sig_lines, x=~x, xend=~xend, y=~y, yend=~yend,
                 line=list(color='blue', dash='dot'), showlegend=F,
                 name='acf')  %>% layout(plot_bgcolor='#e5ecf6') %>%
    layout(annotations=list(
      list(x=max(acfs$lag)/2, y=1.1, text='Autocorrelaciones', showarrow=F, yref='paper')
    ))
  
  subfig <- subplot(acf_plot, hist_residuals, margin=0.0)
  
  fig <- subplot(seq_plot, subfig, nrows=2, margin=0.07) %>%
    layout(showlegend=FALSE, showlegend2=FALSE)
  
  return(fig)
  
}


plot_prewhiten <- function(prewhiten_object, alpha=0.05) {
    lags <- c(prewhiten_object$ccf$lag)
    values <- c(prewhiten_object$ccf$acf)
    n <- prewhiten_object$ccf$n.used
    
    stat <- qnorm(1-alpha/2)/sqrt(n)
    
    sig_lines <- data.frame(x=rep(min(lags), 2), xend=rep(max(lags), 2),
                            y=c(stat, -stat), yend=c(stat, -stat))
    
    corr_plot <- plot_ly() %>%
        add_bars(x=lags, y=values, type='bar', width=(max(lags)-min(lags))/100,
                 marker=list(color='grey'), name='acf') %>%
        add_segments(data=sig_lines, x=~x, xend=~xend, y=~y, yend=~yend,
                     line=list(color='blue', dash='dot'), showlegend=F) %>%
        layout(plot_bgcolor='#e5ecf6', title='Gráfico de correlaciones', 
               xaxis=list(title='lags'), yaxis=list(title='correlations'))
    return(corr_plot)
}

default_colors <- c(
    '#1f77b4',      # or rgb(31, 119, 180)  // muted blue
    '#ff7f0e',      # or rgb(255, 127, 14)  // safety orange
    '#2ca02c',      # or rgb(44, 160, 44)   // cooked asparagus green
    '#d62728',      # or rgb(214, 39, 40)   // brick red
    '#9467bd',      # or rgb(148, 103, 189) // muted purple
    '#8c564b',      # or rgb(140, 86, 75)   // chestnut brown
    '#e377c2',      # or rgb(227, 119, 194)
    '#7f7f7f',      # or rgb(127, 127, 127)
    '#bcbd22',      # or rgb(188, 189, 34)
    '#17becf'       # or rgb(23, 190, 207)
)

pplot <- function(x, y=NULL, title="", color=default_colors[1]) {
  
  if (is.null(y)) {
    y <- x
    x <- 1:length(y)
  }
  
  fig <- plot_ly(type='scatter', mode='lines')%>%
    add_trace(x=x, y=y, line=list(color=color), name=title, showlegend=F) %>%
    layout(xaxis=list(zeroline=F), yaxis=list(zeroline=F),
           plot_bgcolor='#e5ecf6', separators='.')
  
  return(fig)
}

