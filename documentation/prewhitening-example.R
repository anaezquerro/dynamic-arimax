library(plotly)
library(TSA)
library(fpp2)
library(tseries)

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

alpha <- 0.05


# Generar una figura de ejemplo del prewhitening

load("data/patatas.dat")

Y <- ts(patatas[,1])     # log(sales)
X <- ts(patatas[,2])     # price

fig_sales <- plot_ly(type='scatter', mode='lines') %>% 
    add_trace(y=Y, name='log(sales)', line=list(color=default_colors[1]), showlegend=F) %>%
    layout(
        annotations=list(
            x=length(Y)/2, y=1.1, text='log(sales) evolution', font=list(size=16), showarrow=F, yref='paper'
        )
    )

fig_price <- plot_ly(type='scatter', mode='lines') %>% 
    add_trace(y=X, name='price', line=list(color=default_colors[2]), showlegend=F) %>%
    layout(
        annotations=list(
            x=length(X)/2, y=1.05, text='price evolution', font=list(size=16), showarrow=F, yref='paper'
        )
    )


fig_series <- subplot(fig_sales, fig_price, nrows=2, margin=0.04)


# cross-correlation function
ccfs_obj <- ccf(X, Y, plot=F)
ccfs <- ccfs_obj$acf
lags <- ccfs_obj$lag
stat <- qnorm(1-alpha/2)/sqrt(ccfs_obj$n.used)
sig_lines <- data.frame(x=rep(min(lags), 2), xend=rep(max(lags), 2),
                        y=c(stat, -stat), yend=c(stat, -stat))


fig_original <- plot_ly() %>% 
    add_bars(x=c(lags), y=c(ccfs), type='bar', 
             width=(max(lags) - min(lags))/200, 
             marker=list(color='grey'), showlegend=F) %>%
    add_segments(data=sig_lines, x=~x, xend=~xend, y=~y, yend=~yend,
                 line=list(color='blue', dash='dot'), showlegend=F,
                 name='ccf')  %>% layout(plot_bgcolor='#e5ecf6') %>%
    layout(annotations=list(
        list(x=0, y=1.1, font=list(size=16),
             text='CCF with original series', showarrow=F, yref='paper')

))

series_prewhiten <- TSA::prewhiten(diff(X), diff(Y), plot=F)$ccf
ccfs_prewhiten <- series_prewhiten$acf
lags_prewhiten <- series_prewhiten$lag
stat_prewhiten <- qnorm(1-alpha/2)/sqrt(series_prewhiten$n.used)
sig_lines_prewhiten <- data.frame(x=rep(min(lags_prewhiten), 2), xend=rep(max(lags_prewhiten), 2),
                        y=c(stat_prewhiten, -stat_prewhiten), yend=c(stat_prewhiten, -stat_prewhiten))

fig_prewhiten <- plot_ly() %>% 
    add_bars(x=c(lags_prewhiten), y=c(ccfs_prewhiten), type='bar', 
             width=(max(lags_prewhiten) - min(lags))/200, 
             marker=list(color='grey'), showlegend=F) %>%
    add_segments(data=sig_lines_prewhiten, x=~x, xend=~xend, y=~y, yend=~yend,
                 line=list(color='blue', dash='dot'), showlegend=F,
                 name='ccf')  %>% layout(plot_bgcolor='#e5ecf6') %>%
    layout(annotations=list(
        list(x=0, y=1.1, font=list(size=16),
             text='CCF with <b>prewhitened</b> series', showarrow=F, yref='paper', align='center')
        
))


fig_ccfs <- subplot(fig_original, fig_prewhiten, nrows=2, margin=0.07)

fig <- subplot(fig_series, fig_ccfs, nrows=1, widths=c(0.6, 0.4), margin=0.03)
save_image(fig, file='documentation/figures/example prewhitening.pdf', height=700, width=1300)
fig




