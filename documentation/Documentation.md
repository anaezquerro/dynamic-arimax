# Function for automatic fitting ARIMA or ARIMAX models (`auto.fit.arima`)

**Description**: Obtains a valid fitted ARIMA or ARIMAX model that
optimizes a selected information criterion and satisfies that:

1.  All coefficients are statistically significtive.
2.  Model residuals have zero mean and are independent.

**Returns**:

1.  Fitted ARIMA or ARIMAX model (`Arima` object) if it exists and it
    can be optimized.
2.  `NA` in case it does not exists and/or cannnot be optimized.

``` r
auto.fit.arima(serie, xregs = NULL,  ic = c("aicc", "aic", "bic"), d = NA, 
               D = NA, alpha = 0.05, show_info = TRUE, plot_result = FALSE)
```

**Arguments**:

-   `serie` \[`ts`\]: Univariate time series to fit an ARIMA model with.
    If `xregs` is not `NULL` it will act as the dependent variable.
-   `xregs` \[`mts`\]: Matrix of time series that act as covariates in
    an ARIMAX model. By default, `NULL`, meaning that an ARIMA will be
    fitted with `serie`.
-   `ic` \[`character`\]: Information criterion to be used in ARIMA
    orders’ selection. Options available are the same as in
    `forecast::auto.arima()` function: `ic`, `bic` or `aicc`.
-   `d` \[`numeric`\]: Value of the diferentiation order. By default
    `NA`, thus there is no restriction about the *d* order value.
-   `D` \[`numeric`\]: Value fo the seasonal diferentiation order. By
    default `NA`, thus there is no restriction about the *D* order
    value.
-   `alpha` \[`numeric`\]: Significance level of hypothesis tests used
    for checking independence, zero mean and normality of residuals, and
    significance of estimated coefficients.
-   `show_info` \[`boolean`\]: Displaying or not displaying the
    historical of fitted models.
-   `plot_result` \[`boolean`\]: Returning or not with the fitted model
    a sequential plot of time series and residuals.

**Observations**:

-   To check residuals independence the Ljung-Box (`Box.test()`) is
    used. The number of lags is selected from the length of the time
    serie and testing the stationarity of it (see function
    `ljunbox_lag()`).
-   To check null mean of residuals the `t.test()` is used.
-   To check normality of residuals we use the Jarque-Bera
    (`jarque.bera.test()`) and the Shapiro-Wilks tests
    (`shapiro.test()`).
-   Considered models will always have a regular differentiation order
    less or equal to 3 (*d* ≤ 3) and a seasonal differentiation order
    less or equal to 2 (*D* ≤ 2).

**Example**: Influenza evolution in Catalonia.

``` r
dat <- read.csv("data/evolucion_gripe_covid.csv")               # read data
sdgripal <- ts(dat$sdgripal, start=c(2020, 40), frequency=52)   # convert to ts
model_gripal <- auto.fit.arima(sdgripal, plot_result = TRUE)
```

    --------------------------------------------------------------------------------------
    Series: serie 
    ARIMA(1,1,0) 

    Coefficients:
             ar1
          0.2452
    s.e.  0.1089

    sigma^2 = 46663:  log likelihood = -529.48
    AIC=1062.96   AICc=1063.12   BIC=1067.68
    --------------------------------------------------------------------------------------
    Normality hypothesis is rejected
    Model is valid but forecasting asuming normality is not available
    --------------------------------------------------------------------------------------
    --------------------------------------------------------------------------------------
    |                                    FINAL MODEL                                     |
    --------------------------------------------------------------------------------------
    Series: serie 
    ARIMA(1,1,0) 

    Coefficients:
             ar1
          0.2452
    s.e.  0.1089

    sigma^2 = 46663:  log likelihood = -529.48
    AIC=1062.96   AICc=1063.12   BIC=1067.68

``` r
# display series sequential plot
display(model_gripal$fig_serie, "serie gripe", width=1000, height=800)
```

``` r
# display residuals sequential plot
display(model_gripal$fig_residuals, "residuals gripe", width=1000, height=800)
```

**Example**: [Monthly level of Co2 measured in the Observatory of Manua
Loa (Hawaii)](ftp://ftp.cmdl.noaa.gov/ccg/co2/trends/co2_mm_mlo.txt).
Time series starts in March 1958.

``` r
# read and parse data
co2 <- ts(scan('data/co2MaunaLoa.dat'), start=c(1958, 3), frequency=12)

# apply the automatic fitting
model_co2 <- auto.fit.arima(co2, ic="aicc", plot_result=TRUE)
```

    --------------------------------------------------------------------------------------
    Series: serie 
    ARIMA(1,1,1)(0,1,1)[12] 

    Coefficients:
             ar1      ma1     sma1
          0.1645  -0.5210  -0.8684
    s.e.  0.1048   0.0909   0.0208

    sigma^2 = 0.09136:  log likelihood = -135.78
    AIC=279.57   AICc=279.63   BIC=297.23
    --------------------------------------------------------------------------------------
    Removing parameter ar1 since it is not significative
    --------------------------------------------------------------------------------------
    Series: serie 
    ARIMA(0,1,1)(0,1,1)[12] 

    Coefficients:
              ma1     sma1
          -0.3783  -0.8684
    s.e.   0.0415   0.0209

    sigma^2 = 0.09155:  log likelihood = -136.93
    AIC=279.87   AICc=279.91   BIC=293.11
    --------------------------------------------------------------------------------------
    Normality hypothesis is rejected
    Model is valid but forecasting asuming normality is not available
    --------------------------------------------------------------------------------------
    --------------------------------------------------------------------------------------
    |                                    FINAL MODEL                                     |
    --------------------------------------------------------------------------------------
    Series: serie 
    ARIMA(0,1,1)(0,1,1)[12] 

    Coefficients:
              ma1     sma1
          -0.3783  -0.8684
    s.e.   0.0415   0.0209

    sigma^2 = 0.09155:  log likelihood = -136.93
    AIC=279.87   AICc=279.91   BIC=293.11

``` r
# display sequential plot of the time series
display(model_co2$fig_serie, "serie co2", width=1000, height=800)
```

``` r
# display sequential plot of the residuals
display(model_co2$fig_residuals, "residuals co2", width=1000, height=800)
```

# Function for automatic covariates selection and lags in Dynamic Regression models (`drm.select`)

**Description**: Implementation of our new approach of covariates
selection in dynamic regression models where each covariate is
introduced lagged *k* ≥ 0 moments. For more information about this
method, check the paper attached with this file.

**Returns**:

-   \[`Arima`\] Fitted ARIMAX model with significant covariates
    selected.
-   \[`Arima`\] Fitted ARIMA model with the dependent variable in case
    no covariate is selected.

``` r
drm.select(serie, xregs, ic=c('aic', 'bic', 'aicc', alpha=0.05, 
                              st_method='auto.arima', show_info=T, ndiff)
```

**Arguments**:

-   `serie` \[`ts`\]: Univariate time series which acts as the dependent
    variable.
-   `xregs` \[`mts`\]: Set of covariate candidates to model the behavior
    of `serie`.
-   `ic` \[`character`\]: Information criterion used for covariate
    selection and model comparisons. The possibilities are: `aic`, `bic`
    and `aicc`. By default, th AICc is used.
-   `alpha` \[`numeric`\]: Signficance level for hypothesis tests of
    stationary and coefficients significance. By default it is set to
    0.05.
-   `st_method` \[`character`\]: Method used for checking stationary of
    a time series. If it is `auto.arima`, the function
    `forecast::auto.arima()` is used and the differentiation order is
    checked. If it is `adf.test`, the Dickey-Fuller test is used.
-   `show_info` \[`boolean`\]: Displaying or not displaying the
    historical of added covariates.
-   `ndiff` \[`numeric`\]: Internal argument (**do not use**) to apply
    regular differentiations to all data when no ARIMAX model can be
    fitted with stationary errors.

*Note*: Information about inidividual model fitting of each covariate
will not be displayed.

**Example**: Logarithm of weekly sales and *Bluebird* chips price in New
Zeland. The observation period is 104 weeks long (from 20/09/1988 to
10/09/2000).

``` r
load("data/patatas.dat")                    # read data
sales <- patatas[,1]
xregs <- ts(matrix(patatas[,2]))
colnames(xregs) <- c('price')

chips_model <- drm.select(sales, xregs)     # apply selection method
```

    Covariate price has been tested [ic=-71.4371307570267, lag=0]
    Covariate price has been added [aicc=-71.4371307570267, lag=0]
    Series: serie 
    Regression with ARIMA(0,0,4) errors 

    Coefficients:
          ma1     ma2  ma3     ma4  intercept     xreg
            0  0.2884    0  0.5416    15.8559  -2.4682
    s.e.    0  0.0794    0  0.1167     0.1909   0.1100

    sigma^2 = 0.02728:  log likelihood = 41.02
    AIC=-72.05   AICc=-71.44   BIC=-58.83
    --------------------------------------------------------------------------------------
    No more variables will be added
    --------------------------------------------------------------------------------------
    |               Historical of added covariates to the model (ndiff=0)                |
    --------------------------------------------------------------------------------------
       var lag                ic
     price   0 -71.4371307570267
    --------------------------------------------------------------------------------------
    Series: serie 
    Regression with ARIMA(0,0,4) errors 

    Coefficients:
          ma1     ma2  ma3     ma4  intercept     xreg
            0  0.2884    0  0.5416    15.8559  -2.4682
    s.e.    0  0.0794    0  0.1167     0.1909   0.1100

    sigma^2 = 0.02728:  log likelihood = 41.02
    AIC=-72.05   AICc=-71.44   BIC=-58.83

**Example**: Time series about Microsoft *stock*.

``` r
microsoft <- read.csv('data/microsoft-stock.csv')     # read data
close_price <- ts(microsoft$Close)                    # dependent variable
xregs <- ts(microsoft[, c('Open', 'High', 'Low', 'Volume')])

# fit model
fitted_model <- drm.select(close_price, xregs)
```

    Covariate Open has been tested [ic=5923.03932105522, lag=0]
    The optimizer could not fit a model for High
    The optimizer could not fit a model for Low
    The optimizer could not fit a model for Volume
    Covariate Open has been added [aicc=5923.03932105522, lag=0]
    Series: serie 
    Regression with ARIMA(0,0,0) errors 

    Coefficients:
            xreg
          1.0002
    s.e.  0.0004

    sigma^2 = 2.953:  log likelihood = -2959.52
    AIC=5923.03   AICc=5923.04   BIC=5933.67
    --------------------------------------------------------------------------------------
    Covariate High has been tested [ic=5920.18751985925, lag=-1]
    Covariate Low has been tested [ic=5923.03932105528, lag=-1]
    The optimizer could not fit a model for Volume
    Covariate High has been added [aicc=5920.18751985925, lag=-1]
    Series: serie 
    Regression with ARIMA(0,0,0) errors 

    Coefficients:
            Open    High
          0.9476  0.0522
    s.e.  0.0239  0.0237

    sigma^2 = 2.945:  log likelihood = -2957.09
    AIC=5920.17   AICc=5920.19   BIC=5936.13
    --------------------------------------------------------------------------------------
    Covariate Low has been tested [ic=5920.18757955067, lag=-1]
    The optimizer could not fit a model for Volume
    No more variables will be added
    --------------------------------------------------------------------------------------
    |               Historical of added covariates to the model (ndiff=0)                |
    --------------------------------------------------------------------------------------
      var lag               ic
     Open   0 5923.03932105522
     High  -1 5920.18751985925
    --------------------------------------------------------------------------------------
    Series: serie 
    Regression with ARIMA(0,0,0) errors 

    Coefficients:
            Open    High
          0.9476  0.0522
    s.e.  0.0239  0.0237

    sigma^2 = 2.945:  log likelihood = -2957.09
    AIC=5920.17   AICc=5920.19   BIC=5936.13

**Example**: Exitus cases in Spain due to COVID19 virus, considering the
following covariates:

-   Confirmed cases and recovered cases in Spain.
-   Confirmed, exitus and recovered cases in France.
-   Confirmed, exitus and recovered cases in Portugal.
-   Confirmed, exitus and recovered cases in England.

``` r
# read data
confirmed <- read.csv("data/covid-global-confirmed-bycountry.csv")
deaths <- read.csv("data/covid-global-deaths-bycountry.csv")
recovered <- read.csv("data/covid-global-recovered-bycountry.csv")

# Spain data
deaths_spain <- ts(deaths$Spain, frequency=7)
confirmed_spain <- ts(confirmed$Spain, frequency=7)
recovered_spain <- ts(recovered$Spain, frequency=7)

# France data
deaths_france <- ts(deaths$France, frequency=7)
confirmed_france <- ts(confirmed$France, frequency=7)
recovered_france <- ts(recovered$France, frequency=7)

# England data
deaths_england <- ts(deaths$United.Kingdom, frequency=7)
confirmed_england <- ts(confirmed$United.Kingdom, frequency=7)
recovered_england <- ts(recovered$United.Kingdom, frequency=7)

# Portugal data
deaths_portugal <- ts(deaths$Portugal, frequency=7)
confirmed_portugal <- ts(confirmed$Portugal, frequency=7)
recovered_portugal <- ts(recovered$Portugal, frequency=7)

xregs <- ts(cbind(confirmed_spain, recovered_spain, 
                  deaths_france, confirmed_france, recovered_france,
                  deaths_england, confirmed_england, recovered_england,
                  deaths_portugal, confirmed_portugal, recovered_portugal),
            frequency=7)

# fit model
fitted_model <- drm.select(deaths_spain, xregs)
```

    Covariate confirmed_spain has been tested [ic=5596.64110211836, lag=0]
    Covariate recovered_spain has been tested [ic=5620.75355779798, lag=-7]
    Covariate deaths_france has been tested [ic=5644.19614523882, lag=-8]
    Covariate confirmed_france has been tested [ic=5622.6940170394, lag=0]
    Covariate recovered_france has been tested [ic=5629.78160439558, lag=0]
    Covariate deaths_england has been tested [ic=5659.55354812883, lag=-15]
    Covariate confirmed_england has been tested [ic=5659.55354812883, lag=-5]
    Significative correlation with lag<=0 could not be found for recovered_england
    Covariate deaths_portugal has been tested [ic=5659.55354812883, lag=-21]
    Significative correlation with lag<=0 could not be found for confirmed_portugal
    Covariate recovered_portugal has been tested [ic=5630.43143861607, lag=-1]
    Covariate confirmed_spain has been added [aicc=5596.64110211836, lag=0]
    Series: serie 
    Regression with ARIMA(0,2,1)(1,0,1)[7] errors 

    Coefficients:
              ma1    sar1     sma1    xreg
          -0.7846  0.9035  -0.7796  0.0074
    s.e.   0.0292  0.0770   0.1169  0.0011

    sigma^2 = 33006:  log likelihood = -2793.25
    AIC=5596.5   AICc=5596.64   BIC=5616.72
    --------------------------------------------------------------------------------------
    Covariate recovered_spain has been tested [ic=5588.82493175822, lag=-7]
    Covariate deaths_france has been tested [ic=5604.70104752933, lag=-8]
    Covariate confirmed_france has been tested [ic=5582.48568718243, lag=0]
    The optimizer could not fit a model for recovered_france
    Covariate deaths_england has been tested [ic=5580.26620297165, lag=-23]
    Significative correlation with lag<=0 could not be found for confirmed_england
    Significative correlation with lag<=0 could not be found for recovered_england
    Covariate deaths_portugal has been tested [ic=5596.64109720114, lag=-21]
    Covariate confirmed_portugal has been tested [ic=5596.64115395815, lag=0]
    Covariate recovered_portugal has been tested [ic=5575.32330716563, lag=-1]
    Covariate recovered_portugal has been added [aicc=5575.32330716563, lag=-1]
    Series: serie 
    Regression with ARIMA(0,2,1)(1,0,1)[7] errors 

    Coefficients:
              ma1    sar1     sma1  confirmed_spain  recovered_portugal
          -0.7522  0.9283  -0.8180           0.0067             -0.0319
    s.e.   0.0307  0.0561   0.0925           0.0010              0.0066

    sigma^2 = 31288:  log likelihood = -2781.56
    AIC=5575.12   AICc=5575.32   BIC=5599.39
    --------------------------------------------------------------------------------------
    The optimizer could not fit a model for recovered_spain
    The optimizer could not fit a model for deaths_france
    Covariate confirmed_france has been tested [ic=5585.0827651252, lag=-9]
    The optimizer could not fit a model for recovered_france
    The optimizer could not fit a model for deaths_england
    Significative correlation with lag<=0 could not be found for confirmed_england
    Significative correlation with lag<=0 could not be found for recovered_england
    Covariate deaths_portugal has been tested [ic=5575.32353074263, lag=-21]
    Covariate confirmed_portugal has been tested [ic=5575.3235234709, lag=0]
    No more variables will be added
    The global model does not have stationary errors
    Trying to adjust a model that do have stationary errors
    No valid model with stationary errors could be optimized
    Applying regular differentiation (ndiff=1) and calling again the function
    --------------------------------------------------------------------------------------
    --------------------------------------------------------------------------------------
    Covariate confirmed_spain has been tested [ic=5596.64110496284, lag=0]
    Covariate recovered_spain has been tested [ic=5620.75356142264, lag=-7]
    Covariate deaths_france has been tested [ic=5644.19614774329, lag=-8]
    Covariate confirmed_france has been tested [ic=5646.00774268838, lag=0]
    Covariate recovered_france has been tested [ic=5629.781607806, lag=0]
    Covariate deaths_england has been tested [ic=5644.19614774329, lag=-15]
    Covariate confirmed_england has been tested [ic=5644.19614774329, lag=-5]
    Significative correlation with lag<=0 could not be found for recovered_england
    Covariate deaths_portugal has been tested [ic=5644.19614774329, lag=-21]
    Significative correlation with lag<=0 could not be found for confirmed_portugal
    Covariate recovered_portugal has been tested [ic=5652.5499950804, lag=-1]
    Covariate confirmed_spain has been added [aicc=5596.64110496284, lag=0]
    Series: serie 
    Regression with ARIMA(0,1,1)(1,0,1)[7] errors 

    Coefficients:
              ma1    sar1     sma1    xreg
          -0.7846  0.9035  -0.7796  0.0074
    s.e.   0.0292  0.0770   0.1169  0.0011

    sigma^2 = 33006:  log likelihood = -2793.25
    AIC=5596.5   AICc=5596.64   BIC=5616.72
    --------------------------------------------------------------------------------------
    Covariate recovered_spain has been tested [ic=5588.82508685332, lag=-7]
    Covariate deaths_france has been tested [ic=5596.64108624086, lag=-8]
    Covariate confirmed_france has been tested [ic=5582.4857812052, lag=0]
    Covariate recovered_france has been tested [ic=5593.9519575025, lag=-18]
    Covariate deaths_england has been tested [ic=5572.16506296523, lag=-23]
    Significative correlation with lag<=0 could not be found for confirmed_england
    Significative correlation with lag<=0 could not be found for recovered_england
    Covariate deaths_portugal has been tested [ic=5596.64110005419, lag=-21]
    Covariate confirmed_portugal has been tested [ic=5596.64115681075, lag=0]
    Covariate recovered_portugal has been tested [ic=5575.32307482294, lag=-1]
    Covariate deaths_england has been added [aicc=5572.16506296523, lag=-23]
    Series: serie 
    Regression with ARIMA(0,1,1)(1,0,1)[7] errors 

    Coefficients:
              ma1    sar1     sma1  confirmed_spain  deaths_england
          -0.7846  0.9039  -0.7801           0.0074               0
    s.e.   0.0292  0.0772   0.1174           0.0011               0

    sigma^2 = 33164:  log likelihood = -2781.01
    AIC=5572.02   AICc=5572.17   BIC=5592.22
    --------------------------------------------------------------------------------------
    Covariate recovered_spain has been tested [ic=5564.47502910463, lag=-7]
    Covariate deaths_france has been tested [ic=5572.165048401, lag=-8]
    Covariate confirmed_france has been tested [ic=5558.08935421137, lag=0]
    Covariate recovered_france has been tested [ic=5569.49147206183, lag=-18]
    Significative correlation with lag<=0 could not be found for confirmed_england
    Significative correlation with lag<=0 could not be found for recovered_england
    Covariate deaths_portugal has been tested [ic=5572.16510073314, lag=-21]
    Covariate confirmed_portugal has been tested [ic=5572.16533921363, lag=0]
    Covariate recovered_portugal has been tested [ic=5550.9630503792, lag=-1]
    Covariate recovered_portugal has been added [aicc=5550.9630503792, lag=-1]
    Series: serie 
    Regression with ARIMA(0,1,1)(1,0,1)[7] errors 

    Coefficients:
              ma1    sar1     sma1  confirmed_spain  deaths_england
          -0.7521  0.9290  -0.8191           0.0067               0
    s.e.   0.0308  0.0561   0.0927           0.0010               0
          recovered_portugal
                     -0.0320
    s.e.              0.0066

    sigma^2 = 31438:  log likelihood = -2769.38
    AIC=5550.76   AICc=5550.96   BIC=5575
    --------------------------------------------------------------------------------------
    The optimizer could not fit a model for recovered_spain
    Covariate deaths_france has been tested [ic=5550.96307921146, lag=-8]
    Covariate confirmed_france has been tested [ic=5558.87899146256, lag=-9]
    Covariate recovered_france has been tested [ic=5540.65497098646, lag=0]
    Significative correlation with lag<=0 could not be found for confirmed_england
    Significative correlation with lag<=0 could not be found for recovered_england
    Covariate deaths_portugal has been tested [ic=5550.96396788424, lag=-21]
    Covariate confirmed_portugal has been tested [ic=5550.96301198882, lag=0]
    Covariate recovered_france has been added [aicc=5540.65497098646, lag=0]
    Series: serie 
    Regression with ARIMA(0,1,1)(1,0,1)[7] errors 

    Coefficients:
              ma1    sar1     sma1  confirmed_spain  deaths_england
          -0.7652  0.9102  -0.8201           0.0068               0
    s.e.   0.0310  0.0951   0.1359           0.0010               0
          recovered_portugal  recovered_france
                     -0.0314            0.0537
    s.e.              0.0065            0.0153

    sigma^2 = 30647:  log likelihood = -2763.19
    AIC=5540.38   AICc=5540.65   BIC=5568.66
    --------------------------------------------------------------------------------------
    Covariate recovered_spain has been tested [ic=5526.18290825366, lag=-7]
    Covariate deaths_france has been tested [ic=5540.6551615037, lag=-8]
    Covariate confirmed_france has been tested [ic=5522.69957319275, lag=0]
    Significative correlation with lag<=0 could not be found for confirmed_england
    Significative correlation with lag<=0 could not be found for recovered_england
    Covariate deaths_portugal has been tested [ic=5540.65503851333, lag=-21]
    Covariate confirmed_portugal has been tested [ic=5540.65497703673, lag=0]
    Covariate confirmed_france has been added [aicc=5522.69957319275, lag=0]
    Series: serie 
    Regression with ARIMA(1,1,2)(1,0,1)[7] errors 

    Coefficients:
              ar1  ma1      ma2    sar1     sma1  confirmed_spain  deaths_england
          -0.7877    0  -0.5540  0.9241  -0.8338           0.0065               0
    s.e.   0.0399    0   0.0514  0.0765   0.1149           0.0010               0
          recovered_portugal  recovered_france  confirmed_france
                     -0.0296            0.0635           -0.0034
    s.e.              0.0063            0.0149            0.0007

    sigma^2 = 29201:  log likelihood = -2752.13
    AIC=5522.26   AICc=5522.7   BIC=5558.62
    --------------------------------------------------------------------------------------
    Covariate recovered_spain has been tested [ic=5518.39734277427, lag=-7]
    Covariate deaths_france has been tested [ic=5522.47431256101, lag=-8]
    Significative correlation with lag<=0 could not be found for confirmed_england
    Significative correlation with lag<=0 could not be found for recovered_england
    Covariate deaths_portugal has been tested [ic=5522.69927751652, lag=-21]
    Covariate confirmed_portugal has been tested [ic=5504.16966367854, lag=-13]
    Covariate confirmed_portugal has been added [aicc=5504.16966367854, lag=-13]
    Series: serie 
    Regression with ARIMA(0,1,1)(1,0,1)[7] errors 

    Coefficients:
              ma1    sar1     sma1  confirmed_spain  deaths_england
          -0.7838  0.8687  -0.7895           0.0061               0
    s.e.   0.0316  0.1739   0.2167           0.0009               0
          recovered_portugal  recovered_france  confirmed_france
                     -0.0376            0.0626           -0.0034
    s.e.              0.0063            0.0145            0.0007
          confirmed_portugal
                      0.0392
    s.e.              0.0084

    sigma^2 = 27983:  log likelihood = -2742.87
    AIC=5503.73   AICc=5504.17   BIC=5540.09
    --------------------------------------------------------------------------------------
    Covariate recovered_spain has been tested [ic=5500.57302699531, lag=-7]
    Covariate deaths_france has been tested [ic=5504.16941132407, lag=-8]
    Significative correlation with lag<=0 could not be found for confirmed_england
    Significative correlation with lag<=0 could not be found for recovered_england
    Covariate deaths_portugal has been tested [ic=5504.1695768999, lag=-21]
    Covariate recovered_spain has been added [aicc=5500.57302699531, lag=-7]
    Series: serie 
    Regression with ARIMA(1,0,1)(1,0,1)[7] errors 

    Coefficients:
             ar1      ma1    sar1     sma1  confirmed_spain  deaths_england
          0.9723  -0.7506  0.6956  -0.5509           0.0064               0
    s.e.  0.0132   0.0370  0.1895   0.2219           0.0009               0
          recovered_portugal  recovered_france  confirmed_france
                     -0.0337            0.0646           -0.0033
    s.e.              0.0063            0.0142            0.0007
          confirmed_portugal  recovered_spain
                      0.0394          -0.0670
    s.e.              0.0085           0.0176

    sigma^2 = 26721:  log likelihood = -2738.96
    AIC=5499.93   AICc=5500.57   BIC=5544.4
    --------------------------------------------------------------------------------------
    Covariate deaths_france has been tested [ic=5500.57309674975, lag=-8]
    Significative correlation with lag<=0 could not be found for confirmed_england
    Significative correlation with lag<=0 could not be found for recovered_england
    Covariate deaths_portugal has been tested [ic=5490.70642412397, lag=-21]
    Covariate deaths_portugal has been added [aicc=5490.70642412397, lag=-21]
    Series: serie 
    Regression with ARIMA(0,1,1)(1,0,0)[7] errors 

    Coefficients:
              ma1    sar1  confirmed_spain  deaths_england  recovered_portugal
          -0.7751  0.1700           0.0066               0             -0.0325
    s.e.   0.0315  0.0516           0.0008               0              0.0061
          recovered_france  confirmed_france  confirmed_portugal  recovered_spain
                    0.0672           -0.0033              0.0402          -0.0736
    s.e.            0.0140            0.0007              0.0082           0.0174
          deaths_portugal
                        0
    s.e.                0

    sigma^2 = 27120:  log likelihood = -2736.13
    AIC=5490.27   AICc=5490.71   BIC=5526.63
    --------------------------------------------------------------------------------------
    Covariate deaths_france has been tested [ic=5490.70642234221, lag=-8]
    Significative correlation with lag<=0 could not be found for confirmed_england
    Significative correlation with lag<=0 could not be found for recovered_england
    Covariate deaths_france has been added [aicc=5490.70642234221, lag=-8]
    Series: serie 
    Regression with ARIMA(0,1,1)(1,0,0)[7] errors 

    Coefficients:
              ma1    sar1  confirmed_spain  deaths_england  recovered_portugal
          -0.7751  0.1700           0.0066               0             -0.0325
    s.e.   0.0315  0.0516           0.0008               0              0.0061
          recovered_france  confirmed_france  confirmed_portugal  recovered_spain
                    0.0672           -0.0033              0.0402          -0.0736
    s.e.            0.0140            0.0007              0.0082           0.0174
          deaths_portugal  deaths_france
                        0              0
    s.e.                0              0

    sigma^2 = 27120:  log likelihood = -2736.13
    AIC=5490.27   AICc=5490.71   BIC=5526.63
    --------------------------------------------------------------------------------------
    Significative correlation with lag<=0 could not be found for confirmed_england
    Significative correlation with lag<=0 could not be found for recovered_england
    No more variables will be added
    The global model does not have stationary errors
    Trying to adjust a model that do have stationary errors
    --------------------------------------------------------------------------------------
    |               Historical of added covariates to the model (ndiff=1)                |
    --------------------------------------------------------------------------------------
                    var lag               ic
        confirmed_spain   0 5596.64110496284
         deaths_england -23 5572.16506296523
     recovered_portugal  -1  5550.9630503792
       recovered_france   0 5540.65497098646
       confirmed_france   0 5522.69957319275
     confirmed_portugal -13 5504.16966367854
        recovered_spain  -7 5500.57302699531
        deaths_portugal -21 5490.70642412397
          deaths_france  -8 5490.70642234221
    --------------------------------------------------------------------------------------
    Series: serie 
    Regression with ARIMA(1,0,1)(1,0,1)[7] errors 

    Coefficients:
             ar1      ma1    sar1     sma1  confirmed_spain  deaths_england
          0.9724  -0.7508  0.6958  -0.5512           0.0064               0
    s.e.  0.0131   0.0370  0.1892   0.2215           0.0009               0
          recovered_portugal  recovered_france  confirmed_france
                     -0.0337            0.0646           -0.0033
    s.e.              0.0063            0.0142            0.0007
          confirmed_portugal  recovered_spain  deaths_portugal  deaths_france
                      0.0395          -0.0669                0              0
    s.e.              0.0085           0.0176                0              0

    sigma^2 = 26721:  log likelihood = -2738.96
    AIC=5499.93   AICc=5500.57   BIC=5544.4

# Point predictions at horizon *h* and prediction intervals (`forecast_model`)

Once a dynamic regression model has been fitted, we can make point
predictions at horizon *h* with their respective prediction intervals
given the `Arima` object from `drm.select()`.

In order to make point predictions at horizon *h* from a DRM it is
necessary to have covariates’ predictions at the same horizon *h*. These
value can be obtained from:

-   If the lag *r* of the covariate is grater or equal *h*, the
    *r**e**m**a**i**n**i**n**g* values of the original time serie.
-   If the lag *r* of the covariate is less than *h*, by fitting an
    ARIMA model with the covariate observations and predicting *h* − *r*
    values.

**Description**: Given an horizon *h*, compute point predictions and
prediction intervals based in a DRM computed by the function
`drm.select()`.

**Returns**: Point predictions and prediction intervals in original
scale.

``` r
forecast_model(serie, xregs, fitted_model, h, mode='bootstrap', levels=c(80, 90))
```

**Argumentos**:

-   `serie` \[`ts`\]: Univariate time series which acts as the dependent
    variable. It must be the same variable used in `drm.select()`.
-   `xregs` \[`mts`\]: Matrix time series with the selected covariates.
    It must be the same matrix used in `drm.select()`.
-   `fitted_model` \[`Arima`\]: Fitted ARIMAX model obtained.
-   `h` \[`numeric`\]: Horizon value for forecasting.
-   `mode` \[`character`\]: Type of prediction intervals.: based on
    residuals normality (`norm`) or via bootstrap (`bootstrap`). By
    default `bootstrap` .
-   `levels` \[`vector`\]: Numerical vector of the confidence intervals.

**Example**:

``` r
# read data and fit model
load("data/patatas.dat")
sales <- patatas[,1]   
xregs <- ts(matrix(patatas[,2]))
colnames(xregs) <- c('price')
chips_model <- drm.select(sales, xregs)
```

    Covariate price has been tested [ic=-71.4371307570267, lag=0]
    Covariate price has been added [aicc=-71.4371307570267, lag=0]
    Series: serie 
    Regression with ARIMA(0,0,4) errors 

    Coefficients:
          ma1     ma2  ma3     ma4  intercept     xreg
            0  0.2884    0  0.5416    15.8559  -2.4682
    s.e.    0  0.0794    0  0.1167     0.1909   0.1100

    sigma^2 = 0.02728:  log likelihood = 41.02
    AIC=-72.05   AICc=-71.44   BIC=-58.83
    --------------------------------------------------------------------------------------
    No more variables will be added
    --------------------------------------------------------------------------------------
    |               Historical of added covariates to the model (ndiff=0)                |
    --------------------------------------------------------------------------------------
       var lag                ic
     price   0 -71.4371307570267
    --------------------------------------------------------------------------------------
    Series: serie 
    Regression with ARIMA(0,0,4) errors 

    Coefficients:
          ma1     ma2  ma3     ma4  intercept     xreg
            0  0.2884    0  0.5416    15.8559  -2.4682
    s.e.    0  0.0794    0  0.1167     0.1909   0.1100

    sigma^2 = 0.02728:  log likelihood = 41.02
    AIC=-72.05   AICc=-71.44   BIC=-58.83

``` r
# point predictons
preds <- forecast_model(sales, xregs, chips_model, h=10, mode='bootstrap')
display(plot_forecast(preds), name='preds_chips')
```
