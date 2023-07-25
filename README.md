# Welcome to <span style="color:green">`dynamic-arimax`</span> repository!

This is an [R](https://www.r-project.org/) implementation of a new covariates selection method in Dynamic Regression Models. This proposal was published in: 

- [XoveTIC 2022 Conference](https://xovetic.citic.udc.es/).
- [IX Conference of R-Users in Galicia](https://www.r-users.gal/).

The PDF file of both proceedings are attached in this repository in folder `proceedings/`. We recommend reading at least one of the documents to use this code to deeply understand the mathematics of our method.

For any suggestion or issue with our code, please [contact us](mailto:ana.ezquerro@udc.es) in order to solve it and improve our implementation. 

## Table of Contents 

1. [Installation](#installation)
1. [Structure of the module](#structure-of-the-module)
2. [Documentation and examples](#documentation-and-examples)
3. [Table of results](#table-of-results)
4. [Acknowledgements](#acknowledgements)

## Installation

In order to use this implementation and run all files, the following prerequisites are needed:

- [R](https://www.r-project.org/) and [RStudio IDE](https://www.rstudio.com/products/rstudio/download/) for automatic downloading of the necessary packages.
- R packages: [`fpp2`](https://cran.r-project.org/web/packages/fpp2/index.html), [`tseries`](https://cran.r-project.org/web/packages/tseries/index.html), [`TSA`](https://cran.r-project.org/web/packages/TSA/index.html), [`seastests`](https://cran.r-project.org/web/packages/seastests/index.html), [`forecast`](https://cran.r-project.org/web/packages/forecast/index.html), [`plotly`](https://plotly.com/r/), [`prettydoc`](https://prettydoc.statr.me/), [`stringi`](https://cran.r-project.org/web/packages/stringi/index.html), [`stringr`](https://cran.r-project.org/web/packages/stringr/index.html), [`polynom`](https://cran.r-project.org/web/packages/polynom/index.html), [`parallel`](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf).
- [rtools](https://cran.r-project.org/bin/windows/Rtools/).

Once [R](https://www.r-project.org/),  [RStudio IDE](https://www.rstudio.com/products/rstudio/download/) and [rtools](https://cran.r-project.org/bin/windows/Rtools/) have been installed, you can run [`installation.R`](installation.R) to automatically install all R-packages needed and [`test-installation.R`](test-installation.R) to check if all libraries have been correctly added.

## Structure of the module

- [`auto-fit.R`](auto-fit.R): Implementation of automatic fitting in ARIMA or ARIMAX models. It uses the [`forecast::auto.arima()`](https://www.rdocumentation.org/packages/forecast/versions/8.17.0) function and iteratively removes non-significative coefficients.
- [`auto-select.R`](auto-select.R): Implementation of the covariates selection method and their respective correlation lags. For more information we suggest reading the proceedings attached to the repository.
- [`forecasting.R`](forecasting.R): Implementation of dynamic regression models forecasting once a model has been fitted with the selection function.
- [`plot-tools.R`](plot-tools.R): Script for fancy Plotly graphics.
- [`documentation/`](documentation/): Folder where documentation examples are provided.
- [`data/`](data/): Datasets needed to run  examples in [`EXAMPLES.md`](EXAMPLES.md).


## Documentation and examples 

Consult a detailed documentation of the code and examples of use in [`Documentation.html`](documentation/Documentation.html) file.





## Simulation results 

We present in the following tables some metrics obtained via simulating $M=100$ scenarios where a dependent variable was artificially constructed via some randomly generated time series (modelable by an ARIMA model). Specifically, in each scenario:

1. Seven time series were randomly generated. Six of them were used as the set of covariates: $\mathcal{X} = \{X_t^{(1)}, ..., X_t^{(6)}\}$ and the remaining as the residuals of the model $\eta_t$.
2. Random lags $r_i \in[0, 6]$ for $i=1...6$, where selected for each covariate as well as regression coefficients $\beta_0,...,\beta_3$.
3. The dependent variable $Y_t$ was constructed via the DR model formula:

$$ Y_t = \beta_0 + \beta_1 X_{t-r_1}^{(1)} + \beta_2 X_{t-r_2}^{(2)} + \beta_3 X_{t-r_3}^{(3)} + \eta_t$$

We tested our selection method with different configurations:
- With different stationary tests (via `auto.arima()` function or Dickey-Fuller test).
- With different information criterions (AIC, BIC or AICc).

For more detailed information about the simulation procedure, please read the proceeding of the repository.

### Results when $\eta_t \sim \text{ARMA}(p,q)$ is stationary

- Percentage of correctly added covariates to the model (*true positive*):

|                |      AIC      |  BIC          |   AICc        |
|:--------------:|:-------------:|:-------------:|:-------------:|
| **adf.test**   | 97.66%        | 97.66%        | 97.66%        |
| **auto.arima** | 98.33%        | 98.33%        | 98.33%        |

- Percentage of incorrectly added covariates to the model (*false positive*):

|                |      AIC      |  BIC          |   AICc        |
|:--------------:|:-------------:|:-------------:|:-------------:|
| **adf.test**   | 3.66%         | 1.33%         | 3.66%         |
| **auto.arima** | 3.66%         | 1.33%         | 3.66%         |

- Percentage of correctly **not** added covariates to the model (*true negative*):

|                |      AIC      |  BIC          |   AICc        |
|:--------------:|:-------------:|:-------------:|:-------------:|
| **adf.test**   | 96.33%        | 98.66%        | 96.33%        |
| **auto.arima** | 96.33%        | 98.66%        | 96.33%        |

- Percentage of incorrectly **not** added covariates to the model (*false negative*):


|                |      AIC      |  BIC          |   AICc        |
|:--------------:|:-------------:|:-------------:|:-------------:|
| **adf.test**   | 2.33%         | 2.33%         | 2.33%         |
| **auto.arima** | 1.66%         | 1.66%         | 1.66%         |


### Results when $\eta_t \sim \text{ARIMA}(p,d,q)$ is non-stationary

- Percentage of correctly added covariates to the model (*true positive*):

|                |      AIC      |  BIC          |   AICc        |
|:--------------:|:-------------:|:-------------:|:-------------:|
| **adf.test**   | 93.33%        | 93.33%        | 93.33%        |
| **auto.arima** | 94.33%        | 94.66%        | 95.33%        |

- Percentage of incorrectly added covariates to the model (*false positive*):

|                |      AIC      |  BIC          |   AICc        |
|:--------------:|:-------------:|:-------------:|:-------------:|
| **adf.test**   | 4.33%         | 0.30%         | 4.33%         |
| **auto.arima** | 5.00%         | 1.33%         | 5.00%         |

- Percentage of correctly **not** added covariates to the model (*true negative*):

|                |      AIC      |  BIC          |   AICc        |
|:--------------:|:-------------:|:-------------:|:-------------:|
| **adf.test**   | 95.00%        | 98.66%        | 95.00%        |
| **auto.arima** | 94.66%        | 99.66%        | 95.66%        |

- Percentage of incorrectly **not** added covariates to the model (*false negative*):


|                |      AIC      |  BIC          |   AICc        |
|:--------------:|:-------------:|:-------------:|:-------------:|
| **adf.test**   | 6.66%         | 6.66%         | 6.66%         |
| **auto.arima** | 4.66%         | 5.33%         | 4.66%         |

## Acknowledgements

To Banco Santander for the scholarships offered in 2021/2022, which helped the investigation of this proposal, and to [MODES investigation group](https://dm.udc.es/modes/) of [University of A Coruña](https://www.udc.es/).