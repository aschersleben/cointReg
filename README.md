# cointReg

Parameter Estimation and Inference in a Cointegrating Regression

[![Travis-CI Build Status](https://travis-ci.org/aschersleben/cointReg.svg?branch=master)](https://travis-ci.org/aschersleben/cointReg)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/cointReg)](https://cran.r-project.org/web/packages/cointReg)

* Installation:
```r
devtools::install_github("aschersleben/cointReg", build_vignettes = TRUE)
library("cointReg")
```

* Simple example:
```r
set.seed(1909)
x <- cumsum(rnorm(200, mean = 0, sd = 0.1)) + 10
y <- x + rnorm(200, sd = 0.4) + 2
deter <- rep(1, 200)
test <- cointRegFM(x = x, y = y, deter = deter)
print(test)
plot(test)
```

* Package vignette: Provides further examples and explanations.
```r
vignette("cointReg")
```

* Package help page: Overview of all available functions:
```r
package?cointReg
```

* [Package issues](https://github.com/aschersleben/cointReg/issues)



## Theoretical background

Cointegration methods are widely used in empirical macroeconomics and empirical finance. It is well known that in a cointegrating regression the ordinary least squares (OLS) estimator of the parameters is super-consistent, i.e. converges at rate equal to the sample size T. When the regressors are endogenous, the limiting distribution of the OLS estimator is contaminated by so-called second order bias terms, see e.g. Phillips and Hansen (1990). The presence of these bias terms renders inference difficult.

Consequently, several modifications to OLS that lead to zero mean Gaussian mixture limiting distributions have been proposed, which in turn make standard asymptotic inference feasible. These methods include the fully modified OLS (FM-OLS) approach of Phillips and Hansen (1990), the dynamic OLS (D-OLS) approach of Phillips and Loretan (1991), Saikkonen (1991) and Stock and Watson (1993) and the new estimation approach called integrated modified OLS (IM-OLS) of Vogelsang and Wagner (2014).

The latter is based on an augmented partial sum (integration) transformation of the regression model. IM-OLS is similar in spirit to the FM- and D-OLS approaches, with the key difference that it does not require estimation of long run variance matrices and avoids the need to choose tuning parameters (kernels, bandwidths, lags). However, inference does require that a long run variance be scaled out.

This package provides functions for the parameter estimation and inference with all three modified OLS approaches. That includes the automatic bandwidth selection approaches of Andrews (1991) and of Newey and West (1994) as well as the calculation of the long run variance.



## Consistent Monitoring of Stationarity and Cointegrating Relationships

You may also be interested in the `cointmonitoR` package: We propose a consistent monitoring procedure to detect a structural change from a cointegrating relationship to a spurious relationship. The procedure is based on residuals from modified least squares estimation, using either Fully Modified, Dynamic or Integrated Modified OLS. Our R package makes use of the functions in `cointReg`. See its [GitHub page](https://github.com/aschersleben/cointmonitoR) for further information.
