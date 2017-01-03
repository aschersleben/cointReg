# cointReg

Parameter Estimation and Inference in a Cointegrating Regression

[![Travis-CI Build Status](https://travis-ci.org/aschersleben/cointReg.svg?branch=master)](https://travis-ci.org/aschersleben/cointReg)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/aschersleben/cointReg?branch=master&svg=true)](https://ci.appveyor.com/project/aschersleben/cointReg)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/cointReg)](https://cran.r-project.org/package=cointReg)

* Installation:
```r
install.packages("cointReg")
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

* [Package NEWS](https://github.com/aschersleben/cointReg/blob/master/inst/NEWS.md)

* Installation of the current development version from GitHub:
```r
devtools::install_github("aschersleben/cointReg", build_vignettes = TRUE)
```


## Theoretical background

Cointegration methods are widely used in empirical macroeconomics and empirical finance. It is well known that in a cointegrating regression the ordinary least squares (OLS) estimator of the parameters is super-consistent, i.e. converges at rate equal to the sample size $T$. When the regressors are endogenous, the limiting distribution of the OLS estimator is contaminated by so-called second order bias terms, see e.g. [Phillips and Hansen (1990)](http://dx.doi.org/10.2307/2297545). The presence of these bias terms renders inference difficult.

Consequently, several modifications to OLS that lead to zero mean Gaussian mixture limiting distributions have been proposed, which in turn make standard asymptotic inference feasible. These methods include the fully modified OLS (FM-OLS) approach of [Phillips and Hansen (1990)](http://dx.doi.org/10.2307/2297545), the dynamic OLS (D-OLS) approach of [Phillips and Loretan (1991)](http://dx.doi.org/10.2307/2298004), [Saikkonen (1991)](http://dx.doi.org/10.1017/S0266466600004217) and [Stock and Watson (1993)](http://dx.doi.org/10.2307/2951763) and the new estimation approach called integrated modified OLS (IM-OLS) of [Vogelsang and Wagner (2014)](http://dx.doi.org/10.1016/j.jeconom.2013.10.015).

The latter is based on an augmented partial sum (integration) transformation of the regression model. IM-OLS is similar in spirit to the FM- and D-OLS approaches, with the key difference that it does not require estimation of long run variance matrices and avoids the need to choose tuning parameters (kernels, bandwidths, lags). However, inference does require that a long run variance be scaled out.

This package provides functions for the parameter estimation and inference with all three modified OLS approaches. That includes the automatic bandwidth selection approaches of [Andrews (1991)](http://dx.doi.org/10.2307/2938229) and of [Newey and West (1994)](http://dx.doi.org/10.2307/2297912) as well as the calculation of the long run variance.


## References

* Andrews, D.W.K. (1991): "Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation," _Econometrica_, 59, 817--854, [DOI:10.2307/2938229](http://dx.doi.org/10.2307/2938229).

* Newey, W.K. and K.D. West (1994): "Automatic Lag Selection in Covariance Matrix Estimation", _Review of Economic Studies_, 61, 631--653, [DOI:10.2307/2297912](http://dx.doi.org/10.2307/2297912).

* Phillips, P.C.B. and B. Hansen (1990): "Statistical Inference in Instrumental Variables Regression with I(1) Processes," _Review of Economic Studies_, 57, 99--125, [DOI:10.2307/2297545](http://dx.doi.org/10.2307/2297545).

* Phillips, P.C.B. and M. Loretan (1991): "Estimating Long Run Economic Equilibria," _Review of Economic Studies_, 58, 407--436, [DOI:10.2307/2298004](http://dx.doi.org/10.2307/2298004).

* Saikkonen, P. (1991): "Asymptotically Efficient Estimation of Cointegrating Regressions," _Econometric Theory_, 7, 1--21, [DOI:10.1017/S0266466600004217](http://dx.doi.org/10.1017/S0266466600004217).

* Stock, J.H. and M.W. Watson (1993): "A Simple Estimator of Cointegrating Vectors in Higher Order Integrated Systems," _Econometrica_, 61, 783--820, [DOI:10.2307/2951763](http://dx.doi.org/10.2307/2951763).

* Vogelsang, T.J. and M. Wagner (2014): "Integrated Modified OLS Estimation and Fixed-b Inference for Cointegrating Regressions," _Journal of Econometrics_, 148, 741--760, [DOI:10.1016/j.jeconom.2013.10.015](http://dx.doi.org/10.1016/j.jeconom.2013.10.015).


## Consistent Monitoring of Stationarity and Cointegrating Relationships

You may also be interested in the `cointmonitoR` package: We propose a consistent monitoring procedure to detect a structural change from a cointegrating relationship to a spurious relationship. The procedure is based on residuals from modified least squares estimation, using either Fully Modified, Dynamic or Integrated Modified OLS. Our R package makes use of the functions in `cointReg`. See its [GitHub page](https://github.com/aschersleben/cointmonitoR) for further information.
