#' The cointReg package
#'
#' Parameter Estimation and Inference in a Cointegrating Regression
#'
#' See the vignette:\cr
#' \code{vignette("cointReg")}
#'
#' See the DESCRIPTION:\cr
#' \code{help(package = cointReg)}
#'
#' See the README:\cr
#' \url{https://github.com/aschersleben/cointReg/blob/master/README.md}
#'
#' Open the package documentation page:\cr
#' \code{package?cointReg}
#'
#' Further information and bug reporting:\cr
#' \url{https://github.com/aschersleben/cointReg}
#'
#' @section Functions:
#' \itemize{
#'   \item \code{\link{cointReg}(method = c("FM", "D", "IM"), ...)}\cr
#'         General function to estimate parameters of
#'         the given model. Three methods are possible;
#'         they can be choosen directly by using one of the following
#'         functions:
#'   \itemize{
#'     \item \code{\link{cointRegFM}}: Fully Modified OLS
#'     \item \code{\link{cointRegD}}: Dynamic OLS
#'     \item \code{\link{cointRegIM}}: Integrated Modified OLS
#'   }
#'   \item \code{\link[=print.cointReg]{print}}\cr
#'         Print clear results.
#'   \item \code{\link[=plot.cointReg]{plot}}\cr
#'         Plot the residuals of a \code{cointReg} model.
#'   \item Helper functions:
#'   \itemize{
#'     \item Checking inputs and arguments:\cr
#'           \code{\link{checkObject}}, \code{\link{checkVars}}
#'     \item Calculation of bandwidth and long run variance:\cr
#'           \code{\link{getBandwidth}}, \code{\link{getBandwidthAnd}},
#'           \code{\link{getBandwidthNW}}\cr
#'           \code{\link{getLongRunVar}}, \code{\link{getLongRunWeights}}
#'     \item Additional D-OLS functions:\cr
#'           \code{\link{getLeadLag}}, \code{\link{makeLeadLagMatrix}},
#'           \code{\link{getModD}},  \code{\link{checkDoptions}}
#'   }
#' }
#'
#' @docType package
#' @name cointReg-package
NULL

#' Estimation and Inference for cointegrating regressions
#'
#' Computes either the Phillips and Hansen (1990) Fully Modified OLS estimator,
#' or the Saikkonen (1990) Dynamic OLS estimator, or the Vogelsang and Wagner
#' (2014) Integrated Modified OLS estimator.
#'
#' @param method [\code{character(1)}]\cr
#'   Select the method for the estimation of your cointegration model:
#'   \itemize{
#'     \item \code{"FM"}: FM-OLS (default), see details at \link{cointRegFM}
#'     \item \code{"D"}: D-OLS, see details at \link{cointRegD}
#'     \item \code{"IM"}: IM-OLS, see details at \link{cointRegIM}
#'   }
#'
#' @param x [\code{numeric} | \code{matrix} | \code{data.frame}]\cr
#'   RHS variables on which to apply the model estimation.
#'
#' @param y [\code{numeric} | \code{matrix} | \code{data.frame}]\cr
#'   LHS variable(s) on which to apply the model estimation.
#'   Usually one-dimensional, but a \code{matrix} or
#'   \code{data.frame} with more than one column is also possible
#'   (only for FM-OLS).
#'
#' @param ... [\code{any}]\cr Arguments passed to the corresponding \code{cointReg} function,
#'   like:
#'   \itemize{
#'     \item \code{x}, \code{y}, \code{deter}: data to include in the model
#'     \item \code{kernel}, \code{bandwidth}: parameters for calculating the
#'           long-run variance
#'     \item \code{n.lead}, \code{n.lag}, \code{kmax}, \code{info.crit}:
#'           D-OLS specific arguments.
#'     \item \code{selector}, \code{t.test}: IM-OLS specific arguments.
#'     \item \code{check}: Wheather to check (and if necessary convert)
#'           the arguments. See \code{\link{checkVars}} for further information.
#'   }
#'
#' @return [\code{cointReg}] object.
#'
#' @family cointReg
#'
#' @references
#'   \itemize{
#'     \item Phillips, P.C.B. and B. Hansen (1990): "Statistical Inference in
#'           Instrumental Variables Regression with I(1) Processes,"
#'           \emph{Review of Economic Studies}, 57, 99--125,
#'           \href{http://dx.doi.org/10.2307/2297545}{DOI:10.2307/2297545}.
#'     \item Phillips, P.C.B. and M. Loretan (1991): "Estimating Long Run
#'           Economic Equilibria," \emph{Review of Economic Studies}, 58,
#'           407--436,
#'           \href{http://dx.doi.org/10.2307/2298004}{DOI:10.2307/2298004}.
#'     \item Saikkonen, P. (1991): "Asymptotically Efficient Estimation of
#'           Cointegrating Regressions," \emph{Econometric Theory}, 7, 1--21,
#'           \href{http://dx.doi.org/10.1017/S0266466600004217}{DOI:10.1017/S0266466600004217}.
#'     \item Stock, J.H. and M.W. Watson (1993): "A Simple Estimator of
#'           Cointegrating Vectors in Higher Order Integrated Systems,"
#'           \emph{Econometrica}, 61, 783--820,
#'           \href{http://dx.doi.org/10.2307/2951763}{DOI:10.2307/2951763}.
#'     \item Vogelsang, T.J. and M. Wagner (2014): "Integrated Modified OLS
#'           Estimation and Fixed-b Inference for Cointegrating Regressions,"
#'           \emph{Journal of Econometrics}, 148, 741--760,
#'           \href{http://dx.doi.org/10.1016/j.jeconom.2013.10.015}{DOI:10.1016/j.jeconom.2013.10.015}.
#'   }
#'
#' @examples
#' set.seed(1909)
#' x1 = cumsum(rnorm(100, mean = 0.05, sd = 0.1))
#' x2 = cumsum(rnorm(100, sd = 0.1)) + 1
#' x3 = cumsum(rnorm(100, sd = 0.2)) + 2
#' x = cbind(x1, x2, x3)
#' y = x1 + x2 + x3 + rnorm(100, sd = 0.2) + 1
#' deter = cbind(level = 1, trend = 1:100)
#' cointReg("FM", x = x, y = y, deter = deter, kernel = "ba",
#'          bandwidth = "and")
#'
#' # Compare the results of all three models:
#' res = sapply(c("FM", "D", "IM"), cointReg, x = x, y = y, deter = deter)
#' do.call(cbind, lapply(res, "[[", "theta"))
#'
#' @export

cointReg <- function(method = c("FM", "D", "IM"), x, y, ...) {
  method <- match.arg(method)
  switch(method,
         FM = cointRegFM(x, y, ...),
         D = cointRegD(x, y, ...),
         IM = cointRegIM(x, y, ...))
}
