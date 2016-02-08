#' Plot Method for Cointegration Models (Modified OLS).
#'
#' Plotting objects of class \code{"cointReg"}. Currently, only the residuals
#' will be plotted.
#'
#' @param x [\code{cointReg}]\cr
#'   Object of class \code{"cointReg"}, i.e. the result of
#'   \code{\link{cointRegFM}}, \code{\link{cointRegD}},
#'   or \code{\link{cointRegIM}}.
#'
#' @param type [\code{character}]\cr
#'   Plot type (from \code{\link{plot}}). Default is \code{"l"}.
#'
#' @param main,xlab,ylab [\code{character}]\cr
#'   Title and axis titles (from \code{\link{plot}}). Default values will be
#'   generated from the contents of \code{x}.
#'
#' @param axes [\code{logical}]\cr
#'   Whether to add axes (from \code{\link{plot}}) to the plot.
#'
#' @param ... [\code{any}]\cr
#'   Further arguments passed to \code{\link{plot}}.
#'
#' @family cointReg
#'
#' @examples
#' set.seed(42)
#' x = data.frame(x1 = cumsum(rnorm(200)), x2 = cumsum(rnorm(200)))
#' eps1 = rnorm(200, sd = 2)
#' y = x$x1 - x$x2 + 10 + eps1
#' deter = cbind(level = rep(1, 200))
#' test = cointRegFM(x = x, y = y, deter = deter)
#' plot(test)
#'
#' @export

plot.cointReg <- function(x, type, main, xlab, ylab, axes = TRUE, ...) {

  if (!(class(x) == "cointReg"))
    stop("Argument x must be of type \"cointReg\".")

  if (missing(type))
    type <- "l"

  if (missing(main))
    main <- paste("Residuals of Cointegrating Regression",
                  paste0("(", x$mod,"-OLS)"))

  if (missing(xlab))
    xlab <- "Observation Number"

  if (missing(ylab))
    ylab <- paste("Residuals")

  plot.default(x$residuals, type = type, main = main, xlab = xlab,
               ylab = ylab, axes = axes, ...)
}
