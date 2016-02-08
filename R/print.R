#' Print Method for Cointegration Models (Modified OLS).
#'
#' Printing objects of class \code{"cointReg"}.
#'
#' @param x [\code{cointReg}]\cr
#'   Object of class \code{"cointReg"}, i.e. the result of
#'   \code{cointRegFM}, \code{cointRegD} or \code{cointRegIM}.
#' @param ... ignored
#' @param digits [\code{numeric}]\cr
#'   Number of significant digits to be used.
#' @param all.coeffs [\code{logical}]\cr
#'   Whether to show all coefficients (i. e. the "real" regressors AND the
#'   auxiliary regressors). Default is \code{FALSE}.
#'
#' @return
#'   The invisible \code{x} object.
#'
#' @family cointReg
#'
#' @examples
#' set.seed(42)
#' x = data.frame(x1 = cumsum(rnorm(200)), x2 = cumsum(rnorm(200)))
#' eps1 = rnorm(200, sd = 2)
#' y = x$x1 - x$x2 + 10 + eps1
#' deter = cbind(level = rep(1, 200))
#'
#' test.fm = cointRegFM(x = x, y = y, deter = deter)
#' print(test.fm)
#'
#' test.d = cointRegD(x = x, y = y, deter = deter)
#' print(test.d)
#'
#' test.im2 = cointRegIM(x = x, y = y, deter = deter)
#' print(test.im2)
#'
#' @export



print.cointReg <- function(x, ..., digits = getOption("digits"),
                           all.coeffs = FALSE) {

  if (!(class(x) == "cointReg"))
    stop("Argument x must be of type \"cointReg\".")

  cat(paste0("\n### ", x$mod, "-OLS model ###\n\n"))

  cat("Model:\t")
  cat(paste0("\t",  x$name, "\n\n"))

  cat("Parameters:")
  cat(paste0("\t", "Kernel = \"", x$kernel, "\""),
      " //  Bandwidth =", format(x$bandwidth$number, digits = digits),
      paste0("(\"", x$bandwidth$name, "\")\n\n"))
  if(!is.null(x$lead.lag)) {
    cat(paste0("\t\t", "Leads ="), x$lead.lag$n.lead,
        "/ Lags =", x$lead.lag$n.lag,
        ifelse(!is.null(x$lead.lag$kmax),
               paste("/ kmax =", x$lead.lag$kmax), "(set manually)"),
        ifelse(!is.null(x$lead.lag$info.crit),
               paste("/ IC =", x$lead.lag$info.crit, "\n\n"), "\n\n"))
  }

  cat("Coefficients:\n")
  if (!is.null(x$t.theta)) {
    if (!all.coeffs | x$mod == "FM") {
      estim <- x$theta
      sd <- x$sd.theta
      tval <- x$t.theta
      pval <- x$p.theta
    } else {
      estim <- x$theta.all
      sd <- x$sd.theta.all
      tval <- x$t.theta.all
      pval <- x$p.theta.all
    }
    tab <- cbind(Estimate = estim, Std.Err = sd, "t value" = tval,
                 "Pr(|t|>0)" = pval)
    stats::printCoefmat(tab)
  } else {
    print.default(format(drop(x$theta), digits = digits), print.gap = 2L,
                  quote = FALSE)
  }

  return(invisible(x))
}
