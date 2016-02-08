#' Fully Modified OLS
#'
#' Computes the Phillips and Hansen (1990) Fully Modified OLS estimator.
#'
#' @param x [\code{numeric} | \code{matrix} | \code{data.frame}]\cr
#'   RHS variables on which to apply the FM-OLS estimation (see Details).
#'
#' @param y [\code{numeric} | \code{matrix} | \code{data.frame}]\cr
#'   LHS variable(s) on which to apply the FM-OLS estimation (see Details).
#'   Usually one-dimensional, but a \code{matrix} or
#'   \code{data.frame} with more than one column is also possible.
#'
#' @param deter [\code{numeric} | \code{matrix} | \code{data.frame} |
#'               \code{NULL}]\cr
#'   Deterministic variable to include in the equation (see Details). If it's
#'   \code{NULL} or missing, no deterministic variable is included in the model.
#'
#' @param kernel [\code{character(1)}]\cr
#'   The kernel function to use for calculating the long-run variance.
#'   Default is Bartlett kernel (\code{"ba"}), see Details for alternatives.
#'
#' @param bandwidth [\code{character(1)} | \code{integer(1)}]\cr
#'   The bandwidth to use for calculating the long-run variance.
#'   Default is Andrews (1991) (\code{"and"}), an alternative is Newey West
#'   (1994) (\code{"nw"}).
#'
#' @param demeaning [\code{logical}]\cr
#'   Demeaning of residuals in \code{\link{getLongRunVar}}.
#'   Default is \code{FALSE}.
#'
#' @param check [\code{logical}]\cr
#'   Wheather to check (and if necessary convert) the arguments.
#'   See \code{\link{checkVars}} for further information.
#'
#' @param ... Arguments passed to \code{\link{getBandwidthNW}}.
#'
#' @details
#' The equation for which the FM-OLS estimator is calculated:
#' \deqn{y = \delta \cdot D + \beta \cdot x + u}{y = \delta * D + \beta * x + u}
#' with \eqn{D} as the deterministics matrix.
#' Then \eqn{\theta = (\delta', \beta')'} is the full parameter vector.
#'
#' The calculation of t-values and the variance-covariance matrix is only
#' possible, if \code{y} is one-dimensional.
#'
#' @return [\code{cointReg}]. List with components:
#' \describe{
#'   \item{\code{delta} [\code{numeric} | \code{matrix}]}{
#'     coefficients as vector / matrix}
#'
#'   \item{\code{beta} [\code{numeric} | \code{matrix}]}{
#'     coefficients as vector / matrix}
#'
#'   \item{\code{theta} [\code{numeric} | \code{matrix}]}{
#'     combined coefficients of
#'     \code{beta} and \code{delta} as vector / matrix}
#'
#'   \item{\code{sd.theta} [\code{numeric}]}{
#'     standard errors for \code{theta}}
#'
#'   \item{\code{t.theta} [\code{numeric}]}{
#'     t-values for \code{theta}}
#'
#'   \item{\code{p.theta} [\code{numeric}]}{
#'     p-values for \code{theta}}
#'
#'   \item{\code{residuals} [\code{numeric}]}{
#'     FM-OLS residuals (first value is always missing)}
#'
#'   \item{\code{omega.u.v} [\code{numeric}]}{
#'     conditional long-run variance based on OLS residuals.}
#'
#'   \item{\code{varmat} [\code{matrix}]}{
#'     variance-covariance matrix}
#'
#'   \item{\code{Omega} [\code{list}]}{
#'     the whole long-run variance matrix and parts of it}
#'
#'   \item{\code{beta.OLS} [\code{numeric} | \code{matrix}]}{
#'     OLS coefficients as vector / matrix}
#'
#'   \item{\code{delta.OLS} [\code{numeric} | \code{matrix}]}{
#'     OLS coefficients as vector / matrix}
#'
#'   \item{\code{u.OLS} [\code{numeric}]}{
#'     OLS residuals}
#'
#'   \item{\code{bandwidth} [\code{list}]}{
#'     \code{number} and \code{name} of bandwidth}
#'
#'   \item{\code{kernel} [\code{character}]}{
#'     abbr. name of kernel type}
#' }
#'
#' @family cointReg
#'
#' @examples
#' set.seed(1909)
#' x1 = cumsum(rnorm(100, mean = 0.05, sd = 0.1))
#' x2 = cumsum(rnorm(100, sd = 0.1)) + 1
#' x3 = cumsum(rnorm(100, sd = 0.2)) + 2
#' x = cbind(x1, x2, x3)
#' y = x1 + x2 + x3 + rnorm(100, sd = 0.2) + 1
#' deter = cbind(level = 1, trend = 1:100)
#' test = cointRegFM(x, y, deter, kernel = "ba", bandwidth = "and")
#' print(test)
#'
#' @export


cointRegFM <- function(x, y, deter, kernel = c("ba", "pa", "qs", "tr"),
                       bandwidth = c("and", "nw"), demeaning = FALSE,
                       check = TRUE, ...) {

  y.name <- deparse(substitute(y))
  x.name <- deparse(substitute(x))
  d.name <- deparse(substitute(deter))
  mod.name <- paste0(y.name, " ~ ",
                     ifelse(missing(deter) || is.null(deter), "",
                            paste0(d.name, " + ")), x.name)

  if (check) {
    env <- environment()
    checkVars(kernel = kernel, bandwidth = bandwidth,
              demeaning = demeaning, .env = env)
    y <- checkObject(y.fm = y)
    x <- checkObject(x.coint = x)
    if (missing(deter) || is.null(deter)) {
      deter <- matrix(nrow = nrow(x), ncol = 0)
    } else {
      deter <- checkObject(deter = deter)
      if (is.null(colnames(deter))) {
        colnames(deter) <- make.unique(rep(d.name, ncol(deter)))
      }
    }
    if (is.null(colnames(x))) {
      colnames(x) <- make.unique(rep(x.name, ncol(x)))
    }
  }

  x.T <- nrow(x)
  x.k <- ncol(x)
  y.k <- ncol(y)
  d.k <- ncol(deter)

  Z <- cbind(deter, x)

  mod.ols <- lm(y ~ Z - 1)
  theta.ols <- t(mod.ols$coefficients)
  delta.ols <- theta.ols[, 0:d.k, drop = FALSE]
  beta.ols <- theta.ols[, (d.k + 1):(d.k + x.k), drop = FALSE]

  u.ols <- mod.ols$residuals

  if (y.k > 1) {
    u.4var <- u.ols[-1, ]
  } else {
    u.4var <- u.ols[-1]
  }
  x.delta <- colDiffs(x)
  u <- cbind(u.4var, x.delta)

  if (!is.numeric(bandwidth)) {
    bw <- getBandwidth(u, bandwidth = bandwidth, kernel = kernel,
                       check = FALSE, ...)
    band <- switch(bandwidth, and = "Andrews", nw = "Newey-West")
  } else {
    bw <- bandwidth
    band <- "set by user"
  }

  lrvar <- getLongRunVar(u, kernel = kernel, bandwidth = bw,
                         demeaning = demeaning, check = FALSE)
  tmp <- lapply(lrvar, function(x) {
    out <- list()
    out[["all"]] <- x
    out[["uu"]] <- x[1:y.k, 1:y.k, drop = FALSE]
    out[["uv"]] <- x[1:y.k, (y.k + 1):(y.k + x.k), drop = FALSE]
    out[["vu"]] <- t(out[["uv"]])
    out[["vv"]] <- x[(y.k + 1):(y.k + x.k), (y.k + 1):(y.k + x.k), drop = FALSE]
    return(out)
  })

  Omega <- tmp[[1]]
  Delta <- tmp[[2]]

  Omegavv.inv <- trySolve(Omega$vv)
  Omegavv.inv.vu <- Omegavv.inv %*% Omega$vu

  Omega.u.v <- Omega$uu - (Omega$uv %*% Omegavv.inv.vu)
  Delta.vuplus <- Delta$vu - (Delta$vv %*% Omegavv.inv.vu)

  y.plus <- y[-1, , drop = FALSE] - (x.delta %*% Omegavv.inv.vu)
  Zfm <- Z[-1, , drop = FALSE]
  Zfm2s <- trySolve(t(Zfm) %*% Zfm)

  numerat <- t(y.plus) %*% Zfm - cbind(matrix(0, nrow = y.k, ncol = d.k),
                                       x.T * t(Delta.vuplus))
  theta.fm <- drop(Zfm2s %*% t(numerat))
  delta.fm <- theta.fm[0:d.k]
  beta.fm <- theta.fm[(d.k + 1):(d.k + x.k)]

  u.plus <- c(NA, (y - Z %*% theta.fm)[-1, ])

  if (y.k == 1) {
    varmat <- Omega.u.v[1, 1] * Zfm2s
    sd.theta <- sqrt(diag(varmat))
    t.theta <- theta.fm / sd.theta
    df <- x.T - x.k - d.k
    p.theta <- 2 * pt(-abs(t.theta), df = df)
  } else {
    varmat <- NULL
    sd.theta <- NULL
    t.theta <- NULL
    p.theta <- NULL
  }

  out <- list(delta = delta.fm, beta = beta.fm, theta = theta.fm,
              sd.theta = sd.theta, t.theta = t.theta, p.theta = p.theta,
              residuals = u.plus, omega.u.v = Omega.u.v, varmat = varmat,
              Omega = Omega, beta.OLS = beta.ols, delta.OLS = delta.ols,
              u.OLS = u.ols, bandwidth = list(name = band, number = bw),
              kernel = kernel, mod = "FM", name = mod.name)
  class(out) <- "cointReg"

  return(out)
}
