#' Integrated Modified OLS
#'
#' Computes the Vogelsang and Wagner (2014) Integrated Modified OLS estimator.
#'
#' @param x [\code{numeric} | \code{matrix} | \code{data.frame}]\cr
#'   RHS variables on which to apply the IM-OLS estimation (see Details).
#'
#' @param y [\code{numeric} | \code{matrix} | \code{data.frame}]\cr
#'   LHS variable(s) on which to apply the IM-OLS estimation (see Details).
#'   Has to be one-dimensional. If \code{matrix}, it may
#'   have only one row or column, if \code{data.frame} just one column.
#'
#' @param selector [\code{numeric}]\cr
#'   Choose the regression type: \code{1}, \code{2}, or \code{c(1, 2)}
#'   (see Details). Default is \code{1}.
#'
#' @param t.test [\code{logical}]\cr
#'   Wheather to calculate t-values for the coefficients of the first
#'   regression. Default is \code{TRUE}. Attention: Needs more calculation
#'   time, because an additional FM-OLS model has to be fitted to get the
#'   long-run variance.
#'
#' @inheritParams cointRegFM
#'
#' @details
#' The equation for which the IM-OLS estimator is calculated (type 1):
#' \deqn{S_y = \delta \cdot S_{D} + \beta \cdot S_{x} + \gamma \cdot x + u}{
#' S[y] = \delta * S[D] + \beta * S[x] + \gamma * x + u}
#' where \eqn{S_y}{S[y]}, \eqn{S_x}{S[x]} and \eqn{S_D}{S[D]} are the cumulated
#' sums of \eqn{y}, \eqn{x} and \eqn{D} (with \eqn{D} as the deterministics
#' matrix).
#' Then \eqn{\theta = (\delta', \beta', \gamma')'} is the full parameter vector.
#'
#' The equation for which the IM-OLS estimator is calculated (type 2):
#' \deqn{S_y = \delta \cdot S_D + \beta \cdot S_x + \gamma \cdot x +
#' \lambda \cdot Z + u}{S[y] = \delta * S[D] + \beta * S[x] + \gamma * x +
#' \lambda * Z + u}
#' where \eqn{S_y}{S[y]}, \eqn{S_x}{S[x]} and \eqn{S_D}{S[D]} are the cumulated
#' sums of \eqn{y}, \eqn{x} and \eqn{D} (with \eqn{D} as the deterministics
#' matrix) and \eqn{Z} as defined in equation (19) in Vogelsang and Wagner
#' (2015).
#' Then \eqn{\theta = (\delta', \beta', \gamma', \lambda')'} is the full
#' parameter vector.
#'
#' @return [\code{cointReg}]. List with components:
#' \describe{
#'   \item{\code{delta} [\code{numeric}]}{
#'     coefficients of the deterministics (cumulative sum \eqn{S_{deter}})}
#'
#'   \item{\code{beta} [\code{numeric}]}{
#'     coefficients of the regressors (cumulative sum \eqn{S_{x}})}
#'
#'   \item{\code{gamma} [\code{numeric}]}{
#'     coefficients of the regressors (original regressors \eqn{x})}
#'
#'   \item{\code{theta} [\code{numeric}]}{
#'     combined coefficients of \code{beta}, \code{delta}}
#'
#'   \item{\code{sd.theta} [\code{numeric}]}{
#'     standard errors for the \code{theta} coefficients}
#'
#'   \item{\code{t.theta} [\code{numeric}]}{
#'     t-values for the \code{theta} coefficients}
#'
#'   \item{\code{p.theta} [\code{numeric}]}{
#'     p-values for the \code{theta} coefficients}
#'
#'   \item{\code{theta.all} [\code{numeric}]}{
#'     combined coefficients of \code{beta}, \code{delta}, \code{gamma}}
#'
#'   \item{\code{residuals} [\code{numeric}]}{
#'     IM-OLS residuals. Attention: These are the first differences of
#'     \eqn{S_u} -- the original residuals are stored in \code{u.plus}.}
#'
#'   \item{\code{u.plus} [\code{numeric}]}{
#'     IM-OLS residuals, not differenced. See \code{residuals} above.}
#'
#'   \item{\code{omega.u.v} [\code{numeric}]}{
#'     conditional long-run variance based on OLS residuals, via
#'     \code{cointRegFM} (in case of argument \code{t.test} is \code{TRUE})
#'     or \code{NULL}}
#'
#'   \item{\code{varmat} [\code{matrix}]}{
#'     variance-covariance matrix}
#'
#'   \item{\code{Omega} [\code{matrix}]}{
#'     \code{NULL} (no long-run variance matrix for this regression type)}
#'
#'   \item{\code{bandwidth} [\code{list}]}{
#'     \code{number} and \code{name} of bandwidth if \code{t.test = TRUE}}
#'
#'   \item{\code{kernel} [\code{character}]}{
#'     abbr. name of kernel type if \code{t.test = TRUE}}
#'
#'   \item{\code{delta2} [\code{numeric}]}{
#'     coefficients of the deterministics (cumulative sum \eqn{S_{deter}})
#'     for regression type 2}
#'
#'   \item{\code{beta2} [\code{numeric}]}{
#'     coefficients of the regressors (cumulative sum \eqn{S_{x}})
#'     for regression type 2}
#'
#'   \item{\code{gamma2} [\code{numeric}]}{
#'     coefficients of the regressors (original regressors \eqn{x})
#'     for regression type 2}
#'
#'   \item{\code{lambda2} [\code{numeric}]}{
#'     coefficients of the Z regressors for regression type 2}
#'
#'   \item{\code{theta2} [\code{numeric}]}{
#'     combined coefficients of \code{beta2}, \code{delta2}, \code{gamma2} and
#'     \code{lambda2} for regression type 2}
#'
#'   \item{\code{u.plus2} [\code{numeric}]}{
#'     IM-OLS residuals for regression type 2}
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
#' test = cointRegIM(x, y, deter, selector = c(1, 2), t.test = TRUE,
#'                     kernel = "ba", bandwidth = "and")
#' print(test)
#'
#' @export


cointRegIM <- function(x, y, deter, selector = 1, t.test = TRUE,
                       kernel = c("ba", "pa", "qs", "tr"),
                       bandwidth = c("and", "nw"), check = TRUE, ...) {

  y.name <- deparse(substitute(y))
  x.name <- deparse(substitute(x))
  d.name <- deparse(substitute(deter))
  mod.name <- paste0(y.name, " ~ ",
                     ifelse(missing(deter) || is.null(deter), "",
                            paste0(d.name, " + ")), x.name)

  if (check) {
    env <- environment()
    checkVars(y = y, selector = selector, t.test = t.test,
              kernel = kernel, bandwidth = bandwidth, .env = env)
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

  y.T <- nrow(y)
  x.k <- ncol(x)
  d.k <- ncol(deter)

  Sy <- matrixStats::colCumsums(y)
  Sx <- matrixStats::colCumsums(x)
  if (d.k == 0) {
    Sd <- deter
  } else {
    Sd <- matrixStats::colCumsums(deter)
    if (!is.null(colnames(deter))) {
      colnames(Sd) <- paste0("S.", colnames(deter))
    } else {
      colnames(Sd) <- make.unique(rep("S.d", d.k))
    }
  }

  if (!is.null(colnames(y))) {
    colnames(Sy) <- paste0("S.", colnames(y))
  } else {
    colnames(Sy) <- "S.y"
  }
  if (!is.null(colnames(x))) {
    colnames(Sx) <- paste0("S.", colnames(x))
  } else {
    colnames(Sx) <- make.unique(rep("S.x", x.k))
  }

  X <- data.frame(Sd, Sx, x)
  Xmat <- as.matrix(X)

  # Regression 1
  if (1 %in% selector) {
    im.mod1 <- lm(Sy ~ . - 1, data = X)
    theta1 <- im.mod1$coeff
    theta <- theta1[1:(d.k + x.k)]
    u1 <- im.mod1$resid
    delta1 <- theta1[0:d.k]
    beta1 <- theta1[(d.k + 1):(d.k + x.k)]
    gamma1 <- theta1[(d.k + x.k + 1):(d.k + 2 * x.k)]
  } else {
    delta1 <- beta1 <- gamma1 <- theta1 <- NULL
    u1 <- NULL
  }

  # Regression 2
  if (2 %in% selector) {
    Sxrow <- nrow(Sx)
    InvBase <- Sx[Sxrow:1, , drop = FALSE]
    InvBaseCum <- matrixStats::colCumsums(InvBase)
    RevBase <- InvBaseCum[Sxrow:1, , drop = FALSE]
    Zreg <- matrixStats::colCumsums(RevBase)
    colnames(Zreg) <- make.unique(rep("Z", x.k))

    X2 <- data.frame(X, Zreg)
    im.mod2 <- lm(Sy ~ . - 1, data = X2)
    theta2 <- im.mod2$coeff
    u2 <- im.mod2$resid

    delta2 <- theta2[0:d.k]
    beta2 <- theta2[(d.k + 1):(d.k + x.k)]
    gamma2 <- theta2[(d.k + x.k + 1):(d.k + 2 * x.k)]
    lambda2 <- theta2[(d.k + 2 * x.k + 1):(d.k + 3 * x.k)]
  } else {
    delta2 <- beta2 <- gamma2 <- lambda2 <- theta2 <- NULL
    u2 <- NULL
  }

  # VCV estimation
  S <- matrixStats::colCumsums(Xmat)
  DS <- rbind(S[y.T, ], t(S[y.T, ] - t(S[1:(y.T - 1), ])))
  Ha <- trySolve(t(Xmat) %*% Xmat)
  H <- Ha %*% t(DS)
  V <- H %*% t(H)

  # t tests
  if (t.test & (1 %in% selector)) {
    mod.FM <- cointRegFM(x = x, y = y, deter = deter, bandwidth = bandwidth,
                           kernel = kernel, demeaning = FALSE, check = FALSE,
                           ...)
    bw <- mod.FM$bandwidth$number
    band <- mod.FM$bandwidth$name
    omega.u.v <- as.numeric(mod.FM$omega.u.v)
    sd.theta1 <- sqrt(omega.u.v * diag(V))
    sd.theta <- sd.theta1[1:(d.k + x.k)]
    t.theta1 <- as.numeric(theta1 / sd.theta1)
    t.theta <- t.theta1[1:(d.k + x.k)]
    df <- y.T - x.k - d.k
    p.theta1 <- 2 * pt(-abs(t.theta1), df = df)
    p.theta <- p.theta1[1:(d.k + x.k)]
  } else {
    bw <- NULL
    band <- NULL
    omega.u.v <- NULL
    sd.theta <- NULL
    sd.theta1 <- NULL
    t.theta <- NULL
    t.theta1 <- NULL
    p.theta <- NULL
    p.theta1 <- NULL
  }

  u.diff <- c(NA, diff(u1))

  out <- list(delta = delta1, beta = beta1, gamma = gamma1, theta = theta,
              sd.theta = sd.theta, t.theta = t.theta, p.theta = p.theta,
              theta.all = theta1, residuals = u.diff, u.plus = u1,
              omega.u.v = omega.u.v, varmat = V, Omega = NULL,
              bandwidth = list(name = band, number = bw), kernel = kernel,
              delta2 = delta2, beta2 = beta2, gamma2 = gamma2,
              lambda2 = lambda2, theta2 = theta2, u.plus2 = u2,
              mod = "IM", name = mod.name)

  class(out) <- "cointReg"

  return(out)
}
