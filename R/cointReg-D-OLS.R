#' Dynamic OLS
#'
#' Computes the Saikkonen (1990) Dynamic OLS estimator.
#'
#' @param x [\code{numeric} | \code{matrix} | \code{data.frame}]\cr
#'   RHS variables on which to apply the D-OLS estimation (see Details).
#'
#' @param y [\code{numeric} | \code{matrix} | \code{data.frame}]\cr
#'   LHS variable(s) on which to apply the D-OLS estimation (see Details).
#'   Has to be one-dimensional. If \code{matrix}, it may
#'   have only one row or column, if \code{data.frame} just one column.
#'
#' @param n.lead,n.lag [\code{numeric(1)} | \code{NULL}]\cr
#'   Numbers of Leads and Lags (see Details). Default is \code{NULL}.
#'
#' @param kmax [\code{character(1)}]\cr
#'   Maximal value for lags and leads if generated automatically (see Details).
#'   Default is \code{"k4"}.
#'
#' @param info.crit [\code{character(1)}]\cr
#'   Information criterion to use for the automatical calculation of lags and
#'   leads. Default is \code{"AIC"}.
#'
#' @inheritParams cointRegFM
#'
#' @details
#' The equation for which the FM-OLS estimator is calculated:
#' \deqn{y = \delta \cdot D + \beta \cdot x + u}{y = \delta * D + \beta * x + u}
#' with \eqn{D} as the deterministics matrix.
#' Then \eqn{\theta = (\delta', \beta')'} is the full parameter vector.
#'
#' Information about the D-OLS specific arguments:
#' \describe{
#'   \item{\code{n.lag}, \code{n.lead}}{A positive number to set the number
#'   of lags and leads. If at least one of them is equal to \code{NULL}
#'   (default), the function \code{\link{getLeadLag}} will be used to
#'   calculate them automatically (see Choi and Kurozumi (2012) for details).
#'   In that case, the following two arguments are needed.}
#'
#'   \item{\code{kmax}}{Maximal value for lags and leads, when they are
#'   calculated automatically. If "k4", then the maximum is equal to
#'   \code{floor(4 * (x.T/100)^(1/4))}, else it's
#'   \code{floor(12 * (x.T/100)^(1/4))} with \code{x.T} is equal
#'   to the data's length. One of \code{"k4"} or \code{"k12"}.
#'   Default is \code{"k4"}.}
#'
#'   \item{\code{info.crit}}{Information criterion to use for the automatical
#'   calculation of lags and leads. One of \code{"AIC"} or \code{"BIC"}.
#'   Default is \code{"AIC"}.}
#' }
#'
#'
#' @return [\code{cointReg}]. List with components:
#' \describe{
#'   \item{\code{beta} [\code{numeric}]}{
#'     coefficients of the regressors}
#'
#'   \item{\code{delta} [\code{numeric}]}{
#'     coefficients of the deterministics}
#'
#'   \item{\code{theta} [\code{numeric}]}{
#'     combined coefficients of \code{beta} and \code{delta}}
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
#'   \item{\code{theta.all} [\code{numeric}]}{
#'     combined coefficients of \code{beta}, \code{delta} and the auxiliary
#'     leads-and-lags regressors}
#'
#'   \item{\code{residuals} [\code{numeric}]}{
#'     D-OLS residuals (length depends on leads and lags)}
#'
#'   \item{\code{omega.u.v} [\code{numeric}]}{
#'     conditional long-run variance based on OLS residuals}
#'
#'   \item{\code{varmat} [\code{matrix}]}{
#'     variance-covariance matrix}
#'
#'   \item{\code{Omega} [\code{list}]}{
#'     the whole long-run variance matrix and parts of it}
#'
#'   \item{\code{bandwidth} [\code{list}]}{
#'     number and name of the calculated bandwidth}
#'
#'   \item{\code{kernel} [\code{character}]}{
#'     abbr. name of kernel type}
#'
#'   \item{\code{lead.lag} [\code{list}]}{
#'     leads-and-lags parameters}
#' }
#'
#' @family D-OLS
#' @family cointReg
#'
#' @references
#'   \itemize{
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
#'   }
#'
#' @examples
#' set.seed(1909)
#' x1 <- cumsum(rnorm(100, mean = 0.05, sd = 0.1))
#' x2 <- cumsum(rnorm(100, sd = 0.1)) + 1
#' x3 <- cumsum(rnorm(100, sd = 0.2)) + 2
#' x <- cbind(x1, x2, x3)
#' y <- x1 + x2 + x3 + rnorm(100, sd = 0.2) + 1
#' deter <- cbind(level = 1, trend = 1:100)
#' test <- cointRegD(x, y, deter, n.lead = 2, n.lag = 2,
#'                     kernel = "ba", bandwidth = "and")
#' print(test)
#' test2 <- cointRegD(x, y, deter, kmax = "k4", info.crit = "BIC",
#'                      kernel = "ba", bandwidth = "and")
#' print(test2)
#'
#' @export


cointRegD <- function(x, y, deter, kernel = c("ba", "pa", "qs", "tr"),
                      bandwidth = c("and", "nw"), n.lead = NULL, n.lag = NULL,
                      kmax = c("k4", "k12"), info.crit = c("AIC", "BIC"),
                      demeaning = FALSE, check = TRUE, ...) {

  y.name <- deparse(substitute(y))
  x.name <- deparse(substitute(x))
  d.name <- deparse(substitute(deter))
  mod.name <- paste0(y.name, " ~ ",
                     ifelse(missing(deter) || is.null(deter), "",
                            paste0(d.name, " + ")), x.name)

  if (check) {
    env <- environment()
    checkVars(y = y, kernel = kernel, bandwidth = bandwidth,
              demeaning = demeaning, .env = env)
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
  D.opts <- checkDoptions(n.lead = n.lead, n.lag = n.lag, kmax = kmax,
                          info.crit = info.crit)

  x.T <- nrow(x)
  x.k <- ncol(x)
  y.k <- ncol(y)
  d.k <- ncol(deter)

  if (is.null(D.opts$n.lag)) {
    kmax <- switch(D.opts$kmax, NULL,
                   k4 = floor(4 * (x.T / 100)^(1/4)),
                   k12 = floor(12 * (x.T / 100)^(1/4)))
    LagLead <- getLeadLag(x = x, y = y, deter = deter, max.lag = kmax,
                          max.lead = kmax, ic = D.opts$info.crit,
                          symmet = FALSE, check = check)
    n.lag <- D.opts$n.lag <- LagLead["n.lag"]
    n.lead <- D.opts$n.lead <- LagLead["n.lead"]
  } else {
    n.lag <- D.opts$n.lag
    n.lead <- D.opts$n.lead
  }

  d.mod <- getModD(x = x, y = y, deter = deter, n.lag = n.lag, n.lead = n.lead,
                   check = FALSE)

  theta.d.all <- d.mod$coefficients
  theta.d <- theta.d.all[1:(d.k + x.k)]
  delta.d <- theta.d[0:d.k]
  beta.d <- theta.d[(d.k + 1):(d.k + x.k)]
  names(beta.d) <- names(delta.d) <- names(theta.d) <- NULL
  u.dols <- d.mod$residuals
  u <- matrix(u.dols)

  if(!is.numeric(bandwidth)) {
    bw <- getBandwidth(u, bandwidth = bandwidth, kernel, check = FALSE, ...)
    band <- switch(bandwidth, and = "Andrews", nw = "Newey-West")
  } else {
    bw <- bandwidth
    band <- "set by user"
  }

  lrvar <- getLongRunVar(u, bandwidth = bw, kernel = kernel,
                         demeaning = demeaning, check = FALSE)
  Omega.dols <- lrvar[[1]]

  if (n.lag + n.lead == 0) {
    u.4var <- u[-1, , drop = FALSE]
    u.plus <- u.dols
  } else {
    u.4var <- u.dols
    u.plus <- c(NA, rep(NA, n.lag), u.dols, c(rep(NA, n.lead)))
  }
  u.omega <- cbind(u.4var,
                   colDiffs(x)[(n.lag + 1):(x.T - 1 - n.lead), , drop = FALSE])
  lrvar <- getLongRunVar(u.omega, bandwidth = bw, kernel = kernel,
                         demeaning = demeaning, check = FALSE)
  tmp <- lapply(lrvar, function(x) {
    out <- list()
    out[["all"]] <- x
    out[["uu"]] <- x[1:y.k, 1:y.k, drop = FALSE]
    out[["uv"]] <- x[1:y.k, (y.k + 1):(y.k + x.k), drop = FALSE]
    out[["vu"]] <- x[(y.k + 1):(y.k + x.k), 1:y.k, drop = FALSE]
    out[["vv"]] <- x[(y.k + 1):(y.k + x.k), (y.k + 1):(y.k + x.k), drop = FALSE]
    return(out)
  })
  Omega <- tmp[[1]]
  Omegavv.inv <- trySolve(Omega$vv)
  Omega.u.v <- Omega$uu - (Omega$uv %*% Omegavv.inv %*% Omega$vu)

  if (y.k == 1) {
    all.trunc <- d.mod$aux$all.trunc
    alltrunc.inv <- trySolve(t(all.trunc) %*% all.trunc)
    varmat <- Omega.dols[1, 1] * alltrunc.inv
    sd.theta.all <- sqrt(diag(varmat))
    sd.theta <- sd.theta.all[1:(d.k + x.k)]
    t.theta.all <- theta.d.all / sd.theta.all
    t.theta <- t.theta.all[1:(d.k + x.k)]
    df <- x.T - x.k - d.k
    p.theta.all <- 2 * pt(-abs(t.theta.all), df = df)
    p.theta <- p.theta.all[1:(d.k + x.k)]
  } else {
    varmat <- NULL
    sd.theta <- NULL
    t.theta <- NULL
    p.theta <- NULL
  }

  out <- list(beta = beta.d, delta = delta.d, theta = theta.d,
              sd.theta = sd.theta, t.theta = t.theta, p.theta = p.theta,
              theta.all = theta.d.all,
              residuals = u.plus, omega.u.v = Omega.u.v, varmat = varmat,
              Omega = Omega, bandwidth = list(name = band, number = bw),
              kernel = kernel, lead.lag = D.opts, mod = "D", name = mod.name)

  class(out) <- "cointReg"

  return(out)
}



#' Get D OLS model.
#'
#' Generates an \code{lm} model for the Dynamic OLS estimator.
#'
#' @param x [\code{matrix}]\cr
#'   RHS variables of the D OLS estimation.
#' @param y [\code{matrix}]\cr
#'   LHS variable(s) of the D OLS estimation.
#' @param deter [\code{matrix}]\cr
#'   Matrix of deterministic variables to include in the D OLS model.
#' @param n.lag,n.lead [\code{numeric(1)}]\cr
#'   Number of lags and leads, have to be non-negative integer values.
#' @inheritParams cointRegFM
#'
#' @return [\code{lm}]. An \code{\link{lm}} object, containing an additional
#'   list element (\code{aux}) with D-OLS specific objects:
#' \describe{
#'   \item{\code{Z} [\code{matrix}]}{
#'     jointed matrix of deterministics and x}
#'
#'   \item{\code{x.delta} [\code{matrix}]}{
#'     differences of x}
#'
#'   \item{\code{dx.all} [\code{matrix}]}{
#'     leads-and-lags matrix}
#'
#'   \item{\code{all.trunc} [\code{matrix}]}{
#'     truncated version of jointed matrix of \code{Z} and \code{dx.all}}
#'
#'   \item{\code{y.trunc} [\code{matrix}]}{
#'     truncated version of \code{y}}
#'
#' }
#'
#' @family D-OLS
#'
#' @examples
#' set.seed(1909)
#' y <- matrix(cumsum(rnorm(100)), ncol = 1)
#' x <- matrix(rep(y, 4) + rnorm(400, mean = 3, sd = 2), ncol = 4)
#' deter <- cbind(1, 1:100)
#' cointReg:::getModD(x = x, y = y, deter = deter, n.lag = 2, n.lead = 3)

getModD <- function(x, y, deter, n.lag, n.lead, check = FALSE) {

  x.name <- deparse(substitute(x))
  d.name <- deparse(substitute(deter))

  if (check) {
    y <- checkObject(y = y)
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

  Z <- cbind(deter, x)
  x.delta <- colDiffs(x)

  if(n.lag + n.lead == 0) {
    all.trunc <- Z
    y.trunc <- y
    T1 <- x.T - 1
    dx.all <- NULL
  } else {
    dx.all <- makeLeadLagMatrix(x.delta, n.lag, n.lead)
    Zs <- Z[-1, , drop = FALSE]
    all.untrunc <- cbind(Zs, dx.all)
    T1 <- nrow(all.untrunc)
    all.trunc <- all.untrunc[(n.lag + 1):(T1 - n.lead), , drop = FALSE]
    ys <- y[-1, , drop = FALSE]
    T2 <- nrow(ys)
    y.trunc <- ys[(n.lag + 1):(T2 - n.lead), , drop = FALSE]
  }

  dols.mod <- lm(y.trunc ~ all.trunc - 1)
  dols.mod$aux <- list(Z = Z, x.delta = x.delta, dx.all = dx.all,
                       all.trunc = all.trunc, y.trunc = y.trunc)
  return(dols.mod)
}



#' Leads-and-Lags Matrix
#'
#' Generates leads-and-lags matrix for the Dynamic OLS estimator.
#'
#' @param x [\code{matrix}]\cr
#'   Matrix for which to generate the leads-and-lags matrix.
#'
#' @param n.lag,n.lead [\code{numeric(1)}]\cr
#'   Number of lags and leads, have to be non-negative integer values.
#'   If greater than \code{nrow(x)}, produces 0-rows.
#'
#' @return [\code{matrix}]. Leads-and-lags matrix.
#'
#' @family D-OLS
#'
#' @examples
#' x <- matrix(1:20, 2, byrow = TRUE)
#' cointReg:::makeLeadLagMatrix(x = x, n.lag = 2, n.lead = 3)

makeLeadLagMatrix <- function(x, n.lag, n.lead) {

  n <- nrow(x)
  k <- ncol(x)
  x.names <- paste0("xD", 1:k)

  if (n.lag != 0) {
    lag.mat <- do.call(cbind, sapply(1:n.lag, function(i) {
      rbind(matrix(0, i, k), x)[ -((n + 1):(n + i)),]
    }, simplify = FALSE))
    lag.names <- rep(paste0("lag", 1:n.lag), each = k)
    lag.names <- paste(x.names, lag.names, sep = ".")
  } else {
    lag.mat <- NULL
    lag.names <- NULL
  }

  if (n.lead != 0) {
    lead.mat <- do.call(cbind, sapply(1:n.lead, function(i) {
      rbind(x, matrix(0, i, k))[-(1:i), ]
    }, simplify = FALSE))
    lead.names <- rep(paste0("lead", 1:n.lead), each = k)
    lead.names <- paste(x.names, lead.names, sep = ".")
  } else {
    lead.mat <- NULL
    lead.names <- NULL
  }

  DeltaMat <- cbind(x, lag.mat, lead.mat)
  colnames(DeltaMat) <- c(x.names, lag.names, lead.names)

  return(DeltaMat)
}



#' Leads and Lags
#'
#' Generates "optimal" numbers of leads and lags for the Dynamic OLS estimator.
#'
#' @param max.lead,max.lag [\code{numeric(1)}]\cr
#'   Maximal numbers of leads and lags, have to be non-negative integer values.
#'
#' @param ic [\code{character(1)}]\cr
#'   Information criterion (one of \code{"AIC"} or \code{"BIC"}).
#'
#' @param symmet [\code{logical(1)}]\cr
#'   If \code{TRUE}, only looks for equal leads and lags.
#'
#' @inheritParams cointRegD
#'
#' @return [\code{numeric(2)}]. "Optimal" numbers of leads and lags.
#'
#' @family D-OLS
#'
#' @examples
#' set.seed(1909)
#' y <- matrix(cumsum(rnorm(100)), ncol = 1)
#' x <- matrix(rep(y, 4) + rnorm(400, mean = 3, sd = 2), ncol = 4)
#' deter <- cbind(1, 1:100)
#' cointReg:::getLeadLag(x = x, y = y, deter = deter, max.lag = 5,
#'                       max.lead = 5, ic = "AIC", symmet = FALSE)

getLeadLag <- function(x, y, deter, max.lag, max.lead, ic = c("AIC", "BIC"),
                       symmet = FALSE, check = FALSE) {

  y.T <- nrow(y)
  x.k <- ncol(x)

  if (symmet) {
    min.ll <- min(max.lag, max.lead)
    leadlag.table <- data.frame(n.lag = 0:min.ll, n.lead = 0:min.ll)
  } else {
    leadlag.table <- data.frame(expand.grid(n.lag = 0:max.lag,
                                            n.lead = 0:max.lead))
  }

  tmp <- apply(leadlag.table, 1, function(llt) {
    n.lag <- llt[1]
    n.lead <- llt[2]
    if ((n.lag + n.lead) == 0) {
      act.samp <- y.T
    } else {
      act.samp <- y.T - 1 - n.lag - n.lead
    }
    if (act.samp <= 0) {
      out <- Inf
    } else {
      d.mod <- getModD(x = x, y = y, deter = deter, n.lag = n.lag,
                       n.lead = n.lead, check = check)
      u.dols <- d.mod$residuals
      SSR <- sum(u.dols^2)
      if (ic == "AIC") {
        out <- act.samp * log(SSR / act.samp) +
          2 * (x.k * (n.lag + n.lead + 2) + 2)
      } else {
        out <- act.samp * log(SSR / act.samp) + log(act.samp) *
          (x.k * (n.lag + n.lead + 2) + 2)
      }
    }
    names(out) <- ic
    return(out)
  })

  res <- cbind(leadlag.table, ic = tmp)

  out <- as.numeric(res[which.min(res[, "ic"]), c("n.lag", "n.lead")])
  names(out) <- c("n.lag", "n.lead")
  return(out)
}
