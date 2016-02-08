#' Long-Run Variance
#'
#' This function computes the long-run variance Omega,
#' the one sided long-run variance Delta (starting with lag 0)
#' and the variance Sigma from an input matrix of residuals.
#'
#' @param u [\code{numeric} | \code{matrix}]\cr
#'   Data on which to apply the calculation of the long-run variance.
#'
#' @param bandwidth [\code{numeric(1)}]\cr
#'   The bandwidth to use for calculating the long-run variance
#'   as a positive intergerish value.
#'
#' @param demeaning [\code{logical}]\cr
#'   Demeaning of the data before the calculation (default is \code{FALSE}).
#'
#' @inheritParams getBandwidth
#' @inheritParams cointRegFM
#'
#' @details
#' The bandwidth can be one of the following:
#' \itemize{
#'   \item \code{"ba"}: Bartlett kernel
#'   \item \code{"bo"}: Bohmann kernel
#'   \item \code{"da"}: Daniell kernel
#'   \item \code{"pa"}: Parzen kernel
#'   \item \code{"qs"}: Quadratic Spectral kernel
#'   \item \code{"tr"}: Truncated kernel
#' }
#'
#' @return [\code{list}] with components:
#' \describe{
#'   \item{\code{Omega} [\code{matrix}]}{
#'     Long-run variance matrix}
#'
#'   \item{\code{Delta} [\code{matrix}]}{
#'     One-sided long-run variance matrix}
#'
#'   \item{\code{Sigma} [\code{matrix}]}{
#'     Variance matrix}
#' }
#'
#' @seealso \code{\link{getBandwidth}}
#'
#' @examples
#' set.seed(1909)
#' x <- rnorm(100)
#' band <- getBandwidthAnd(x, kernel = "ba")
#' getLongRunVar(x, kernel = "ba", bandwidth = band)
#' # shorter:
#' getLongRunVar(x, kernel = "ba", bandwidth = "and")
#'
#' x2 <- arima.sim(model = list(ar = c(0.7, 0.2)), innov = x, n = 100)
#' x2 <- cbind(a = x2, b = x2 + rnorm(100))
#' getLongRunVar(x2, kernel = "ba", bandwidth = "nw")
#'
#' @export

getLongRunVar <- function(u, bandwidth = c("and", "nw"),
                          kernel = c("ba", "bo", "da", "pa", "qs", "tr"),
                          demeaning = FALSE, check = TRUE, ...) {

  # Check arguments
  if (check) {
    env <- environment()
    checkVars(kernel = kernel, bandwidth = bandwidth, demeaning = demeaning,
              .env = env)
    u <- checkObject(x.coint = u)
  }

  u.T <- nrow(u)

  if(demeaning) {
    u <- apply(u, 2, function(x) x - mean(x))
  }

  if(!is.numeric(bandwidth)) {
    bw <- getBandwidth(u, bandwidth = bandwidth, kernel = kernel,
                       check = FALSE, ...)
  } else {
    bw <- bandwidth
  }

  weights <- getLongRunWeights(n = u.T, kernel = kernel, bandwidth = bw)
  w <- weights[[1]]
  upper <- weights[[2]]

  Omega <- Delta <- Sigma <- t(u) %*% u / u.T

  for (j in 1:upper) {
    matj <- t(u[(j + 1):u.T, , drop = FALSE]) %*%
      u[1:(u.T - j), , drop = FALSE] / u.T
    Omega <- Omega + w[j] * (matj + t(matj))
    Delta <- Delta + w[j] * matj
  }

  return(list(Omega = Omega, Delta = Delta, Sigma = Sigma))
}



#' Weights for Long-Run Variance
#'
#' Compute the weights corresponding to some kernel funtions.
#'
#' @param n [\code{numeric(1)}]\cr
#'   Length of weights' vector.
#'
#' @param bandwidth [\code{numeric(1)}]\cr
#'   The bandwidth (as number).
#'
#' @param kernel [\code{character(1)}]\cr
#'   The kernel function (see \code{\link{getLongRunVar}} for possible values).
#'
#' @return [\code{list}] with components:
#' \describe{
#'   \item{\code{w} [\code{numeric}]}{
#'     Vector of weights}
#'
#'   \item{\code{upper} [\code{numeric(1)}]}{
#'     Index to largest non-zero entry in w}
#' }
#'
#' @examples
#' lrw.ba = cointReg:::getLongRunWeights(100, kernel = "ba", bandwidth = 25)
#' plot(lrw.ba$w)
#'
#' @seealso \code{\link{getLongRunVar}}
#'

getLongRunWeights <- function(n, bandwidth, kernel) {

  w <- numeric(n - 1)
  bw <- bandwidth

  if (kernel == "tr") {
    w <- w + 1
    upper <- min(bw, n - 1)
  }

  else if (kernel == "ba") {
    upper <- ceiling(bw) - 1
    if (upper > 0) {
      j <- 1:upper
    } else {
      j <- 1
    }
    w[j] <- 1 - j / bw
  }

  else if (kernel == "pa") {
    upper1 <- floor(bw / 2)
    if (upper1 > 0) {
      j <- 1:upper1
    } else {
      j <- 1
    }
    jj <- j / bw
    w[j] <- 1 - 6 * jj^2 + 6 * jj^3
    j2 <- (floor(bw / 2) + 1):bw
    jj2 <- j2 / bw
    w[j2] <- 2 * (1 - jj2)^3
    upper <- ceiling(bw) - 1
  }

  else if (kernel == "bo") {
    upper <- ceiling(bw) - 1
    if (upper > 0) {
      j <- 1:upper
    } else {
      j <- 1
    }
    jj <- j / bw
    w[j] <- (1 - jj) * cos(pi * jj) + sin(pi * jj) / pi
  }

  else if (kernel == "da") {
    upper <- n - 1
    j <- 1:upper
    w[j] <- sin(pi * j / bw) / (pi * j / bw)
  }

  else if (kernel == "qs") {
    sc <- 1.2 * pi
    upper <- n - 1
    j <- 1:upper
    jj <- j / bw
    w[j] <- 25 / (12 * pi^2 * jj^2) * (sin(sc * jj) / (sc * jj) - cos(sc * jj))
  }

  if (upper <= 0) upper <- 1

  return(list(w = w, upper = upper))
}
