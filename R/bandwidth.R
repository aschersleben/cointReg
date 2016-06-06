#' Automatic Bandwidth Selection
#'
#' Automatic bandwidth selection of Andrews (1991) and of Newey and West (1994).
#'
#' @param u [\code{numeric}]\cr
#'   Data on which to apply the bandwidth selction.
#'
#' @param bandwidth [\code{character(1)}]\cr
#'   The bandwidth selection method to use. Default is Andrews (1991) (\code{"and"}),
#'   an alternative is Newey West (1994) (\code{"nw"}).
#'
#' @param kernel [\code{character(1)}]\cr
#'   The kernel function to use for selecting the bandwidth.
#'   Default is Bartlett kernel (\code{"ba"}), see Details for alternatives.
#'
#' @param inter [\code{logical}]\cr
#'   The first column will be ignored, if \code{TRUE} (intercept).
#'   Default is \code{FALSE}.
#'
#' @param u.weights [\code{numeric}]\cr
#'   How to weight the columns of \code{u}.
#'   If \code{NULL} (default), uses identical weights for all columns.
#'
#' @inheritParams cointRegFM
#'
#' @param ... Arguments passed to \code{getBandwidthNW}.
#'
#' @details
#' For Andrews (1991), the AR(1) individual version is implemented.
#'
#' The kernel that is used for calculating the long-run variance can be
#' one of the following:
#' \itemize{
#'   \item \code{"ba"}: Bartlett kernel
#'   \item \code{"pa"}: Parzen kernel
#'   \item \code{"qs"}: Quadratic Spectral kernel
#'   \item \code{"th"}: Tukey-Hanning kernel (only if \code{bandwidth = "and"})
#'   \item \code{"tr"}: Truncated kernel (only if \code{bandwidth = "and"})
#' }
#'
#' @return [\code{numeric(1)}]. Bandwidth
#'
#' @seealso \code{\link{getLongRunVar}}
#' @family bandwidth
#'
#' @references
#'   \itemize{
#'     \item Andrews, D.W.K. (1991): "Heteroskedasticity and Autocorrelation
#'     Consistent Covariance Matrix Estimation," \emph{Econometrica}, 59,
#'     817--854, \href{http://dx.doi.org/10.2307/2938229}{DOI:10.2307/2938229}.
#'     \item Newey, W.K. and K.D. West (1994): "Automatic Lag Selection in
#'     Covariance Matrix Estimation", \emph{Review of Economic Studies}, 61,
#'     631--653, \href{http://dx.doi.org/10.2307/2297912}{DOI:10.2307/2297912}.
#'   }
#'
#' @examples
#' set.seed(1909)
#' x <- rnorm(100)
#' getBandwidth(x, kernel = "ba")
#' getBandwidth(x, bandwidth = "nw", kernel = "ba")
#'
#' x2 <- arima.sim(model = list(ar = c(0.7, 0.2)), innov = x, n = 100)
#' getBandwidth(x2, kernel = "qs")
#' getBandwidth(x2, bandwidth = "nw", kernel = "qs")
#'
#' @export


getBandwidth <- function(u, bandwidth = c("and", "nw"), kernel, ...,
                         check = TRUE) {
  bandwidth <- match.arg(bandwidth)
  if (bandwidth == "and") {
    return(getBandwidthAnd(u, kernel = kernel, check = check))
  } else {
    return(getBandwidthNW(u, kernel = kernel, ..., check = check))
  }
}

#' @describeIn getBandwidth
#' Automatic bandwidth selection of Andrews (1991).
#' @export
getBandwidthAnd <- function(u, kernel = c("ba", "pa", "qs", "th", "tr"),
                            check = TRUE) {

  # Check arguments
  if (check) {
    kernel <- match.arg(kernel)
    u <- checkObject(x.coint = u)
  }

  u.T <- nrow(u)
  u.k <- ncol(u)

  rhovec <- sigma2vec <- numeric(u.k)

  for (j in 1:u.k) {
    mod <- lm(u[-1, j] ~ u[-u.T, j] - 1)
    rhovec[j] <- mod$coeff
    sigma2vec[j] <- (1 / u.T) * sum(mod$resid^2)
  }

  denom <- sum(sigma2vec^2 / (1 - rhovec)^4)

  if (kernel == "ba") {
    numer <- sum(4 * rhovec^2 * sigma2vec^2 / ((1 - rhovec)^6 * (1 + rhovec)^2))
  } else {
    numer <- sum(4 * rhovec^2 * sigma2vec^2 / (1 - rhovec)^8)
  }

  a <- numer / denom

  bandwidth <- switch(kernel,
                      ba = 1.1447 * (a * u.T)^(1 / 3),
                      pa = 2.6614 * (a * u.T)^(1 / 5),
                      qs = 1.3221 * (a * u.T)^(1 / 5),
                      th = 1.7462 * (a * u.T)^(1 / 5),
                      tr = 0.661 * (a * u.T)^(1 / 5))

  if(bandwidth > (u.T - 1)) bandwidth <- u.T - 1

  return(bandwidth)
}


#' @describeIn getBandwidth
#' Automatic bandwidth selection of Newey and West (1994).
#' @export
getBandwidthNW <- function(u, kernel = c("ba", "pa", "qs"), inter = FALSE,
                           u.weights = NULL, check = TRUE) {

  # Check arguments
  if (check) {
    kernel <- match.arg(kernel)
    u <- checkObject(x.coint = u)
  }

  u.T <- nrow(u)
  u.k <- ncol(u)

  if (is.null(u.weights)) {
    u.weights <- rep(1, u.k)
    if (inter) u.weights[1] <- 0
  }

  npower <- switch(kernel, ba = 2 / 9, pa = 4 / 25, qs = 2 / 25)
  n <- floor(4 * (u.T / 100)^npower)

  umatw <- as.numeric(u %*% u.weights)

  sigma <- sapply(0:n, function(j) {
    sum(umatw[(j + 1):u.T] * umatw[1:(u.T - j)]) / u.T
  })

  s0 <- sigma[1] + 2 * sum(sigma[-1])

  if(kernel == "ba") {
    s1 <- 2 * sum(1:n * sigma[-1])
    q <- 1
  } else {
    s2 <- 2 * sum((1:n)^2 * sigma[-1])
    q <- 2
  }

  Tpower <- 1 / (2 * q + 1)

  gamma <- switch(kernel,
                  ba = 1.1447 * ((s1 / s0)^2)^Tpower,
                  pa = 2.6614 * ((s2 / s0)^2)^Tpower,
                  qs = 1.3221 * ((s2 / s0)^2)^Tpower)

  bandwidth <- gamma * u.T^Tpower

  return(bandwidth)
}
