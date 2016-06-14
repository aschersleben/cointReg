#' @importFrom graphics plot.default
#' @importFrom stats lm
#' @importFrom stats pt
#' @import checkmate
#' @importFrom matrixStats colCumsums
#' @importFrom matrixStats colDiffs
#' @importFrom MASS ginv

trySolve <- function(x) {
  out <- try(solve(x), silent = TRUE)
  if (inherits(out, "try-error")) {
    out <- MASS::ginv(x)
  }
  return(out)
}

.onAttach <- function(libname, pkgname) {
  psm = paste("cointReg", paste0("(v", utils::packageVersion("cointReg"), "):"),
              "Parameter Estimation and Inference in a",
              "Cointegrating Regression.")
  packageStartupMessage(psm)
}
