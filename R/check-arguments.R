#' Multiple variable checks for certain values.
#'
#' Checking the arguments and convert them for internal use, if necessary.
#'
#' @param ... [\code{any}]\cr
#'   Variables to check, see details.
#'
#' @inheritParams checkObject
#'
#' @return
#'   The checked and converted arguments are assigned to
#'   the given environment (\code{.env}) or invisibly returned as a \code{list}.
#'
#' @details
#' See \code{\link{checkObject}} for details.
#'
#' @family check
#'
#' @examples
#' env = environment()
#' x.data = data.frame(a = 1:10, b = 20:11)
#' y.data = 1:10
#' checkVars(x.coint = x.data, y = y.data, .env = env)
#' x.coint
#' y
#'
#' test = checkVars(x.coint = x.data, y = y.data, out = "return")
#' str(test)
#'
#' # If the variables already have the "right" name,
#' # there's no need to do something like
#' # checkVars(kernel = kernel, bandwidth = bandwidth) -
#' # just call checkVars without specifing the arguments:
#' kernel = "q"
#' bandwidth = "a"
#' (checkVars(kernel, bandwidth, out = "return"))
#'
#' @export

checkVars <- function(..., out = "assign", .env) {
  arg.list <- list(...)
  arg.names <- names(arg.list)
  arg.names2 <- as.character(eval(substitute(alist(...))))
  if(is.null(arg.names)) {
    arg.names <- arg.names2
  } else {
    arg.noname <- which(arg.names == "")
    arg.names[arg.noname] <- arg.names2[arg.noname]
  }
  names(arg.list) <- arg.names
  if ("assign" %in% out) {
    if (missing(.env) || !is.environment(.env)) {
      stop("argument \".env\" is missing or wrong type: Needs an environment.")
    } else {
      arg.list2 <- list(out = out, .env = .env)
    }
  } else {
    arg.list2 <- list(out = out)
  }
  tmp <- mapply(checkObject, obj = arg.list, obj.name = arg.names,
                MoreArgs = arg.list2, SIMPLIFY = FALSE)
  if("return" %in% out) {
    return(invisible(tmp))
  }
}


#' Variable check for single objects.
#'
#' Checking the variable and convert it for internal use, if necessary.
#' (Also used by the \code{cointmonitoR} package.)
#'
#' @param obj [\code{any}]\cr
#'   Variable or value to check and convert.
#'
#' @param obj.name [\code{character(1)}]\cr
#'   Name of the object to check.
#'   If \code{missing}, the name of \code{obj} has to be one of the possible
#'   names (see details).
#'
#' @param ... [\code{any}]\cr
#'   An alternative to the use of the \code{obj} and \code{obj.name} arguments
#'   is to directly give the name and the variable to be checked via
#'   \code{name = variable} arguments (see examples). In the case of more than
#'   one \code{...} argument, \code{\link{checkVars}} will be called internally.
#'
#' @param out [\code{character}]\cr
#'   Whether to \code{"return"} or to \code{"assign"} the checked
#'   (and converted) object. Also possible: \code{c("return", "assign")}.
#'
#' @param .env [\code{environment}]\cr
#'   Environment to which to assign the converted \code{obj} (usually the
#'   same on that contains \code{obj}, if it's a variable).\cr
#'   Required, if argument \code{out} contains \code{"assign"}.
#'
#' @details
#' Possible values of \code{obj.name} to check:
#' \describe{
#'   \item{\code{"y"}, \code{"x.stat"}:}{
#'   Of type \code{numeric}, \code{matrix} or \code{data.frame}.
#'   Only the first row/column will be used.\cr
#'   Converted to object: column matrix}
#'
#'   \item{\code{"y.fm"}, \code{"x.coint"}, \code{"deter"}:}{
#'   Of type \code{numeric}, \code{matrix} or \code{data.frame}.\cr
#'   Converted to object: column matrix}
#'
#'   \item{\code{"m"}:}{
#'   Of type \code{numeric(1)}, has to be greater than 0.}
#'
#'   \item{\code{"model"}:}{
#'   One of \code{c("FM", "D", "IM")}.}
#'
#'   \item{\code{"trend"}:}{
#'   One of \code{c("level", "trend")}.}
#'
#'   \item{\code{"signif.level"}:}{
#'   Of type \code{numeric(1)}, has to be in the interval [0.01, 0.1].}
#'
#'   \item{\code{"return.stats"}, \code{"return.input"},
#'         \code{"demeaning"}, \code{"t.test"}:}{
#'   Converted to object: \code{logical(1)}.}
#'
#'   \item{\code{"kernel"}:}{
#'   One of \code{c("ba", "bo", "da", "pa", "qs", "th", "tr")}.}
#'
#'   \item{\code{"bandwidth"}:}{
#'   One of \code{c("and", "nw")}.}
#'
#'   \item{\code{"selector"}:}{
#'   One or both \code{c(1, 2)}.}
#'
#' }
#'
#' @return
#'   The checked and converted argument is assigned to
#'   the given environment (\code{.env}) and/or returned (depending on the
#'   argument \code{out}).
#'
#' @family check
#'
#' @examples
#' x = matrix(1:20, nrow = 2)
#' x2 = checkObject(x, "x.coint")
#' x2
#'
#' env = environment()
#' y = 1:10
#' checkObject(y, out = "assign", .env = env)
#' y
#'
#' # example for the use of the ... argument:
#' det = rbind(1, 1:10)
#' x3 = sin(10:20)
#' det2 = checkObject(deter = det)
#' det2
#' (checkObject(deter = det, x.stat = x3))
#'
#' @export

checkObject <- function(obj, obj.name, ..., out = "return", .env) {

  arg.list <- list(...)
  arg.names <- names(arg.list)
  if (length(arg.list) > 1) {
    return(checkVars(..., out = out, .env = .env))
  }
  if (length(arg.list) == 1) {
    obj <- arg.list[[1]]
    if (arg.names == "") {
      obj.name <- as.character(eval(substitute(alist(...))))
    } else {
      obj.name <- arg.names[1]
    }
  }

  if (missing(obj.name)) {
    obj.name <- deparse(substitute(obj))
  }

  out <- match.arg(out, c("return", "assign"), several.ok = TRUE)

  if (("assign" %in% out) & (missing(.env) || !is.environment(.env))) {
    stop("argument \".env\" is missing or wrong type: Needs an environment.")
  }

  bk <- c("ba", "bo", "da", "pa", "qs", "th", "tr")
  bd <- c("and", "nw")
  mod <- c("FM", "D", "IM")
  trend <- c("level", "trend")
  sel <- c(1, 2)

  ### check "obj.name" for possible arguments

  poss.args <- c("y", "y.fm", "x.stat", "x.coint", "m", "model", "trend",
                 "signif.level", "return.stats", "return.input", "deter",
                 "kernel", "bandwidth", "demeaning", "t.test", "selector")
  obj.name <- match.arg(obj.name, choices = poss.args)

  ### check "x.coint", "deter"

  if (testChoice(obj.name, c("y.fm", "x.coint", "deter"))) {
    assert(checkNumeric(obj), checkMatrix(obj), checkDataFrame(obj),
           .var.name = obj.name)
    if (testMatrix(obj) || testDataFrame(obj)) {
      if (nrow(obj) < ncol(obj)) {
        tmp <- t(as.matrix(obj))
      } else {
        tmp <- as.matrix(obj)
      }
    } else {
      tmp <- matrix(obj, ncol = 1, dimnames = list(NULL, obj.name))
    }
  }

  ### check "x.stat", "y"

  if (testChoice(obj.name, c("x.stat", "y"))) {
    assert(checkNumeric(obj), checkMatrix(obj), checkDataFrame(obj),
           .var.name = obj.name)
    if (testMatrix(obj) || testDataFrame(obj)) {
      if (nrow(obj) < ncol(obj)) {
        tmp <- t(as.matrix(obj[1, , drop = FALSE]))
      } else {
        tmp <- as.matrix(obj[, 1, drop = FALSE])
      }
      if (nrow(obj) > 1 & ncol(obj) > 1) {
        if (nrow(obj) < ncol(obj)) {
          what <- "rows"
          hmany <- nrow(obj)
        } else {
          what <- "columns"
          hmany <- ncol(obj)
        }
        warning(obj.name, " has to many ", what, " (", hmany,
                ", but may have 1). Only the first one will be used.",
                call. = FALSE)
      }
    } else if (testNumeric(obj)) {
      tmp <- matrix(obj, ncol = 1, dimnames = list(NULL, obj.name))
    }
  }

  ### check "m"

  if (obj.name == "m") {
    assertNumber(obj, lower = 0)
    tmp <- obj
  }

  ### check "model"

  if (obj.name == "model") {
    tmp <- match.arg(obj, mod)
  }

  ### check "trend"

  if (obj.name == "trend") {
    tmp <- match.arg(tolower(obj[1]), trend)
  }

  ### check "signif.level"

  if (obj.name == "signif.level") {
    assertNumber(obj, lower = 0.01, upper = 0.1)
    tmp <- obj
  }

  ### check "kernel"

  if (obj.name == "kernel") {
    tmp <- match.arg(obj[1], bk)
  }

  ### check "bandwidth"

  if (obj.name == "bandwidth") {
    if (is.character(obj)) {
      tmp <- match.arg(tolower(obj[1]), bd)
    } else {
      assertNumber(obj, lower = 0, finite = TRUE)
      tmp <- obj
    }
  }

  ### check "demeaning", "return.stats", "return.input", "t.test"

  if (obj.name == "demeaning" || obj.name == "return.stats" ||
      obj.name == "return.input" || obj.name == "t.test") {
    tmp <- as.logical(obj)
    assert(checkFlag(tmp), .var.name = obj.name)
  }

  ### check "selector"

  if (obj.name == "selector") {
    tmp <- as.numeric(obj)
    assertSubset(tmp, sel)
  }

  if ("assign" %in% out) {
    assign(obj.name, value = tmp, envir = .env)
  }

  if ("return" %in% out) {
    return(invisible(tmp))
  }

}



#' Check list D.options.
#'
#' Checking the list D.options, that is an argument of
#' \code{\link{cointRegD}}.
#'
#' @param n.lag,n.lead [\code{NULL} | \code{numeric(1)}]\cr
#'   Have to be "integerish" and > 0.
#' @param kmax [\code{NULL} | \code{character(1)}]\cr
#'   One of \code{"k4"} or \code{"k12"}.
#' @param info.crit [\code{NULL} | \code{character(1)}]\cr
#'   One of \code{"AIC"} or \code{"BIC"}.
#'
#' @return [\code{list}]. List with the checked and (if necessary)
#'   converted arguments.
#'
#'   If one of \code{n.lag} and \code{n.lead} is
#'   \code{NULL}, only \code{kmax} and \code{info.crit} will be not \code{NULL}.
#'
#' @family check
#'
#' @examples
#' checkDoptions(n.lag = 3, n.lead = 4)
#' checkDoptions(info.crit = "BIC")
#' checkDoptions()
#'
#' # It's not sufficient to include only one of "n.lag" and "n.lead":
#' checkDoptions(n.lag = 2)
#'
#' @export

checkDoptions <- function(n.lag = NULL, n.lead = NULL, kmax = c("k4", "k12"),
                          info.crit = c("AIC", "BIC")) {

  assert(checkNull(n.lag), checkNumber(n.lag, lower = 0))

  assert(checkNull(n.lead), checkNumber(n.lead, lower = 0))

  if (testNumber(n.lead) && testNumber(n.lag)) {
    n.lag <- as.integer(n.lag)
    n.lead <- as.integer(n.lead)
    kmax <- info.crit <- NULL
  } else {
    n.lag <- n.lead <- NULL
    kmax <- match.arg(kmax)
    info.crit <- match.arg(info.crit)
  }

  return(list(n.lag = n.lag, n.lead = n.lead, kmax = kmax,
              info.crit = info.crit))
}
