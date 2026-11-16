#' The spread formula
#'
#' `~ 0 + RHS`, which has no intercept, or `~ 1 + RHS` with `intercept = TRUE`
#' (`fireSenseUtils::spreadInterceptTxt` is its coefficient in `fireSense_spreadFit`). Without an
#' intercept the string is what it always was, so cache keys that hold it do not move.
#'
#' @param RHS character; the covariate terms joined with `" + "`. Empty (`character()`) is the null
#'   model, `~ 1`, which needs `intercept = TRUE`.
#' @param intercept logical; `TRUE` for a formula with an intercept.
#' @return a character string.
spreadFormula <- function(RHS, intercept = FALSE) {
  if (!length(RHS)) {
    stopifnot(isTRUE(intercept))
    return("~ 1")
  }
  paste0("~ ", if (isTRUE(intercept)) "1" else "0", " + ", RHS)
}

#' Check the `spreadCovariates` parameter
#'
#' `"all"` (the formula has the covariates) or `"none"` (the null model, `~ 1`), which has to have an
#' intercept: without one it would have no coefficient at all.
#'
#' @param spreadCovariates the module's `spreadCovariates`.
#' @param spreadIntercept the module's `spreadIntercept`.
#' @return `NULL`, invisibly; stops on a value that is not allowed.
checkSpreadCovariates <- function(spreadCovariates, spreadIntercept) {
  if (!is.character(spreadCovariates) || length(spreadCovariates) != 1L || !spreadCovariates %in% c("all", "none"))
    stop("fireSense_dataPrepFit: parameter 'spreadCovariates' must be \"all\" or \"none\"; got: ",
         paste(format(spreadCovariates), collapse = ", "))
  if (identical(spreadCovariates, "none") && !isTRUE(spreadIntercept))
    stop("fireSense_dataPrepFit: spreadCovariates = \"none\" fits the formula `~ 1`, which needs an intercept: ",
         "set spreadIntercept = TRUE")
  invisible(NULL)
}
