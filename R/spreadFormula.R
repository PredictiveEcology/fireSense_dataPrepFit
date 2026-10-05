#' The spread formula
#'
#' `~ 0 + RHS`, which has no intercept, or `~ 1 + RHS` with `intercept = TRUE`
#' (`fireSenseUtils::spreadInterceptTxt` is its coefficient in `fireSense_spreadFit`). Without an
#' intercept the string is what it always was, so cache keys that hold it do not move.
#'
#' @param RHS character; the covariate terms joined with `" + "`.
#' @param intercept logical; `TRUE` for a formula with an intercept.
#' @return a character string.
spreadFormula <- function(RHS, intercept = FALSE) {
  paste0("~ ", if (isTRUE(intercept)) "1" else "0", " + ", RHS)
}
