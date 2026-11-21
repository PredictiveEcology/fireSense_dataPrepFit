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

## the one value `spreadInteractions` takes: each spread climate variable x each fuel covariate
spreadInteractionClimateFuel <- "climate:fuel"

#' The right-hand side of the spread formula
#'
#' The spread climate variable(s), `youngAge`, then the other covariates (fuel, non-forest), joined with
#' `" + "`; then, with `spreadInteractions = "climate:fuel"`, one interaction term (`climate:fuel`) for each climate
#' variable and each fuel covariate (`fireSenseUtils::spreadFuelCols()`: not `youngAge`, `nfLCC_*` or `treedWetland`),
#' climate first, after all the main effects. With `climateSplit` the climate part of each is the between-year term
#' (`C_yr:fuel`); the within-year term (`C_dev`) does not interact. Without interactions the string is what it always was, so the cache keys
#' that hold it do not move.
#'
#' @param spreadCovariates the module's `spreadCovariates`: `"none"` gives `character()`, the null model.
#' @param spreadVars character; the spread climate variable(s).
#' @param vegColsForRHS character; the other covariates (not `youngAge` or the climate variables).
#' @param spreadInteractions the module's `spreadInteractions`; `NULL` or `character(0)` for none.
#' @param climateSplit the module's `spreadClimateSplit`: `TRUE` replaces each climate variable `C` by `C_yr` and
#'   `C_dev` (`fireSenseUtils::climateSplitCols()`).
#' @return a character string, or `character()` for the null model.
spreadRHS <- function(spreadCovariates, spreadVars, vegColsForRHS, spreadInteractions = character(), climateSplit = FALSE) {
  if (identical(spreadCovariates, "none")) return(character())
  fuel <- if (spreadInteractionClimateFuel %in% spreadInteractions) fireSenseUtils::spreadFuelCols(vegColsForRHS)
  climCols <- if (isTRUE(climateSplit)) fireSenseUtils::climateSplitCols(spreadVars) else spreadVars
  ## a split climate variable interacts through its between-year term only, never `_dev`
  interClim <- if (isTRUE(climateSplit)) paste0(spreadVars, fireSenseUtils::climateSplitYrSuffix) else spreadVars
  paste(c(paste0(climCols, collapse = " + "), fireSenseUtils::youngAgeTxt,
          paste0(vegColsForRHS, collapse = " + "),
          paste(rep(interClim, each = length(fuel)), rep(fuel, times = length(interClim)), sep = ":")),
        collapse = " + ")
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

#' Check the `spreadInteractions` parameter
#'
#' `character(0)` (none; `NULL` too) or `"climate:fuel"` (an interaction of each spread climate variable with each fuel
#' covariate). Interactions are terms of the covariates, so they need covariates: `spreadCovariates = "none"` (the
#' null model) has none.
#'
#' @param spreadInteractions the module's `spreadInteractions`.
#' @param spreadCovariates the module's `spreadCovariates`.
#' @return `NULL`, invisibly; stops on a value that is not allowed.
checkSpreadInteractions <- function(spreadInteractions, spreadCovariates) {
  if (length(spreadInteractions) && (!is.character(spreadInteractions) || !all(spreadInteractions %in% spreadInteractionClimateFuel)))
    stop("fireSense_dataPrepFit: parameter 'spreadInteractions' must be character(0) or \"", spreadInteractionClimateFuel, "\"; got: ",
         paste(format(spreadInteractions), collapse = ", "))
  if (length(spreadInteractions) && identical(spreadCovariates, "none"))
    stop("fireSense_dataPrepFit: spreadInteractions = \"", spreadInteractionClimateFuel, "\" needs covariates: ",
         "the null model (spreadCovariates = \"none\") has none")
  invisible(NULL)
}
