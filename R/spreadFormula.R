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

## the values `spreadCovariates` takes: "all" (the default), "none" (the null model), and "fuel" or "climate" (one group of
## "all"'s covariates, for an ANOVA of what drives spread-fit quality)
spreadCovariatesValues <- c("all", "none", "fuel", "climate")

## the one value `spreadInteractions` takes: each spread climate variable x each fuel covariate
spreadInteractionClimateFuel <- "climate:fuel"

#' The right-hand side of the spread formula
#'
#' The spread climate variable(s), `youngAge`, then the other covariates (fuel, non-forest), joined with
#' `" + "`; then, with `spreadInteractions = "climate:fuel"`, one interaction term (`climate:fuel`) for each climate
#' variable and each fuel covariate (`fireSenseUtils::spreadFuelCols()`: not `youngAge`, `nfLCC_*` or `treedWetland`),
#' climate first, after all the main effects. Without interactions the string is what it always was, so the cache keys
#' that hold it do not move.
#'
#' @param spreadCovariates the module's `spreadCovariates`: `"none"` gives `character()`, the null model; `"climate"`
#'   keeps the spread climate variable(s) only; `"fuel"` keeps everything else (`youngAge`, fuel, non-forest: the
#'   vegetation state), so each of the two is a subset of `"all"`'s terms.
#' @param spreadVars character; the spread climate variable(s).
#' @param vegColsForRHS character; the other covariates (not `youngAge` or the climate variables).
#' @param spreadInteractions the module's `spreadInteractions`; `NULL` or `character(0)` for none.
#' @return a character string, or `character()` for the null model.
spreadRHS <- function(spreadCovariates, spreadVars, vegColsForRHS, spreadInteractions = character()) {
  if (identical(spreadCovariates, "none")) return(character())
  fuel <- if (spreadInteractionClimateFuel %in% spreadInteractions) fireSenseUtils::spreadFuelCols(vegColsForRHS)
  climateTerms <- if (!identical(spreadCovariates, "fuel")) paste0(spreadVars, collapse = " + ")
  fuelTerms <- if (!identical(spreadCovariates, "climate"))
    c(fireSenseUtils::youngAgeTxt, paste0(vegColsForRHS, collapse = " + "))
  paste(c(climateTerms, fuelTerms,
          paste(rep(spreadVars, each = length(fuel)), rep(fuel, times = length(spreadVars)), sep = ":")),
        collapse = " + ")
}

#' Check the `spreadCovariates` parameter
#'
#' `"all"` (the formula has the covariates), `"none"` (the null model, `~ 1`), `"climate"` (the climate variable(s)
#' only) or `"fuel"` (`youngAge`, fuel and non-forest only). `"none"` has to have an intercept: without one it would
#' have no coefficient at all.
#'
#' @param spreadCovariates the module's `spreadCovariates`.
#' @param spreadIntercept the module's `spreadIntercept`.
#' @return `NULL`, invisibly; stops on a value that is not allowed.
checkSpreadCovariates <- function(spreadCovariates, spreadIntercept) {
  if (!is.character(spreadCovariates) || length(spreadCovariates) != 1L || !spreadCovariates %in% spreadCovariatesValues)
    stop("fireSense_dataPrepFit: parameter 'spreadCovariates' must be one of ",
         paste0("\"", spreadCovariatesValues, "\"", collapse = ", "), "; got: ",
         paste(format(spreadCovariates), collapse = ", "))
  if (identical(spreadCovariates, "none") && !isTRUE(spreadIntercept))
    stop("fireSense_dataPrepFit: spreadCovariates = \"none\" fits the formula `~ 1`, which needs an intercept: ",
         "set spreadIntercept = TRUE")
  invisible(NULL)
}

#' Check the `spreadInteractions` parameter
#'
#' `character(0)` (none; `NULL` too) or `"climate:fuel"` (an interaction of each spread climate variable with each fuel
#' covariate). Interactions are terms of the covariates, so they need both climate and fuel: `spreadCovariates = "none"` (the
#' null model), `"fuel"` and `"climate"` lack one or both.
#'
#' @param spreadInteractions the module's `spreadInteractions`.
#' @param spreadCovariates the module's `spreadCovariates`.
#' @return `NULL`, invisibly; stops on a value that is not allowed.
checkSpreadInteractions <- function(spreadInteractions, spreadCovariates) {
  if (length(spreadInteractions) && (!is.character(spreadInteractions) || !all(spreadInteractions %in% spreadInteractionClimateFuel)))
    stop("fireSense_dataPrepFit: parameter 'spreadInteractions' must be character(0) or \"", spreadInteractionClimateFuel, "\"; got: ",
         paste(format(spreadInteractions), collapse = ", "))
  if (length(spreadInteractions) && spreadCovariates %in% c("none", "fuel", "climate"))
    stop("fireSense_dataPrepFit: spreadInteractions = \"", spreadInteractionClimateFuel, "\" needs covariates (both climate and fuel): ",
         "spreadCovariates = \"", spreadCovariates, "\" lacks ",
         switch(spreadCovariates, none = "both (the null model)", fuel = "climate", climate = "fuel"))
  invisible(NULL)
}
