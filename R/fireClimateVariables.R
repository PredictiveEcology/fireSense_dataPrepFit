## Which climate variables the fire models use, and the climate layers that requires.
##
## `climateVariablesForFire` names the variables for each process: IgnitionFit (xgboost) takes all of
## `ignition`, as it copes with correlated covariates; SpreadFit takes `spread`. Unless a user supplies
## them, this module sets both it and `climateVariables` (what canClimateData prepares) from the
## default below, so the two always agree and a user need not choose. canClimateData leaves
## `climateVariables` to this module whenever both are loaded (this module lists it as an output), and
## keeps its own default otherwise.

#' Default climate variables for the fire models
#'
#' ClimateNA names; `"cumMDC"` is derived by `climateData::calcCumMDC()`. `spread = "auto"` picks,
#' per study area, the `ignition` variable that best tracks annual area burned
#' ([selectSpreadClimateVariable()]).
#' @keywords internal
defaultClimateVariablesForFire <- list(ignition = c("CMD", "cumMDC", "CMD_sm", "CMD_sp"),
                                       spread = "auto")

#' ClimateNA names of climate variables given with or without underscores
#'
#' The fire modules name climate layers without underscores (`"CMDsm"`), ClimateNA with them
#' (`"CMD_sm"`). Both forms map to the ClimateNA name.
#'
#' @param x Character vector of variable names.
#' @return Character vector of ClimateNA names, same length as `x`.
#' @keywords internal
climateNAnames <- function(x) {
  allowed <- getFromNamespace(".allowedClimateVars", ns = asNamespace("climateData"))
  out <- allowed[match(gsub("_", "", x), gsub("_", "", allowed))]
  if (anyNA(out))
    stop("Unknown climate variable(s): ", paste(x[is.na(out)], collapse = ", "),
         ". Use ClimateNA names (e.g., CMD_sm) or climateData's derived variables (e.g., cumMDC).")
  out
}

#' Climate layers needed by the fire models
#'
#' @param climateVariablesForFire List with `ignition` and `spread` character vectors.
#' @param historicalYears Years for fitting (`P(sim)$fireYears`).
#' @param projected Logical; also request projected climate.
#' @param projectedYears Projected years.
#' @return A list for `sim$climateVariables`, as built by `climateData::climateLayers()`.
#' @keywords internal
fireClimateLayers <- function(climateVariablesForFire, historicalYears, projected = TRUE,
                              projectedYears = 2011:2100) {
  vars <- unlist(climateVariablesForFire[c("ignition", "spread")], use.names = FALSE)
  vars <- unique(climateNAnames(setdiff(vars, "auto")))    # "auto" chooses among the ignition variables
  climateData::climateLayers(vars, fun = quote(calcAsIs), historicalYears = historicalYears,
                             projected = projected, projectedYears = projectedYears)
}

#' Pick the spread climate variable that best tracks annual area burned
#'
#' For `climateVariablesForFire$spread = "auto"`. Each candidate's study-area mean per year is ranked
#' against the area burned that year (Spearman); the highest positive correlation wins, so a candidate
#' must increase with dryness. A year-level correlation need not stay exactly right once the fitted
#' model predicts, but a variable that ranks first here has no reason to rank worst there.
#'
#' @param climateRasters Named list of `SpatRaster`s, one per climate variable, layers `year<YYYY>`.
#' @param burnedByYear Named numeric: pixels burned per fire year, names `year<YYYY>`; 0 for fire
#'   years without fires (a dry year with no fire is informative too).
#' @param candidates Names in `climateRasters` to choose from.
#' @return `data.table` with `var`, `rho`, `nYears`, and `chosen` (one `TRUE`).
#' @keywords internal
selectSpreadClimateVariable <- function(climateRasters, burnedByYear, candidates) {
  yrs <- names(burnedByYear)
  scores <- data.table::rbindlist(lapply(candidates, function(v) {
    r <- climateRasters[[v]]
    lyr <- if (is.null(r)) character(0) else intersect(yrs, names(r))
    rho <- if (length(lyr) < 3) NA_real_ else {
      m <- unlist(terra::global(r[[lyr]], "mean", na.rm = TRUE), use.names = FALSE)
      suppressWarnings(stats::cor(m, unname(burnedByYear[lyr]), method = "spearman"))
    }
    data.table::data.table(var = v, rho = rho, nYears = length(lyr))
  }))
  ok <- is.finite(scores$rho) & scores$rho > 0
  best <- if (any(ok)) scores$var[ok][which.max(scores$rho[ok])] else {
    warning("No climate variable correlates positively with annual area burned; using ", candidates[1])
    candidates[1]
  }
  scores[, chosen := var == best]
  scores[]
}
