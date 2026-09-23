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
#' ClimateNA names; `"cumMDC"` is derived by `climateData::calcCumMDC()`.
#' @keywords internal
defaultClimateVariablesForFire <- list(ignition = c("CMD_sm", "cumMDC", "CMD_sp"),
                                       spread = "CMD_sm")

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
  vars <- unique(climateNAnames(unlist(climateVariablesForFire[c("ignition", "spread")], use.names = FALSE)))
  climateData::climateLayers(vars, fun = quote(calcAsIs), historicalYears = historicalYears,
                             projected = projected, projectedYears = projectedYears)
}
