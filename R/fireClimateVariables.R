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

#' Pick the spread climate variable that best separates the bad fire years
#'
#' For `climateVariablesForFire$spread = "auto"`. What matters most is doing well in bad fire years,
#' so each candidate is scored by how well its study-area yearly mean separates the worst years
#' from the rest: the bad years are the top `badYearFraction` of fire years by area burned (and must
#' have burned something), and the score is the AUC, the probability that a bad year is drier on that
#' variable than another year (1 = every bad year drier, 0.5 = chance). The highest AUC wins; the
#' Spearman correlation with area burned breaks ties. If no candidate reaches `minAUC`, climate does not
#' separate this study area's bad years, and the first candidate is used, with a warning.
#'
#' A year-level score need not stay exactly right once the fitted model predicts, but a variable that
#' ranks first here has no reason to rank worst there.
#'
#' @param climateRasters Named list of `SpatRaster`s, one per climate variable, layers `year<YYYY>`.
#' @param burnedByYear Named numeric: pixels burned per fire year, names `year<YYYY>`; 0 for fire
#'   years without fires.
#' @param candidates Names in `climateRasters` to choose from; each must increase with dryness.
#' @param badYearFraction Share of fire years counted as bad (default 0.25).
#' @param minAUC Below this best AUC, fall back to `candidates[1]` (default 0.6).
#' @return `data.table` with `var`, `auc`, `rho`, `nYears`, `nBad`, and `chosen` (one `TRUE`).
#' @keywords internal
selectSpreadClimateVariable <- function(climateRasters, burnedByYear, candidates,
                                        badYearFraction = 0.25, minAUC = 0.6) {
  yrs <- names(burnedByYear)
  scores <- data.table::rbindlist(lapply(candidates, function(v) {
    r <- climateRasters[[v]]
    lyr <- if (is.null(r)) character(0) else intersect(yrs, names(r))
    out <- data.table::data.table(var = v, auc = NA_real_, rho = NA_real_, nYears = length(lyr), nBad = 0L)
    if (length(lyr) < 4) return(out)
    x <- unlist(terra::global(r[[lyr]], "mean", na.rm = TRUE), use.names = FALSE)
    a <- unname(burnedByYear[lyr])
    k <- max(1L, round(badYearFraction * length(a)))
    bad <- a >= sort(a, decreasing = TRUE)[k] & a > 0
    out$nBad <- sum(bad)
    if (sum(bad) >= 2 && sum(!bad) >= 2) {
      d <- outer(x[bad], x[!bad], "-")
      out$auc <- (sum(d > 0) + 0.5 * sum(d == 0)) / length(d)
    }
    out$rho <- suppressWarnings(stats::cor(x, a, method = "spearman"))
    out
  }))
  ok <- is.finite(scores$auc) & scores$auc >= minAUC
  best <- if (any(ok)) {
    s <- scores[ok][order(-auc, -rho)]
    s$var[1]
  } else {
    warning("No climate variable separates the bad fire years (best AUC ",
            round(max(scores$auc, na.rm = TRUE), 2), " < ", minAUC, "); using ", candidates[1])
    candidates[1]
  }
  scores[, chosen := var == best]
  scores[]
}

#' Add the climate variables of existing fits to what is being prepared
#'
#' For prediction from existing fits, possibly of several ELFs that chose different spread variables
#' (`spread = "auto"`): every variable a fit uses must be prepared. Variables already in
#' `climateVariables` keep their definition (and years); only missing ones are added.
#'
#' @param climateVariables The list for `sim$climateVariables` (names like `historical_CMDsm`).
#' @param fitClimVars Climate variables the fits use, ClimateNA names.
#' @inheritParams fireClimateLayers
#' @return `climateVariables`, with any missing variable added.
#' @keywords internal
addFitClimateVariables <- function(climateVariables, fitClimVars, historicalYears, projected = TRUE,
                                   projectedYears = 2011:2100) {
  prepared <- unique(sub("^[^_]+_", "", names(climateVariables)))
  toAdd <- fitClimVars[!gsub("_", "", fitClimVars) %in% prepared]
  if (!length(toAdd)) return(climateVariables)
  c(climateVariables, fireClimateLayers(list(ignition = toAdd), historicalYears = historicalYears,
                                        projected = projected, projectedYears = projectedYears))
}
