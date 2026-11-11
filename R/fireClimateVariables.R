## The spread climate variable chosen for a fit (`spread = "auto"`). The fire models' default climate
## variables and the climate layers they need are in fireSenseUtils (fireClimateVariables.R).

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
#' @param candidates Names in `climateRasters` to choose from; each must increase with dryness, except
#'   precipitation (names starting `PPT`), which is scored as decreasing with dryness.
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
    if (grepl("^PPT", v)) x <- -x            # precipitation: drier means less of it
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


#' The climate variables of `climateVariablesForFire`, with the spread candidates when `spread = "auto"`
#'
#' What `climateVariables` must prepare: the ignition and spread variables, and, when the spread
#' variable is chosen ("auto"), every candidate it chooses among. `ignition` is not changed, so
#' IgnitionFit still uses exactly its own list.
#'
#' @param climateVariablesForFire List with `ignition` and `spread`.
#' @param spreadClimateCandidates Character: the module's `spreadClimateCandidates`.
#' @return `climateVariablesForFire`, with `spread` replaced by the candidates if it was `"auto"`.
#' @keywords internal
withSpreadCandidates <- function(climateVariablesForFire, spreadClimateCandidates) {
  if ("auto" %in% climateVariablesForFire[["spread"]])
    climateVariablesForFire[["spread"]] <- union(setdiff(climateVariablesForFire[["spread"]], "auto"),
                                                 spreadClimateCandidates)
  climateVariablesForFire
}

#' Score the spread climate candidates against the fire years
#'
#' `spread = "auto"`: area burned per fire year is the number of buffer-1 pixels in
#' `fireBufferedListDT`; the candidates, named without underscores as `historicalClimateRasters` is, are
#' scored by [selectSpreadClimateVariable()].
#'
#' @param climateRasters `sim$historicalClimateRasters`.
#' @param fireBufferedListDT `sim$fireBufferedListDT`: per fire year a table with `buffer`.
#' @param fireYears Fire years.
#' @param candidates `spreadClimateCandidates`, with or without underscores.
#' @return As [selectSpreadClimateVariable()].
#' @keywords internal
autoSpreadClimate <- function(climateRasters, fireBufferedListDT, fireYears, candidates) {
  fy <- paste0(fireSenseUtils::yearTxt, fireYears)
  burned <- setNames(numeric(length(fy)), fy)
  b <- vapply(fireBufferedListDT, function(d) as.numeric(sum(d$buffer == 1)), numeric(1))
  b <- b[names(b) %in% fy]
  burned[names(b)] <- b
  selectSpreadClimateVariable(climateRasters, burned, candidates = gsub("_", "", candidates))
}
