#' Stop unless non-forest can be `youngAge`
#'
#' `youngAge` is resolved for each fire year from time since disturbance, which covers forest and
#' non-forest pixels alike, so the mode where non-forest is never young (whose exclusivity needs
#' the non-forest land-cover columns in the annual table) is not available.
#'
#' @param nonForestCanBeYoungAge `P(sim)$nonForestCanBeYoungAge`.
#' @return `NULL`, invisibly.
stopIfNonForestCannotBeYoung <- function(nonForestCanBeYoungAge) {
  if (!isTRUE(nonForestCanBeYoungAge)) {
    stop("nonForestCanBeYoungAge = FALSE is not supported: youngAge is resolved for each fire year ",
         "from time since disturbance, which treats forest and non-forest pixels alike. ",
         "Set nonForestCanBeYoungAge = TRUE.", call. = FALSE)
  }
  invisible(NULL)
}

#' Every fire, as the pixels burned in each year
#'
#' All fires, not only those being fitted, because any of them resets `youngAge`. From
#' `historicalFireRaster` when supplied (as `fireSenseUtils::makeTSD()` does), otherwise from
#' `firePolysForAge`.
#'
#' @param sim a `simList` with `firePolysForAge` or `historicalFireRaster`, and `flammableRTMs`.
#' @return list named by year of integer `pixelID`s, see `fireSenseUtils::firePixelsByYear()`.
allFirePixelsByYear <- function(sim) {
  template <- sim$flammableRTMs[[1]]
  if (!is.null(sim$historicalFireRaster)) {
    terra::compareGeom(sim$historicalFireRaster, template, stopOnError = TRUE)
    return(fireSenseUtils::firePixelsByYear(fireRaster = sim$historicalFireRaster))
  }
  fireSenseUtils::firePixelsByYear(firePolys = sim$firePolysForAge, template = template)
}

#' Add `youngAge` to each fire year's annual table
#'
#' @param annualCovariates list named by data year of lists named `year<year>` of data.tables with
#'   `pixelID`, as `prepare_SpreadFit()` builds them. The tables are changed by reference.
#' @param tsds list named `year<dataYear>` of `SpatRaster`s, time since disturbance at each data year.
#' @param firePixelsByYear see `allFirePixelsByYear()`.
#' @param cutoffForYoungAge numeric.
#' @return `annualCovariates`, each table with a `youngAge` column (0/1), from
#'   `fireSenseUtils::youngAgeAtYear()`.
addAnnualYoungAge <- function(annualCovariates, tsds, firePixelsByYear, cutoffForYoungAge) {
  for (dy in names(annualCovariates)) {
    dataYear <- as.integer(gsub("[^0-9]", "", dy))
    ras <- tsds[[paste0(fireSenseUtils::yearTxt, dataYear)]]
    tsd <- data.table::data.table(pixelID = seq_len(terra::ncell(ras)),
                                  tsd = terra::values(ras, mat = FALSE))
    for (yrChar in names(annualCovariates[[dy]])) {
      x <- annualCovariates[[dy]][[yrChar]]
      data.table::set(x, NULL, fireSenseUtils::youngAgeTxt,
                      fireSenseUtils::youngAgeAtYear(tsd, dataYear = dataYear,
                                                     year = as.integer(gsub("[^0-9]", "", yrChar)),
                                                     firePixelsByYear = firePixelsByYear,
                                                     cutoffForYoungAge = cutoffForYoungAge,
                                                     pixelID = x$pixelID))
    }
  }
  annualCovariates
}
