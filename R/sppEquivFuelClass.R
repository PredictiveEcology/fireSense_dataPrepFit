#' Add the estimated fuel class to `sppEquiv`, keeping every species
#'
#' `assignedFuelClass` has only the species that have cohorts in the study area. `sppEquiv` is shared
#' with `sppNameVector` and `sppColorVect`, which keep all species, so no row may be dropped here:
#' a species without cohorts gets `NA` for its fuel class.
#'
#' @param sppEquiv data.table; replaced by a copy with `fuelClassCol` overwritten (added last).
#' @param assignedFuelClass data.table with the species (`sppEquivCol`) and its fuel class (`fuelClassCol`).
#' @param sppEquivCol,fuelClassCol column names, as the module parameters.
#' @return `sppEquiv` with the same rows, in the same order.
#' @noRd
setSppEquivFuelClass <- function(sppEquiv, assignedFuelClass, sppEquivCol, fuelClassCol) {
  out <- copy(sppEquiv)
  set(out, NULL, fuelClassCol, NULL)
  out[assignedFuelClass, (fuelClassCol) := get(paste0("i.", fuelClassCol)), on = sppEquivCol]
  out
}
