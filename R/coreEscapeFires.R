## The fires an ELF is fitted to are those that start in its core (Eliot, 2026-11-10): the buffer only removes
## edge effects. A fire's start is its polygon centroid moved to the nearest burnable pixel of the fire
## (`harmonizeFireData()`), never the NFDB point; a fire that starts in the buffer is not fitted at all, and is
## not an escape or an ignition of this ELF. The sets themselves come from `fireSenseUtils::escapeFireSet()`.

#' Every escaped fire's start, and whether it is in the core
#'
#' `fireSenseUtils::escapeFireSet()` for each data-year group, with the starts the spread fit was given
#' (`harmonizeFireData()`), so a fire is fitted from the start it is counted at.
#'
#' @param polys,points lists named `year<year>` of the harmonized fire polygons and their start points. With
#'   `points = NULL` the starts are the polygons' centroids moved to the nearest burnable pixel (`idCol` is not needed).
#' @param flammableRTMs list of `SpatRaster`s named `year<dataYear>`.
#' @param allYears list named by data year of the fire year names, as `mod$allYears`.
#' @param coreRaster `rasterToMatchLargeELF`, or `NULL` when the run has no ELF (every fire is core).
#' @param definition from `fireSenseUtils::escapeDefinitionOf()`.
#' @param idCol the fire ID column of `polys` and `points`.
#' @return a `SpatVector` of start points, with `core` (see `fireSenseUtils::fireStarts()`), or `NULL` with no fires.
escapeStartsByDataYear <- function(polys, points, flammableRTMs, allYears, coreRaster, definition, idCol) {
  sets <- Map(yrs = allYears, dy = names(allYears), function(yrs, dy) {
    flammableRTM <- flammableRTMs[[fireSenseUtils::dataYearElement(flammableRTMs, as.integer(dy))]]
    fireSenseUtils::escapeFireSet(polys[names(polys) %in% yrs], flammableRTM, coreRaster = coreRaster, definition = definition,
                                  starts = if (!is.null(points)) points[names(points) %in% yrs], idCol = idCol)
  })
  sets <- Filter(NROW, sets)
  if (!length(sets)) return(NULL)
  do.call(rbind, unname(sets))
}

#' The IDs of the fires that start in the core, by fire year
#'
#' @param fireSet from [escapeStartsByDataYear()].
#' @param idCol the fire ID column of `fireSet`.
#' @return list named `year<year>` of fire IDs; empty if `fireSet` is `NULL`.
coreFireIDs <- function(fireSet, idCol) {
  if (is.null(fireSet)) return(list())
  v <- terra::values(fireSet)
  v <- v[v$core, ]
  split(v[[idCol]], paste0(fireSenseUtils::yearTxt, v$YEAR))
}
