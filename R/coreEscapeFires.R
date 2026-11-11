## The fires an ELF is fitted to are those that start in its core (Eliot, 2026-11-10): the buffer only removes
## edge effects. A fire's start is its polygon centroid moved to the nearest burnable pixel of the fire
## (`harmonizeFireData()`), never the NFDB point; a fire that starts in the buffer is not fitted at all, and is
## not an escape or an ignition of this ELF. The sets themselves come from `fireSenseUtils::escapeFireSet()`.

#' The escape definition of this run
#'
#' @param Par the module's parameters.
#' @return `fireSenseUtils::escapeDefinition()` of `escapeSizeHa` and `fireCauses`.
escapeDefinitionOf <- function(Par) {
  fireSenseUtils::escapeDefinition(Par$escapeSizeHa, Par$fireCauses)
}

#' The fire-size column of fire polygons or points
#'
#' @param fires `sf` or `SpatVector` of fires.
#' @return the first column name ending in `_HA` (it has been `SIZE_HA` and `POLY_HA`).
fireSizeColumn <- function(fires) {
  grep("_HA$", names(fires), value = TRUE)[1]
}

#' The name of the flammable raster of a data-year group
#'
#' @param flammableRTMs list of `SpatRaster`s named for their data year (`year<dataYear>`).
#' @param dataYear the data year, as in `names(mod$allYears)`.
#' @return the name in `flammableRTMs` that contains `dataYear`.
flammableRTMName <- function(flammableRTMs, dataYear) {
  grep(dataYear, names(flammableRTMs), value = TRUE)
}

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
#' @param definition from [escapeDefinitionOf()].
#' @param idCol the fire ID column of `polys` and `points`.
#' @param sizeCol the fire-size column of `polys`.
#' @return a `SpatVector` of start points, with `core` (see `fireSenseUtils::fireStarts()`), or `NULL` with no fires.
escapeStartsByDataYear <- function(polys, points, flammableRTMs, allYears, coreRaster, definition, idCol, sizeCol) {
  sets <- Map(yrs = allYears, dy = names(allYears), function(yrs, dy) {
    fireSenseUtils::escapeFireSet(polys[names(polys) %in% yrs], flammableRTMs[[flammableRTMName(flammableRTMs, dy)]],
                                  coreRaster = coreRaster, definition = definition,
                                  starts = if (!is.null(points)) points[names(points) %in% yrs], idCol = idCol, sizeCol = sizeCol)
  })
  sets <- Filter(NROW, sets)
  if (!length(sets)) return(NULL)
  do.call(rbind, unname(sets))
}

#' The pixels of the core
#'
#' @param rasterToMatch `SpatRaster` whose pixels are returned.
#' @param coreRaster `rasterToMatchLargeELF`: its cells equal to `coreValue` are the core. `NULL`: every pixel.
#' @param coreValue the value of the core, `fireSenseUtils::fireSenseELFCoreValue`.
#' @return integer `pixelID`s (cells of `rasterToMatch`) whose centre is in the core.
corePixelIDs <- function(rasterToMatch, coreRaster, coreValue = fireSenseUtils::fireSenseELFCoreValue) {
  cells <- seq_len(terra::ncell(rasterToMatch))
  if (is.null(coreRaster)) return(cells)
  v <- terra::extract(coreRaster, terra::xyFromCell(rasterToMatch, cells))[[1]]
  cells[!is.na(v) & v %in% coreValue]
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
