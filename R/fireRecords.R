## Canadian fire records for fitting.
##
## The archives are fetched with reproducible::preProcess(), which records where each file came from,
## so reproducible::preProcessCheckURLs() can recheck them against the server when asked. They are read
## with fireregimetools, which harmonises the YEAR and SIZE_HA columns across releases; fires of every
## size are kept (`min_size_ha = 0`).

#' Path to the shapefile in a fire-record archive, downloading it if needed
#'
#' The file name carries the release (e.g. NFDB_point_20260811.shp), so callers key their Cache on
#' it: the NFDB URL stays the same from one release to the next.
#'
#' @param url URL of the zip archive.
#' @param destinationPath directory to download and extract into.
#' @param ... passed to `reproducible::preProcess()`, e.g. `archive`.
#' @return character, path of the `.shp` file.
fireRecordShapefile <- function(url, destinationPath, ...) {
  ## reproducible:: by name: caret (a fireSense_ignitionFit reqdPkg) also has a preProcess(), which
  ## masks this one when attached later. With `fun = NA`, preProcess() returns every file in the
  ## archive, and the CFS archives also carry metadata, hence the grep.
  files <- reproducible::preProcess(url = url, destinationPath = destinationPath, fun = NA, ...)$targetFilePath
  grep("\\.shp$", files, value = TRUE)
}

#' NBAC fire perimeters in a study area, by year
#'
#' @param shp path of the NBAC shapefile.
#' @param years integer years to keep.
#' @param studyArea polygon to crop to; the result is in its CRS.
#' @return list of `SpatVector` polygons, one element per year, named `year<year>`; `NULL` where a
#'   year has no fires. `POLY_HA` is the area (ha) inside `studyArea`.
firePolysByYear <- function(shp, years, studyArea) {
  polys <- fireregimetools::load_nbac_polys(shp, study_area = studyArea, fire_years = years, min_size_ha = 0)
  polys$POLY_HA <- round(terra::expanse(polys, unit = "ha"), 2)
  out <- lapply(years, function(yr) {
    p <- polys[polys$YEAR == yr, ]
    if (nrow(p) > 0) p else NULL
  })
  names(out) <- paste0(fireSenseUtils::yearTxt, years)
  out
}

#' NFDB fire points in a study area, of every size and cause
#'
#' @param shp path of the NFDB point shapefile.
#' @param years integer years to keep.
#' @param studyArea polygon to crop to; the result is in its CRS.
#' @return spatial points, with `YEAR` and `SIZE_HA` harmonized by fireregimetools.
nfdbFirePoints <- function(shp, years, studyArea) {
  fireregimetools::load_nfdb_points(shp, study_area = studyArea, fire_years = years, min_size_ha = 0)
}

#' Rows of `x` whose ID is among the IDs of `y`
#'
#' `x[[IDvar]]` on a `SpatVector` is a one-column `data.frame`, and `%in%` between two `data.frame`s compares
#' each whole column as one element: one `TRUE` or `FALSE` for the lot, `TRUE` only when the IDs are identical
#' and in the same order. `harmonizeFireData()` reorders points (moved centroids go last), so a year whose
#' points came back reordered lost every polygon (ELF 6.6.1: 64 of 127, all escapes of 9 years). This compares
#' the ID values.
#'
#' @param x,y `SpatVector`s with an `IDvar` column.
#' @param IDvar Name of the ID column.
#' @return `x` without the rows whose ID is not in `y`.
withIDsIn <- function(x, y, IDvar) {
  x[terra::values(x)[[IDvar]] %in% terra::values(y)[[IDvar]], ]
}

#' NFDB fire points with one row per fire
#'
#' Some NFDB rows repeat another row's date (`YEAR`, `MONTH`, `DAY`) and location exactly, nearly always under
#' a different `FIRE_ID`: 554 of the 134,007 natural fires of 1986-2025 in release 20260811, 319 of them in BC.
#' Each such group is taken as one fire recorded more than once, and the row with the largest `SIZE_HA` is kept.
#'
#' @param pts `SpatVector` of NFDB points with `YEAR`, `MONTH`, `DAY` and `SIZE_HA`.
#' @return `pts` without the repeated rows, in their original order.
dropRepeatedFirePoints <- function(pts) {
  xy <- terra::crds(pts)
  key <- paste(pts$YEAR, pts$MONTH, pts$DAY, xy[, 1], xy[, 2])
  ord <- order(pts$SIZE_HA, decreasing = TRUE)
  pts[sort(ord[!duplicated(key[ord])]), ]
}
