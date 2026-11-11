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

#' Rows of `x` whose ID is among the IDs of `y`
#'
#' `x[[IDvar]]` on a `SpatVector` is a one-column `data.frame`, and `%in%` between two `data.frame`s compares
#' each whole column as one element: one `TRUE` or `FALSE` for the lot, `TRUE` only when the IDs are identical
#' and in the same order. `harmonizeFireData()` reorders points (moved centroids go last), so a year whose
#' points came back reordered lost every polygon (ELF 6.6.1: 64 of 127, all escapes of 9 years). This compares
#' the ID values.
#'
#' @param x `SpatVector` with an `IDvar` column.
#' @param y `SpatVector` with an `IDvar` column, or a vector of IDs.
#' @param IDvar Name of the ID column.
#' @return `x` without the rows whose ID is not in `y`.
withIDsIn <- function(x, y, IDvar) {
  if (inherits(y, "SpatVector")) y <- terra::values(y)[[IDvar]]
  x[terra::values(x)[[IDvar]] %in% y, ]
}

#' Fire polygons of the given causes
#'
#' The cause is the polygon's `FIRECAUS` (`fireSenseUtils::polygonIgnitions()` reads the same column). A polygon with
#' no value is not of any cause asked for.
#'
#' @param polys `SpatVector` of fire polygons.
#' @param causes `FIRECAUS` values to keep, e.g. from `fireSenseUtils::fireCauseCodes()`.
#' @param causeCol Name of the cause column.
#' @return `polys` without the other rows; `NULL` if none is left.
polysOfCauses <- function(polys, causes, causeCol = "FIRECAUS") {
  if (is.null(polys)) return(NULL)
  if (!causeCol %in% names(polys)) stop("fire polygons have no '", causeCol, "' column: cannot select by cause")
  polys <- polys[terra::values(polys)[[causeCol]] %in% causes, ]
  if (nrow(polys) > 0) polys else NULL
}
