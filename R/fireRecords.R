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
  ## reproducible:: by name: caret (a fireSense_IgnitionFit reqdPkg) also has a preProcess(), which
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
