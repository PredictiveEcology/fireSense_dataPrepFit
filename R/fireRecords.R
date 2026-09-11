## Canadian fire records for fitting.
##
## The archives are fetched with reproducible::preProcess(), which records where each file came from,
## so reproducible::preProcessCheckURLs() can recheck them against the server when asked. They are read
## with fireregimetools, which harmonises the YEAR and SIZE_HA columns across releases; fires of every
## size are kept (`min_size_ha = 0`).

## The shapefile in the fire-record archive at `url`, downloaded into `destinationPath` if needed.
## Its name carries the release (e.g. NFDB_point_20260811.shp), so callers key their Cache on it: the
## NFDB URL stays the same from one release to the next. With `fun = NA`, preProcess() returns every
## file in the archive, and the CFS archives also carry metadata (schema.ini, a pdf, a spreadsheet).
fireRecordShapefile <- function(url, destinationPath, ...) {
  files <- preProcess(url = url, destinationPath = destinationPath, fun = NA, ...)$targetFilePath
  grep("\\.shp$", files, value = TRUE)
}

## NBAC fire perimeters from shapefile `shp` in `studyArea` (in its CRS), as a list with one element
## per year in `years`, named "year<YYYY>" and NULL where that year has no fires. POLY_HA is the area
## inside `studyArea`.
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

## NFDB fire points of every size and cause from shapefile `shp` in `studyArea` (in its CRS) during
## `years`.
nfdbFirePoints <- function(shp, years, studyArea) {
  fireregimetools::load_nfdb_points(shp, study_area = studyArea, fire_years = years, min_size_ha = 0)
}
