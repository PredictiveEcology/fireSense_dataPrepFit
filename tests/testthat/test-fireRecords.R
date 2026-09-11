## Fire records used to come from fireSenseUtils::getFirePoints_NFDB_V2(), whose NFDB URL
## (current_version/NFDB_point.zip) CFS renamed, so it returned HTTP 404 and no newer release could be
## fetched, and from fireSenseUtils::getFirePolygons(). They now come through
## reproducible::preProcess() and fireregimetools' loaders.

## A zipped shapefile beside a non-spatial file, as in the CFS archives.
fireArchive <- function(v, name, dir) {
  src <- file.path(dir, "src")
  dir.create(src, showWarnings = FALSE)
  terra::writeVector(v, file.path(src, paste0(name, ".shp")), overwrite = TRUE)
  writeLines("[fires]", file.path(src, "schema.ini"))
  zipFile <- file.path(dir, paste0(name, ".zip"))
  withr::with_dir(src, utils::zip(zipFile, dir(src), flags = "-q"))
  zipFile
}

crs <- "EPSG:3978"
studyArea <- terra::vect(terra::ext(0, 10000, 0, 10000), crs = crs)
square <- function(x0, y0, side) terra::as.polygons(terra::ext(x0, x0 + side, y0, y0 + side), crs = crs)

test_that("firePolysByYear returns the NBAC perimeters in the study area, split by year", {
  polys <- rbind(square(1000, 1000, 2000),   # 2001
                 square(9000, 5000, 2000),   # 2001, half outside the study area
                 square(3000, 3000, 50),     # 2003, 0.25 ha
                 square(20000, 20000, 1000), # 2003, outside the study area
                 square(5000, 5000, 1000))   # 1999, before the years asked for
  polys$YEAR <- c(2001, 2001, 2003, 2003, 1999)
  polys$NFIREID <- 1:5
  polys$POLY_HA <- polys$ADJ_HA <- c(400, 400, 0.25, 100, 100)
  dir <- withr::local_tempdir()

  out <- firePolysByYear(url = NULL, archive = fireArchive(polys, "NBAC_test", dir), years = 2001:2003,
                         studyArea = studyArea, destinationPath = dir)

  expect_identical(names(out), paste0(fireSenseUtils::yearTxt, 2001:2003))
  expect_null(out[["year2002"]])
  y2001 <- out[["year2001"]]
  expect_equal(sort(y2001$NFIREID), c(1, 2))
  ## POLY_HA is the area inside the study area: the straddling 400 ha fire keeps 200 ha
  ## (terra::expanse() measures on the ellipsoid, hence the tolerance)
  expect_equal(y2001$POLY_HA[order(y2001$NFIREID)], c(400, 200), tolerance = 1e-3)
  ## no size floor: the 0.25 ha fire is kept
  expect_equal(out[["year2003"]]$NFIREID, 3)
  expect_true(terra::same.crs(y2001, studyArea))
})

test_that("nfdbFirePoints keeps fires of every size and cause in the study area and years", {
  pts <- terra::vect(cbind(c(1000, 2000, 3000, 4000, 20000), c(1000, 2000, 3000, 4000, 20000)), crs = crs)
  pts$FIRE_ID <- c("a", "b", "c", "d", "e")
  pts$YEAR <- c(1995, 2002, 2003, 2010, 2002)
  pts$SIZE_HA <- c(0, 0.1, 50, 1, 5)
  pts$CAUSE <- c("L", "N", "H", "L", "L")
  dir <- withr::local_tempdir()

  out <- nfdbFirePoints(url = NULL, archive = fireArchive(pts, "NFDB_point_test", dir), years = 1995:2005,
                        studyArea = studyArea, destinationPath = dir)

  ## "d" is from 2010, "e" is outside the study area
  expect_identical(sort(out$FIRE_ID), c("a", "b", "c"))
  expect_equal(out$SIZE_HA[order(out$FIRE_ID)], c(0, 0.1, 50))
  expect_identical(out$CAUSE[order(out$FIRE_ID)], c("L", "N", "H"))
})
