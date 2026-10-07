## dataPrepBuild() only clipped sim$ignitionFirePoints inside an `if (!same.crs(points, rasterToMatch))`
## branch, and clipped to the rasterToMatch rectangle. When the points already arrived in the RTM's
## CRS -- the normal case -- nothing was clipped at all, so a point outside the studyArea polygon
## survived to the `all ignitionFirePoints are not within studyArea` stopifnot() in
## prepare_IgnitionFit(). ELF 5.1.2: 1 of 3307 NFDB points sat 618 m outside the polygon but inside
## the RTM extent, and failed the run three times.
test_that("clipPointsToStudyArea drops points outside the polygon, same CRS or not", {
  crsSA <- "EPSG:3978"
  sa <- terra::vect("POLYGON ((0 0, 10000 0, 10000 10000, 0 10000, 0 0))", crs = crsSA)

  ## inside, and outside the polygon but inside the raster-aligned extent around it
  pts <- terra::vect(cbind(c(5000, 2000, 12000), c(5000, 8000, 5000)), crs = crsSA)
  pts$id <- c("in1", "in2", "out")

  clipped <- clipPointsToStudyArea(pts, sa)
  expect_equal(nrow(clipped), 2)
  expect_identical(sort(clipped$id), c("in1", "in2"))

  ## and when the points need reprojecting first
  ptsLL <- terra::project(pts, "EPSG:4326")
  clippedLL <- clipPointsToStudyArea(ptsLL, sa)
  expect_equal(nrow(clippedLL), 2)
  expect_true(terra::same.crs(clippedLL, sa))

  ## an sf study area is accepted too (dataPrepBuild() st_union()s an sf studyArea)
  clippedSF <- clipPointsToStudyArea(pts, sf::st_as_sf(sa))
  expect_equal(nrow(clippedSF), 2)
})
