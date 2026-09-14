## prepare_EscapeFit() took the escape points' coordinates with `terra::geom(escapes)[, c("x", "y")]`.
## With exactly one escape the one-row matrix became a plain vector and terra::cellFromXY() stopped:
## "unable to find an inherited method for function 'cellFromXY' for signature
## 'object = "SpatRaster", xy = "numeric"'" (ELF 12.1).
test_that("pointCoords returns a two-column matrix, also for a single point", {
  ras <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10000, ymin = 0, ymax = 10000,
                     crs = "EPSG:3978")
  one <- terra::vect(cbind(4500, 5500), crs = "EPSG:3978")
  two <- terra::vect(cbind(c(4500, 1500), c(5500, 2500)), crs = "EPSG:3978")

  cases <- list(list(one, 1L), list(two, 2L), list(sf::st_as_sf(one), 1L), list(sf::st_as_sf(two), 2L))
  for (case in cases) {
    xy <- pointCoords(case[[1]])
    expect_true(is.matrix(xy))
    expect_identical(ncol(xy), 2L)
    expect_identical(nrow(xy), case[[2]])
  }
  expect_identical(terra::cellFromXY(ras, pointCoords(one)), 45)
  expect_identical(terra::cellFromXY(ras, pointCoords(sf::st_as_sf(one))), 45)
})
