## spreadFitCauses: NULL keeps every fire (unchanged); "Natural" fits the natural-cause polygons only, and the
## start points follow their polygons.

sqWkt <- function(x0, y0, w) sprintf("POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))",
                                     x0, y0, x0 + w, y0, x0 + w, y0 + w, x0, y0 + w, x0, y0)

test_that("polysOfCauses keeps the polygons of the given causes, and none without a value", {
  polys <- terra::vect(sqWkt(c(0, 200, 400, 600), 0, 100), crs = "EPSG:3978")
  polys$FIRECAUS <- c("Natural", "Human", "Unknown", NA)
  expect_identical(polysOfCauses(polys, "Natural")$FIRECAUS, "Natural")
  expect_identical(polysOfCauses(polys, c("Natural", "Human"))$FIRECAUS, c("Natural", "Human"))
  expect_null(polysOfCauses(polys, "Lightning"))
  expect_null(polysOfCauses(NULL, "Natural"))
  polys$FIRECAUS <- NULL
  expect_error(polysOfCauses(polys, "Natural"), "FIRECAUS")
})

test_that("prepare_SpreadFitFire_Vector: NULL spreadFitCauses leaves the fires as they are, 'Natural' drops the rest", {
  polys <- terra::vect(sqWkt(c(0, 2000, 4000), 0, 1000), crs = "EPSG:3978") # 100 ha each
  polys$YEAR <- 2000L
  polys$NFIREID <- c(1, 2, 3)
  polys$POLY_HA <- 100
  polys$FIRECAUS <- c("Natural", "Human", "Natural")
  pts <- terra::centroids(polys)

  run <- function(causes) {
    sim <- new.env()
    sim$spreadFirePolys <- list(year2000 = polys)
    sim$spreadFirePoints <- list(year2000 = pts)
    sim$flammableRTMs <- list(year1985 = terra::rast(nrows = 100, ncols = 100, xmin = 0, xmax = 10000, ymin = 0,
                                                      ymax = 10000, crs = "EPSG:3978"))
    env <- new.env(parent = environment(prepare_SpreadFitFire_Vector))
    env$Par <- list(escapeSizeHa = 50, spreadFitCauses = causes)
    env$P <- function(sim) list(areaMultiplier = 1, minBufferSize = 0)
    env$mod <- new.env()
    env$mod$allYears <- list(year1985 = "year2000")
    env$Cache <- function(x, ...) x
    ## the harmonization is tested elsewhere: pass the fires through
    env$harmonizeFireData <- function(firePolys, spreadFirePoints, ...)
      list(fireBufferedListDT = Map(function(x) 1, firePolys), spreadFirePoints = spreadFirePoints, firePolys = firePolys)
    f <- prepare_SpreadFitFire_Vector
    environment(f) <- env
    suppressMessages(f(sim))
    sim
  }

  all <- run(NULL)
  expect_identical(all$spreadFirePolys$year2000$NFIREID, c(1, 2, 3))
  expect_identical(all$spreadFirePoints$year2000$NFIREID, c(1, 2, 3))

  nat <- run("Natural")
  expect_identical(nat$spreadFirePolys$year2000$NFIREID, c(1, 3))
  expect_identical(nat$spreadFirePoints$year2000$NFIREID, c(1, 3))

  expect_error(run("Lightning"), "spreadFitCauses")
})
