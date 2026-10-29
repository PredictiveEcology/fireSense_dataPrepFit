## fireCauses: NULL is the original behaviour in all three models (ignition points L/N, polygon ignitions
## "Natural", spread fires of every cause); a value switches all three to those cause groups together.

sqWkt <- function(x0, y0, w) sprintf("POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))",
                                     x0, y0, x0 + w, y0, x0 + w, y0 + w, x0, y0 + w, x0, y0)

test_that("pointsOfCauses: natural keeps L and N; adding unknown keeps U, blank and missing", {
  pts <- terra::vect(cbind(1:7, 1), crs = "EPSG:3978")
  pts$CAUSE <- c("L", "N", "H", "U", "", NA, "H-PB")
  expect_identical(pointsOfCauses(pts, "natural")$CAUSE, c("L", "N"))
  expect_equal(NROW(pointsOfCauses(pts, c("natural", "unknown"))), 5)
  expect_identical(pointsOfCauses(pts, c("natural", "unknown"))$CAUSE[1:4], c("L", "N", "U", ""))
})

test_that("polysOfCauses keeps the polygons of the given causes, and none without a value", {
  polys <- terra::vect(sqWkt(c(0, 200, 400, 600), 0, 100), crs = "EPSG:3978")
  polys$FIRECAUS <- c("Natural", "Human", "Undetermined", NA)
  expect_identical(polysOfCauses(polys, "Natural")$FIRECAUS, "Natural")
  expect_identical(polysOfCauses(polys, fireSenseUtils::fireCauseCodes(c("natural", "unknown"), "nbac"))$FIRECAUS,
                   c("Natural", "Undetermined", NA))
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
  polys$FIRECAUS <- c("Natural", "Human", "Undetermined")
  pts <- terra::centroids(polys)

  run <- function(causes) {
    sim <- new.env()
    sim$spreadFirePolys <- list(year2000 = polys)
    sim$spreadFirePoints <- list(year2000 = pts)
    sim$flammableRTMs <- list(year1985 = terra::rast(nrows = 100, ncols = 100, xmin = 0, xmax = 10000, ymin = 0,
                                                      ymax = 10000, crs = "EPSG:3978"))
    env <- new.env(parent = environment(prepare_SpreadFitFire_Vector))
    env$Par <- list(escapeSizeHa = 50, fireCauses = causes)
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

  nat <- run("natural")
  expect_identical(nat$spreadFirePolys$year2000$NFIREID, 1)
  expect_identical(nat$spreadFirePoints$year2000$NFIREID, 1)

  natUnk <- run(c("natural", "unknown"))
  expect_identical(natUnk$spreadFirePolys$year2000$NFIREID, c(1, 3))
  expect_identical(natUnk$spreadFirePoints$year2000$NFIREID, c(1, 3))

  expect_error(run("other"), "fireCauses")
})

test_that("prepare_IgnitionFit: an unknown-cause polygon with no point adds an ignition only when fireCauses has unknown", {
  skip_if_not_installed("fireSenseUtils")
  polys <- terra::vect(sqWkt(c(100, 3100, 6100), 100, 800), crs = "EPSG:3978")
  polys$YEAR <- 2000L
  polys$FIRECAUS <- c("Natural", "Undetermined", "Human")
  polys$NFIREID <- c(1, 2, 3)
  starts <- terra::vect(cbind(c(500, 3500, 6500), 500), crs = "EPSG:3978")
  starts$YEAR <- 2000L
  starts$NFIREID <- c(1, 2, 3)
  pts <- terra::vect(cbind(300, 300), crs = "EPSG:3978") # in polygon 1
  pts$YEAR <- 2000L; pts$CAUSE <- "L"; pts$SIZE_HA <- 10

  run <- function(causes) {
    sim <- new.env()
    sim$ignitionFirePoints <- pts
    sim$spreadFirePolys <- list(year2000 = polys)
    sim$spreadFirePoints <- list(year2000 = starts)
    sim$flammableRTMs <- list(year1985 = terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 10000, ymin = 0,
                                                      ymax = 10000, crs = "EPSG:3978"))
    env <- new.env(parent = environment(prepare_IgnitionFit))
    env$sim <- sim
    env$Par <- list(escapeSizeHa = 50, fireCauses = causes)
    env$mod <- new.env()
    env$mod$studyAreaUnion <- terra::as.polygons(terra::ext(0, 10000, 0, 10000), crs = "EPSG:3978")
    eval(as.list(body(prepare_IgnitionFit))[-1][[1]], env)
    sim$ignitionFirePoints
  }
  expect_equal(NROW(run(NULL)), 1)                          # as before: only "Natural" polygons add an ignition
  expect_equal(NROW(run("natural")), 1)
  expect_equal(NROW(run(c("natural", "unknown"))), 2)       # polygon 2 adds one; the human polygon never does
  expect_equal(unname(terra::crds(run(c("natural", "unknown")))[2, ]), c(3500, 500))
})
