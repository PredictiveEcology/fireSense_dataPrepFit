## Eliot, 2026-11-10: the buffer of an ELF only removes edge effects. A fire belongs to the ELF if its start (the
## polygon centroid moved to the nearest burnable pixel, as harmonizeFireData() gives it) is in the core
## (rasterToMatchLargeELF == 2); a fire that starts in the buffer is not fitted, and its NFDB point is not an
## ignition or an escape of the ELF. The sets come from fireSenseUtils::escapeFireSet() and escapeIgnitionPoints().

## 10 km square of 250 m (6.25 ha) pixels; core (2) in the west half, buffer (1) in the east half
coreFixture <- function() {
  crs <- "EPSG:3978"
  rtm <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 10000, ymin = 0, ymax = 10000, crs = crs, vals = 1)
  core <- rtm
  terra::values(core) <- ifelse(terra::xFromCell(core, seq_len(terra::ncell(core))) < 5000, 2, 1)
  sq <- function(x0, y0, w) sprintf("POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))",
                                    x0, y0, x0 + w, y0, x0 + w, y0 + w, x0, y0 + w, x0, y0)
  ## 100 ha each: fire 1 starts in the core, fire 2 in the buffer, fire 3 is small (20 ha), in the core
  polys <- terra::vect(c(sq(1000, 4000, 1000), sq(7000, 4000, 1000), sq(1000, 7000, 500)), crs = crs)
  polys$YEAR <- 2000L
  polys$NFIREID <- 1:3
  polys$FIRECAUS <- "Natural"
  polys$POLY_HA <- c(100, 100, 25)
  starts <- terra::centroids(polys)
  list(rtm = rtm, core = core, polys = polys, starts = starts, crs = crs)
}

test_that("a fire that starts in the buffer is not fitted: prepare_SpreadFitFire_Vector keeps core starts only", {
  skip_if_not_installed("fireSenseUtils")
  skip_if_not("escapeFireSet" %in% getNamespaceExports("fireSenseUtils"))
  fx <- coreFixture()
  sim <- new.env()
  sim$spreadFirePolys <- list(year2000 = fx$polys)
  sim$spreadFirePoints <- list(year2000 = fx$starts)
  sim$flammableRTMs <- list(year1985 = fx$rtm)
  sim$rasterToMatchLargeELF <- fx$core
  env <- new.env(parent = environment(prepare_SpreadFitFire_Vector))
  env$Par <- list(escapeSizeHa = 50, fireCauses = "natural")
  env$P <- function(sim) list(areaMultiplier = 10, minBufferSize = 100)
  env$Cache <- function(x, ...) x
  env$mod <- new.env()
  env$mod$allYears <- list(`1985` = "year2000")
  f <- prepare_SpreadFitFire_Vector
  environment(f) <- env
  suppressMessages(suppressWarnings(f(sim)))

  ## the two 100 ha fires are escapes; the 25 ha fire is not
  expect_setequal(sim$escapeFireSet$NFIREID, 1:2)
  expect_identical(sim$escapeFireSet$core[order(sim$escapeFireSet$NFIREID)], c(TRUE, FALSE))
  ## only the core one is fitted
  expect_identical(terra::values(sim$spreadFirePoints$year2000)$NFIREID, 1)
  expect_setequal(unique(sim$fireBufferedListDT$year2000$ids), 1)
  ## polygons keep the buffer-start escape, so its NFDB point can be left out of the ignition data
  expect_setequal(terra::values(sim$spreadFirePolys$year2000)$NFIREID, 1:2)
})

test_that("the core ignition training set leaves out the NFDB points of buffer-start fires and of the buffer", {
  skip_if_not_installed("fireSenseUtils")
  skip_if_not("escapeIgnitionPoints" %in% getNamespaceExports("fireSenseUtils"))
  fx <- coreFixture()
  polys <- fx$polys[1:2, ]
  def <- fireSenseUtils::escapeDefinition(50, "natural")
  fireSet <- escapeStartsByDataYear(list(year2000 = polys), list(year2000 = fx$starts[1:2, ]),
                                    list(year1985 = fx$rtm), list(`1985` = "year2000"), fx$core, def,
                                    "NFIREID", "POLY_HA")
  pts <- terra::vect(cbind(c(1100, 7100, 3000, 8500), c(4100, 4100, 8000, 8000)), crs = fx$crs)
  pts$YEAR <- 2000L
  pts$CAUSE <- "L"
  pts$SIZE_HA <- c(100, 100, 2, 2) # in the core fire, in the buffer fire, a small core fire, a small buffer fire

  sim <- new.env()
  sim$ignitionFirePoints <- pts
  sim$spreadFirePolys <- list(year2000 = polys)
  sim$escapeFireSet <- fireSet
  sim$rasterToMatchLargeELF <- fx$core
  sim$flammableRTMs <- list(year1985 = fx$rtm)
  env <- new.env(parent = environment(prepare_IgnitionFit))
  env$sim <- sim
  env$Par <- list(escapeSizeHa = 50, fireCauses = "natural")
  env$mod <- new.env()
  env$mod$studyAreaUnion <- terra::as.polygons(terra::ext(0, 10000, 0, 10000), crs = fx$crs)
  first <- as.list(body(prepare_IgnitionFit))[[2]]  # the block that makes the core ignition points
  expect_match(paste(deparse(first), collapse = " "), "escapeIgnitionPoints", fixed = TRUE)
  eval(first, env)

  got <- sim$ignitionFirePoints
  xy <- terra::crds(got)
  ## the core escape's point is replaced by its fire's start; the small core point stays, non-escaped
  expect_true(all(xy[, 1] < 5000))
  expect_equal(sum(got$escaped), 1)
  expect_equal(NROW(got), 2)
  expect_false(any(xy[, 1] > 5000)) # no buffer point, with or without a fire
})

test_that("corePixelIDs gives the pixels whose centre is in the core, and every pixel without a core raster", {
  fx <- coreFixture()
  ids <- corePixelIDs(fx$rtm, fx$core)
  expect_equal(length(ids), 40 * 20)
  expect_true(all(terra::xFromCell(fx$rtm, ids) < 5000))
  expect_identical(corePixelIDs(fx$rtm, NULL), seq_len(terra::ncell(fx$rtm)))
})

## A run can ask for summer maximum temperature as the spread climate (or an ignition covariate) without any
## code here: ClimateNA's seasonal Tmax_sm is the mean of the monthly Tmax06, Tmax07 and Tmax08, climateData
## serves it as it does CMD_sm (`.allowedClimateVars_S`), and the module reads layers as "Tmaxsm".
## It is a candidate only when asked for: not in the defaults, so not in the "auto" choice.
test_that("Tmax_sm can be requested as a climate variable and is not in the defaults", {
  skip_if_not_installed("fireSenseUtils")
  skip_if_not_installed("climateData")
  expect_false(any(grepl("Tmax", unlist(fireSenseUtils::defaultClimateVariablesForFire))))
  cl <- fireSenseUtils::fireClimateLayers(list(ignition = "CMD_sm", spread = "Tmax_sm"),
                                          historicalYears = 1990:1995, projectedYears = NULL)
  expect_true("historical_Tmaxsm" %in% names(cl))
  expect_identical(cl$historical_Tmaxsm$vars, "historical_Tmax_sm")
  expect_identical(fireSenseUtils::climateNAnames("Tmaxsm"), "Tmax_sm")
})
