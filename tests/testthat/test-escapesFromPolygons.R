## 2026-10-25 (Eliot): "The escape model is about the transition from point to polygon." An escape is an
## ignition that became a fire polygon of at least escapeSizeHa, one per polygon; before, it was an NFDB point
## whose reported SIZE_HA reached escapeSizeHa. Since 2026-11-10 `prepare_IgnitionFit()` marks them (`escaped`,
## `fireSenseUtils::escapeIgnitionPoints()`; see test-coreEscapeFires.R) and `prepare_EscapeFit()` counts them.

test_that("prepare_EscapeFit counts the escaped ignition points, one per escape", {
  skip_if_not_installed("fireSenseUtils")
  ras <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10000, ymin = 0, ymax = 10000,
                     crs = "EPSG:3978", vals = 1)  # coarse cells of 100 ha
  sq <- function(x0, y0, w) sprintf("POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))",
                                    x0, y0, x0 + w, y0, x0 + w, y0 + w, x0, y0 + w, x0, y0)
  polys <- terra::vect(c(sq(100, 100, 800), sq(2100, 100, 400)), crs = "EPSG:3978") # 64 ha and 16 ha
  polys$YEAR <- 2000L
  ## two points in the 64 ha polygon (10 and 30 ha reported), a 60 ha point in the 16 ha polygon,
  ## an 80 ha point in no polygon
  pts <- terra::vect(cbind(c(300, 700, 2300, 5500), c(300, 700, 300, 5500)), crs = "EPSG:3978")
  pts$YEAR <- 2000L
  pts$CAUSE <- "L"
  pts$SIZE_HA <- c(10, 30, 60, 80)
  pts$escaped <- c(TRUE, FALSE, FALSE, FALSE) # one escape: the 64 ha polygon, matched once for its two points
  cells <- terra::cellFromXY(ras, terra::crds(pts))

  sim <- new.env()
  sim$ignitionFirePoints <- pts
  sim$spreadFirePolys <- list(year2000 = polys)
  sim$flammableRTMs <- list(year1985 = terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 10000, ymin = 0,
                                                    ymax = 10000, crs = "EPSG:3978")) # 6.25 ha pixels
  sim$historicalClimateRasters <- list(CMD = ras)
  sim$climateVariablesForFire <- list(ignition = "CMD")
  sim$fireSense_ignitionCovariates <- data.table::data.table(
    pixelID = unique(cells), year = 2000, ignitions = c(2L, 1L, 1L), fuelA = 1, CMD = 1)
  env <- new.env(parent = environment(prepare_EscapeFit))
  env$Par <- list(escapeSizeHa = 50)
  env$P <- function(sim) list(igAggFactor = 1)
  env$Cache <- function(x, ...) x
  f <- prepare_EscapeFit
  environment(f) <- env
  suppressWarnings(f(sim))

  esc <- sim$fireSense_escapeCovariates
  expect_identical(esc[pixelID == cells[1], escapes], 1L)
  expect_identical(sum(esc$escapes), 1L)                    # the point-size rule had 2: the 60 and 80 ha points

  sim$spreadFirePolys <- NULL
  expect_error(f(sim), "needs firePolys")
})
