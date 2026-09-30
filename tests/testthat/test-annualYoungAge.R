## youngAge is resolved for each fire year, not once per data year: a pixel that burned between the
## data year and fire year t is young at t, and the non-annual (per data year) tables carry no youngAge
## and have no fuel zeroed. The per-year calculation is fireSenseUtils::youngAgeAtYear(); these tests
## cover how this module feeds it.

tsdRaster <- function(vals) terra::rast(nrows = 2, ncols = 2, vals = vals)

test_that("addAnnualYoungAge puts youngAge in each fire year's table, from the data year's TSD and later fires", {
  ## data year 2000: pixel 1 is old, pixel 2 is young (TSD 3), pixel 3 old, pixel 4 old
  tsds <- list(year2000 = tsdRaster(c(100, 3, 100, 100)))
  ## pixel 1 burns in 2002, pixel 4 in 2001 (a fire not being fitted: it is in no annual table's buffer)
  fires <- list("2002" = 1L, "2001" = 4L)
  mk <- function() data.table::data.table(pixelID = 1:3, CMD = 1)   # pixel 4 is not in a buffer
  ac <- list(`2000` = list(year2001 = mk(), year2004 = mk(), year2020 = mk()))

  out <- addAnnualYoungAge(ac, tsds, fires, cutoffForYoungAge = 15)[["2000"]]

  expect_identical(out$year2001$youngAge, c(0L, 1L, 0L))  # 2001: pixel 2 still young; nothing else burned yet
  expect_identical(out$year2004$youngAge, c(1L, 1L, 0L))  # 2002 fire made pixel 1 young, though old in 2000
  expect_identical(out$year2020$youngAge, c(0L, 0L, 0L))  # all aged out
  expect_true(all(c("pixelID", "CMD") %in% names(out$year2004)))
})

test_that("addAnnualYoungAge keeps the empty table of a year without fires", {
  tsds <- list(year2000 = tsdRaster(c(100, 3, 100, 100)))
  ac <- list(`2000` = list(year2001 = data.table::data.table(pixelID = 1:3, CMD = 1),
                           year2002 = data.table::data.table(pixelID = integer(0), CMD = numeric(0))))
  out <- addAnnualYoungAge(ac, tsds, list(), 15)[["2000"]]
  expect_identical(nrow(out$year2002), 0L)
  expect_true("youngAge" %in% names(out$year2002))
})

test_that("allFirePixelsByYear uses every fire in firePolysForAge, or historicalFireRaster when supplied", {
  tmpl <- terra::rast(nrows = 2, ncols = 2, vals = 1)
  left <- terra::as.polygons(terra::ext(terra::xmin(tmpl), 0, terra::ymin(tmpl), terra::ymax(tmpl)),
                             crs = terra::crs(tmpl))
  sim <- new.env()
  sim$flammableRTMs <- list(year2000 = tmpl)
  sim$firePolysForAge <- list(year1999 = left, year2000 = NULL)
  sim$historicalFireRaster <- NULL
  expect_identical(allFirePixelsByYear(sim), list(`1999` = c(1L, 3L)))

  sim$historicalFireRaster <- terra::rast(nrows = 2, ncols = 2, vals = c(2001, NA, NA, 2001))
  expect_identical(allFirePixelsByYear(sim), list(`2001` = c(1L, 4L)))
})

test_that("nonForestCanBeYoungAge = FALSE stops; TRUE does not", {
  expect_error(stopIfNonForestCannotBeYoung(FALSE), "not supported")
  expect_silent(stopIfNonForestCannotBeYoung(TRUE))
})

test_that("prepare_SpreadFit and prepare_IgnitionFit take youngAge from the per-year calculation", {
  spread <- paste(deparse(body(prepare_SpreadFit)), collapse = " ")
  ## fuels are built without youngAge, so they are not zeroed at the data year
  expect_match(spread, "youngAge = FALSE", fixed = TRUE)
  expect_match(spread, "addAnnualYoungAge(", fixed = TRUE)
  expect_match(spread, "stopIfNonForestCannotBeYoung(", fixed = TRUE)
  ## the annual tables get youngAge before they are combined
  expect_lt(regexpr("addAnnualYoungAge(", spread, fixed = TRUE),
            regexpr("fireSense_annualSpreadFitCovariates <- do.call", spread, fixed = TRUE))

  ignition <- paste(deparse(body(prepare_IgnitionFit)), collapse = " ")
  expect_match(ignition, "prepare_FuelCovsCoarseByYear", fixed = TRUE)
  expect_false(grepl("prepare_FuelCovsCoarse(", ignition, fixed = TRUE))
  expect_match(ignition, "stopIfNonForestCannotBeYoung(", fixed = TRUE)
})
