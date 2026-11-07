## With wholeLandscapeSpreadCovariates the spread-fit covariate tables have every burnable pixel, made by the
## same functions (annualSpreadCovariates(), nonAnnualSpreadCovariates()) called with all burnable pixels
## instead of the buffer pixels. They must have the buffer tables' names and columns, and the same values
## on the buffer pixels.

wholeFixture <- function() {
  r <- terra::rast(nrows = 4, ncols = 5)
  years <- 2001:2004
  clim <- list(CMD = terra::rast(lapply(years, function(y) { x <- r; terra::values(x) <- seq_len(20) + y; x })))
  names(clim$CMD) <- paste0("year", years)
  tsd <- function(v) { x <- r; terra::values(x) <- v; x }
  tsds <- list(year2000 = tsd(c(rep(100, 19), 3)), year2002 = tsd(rep(100, 20)))
  allYears <- list(`2000` = c("year2001", "year2002"), `2002` = c("year2003", "year2004"))
  burnable <- seq_len(16)                                        # pixels 17:20 cannot burn
  burnablePixels <- setNames(rep(list(burnable), 4), unlist(allYears))
  buffers <- list(year2001 = data.table::data.table(pixelID = c(2L, 3L, 5L)),
                  year2002 = data.table::data.table(pixelID = integer(0)),    # a year without fire
                  year2003 = data.table::data.table(pixelID = c(3L, 16L, 17L)),
                  year2004 = data.table::data.table(pixelID = c(7L, 8L)))
  list(clim = clim, tsds = tsds, allYears = allYears, burnablePixels = burnablePixels, buffers = buffers,
       fires = list(`2003` = 9L))
}

test_that("whole annual tables have the buffer tables' names and rows on the buffer pixels, and every burnable pixel", {
  skip_if_not_installed("fireSenseUtils")
  withr::local_options(reproducible.cachePath = withr::local_tempdir(), reproducible.useMemoise = FALSE)
  f <- wholeFixture()
  run <- function(pixels) annualSpreadCovariates(pixels, f$clim, f$allYears, f$burnablePixels,
                                                 tsds = f$tsds, firePixelsByYear = f$fires, cutoffForYoungAge = 15)
  buffer <- do.call(c, unname(run(f$buffers)))
  whole <- do.call(c, unname(run(lapply(f$burnablePixels, function(p) data.table::data.table(pixelID = p))))) 

  expect_identical(names(whole), names(buffer))
  expect_identical(names(whole), unname(unlist(f$allYears)))
  for (fy in names(buffer)) {
    expect_identical(names(whole[[fy]]), names(buffer[[fy]]))
    expect_setequal(whole[[fy]]$pixelID, f$burnablePixels[[fy]])
    w <- whole[[fy]][match(buffer[[fy]]$pixelID, whole[[fy]]$pixelID), ]
    expect_equal(w$CMD, buffer[[fy]]$CMD)
    expect_equal(w$youngAge, buffer[[fy]]$youngAge)
  }
  expect_identical(nrow(buffer$year2002), 0L)
  expect_false(17L %in% buffer$year2003$pixelID)   # a non-burnable buffer pixel is dropped from both
  expect_identical(whole$year2001$youngAge[whole$year2001$pixelID == 9L], 0L)
  expect_identical(whole$year2004$youngAge[whole$year2004$pixelID == 9L], 1L)  # burned in 2003
})

test_that("whole non-annual tables have the buffer tables' names and rows on the buffer pixels", {
  veg <- data.table::data.table(pixelID = rep(1:6, 2), year = rep(c(2000L, 2002L), each = 6),
                                fuel = c(1, 2, NA, 4, 5, 6, 7, 8, 9, 10, 11, 12), burned = 0L)
  cols <- c("pixelID", "fuel")
  inBuffer <- veg[pixelID %in% c(2L, 3L, 5L)]
  buffer <- Map(dy = c("2000", "2002"), function(dy) nonAnnualSpreadCovariates(inBuffer[year == as.integer(dy)], cols))
  whole <- Map(dy = c("2000", "2002"), function(dy) nonAnnualSpreadCovariates(veg[year == as.integer(dy)], cols))

  expect_identical(names(whole), names(buffer))
  expect_identical(whole[["2000"]]$pixelID, c(1L, 2L, 4L, 5L, 6L))   # the NA fuel row is dropped from both
  expect_identical(whole[["2002"]]$pixelID, 1:6)
  for (dy in names(buffer)) {
    w <- whole[[dy]][match(buffer[[dy]]$pixelID, whole[[dy]]$pixelID), ]
    expect_identical(w$fuel, buffer[[dy]]$fuel)
  }
})

test_that("prepare_SpreadFit picks the buffer or the burnable pixels by the parameter, into the existing outputs", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_false(md$parameters$default[md$parameters$paramName == "wholeLandscapeSpreadCovariates"][[1]])
  expect_false(any(grepl("Whole$", md$outputObjects$objectName)))
  body <- paste(deparse(body(prepare_SpreadFit)), collapse = " ")
  expect_match(body, "wholeLandscapeSpreadCovariates", fixed = TRUE)
  expect_match(body, "pixelsByYear <- if (wholeSpreadCovariates)", fixed = TRUE)
  expect_equal(lengths(regmatches(body, gregexpr("annualSpreadCovariates(", body, fixed = TRUE))), 1L) # one call, whichever pixels
})
