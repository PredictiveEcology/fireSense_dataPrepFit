## `spread = "auto"`: selectSpreadClimateVariable() picks this ELF's spread climate variable. The default
## climate variables and their layers are tested in fireSenseUtils (test-fireClimateVariables.R).

## one layer per year, each a constant field with that year's value
yearlyRaster <- function(values, years) {
  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  s <- terra::rast(lapply(values, function(v) terra::setValues(r, rep(v, 4))))
  names(s) <- paste0(fireSenseUtils::yearTxt, years)
  s
}

test_that("auto picks the variable that separates the bad fire years, not the best overall correlation", {
  years <- 2001:2012
  burned <- setNames(c(1:9, 100, 200, 300), paste0(fireSenseUtils::yearTxt, years))   # 3 bad years
  ordinary <- c(1:9, 5, 5, 5)                     # tracks the quiet years, says nothing about the bad ones
  badYears <- c(9, 1, 8, 2, 7, 3, 6, 4, 5, 10, 11, 12)   # noisy, but driest in exactly the bad years
  rasters <- list(ordinary = yearlyRaster(ordinary, years), badYears = yearlyRaster(badYears, years))
  sel <- selectSpreadClimateVariable(rasters, burned, c("ordinary", "badYears"))
  expect_gt(sel[var == "ordinary", rho], sel[var == "badYears", rho])  # a correlation would pick "ordinary"
  expect_equal(sel[var == "badYears", auc], 1)
  expect_identical(sel[chosen == TRUE, var], "badYears")
  expect_identical(sel$nBad, c(3L, 3L))
})

test_that("a variable that runs against fire, or is missing, is not chosen", {
  years <- 2001:2012
  set.seed(3)
  burned <- setNames(c(0, 5, 40, 12, 0, 90, 30, 3, 60, 8, 0, 25), paste0(fireSenseUtils::yearTxt, years))
  dry <- rank(burned) + stats::rnorm(12, 0, 1)
  rasters <- list(dry = yearlyRaster(dry, years), wet = yearlyRaster(-dry, years))
  sel <- selectSpreadClimateVariable(rasters, burned, c("missing", "wet", "dry"))
  expect_identical(sel[chosen == TRUE, var], "dry")
  expect_lt(sel[var == "wet", auc], 0.5)
  expect_true(is.na(sel[var == "missing", auc]))
  expect_equal(sum(sel$chosen), 1L)
})

test_that("auto falls back to the first candidate when no variable separates the bad years", {
  years <- 2001:2012
  burned <- setNames(c(1:9, 100, 200, 300), paste0(fireSenseUtils::yearTxt, years))
  rasters <- list(a = yearlyRaster(c(12:4, 3, 2, 1), years), b = yearlyRaster(c(9, 1, 8, 2, 7, 3, 6, 4, 5, 5, 5, 5), years))  # AUC 0.5
  expect_warning(sel <- selectSpreadClimateVariable(rasters, burned, c("a", "b")), "using a")
  expect_identical(sel[chosen == TRUE, var], "a")
})
