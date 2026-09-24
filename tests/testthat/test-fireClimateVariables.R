## The fire models' climate variables are set here unless a user supplies them, and canClimateData's
## `climateVariables` follows from them (Eliot, 2026-09-23: "I want the user to NOT need to select").
## Ignition (xgboost) takes several correlated variables; spread takes one.

test_that("the default gives ignition several variables and lets spread choose among them", {
  expect_identical(defaultClimateVariablesForFire$ignition, c("CMD", "cumMDC", "CMD_sm", "CMD_sp"))
  expect_identical(defaultClimateVariablesForFire$spread, "auto")
})

test_that("climate variable names are accepted with or without underscores", {
  expect_identical(climateNAnames(c("CMD_sm", "CMDsm", "cumMDC", "CMDsp")), c("CMD_sm", "CMD_sm", "cumMDC", "CMD_sp"))
  expect_error(climateNAnames(c("CMDsm", "notAVariable")), "notAVariable")
})

test_that("the climate layers cover every variable once, for the fire years", {
  cl <- fireClimateLayers(defaultClimateVariablesForFire, historicalYears = 1985:2024, projected = FALSE)
  expect_setequal(names(cl), c("historical_CMD", "historical_CMDsm", "historical_cumMDC", "historical_CMDsp"))
  expect_identical(cl$historical_CMDsm$.dots$historical_years, 1985:2024)
  expect_identical(cl$historical_CMDsm$fun, quote(calcAsIs))
  ## cumMDC is derived from monthly variables and starts 5 years early (climateData::climateLayers)
  expect_identical(cl$historical_cumMDC$fun, quote(calcCumMDC))
  expect_identical(cl$historical_cumMDC$.dots$historical_years, 1980:2024)
  ## no-underscore names from an older project give the same layers
  cl2 <- fireClimateLayers(list(ignition = c("CMD", "cumMDC", "CMDsm", "CMDsp"), spread = "auto"), 1985:2024, FALSE)
  expect_identical(cl2, cl)
  ## projected layers only when asked
  expect_true(any(grepl("^projected_", names(fireClimateLayers(defaultClimateVariablesForFire, 1985:2024, TRUE, 2025:2044)))))
})

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
