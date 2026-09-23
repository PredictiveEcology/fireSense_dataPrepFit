## The fire models' climate variables are set here unless a user supplies them, and canClimateData's
## `climateVariables` follows from them (Eliot, 2026-09-23: "I want the user to NOT need to select").
## Ignition (xgboost) takes several correlated variables; spread takes one.

test_that("the default gives ignition several variables and spread one", {
  expect_identical(defaultClimateVariablesForFire$ignition, c("CMD_sm", "cumMDC", "CMD_sp"))
  expect_identical(defaultClimateVariablesForFire$spread, "CMD_sm")
})

test_that("climate variable names are accepted with or without underscores", {
  expect_identical(climateNAnames(c("CMD_sm", "CMDsm", "cumMDC", "CMDsp")), c("CMD_sm", "CMD_sm", "cumMDC", "CMD_sp"))
  expect_error(climateNAnames(c("CMDsm", "notAVariable")), "notAVariable")
})

test_that("the climate layers cover every variable once, for the fire years", {
  cl <- fireClimateLayers(defaultClimateVariablesForFire, historicalYears = 1985:2024, projected = FALSE)
  expect_setequal(names(cl), c("historical_CMDsm", "historical_cumMDC", "historical_CMDsp"))
  expect_identical(cl$historical_CMDsm$.dots$historical_years, 1985:2024)
  expect_identical(cl$historical_CMDsm$fun, quote(calcAsIs))
  ## cumMDC is derived from monthly variables and starts 5 years early (climateData::climateLayers)
  expect_identical(cl$historical_cumMDC$fun, quote(calcCumMDC))
  expect_identical(cl$historical_cumMDC$.dots$historical_years, 1980:2024)
  ## no-underscore names from an older project give the same layers
  cl2 <- fireClimateLayers(list(ignition = c("CMDsm", "cumMDC", "CMDsp"), spread = "CMDsm"), 1985:2024, FALSE)
  expect_identical(cl2, cl)
  ## projected layers only when asked
  expect_true(any(grepl("^projected_", names(fireClimateLayers(defaultClimateVariablesForFire, 1985:2024, TRUE, 2025:2044)))))
})
