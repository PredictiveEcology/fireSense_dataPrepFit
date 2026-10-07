## An explicitly requested fire year must never be dropped quietly. The module used to warn and
## truncate P(sim)$fireYears to whatever the climate rasters happened to cover, so a short climate
## raster silently changed the fitting window and the fit still reported success. With the launcher
## pinning one window per campaign, that could leave some ELFs fit on 1985:2024 and others on
## 1985:2022 with nothing in the results to say so.

test_that("checkClimateYears passes when climate covers every requested year", {
  allYearsVect <- paste0(fireSenseUtils::yearTxt, 2002:2010)
  whAvailable <- rep(TRUE, length(allYearsVect))

  expect_silent(checkClimateYears(allYearsVect, whAvailable))
  expect_identical(checkClimateYears(allYearsVect, whAvailable), allYearsVect)
})

test_that("checkClimateYears stops and names every missing year", {
  allYearsVect <- paste0(fireSenseUtils::yearTxt, 1985:2024)
  ## the real case this came from: the climate tile index stopped at 2022
  whAvailable <- as.integer(gsub(fireSenseUtils::yearTxt, "", allYearsVect)) <= 2022

  expect_error(checkClimateYears(allYearsVect, whAvailable), "2023, 2024")
  expect_error(checkClimateYears(allYearsVect, whAvailable), "2 of the 40 requested fireYears")
  ## it must not offer truncation as a silent fallback
  expect_error(checkClimateYears(allYearsVect, whAvailable), "NOT truncated automatically")
})

test_that("checkClimateYears reports missing years from anywhere in the range, not just the end", {
  allYearsVect <- paste0(fireSenseUtils::yearTxt, 2010:2015)
  whAvailable <- !(allYearsVect %in% paste0(fireSenseUtils::yearTxt, c(2012, 2013)))

  err <- tryCatch(checkClimateYears(allYearsVect, whAvailable), error = conditionMessage)
  expect_match(err, "2012, 2013")
  ## the old behaviour kept min:max of what was available, which would have silently kept
  ## 2010:2015 and dropped the hole in the middle without saying so
  expect_no_match(err, "2010:2015")
})
