## The fitting window is 1985 to the latest year every input can supply, and climate is the last to
## reach a year. A hard-coded end (2002:2025 before 1.2.0.9010) ran past climate, so every project had
## to set `fireYears` itself; the first dataYear postdated 1985, so a 1985 start could not be the default.

test_that("the default fire years run from the first data year to the latest climate year", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  defaults <- stats::setNames(md$parameters$default, md$parameters$paramName)
  fireYears <- defaults[["fireYears"]]
  dataYears <- defaults[["dataYears"]]

  expect_identical(min(fireYears), 1985L)
  expect_identical(max(fireYears), climateData::latestHistoricalYear())
  expect_identical(fireYears, min(fireYears):max(fireYears))
  expect_identical(dataYears, c(1985L, 1990L, 2000L, 2010L, 2020L))
  ## the defaults satisfy the module's own constraints: no fire year before the first data year,
  ## and every data year has fire years
  expect_no_error(yearGroups(dataYears, fireYears))
})
