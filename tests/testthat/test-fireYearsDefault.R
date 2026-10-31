## The fitting window is 1985 to the latest year every input can supply, and climate is the last to
## reach a year. A hard-coded end (2002:2025 before 1.2.0.9010) ran past climate, so every project had
## to set `fireYears` itself; the first dataYear postdated 1985, so a 1985 start could not be the default.
## A fire year needs vegetation from an earlier year, so the default starts the year after the first data year.

test_that("the default fire years run from the year after the first data year to the latest climate year", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  defaults <- stats::setNames(md$parameters$default, md$parameters$paramName)
  fireYears <- defaults[["fireYears"]]
  dataYears <- defaults[["dataYears"]]

  expect_identical(min(fireYears), 1986L)
  expect_identical(max(fireYears), climateData::latestHistoricalYear())
  expect_identical(fireYears, min(fireYears):max(fireYears))
  ## every SCANFI V2 year with a fire year after it: 1985, 1990, ..., 2020 (not the old five years)
  expect_identical(dataYears, fireSenseUtils::defaultDataYears())
  expect_identical(dataYears, seq(1985L, max(dataYears), 5L))
  expect_lt(max(dataYears), max(fireYears))
  ## the defaults satisfy the module's own constraints: every fire year after the first data year,
  ## and every data year has fire years
  expect_no_error(yearGroups(dataYears, fireYears))
})

## the default expression of a parameter, as written in the module's `defineParameter()` call
defaultExpr <- function(param) {
  exprs <- parse(file.path(modulePath, moduleName, paste0(moduleName, ".R")), keep.source = FALSE)
  found <- list()
  walk <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], as.name("defineParameter")) && identical(x[[2]], param))
        found[[length(found) + 1L]] <<- x
      for (i in seq_along(x)[-1L]) walk(x[[i]])
    }
  }
  for (e in exprs) walk(e)
  expect_length(found, 1L)
  found[[1]][[4]]
}

test_that("the default is fireSenseUtils::defaultFireYears(), which fireSense_ELFs shares", {
  ## the default is defined once, in fireSenseUtils, so the ELF merging and the fit count the same years
  expect_identical(defaultExpr("fireYears"), quote(fireSenseUtils::defaultFireYears()))
})

test_that("the default dataYears is fireSenseUtils::defaultDataYears(), defined once with the fire years", {
  expect_identical(defaultExpr("dataYears"), quote(fireSenseUtils::defaultDataYears()))
})
