## yearGroups() built the groups with rep(mm[1:2], diff(mm)), which only works for exactly three
## dataYears: four or more stopped with "invalid 'times' argument", so phase 1 could not use
## dataYears = c(1985, 1990, 2000, 2010, 2020).
test_that("yearGroups gives each fire year to the latest data year at or before it", {
  dataYears <- c(1985L, 1990L, 2000L, 2010L, 2020L)
  fireYears <- 1985:2025

  groups <- yearGroups(dataYears, fireYears, minmaxOnly = FALSE)
  expect_identical(names(groups), as.character(dataYears))
  expect_identical(groups[["1985"]], 1985:1989)
  expect_identical(groups[["1990"]], 1990:1999)
  expect_identical(groups[["2020"]], 2020:2025)
  expect_identical(unname(unlist(groups)), fireYears)

  expect_identical(yearGroups(dataYears, fireYears)[["2000"]], c(2000L, 2009L))
})

test_that("yearGroups keeps the grouping for three data years", {
  groups <- yearGroups(c(2000L, 2010L, 2020L), 2002:2025, minmaxOnly = FALSE)
  expect_identical(groups, list(`2000` = 2002:2009, `2010` = 2010:2019, `2020` = 2020:2025))
})

test_that("yearGroups stops on fire years without vegetation data and on data years without fire years", {
  expect_error(yearGroups(c(1985L, 1990L), 1984:1995), "before the first dataYear")
  expect_error(yearGroups(c(2000L, 2010L, 2020L), 2012:2025), "no fireYears")
})
