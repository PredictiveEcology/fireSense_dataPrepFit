## A satellite vegetation year already shows part of that year's burns (6.6.2: pixels burned in 2010
## have mean biomass 64 in SCANFI 2005 and 26 in SCANFI 2010, unburned pixels ~79 in both). A fire year
## therefore uses the latest data year STRICTLY BEFORE it, so the fit does not learn that fire burns where
## the vegetation already looks burned.
test_that("yearGroups gives each fire year the latest data year strictly before it", {
  dataYears <- c(1985L, 1990L, 2000L, 2010L, 2020L)
  fireYears <- 1986:2025

  groups <- yearGroups(dataYears, fireYears, minmaxOnly = FALSE)
  expect_identical(names(groups), as.character(dataYears))
  expect_identical(groups[["1985"]], 1986:1990)
  expect_identical(groups[["1990"]], 1991:2000)
  expect_identical(groups[["2000"]], 2001:2010)
  expect_identical(groups[["2010"]], 2011:2020)
  expect_identical(groups[["2020"]], 2021:2025)
  expect_identical(unname(unlist(groups)), fireYears)

  expect_identical(yearGroups(dataYears, fireYears)[["2000"]], c(2001L, 2010L))
})

test_that("a fire year pairs with the data year before it: 1990 with 1985, 1991 with 1990, 2020 with 2010, 2024 with 2020", {
  dataYears <- c(1985L, 1990L, 2000L, 2010L, 2020L)
  fireYears <- c(1990L, 1991L, 2020L, 2024L)
  expect_identical(vegetationDataYear(dataYears, fireYears), c(1985L, 1990L, 2010L, 2020L))
})

test_that("yearGroups keeps the grouping for three data years", {
  groups <- yearGroups(c(2000L, 2010L, 2020L), 2002:2025, minmaxOnly = FALSE)
  expect_identical(groups, list(`2000` = 2002:2010, `2010` = 2011:2020, `2020` = 2021:2025))
})

test_that("yearGroups stops on fire years not after the first data year and on data years without fire years", {
  expect_error(yearGroups(c(1985L, 1990L), 1984:1995), "after the first of dataYears")
  ## a fire year equal to the first data year is not allowed: that vegetation already shows its fires
  expect_error(yearGroups(c(1985L, 1990L), 1985:1995), "after the first of dataYears.*1985")
  expect_error(yearGroups(c(2000L, 2010L, 2020L), 2012:2025), "no fireYears")
})
