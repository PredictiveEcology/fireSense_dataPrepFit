## A hindcast (fireSense_dataPrepPredict, vegetationYear = 1985) overwrites sim$rstLCC_RTM and
## sim$flammableRTM with 1985 before prepare_SpreadFit() runs. Whether treed wetland is a spread
## covariate (minCovariateProp) and the fuel-class land-cover shares must not follow that: they come
## from the data year the decision year's vegetation uses -- the latest data year strictly before it, the
## rule yearGroups() assigns vegetation to fire years with -- read from the per-year lists.

lccRaster <- function(vals) terra::rast(nrows = 2, ncols = 5, vals = vals)

lccSim <- function() {
  lcc2010 <- lccRaster(c(rep(1L, 8), 2L, 2L))   # class 2 is 20%
  lcc2020 <- lccRaster(c(rep(1L, 9), 2L))       # class 2 is 10%
  flam <- lccRaster(rep(1L, 10))
  list(rstLCCs = list(year2010 = lcc2010, year2020 = lcc2020),
       flammableRTMs = list(year2010 = flam, year2020 = flam))
}

test_that("yearGroups assigns each fire year the data year vegetationDataYear() gives", {
  dataYears <- c(1985L, 1990L, 2000L, 2010L, 2020L)
  groups <- yearGroups(dataYears, 1986:2024, minmaxOnly = FALSE)
  expect_identical(unlist(Map(function(dy, ys) rep(dy, length(ys)), as.integer(names(groups)), groups),
                          use.names = FALSE),
                   vegetationDataYear(dataYears, 1986:2024))
})

test_that("lccShareForYear reads the list element of the data year the year's vegetation uses", {
  sim <- lccSim()
  share2010 <- c(`1` = 0.8, `2` = 0.2)
  share2020 <- c(`1` = 0.9, `2` = 0.1)
  expect_equal(lccShareForYear(sim, 2011L), share2010)
  expect_equal(lccShareForYear(sim, 2015L), share2010) # 2015 uses 2010, not the latest year
  expect_equal(lccShareForYear(sim, 2020L), share2010) # 2020 burns are already in the 2020 vegetation
  expect_equal(lccShareForYear(sim, 2021L), share2020)
  expect_equal(lccShareForYear(sim, 2024L), share2020)
  expect_error(lccShareForYear(sim, 2010L), "before the first dataYear") # no vegetation before 2010

  ## the shared objects a hindcast overwrites with another year change nothing
  sim$rstLCC_RTM <- sim$rstLCCs$year2010
  sim$flammableRTM <- sim$flammableRTMs$year2010
  expect_equal(lccShareForYear(sim, 2024L), share2020)
})

test_that("the fit's structure decisions use the fit's decision year, not the shared year's rasters", {
  for (f in c("dataPrepBuild", "prepare_SpreadFit")) {
    txt <- paste(deparse(body(get(f))), collapse = " ")
    expect_true(grepl("lccShareForYear\\(sim,\\s*fitDecisionYear\\(sim\\)\\)", txt), info = f)
    expect_false(grepl("lccFlammableShare(sim$rstLCC_RTM", txt, fixed = TRUE), info = f)
  }
  expect_true(grepl("max(P(sim)$fireYears)", paste(deparse(body(fitDecisionYear)), collapse = " "), fixed = TRUE))
  ## the dominant/secondary fuel classes are chosen on the decision year's element of the lists, not tail()
  txt <- paste(deparse(body(prepare_SpreadFit)), collapse = " ")
  expect_true(grepl("fitDecisionYearName(sim)", txt, fixed = TRUE))
  expect_false(grepl("tail(sim$cohortDatas", txt, fixed = TRUE))
})

test_that("fitDecisionYearName is the list element of the decision year's data year", {
  ## P(sim)$fireYears is not available outside a simList: the element rule alone
  sim <- lccSim()
  expect_identical(fireSenseUtils::dataYearElement(sim$rstLCCs, 2024L, before = TRUE), "year2020")
  expect_identical(fireSenseUtils::dataYearElement(sim$rstLCCs, 2020L, before = TRUE), "year2010")
  expect_true(grepl("dataYearElement(sim$rstLCCs, fitDecisionYear(sim), before = TRUE)",
                    paste(trimws(deparse(body(fitDecisionYearName))), collapse = " "), fixed = TRUE))
})
