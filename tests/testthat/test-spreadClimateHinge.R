## spreadClimateHinge: an extra covariate C_hinge = pmax(0, C - k) per spread climate variable, k the ELF's median of C over the
## burnable pixels of all the fit years, stored with the covariates for the ledger. FALSE leaves everything as it was.

hingeRun <- function(climateSplit = FALSE, climateHinge = FALSE) {
  f <- wholeFixture()
  withr::local_options(reproducible.cachePath = withr::local_tempdir(), reproducible.useMemoise = FALSE,
                       .local_envir = parent.frame())
  list(f = f, out = annualSpreadCovariates(f$buffers, f$clim, f$allYears, f$burnablePixels, tsds = f$tsds,
                                           firePixelsByYear = f$fires, cutoffForYoungAge = 15,
                                           climateSplit = climateSplit, climateHinge = climateHinge))
}

test_that("with the hinge, each table has C_hinge = pmax(0, C - k), k the median over every burnable pixel of every year", {
  skip_if_not_installed("fireSenseUtils")
  r <- hingeRun(climateHinge = TRUE)
  ## burnable pixels 1:16 in years 2001:2004: CMD is pixel + year
  k <- stats::median(unlist(lapply(2001:2004, function(y) 1:16 + y)))
  expect_identical(unname(unlist(attr(r$out, fireSenseUtils::spreadFitClimateHingeKnotsTxt))), k)
  tabs <- do.call(c, unname(r$out))
  for (x in tabs) {
    expect_true(all(c("CMD", "CMD_hinge") %in% names(x)))
    expect_equal(x$CMD_hinge, pmax(0, x$CMD - k))
  }
  expect_gt(max(unlist(lapply(tabs, function(x) x$CMD_hinge))), 0)        # the fixture reaches above k
  expect_identical(min(unlist(lapply(tabs, function(x) x$CMD_hinge))), 0)
})

test_that("a reused fit's knots are used as they are", {
  skip_if_not_installed("fireSenseUtils")
  r <- hingeRun(climateHinge = list(CMD = 12))
  for (x in do.call(c, unname(r$out))) expect_equal(x$CMD_hinge, pmax(0, x$CMD - 12))
  expect_identical(attr(r$out, fireSenseUtils::spreadFitClimateHingeKnotsTxt), list(CMD = 12))
})

test_that("the hinge is of the pixel's value with the split too", {
  skip_if_not_installed("fireSenseUtils")
  r <- hingeRun(climateSplit = TRUE, climateHinge = list(CMD = 12))
  x <- r$out[[1]][[1]]
  expect_true(all(c("CMD_yr", "CMD_dev", "CMD_hinge") %in% names(x)))
  expect_false("CMD" %in% names(x))
  expect_equal(x$CMD_hinge, pmax(0, x$CMD_yr + x$CMD_dev - 12))
})

test_that("without the hinge the tables and the formula are what they were", {
  skip_if_not_installed("fireSenseUtils")
  plain <- hingeRun()$out
  expect_null(attr(plain, fireSenseUtils::spreadFitClimateHingeKnotsTxt))
  expect_false(any(vapply(do.call(c, unname(plain)), function(x) "CMD_hinge" %in% names(x), logical(1))))
  expect_identical(spreadRHS("all", "CMD", "fuelA"), "CMD + youngAge + fuelA")
  expect_identical(spreadRHS("all", "CMD", "fuelA", climateHinge = FALSE), spreadRHS("all", "CMD", "fuelA"))
})

test_that("the formula has the hinge term after the variable's other climate terms, and it does not interact", {
  expect_identical(spreadRHS("all", "CMD", "fuelA", climateHinge = TRUE), "CMD + CMD_hinge + youngAge + fuelA")
  expect_identical(spreadRHS("all", c("CMD", "Tmax"), "fuelA", climateHinge = TRUE),
                   "CMD + Tmax + CMD_hinge + Tmax_hinge + youngAge + fuelA")
  expect_identical(spreadRHS("all", "CMD", "fuelA", "climate:fuel", climateSplit = TRUE, climateHinge = TRUE),
                   "CMD_yr + CMD_dev + CMD_hinge + youngAge + fuelA + CMD_yr:fuelA")
  expect_identical(spreadRHS("fuel", "CMD", "fuelA", climateHinge = TRUE), spreadRHS("fuel", "CMD", "fuelA"))
})

test_that("the module defines spreadClimateHinge, FALSE by default", {
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("defineModule")),
                parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE))
  params <- eval(as.list(def[[1]][[3]])[["parameters"]], envir = asNamespace("SpaDES.core"))
  i <- match("spreadClimateHinge", params$paramName)
  expect_false(is.na(i))
  expect_identical(params$default[[i]], FALSE)
})

## a reused ledger fit keeps its own hinge terms and knots, whatever the parameter says
hingeLedger <- function(withKnots = TRUE) {
  ledger <- ownAndNeighbourLedger()
  own <- ledger$params[[2]]
  own$CMD_hinge <- 0
  ledger$params[[2]] <- own
  if (withKnots) ledger[[fireSenseUtils::spreadFitClimateHingeKnotsTxt]] <- I(list(NULL, list(CMD = 12)))
  ledger
}
initHinge <- function(ledger, doFit, parHinge) {
  e <- ledgerInitEnv(ledger)
  e$sim$.ELFind <- "4.2.2"
  e$sim$climateVariablesForFire <- list(ignition = c("CMD", "cumMDC", "CMDsm", "CMDsp"), spread = "auto")
  e$Par$doFit <- doFit
  e$Par$spreadClimateHinge <- parHinge
  for (s in initBody()) suppressMessages(eval(s, e))
  e
}

test_that("a reused fit with a hinge term keeps it and its knots even when spreadClimateHinge is FALSE", {
  skip_if_not_installed("fireSenseUtils")
  e <- initHinge(hingeLedger(), doFit = NA, parHinge = FALSE)
  expect_true(e$mod$spreadClimateHinge)
  expect_identical(e$mod$climateHingeKnots, list(CMD = 12))
  expect_true(useClimateHinge(e$mod, e$Par))
})

test_that("a reused fit without a hinge term stays without, also when spreadClimateHinge is TRUE", {
  skip_if_not_installed("fireSenseUtils")
  e <- initHinge(ownAndNeighbourLedger(), doFit = NA, parHinge = TRUE)
  expect_false(e$mod$spreadClimateHinge)
  expect_false(useClimateHinge(e$mod, e$Par))
})

test_that("a refit follows the parameter", {
  skip_if_not_installed("fireSenseUtils")
  e <- initHinge(hingeLedger(), doFit = TRUE, parHinge = FALSE)
  expect_null(e$mod$spreadClimateHinge)
  expect_false(useClimateHinge(e$mod, e$Par))
  expect_true(useClimateHinge(e$mod, list(spreadClimateHinge = TRUE)))
})
