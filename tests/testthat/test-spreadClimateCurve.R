## spreadClimateCurve: an extra covariate C_curve = softplus(C - k, s) per spread climate variable, k the ELF's median and s its sd of C over the
## burnable pixels of all the fit years, stored with the covariates for the ledger. FALSE leaves everything as it was.

curveRun <- function(climateSplit = FALSE, climateCurve = FALSE) {
  f <- wholeFixture()
  withr::local_options(reproducible.cachePath = withr::local_tempdir(), reproducible.useMemoise = FALSE,
                       .local_envir = parent.frame())
  list(f = f, out = annualSpreadCovariates(f$buffers, f$clim, f$allYears, f$burnablePixels, tsds = f$tsds,
                                           firePixelsByYear = f$fires, cutoffForYoungAge = 15,
                                           climateSplit = climateSplit, climateCurve = climateCurve))
}

test_that("with the curve, each table has C_curve = softplus(C - k, s), k the median and s the sd over every burnable pixel of every year", {
  skip_if_not_installed("fireSenseUtils")
  r <- curveRun(climateCurve = TRUE)
  ## burnable pixels 1:16 in years 2001:2004: CMD is pixel + year
  all <- unlist(lapply(2001:2004, function(y) 1:16 + y))
  k <- stats::median(all); s <- stats::sd(all)
  expect_identical(attr(r$out, fireSenseUtils::spreadFitClimateCurveParsTxt), list(CMD = c(k = k, s = s)))
  tabs <- do.call(c, unname(r$out))
  for (x in tabs) {
    expect_true(all(c("CMD", "CMD_curve") %in% names(x)))
    expect_equal(x$CMD_curve, fireSenseUtils::softplus(x$CMD - k, s))
  }
  expect_true(all(unlist(lapply(tabs, function(x) x$CMD_curve)) > 0))      # smooth: never exactly 0
})

test_that("a reused fit's k and s are used as they are", {
  skip_if_not_installed("fireSenseUtils")
  r <- curveRun(climateCurve = list(CMD = c(k = 12, s = 3)))
  for (x in do.call(c, unname(r$out))) expect_equal(x$CMD_curve, fireSenseUtils::softplus(x$CMD - 12, 3))
  expect_identical(attr(r$out, fireSenseUtils::spreadFitClimateCurveParsTxt), list(CMD = c(k = 12, s = 3)))
})

test_that("the curve is of the pixel's value with the split too", {
  skip_if_not_installed("fireSenseUtils")
  r <- curveRun(climateSplit = TRUE, climateCurve = list(CMD = c(k = 12, s = 3)))
  x <- r$out[[1]][[1]]
  expect_true(all(c("CMD_yr", "CMD_dev", "CMD_curve") %in% names(x)))
  expect_false("CMD" %in% names(x))
  expect_equal(x$CMD_curve, fireSenseUtils::softplus(x$CMD_yr + x$CMD_dev - 12, 3))
})

test_that("without the curve the tables and the formula are what they were", {
  skip_if_not_installed("fireSenseUtils")
  plain <- curveRun()$out
  expect_null(attr(plain, fireSenseUtils::spreadFitClimateCurveParsTxt))
  expect_false(any(vapply(do.call(c, unname(plain)), function(x) "CMD_curve" %in% names(x), logical(1))))
  expect_identical(spreadRHS("all", "CMD", "fuelA"), "CMD + youngAge + fuelA")
  expect_identical(spreadRHS("all", "CMD", "fuelA", climateCurve = FALSE), spreadRHS("all", "CMD", "fuelA"))
})

test_that("the formula has the curve term after the variable's other climate terms, and it does not interact", {
  expect_identical(spreadRHS("all", "CMD", "fuelA", climateCurve = TRUE), "CMD + CMD_curve + youngAge + fuelA")
  expect_identical(spreadRHS("all", c("CMD", "Tmax"), "fuelA", climateCurve = TRUE),
                   "CMD + Tmax + CMD_curve + Tmax_curve + youngAge + fuelA")
  expect_identical(spreadRHS("all", "CMD", "fuelA", "climate:fuel", climateSplit = TRUE, climateCurve = TRUE),
                   "CMD_yr + CMD_dev + CMD_curve + youngAge + fuelA + CMD_yr:fuelA")
  expect_identical(spreadRHS("fuel", "CMD", "fuelA", climateCurve = TRUE), spreadRHS("fuel", "CMD", "fuelA"))
})

test_that("the module defines spreadClimateCurve, FALSE by default", {
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("defineModule")),
                parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE))
  params <- eval(as.list(def[[1]][[3]])[["parameters"]], envir = asNamespace("SpaDES.core"))
  i <- match("spreadClimateCurve", params$paramName)
  expect_false(is.na(i))
  expect_identical(params$default[[i]], FALSE)
})

## a reused ledger fit keeps its own curve terms and knots, whatever the parameter says
curveLedger <- function(withPars = TRUE) {
  ledger <- ownAndNeighbourLedger()
  own <- ledger$params[[2]]
  own$CMD_curve <- 0
  ledger$params[[2]] <- own
  if (withPars) ledger[[fireSenseUtils::spreadFitClimateCurveParsTxt]] <- I(list(NULL, list(CMD = c(k = 12, s = 3))))
  ledger
}
initCurve <- function(ledger, doFit, parCurve) {
  e <- ledgerInitEnv(ledger)
  e$sim$.ELFind <- "4.2.2"
  e$sim$climateVariablesForFire <- list(ignition = c("CMD", "cumMDC", "CMDsm", "CMDsp"), spread = "auto")
  e$Par$doFit <- doFit
  e$Par$spreadClimateCurve <- parCurve
  for (s in initBody()) suppressMessages(eval(s, e))
  e
}

test_that("a reused fit with a curve term keeps it and its k and s even when spreadClimateCurve is FALSE", {
  skip_if_not_installed("fireSenseUtils")
  e <- initCurve(curveLedger(), doFit = NA, parCurve = FALSE)
  expect_true(e$mod$spreadClimateCurve)
  expect_identical(e$mod$climateCurvePars, list(CMD = c(k = 12, s = 3)))
  expect_true(useClimateCurve(e$mod, e$Par))
})

test_that("a reused fit without a curve term stays without, also when spreadClimateCurve is TRUE", {
  skip_if_not_installed("fireSenseUtils")
  e <- initCurve(ownAndNeighbourLedger(), doFit = NA, parCurve = TRUE)
  expect_false(e$mod$spreadClimateCurve)
  expect_false(useClimateCurve(e$mod, e$Par))
})

test_that("a refit follows the parameter", {
  skip_if_not_installed("fireSenseUtils")
  e <- initCurve(curveLedger(), doFit = TRUE, parCurve = FALSE)
  expect_null(e$mod$spreadClimateCurve)
  expect_false(useClimateCurve(e$mod, e$Par))
  expect_true(useClimateCurve(e$mod, list(spreadClimateCurve = TRUE)))
})
