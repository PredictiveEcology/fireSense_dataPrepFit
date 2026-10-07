## Empty spread covariates are dropped before the formula is built. Fuel is stored as log biomass, so a
## fuel class with no biomass is on the logMinB floor (3.6) everywhere, not 0: a colSums == 0 test never
## dropped it and it entered the model carrying no information. fireSenseUtils::emptySpreadCovariates()
## knows the floor; this checks the module uses it, and what it does to a covariate table like the module's.

moduleBody <- function(name) {
  defs <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name(name)),
                 parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE))
  stopifnot(length(defs) == 1L)
  body(eval(defs[[1]][[3]]))
}

test_that("prepare_SpreadFit drops empty covariates with emptySpreadCovariates(), not a column sum", {
  skip_if_not_installed("fireSenseUtils")
  body <- moduleBody("prepare_SpreadFit")
  expect_true("emptySpreadCovariates" %in% all.names(body))
  expect_false("colSums" %in% all.names(body))
})

test_that("a fuel column on the logMinB floor everywhere is dropped, as an all-zero indicator is", {
  skip_if_not_installed("fireSenseUtils")
  dt <- data.table::data.table(
    BlkSprc = fireSenseUtils::logMinB(c(0, 0, 0)), WhtSprc = fireSenseUtils::logMinB(c(0, 400, 9000)),
    nfLCC_100 = c(0, 0, 0), nfLCC_50 = c(1, 0, 0))
  expect_setequal(fireSenseUtils::emptySpreadCovariates(dt, names(dt)), c("BlkSprc", "nfLCC_100"))
})
