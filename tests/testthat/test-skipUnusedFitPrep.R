## `skipUnusedFitPrep`: a job that only predicts from the ledger's fits makes only the spread FIRE data
## (`prepSpreadFitFireData`), not the spread covariates. The ignition and escape data are made as in a full run,
## from the same fire data, so `ignitionFirePoints` and `fireSense_escapeCovariates` are the same.

bothFits <- c("fireSense_ignitionFit", "fireSense_spreadFit")
all4 <- c("prepSpreadFitData", "prepIgnitionFitData", "prepEscapeFitData")
fireOnly <- c("prepSpreadFitFireData", "prepIgnitionFitData", "prepEscapeFitData")

## a simList with this module current, `doFit` as given, and optionally a spreadFit with its own parameters
prepSim <- function(doFit = FALSE, heldOutFold = NA, spreadFit = NULL, skip = TRUE) {
  sim <- suppressMessages(SpaDES.core::simInit())
  sim@params <- list(fireSense_dataPrepFit = list(doFit = doFit, heldOutFold = heldOutFold,
                                                  skipUnusedFitPrep = skip))
  if (!is.null(spreadFit)) sim@params$fireSense_spreadFit <- c(list(skipUnusedFitPrep = skip, doFit = doFit), spreadFit)
  sim@current <- data.table::data.table(moduleName = "fireSense_dataPrepFit")
  sim
}
scheduled <- function(sim) prepEventsToSchedule(bothFits, fitPrepNeeded(sim))

test_that("by default every prep event is scheduled, as before", {
  expect_identical(scheduled(prepSim(FALSE, skip = FALSE)), all4)
})

test_that("skipUnusedFitPrep: a predict-only job makes the fire data, not the spread covariates", {
  expect_identical(scheduled(prepSim(FALSE)), fireOnly)
  expect_identical(scheduled(prepSim(FALSE, spreadFit = list(mode = "fit", heldOutFold = NA))), fireOnly)
})

test_that("skipUnusedFitPrep: a job that may fit, validate, profile or hold out a fold prepares everything", {
  expect_identical(scheduled(prepSim(TRUE)), all4)
  expect_identical(scheduled(prepSim(NA)), all4)                 # the ledgers decide: not known here
  expect_identical(scheduled(prepSim(FALSE, heldOutFold = 2L)), all4)
  expect_identical(scheduled(prepSim(FALSE, spreadFit = list(mode = c("fit", "validate")))), all4)
  expect_identical(scheduled(prepSim(FALSE, spreadFit = list(mode = "profile"))), all4)
})

test_that("a project that prepares only the spread fit schedules nothing when it only predicts", {
  expect_length(prepEventsToSchedule("fireSense_spreadFit", FALSE), 0L)
})

test_that("the modules must agree on skipUnusedFitPrep and doFit", {
  sim <- prepSim(FALSE, spreadFit = list())
  sim@params$fireSense_spreadFit$skipUnusedFitPrep <- FALSE
  expect_error(fitPrepNeeded(sim), "skipUnusedFitPrep")
  sim <- prepSim(FALSE, spreadFit = list())
  sim@params$fireSense_spreadFit$doFit <- NA
  expect_error(fitPrepNeeded(sim), "doFit")
})

## ignitionFirePoints and the escape covariates read the fire data (spreadFirePolys, escapeFireSet) and nothing else
## that only the spread covariates make; the fire half is one function that a full run calls too
test_that("the fire data are made by one function that the full spread prep calls", {
  defs <- moduleDefs()
  expect_true("prepare_SpreadFitFires" %in% all.names(defs[["prepare_SpreadFit"]]))
  expect_false(any(c("prepare_SpreadFitFire_Raster", "prepare_SpreadFitFire_Vector") %in%
                     all.names(defs[["prepare_SpreadFit"]])))
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  made <- md$outputObjects$objectName[vapply(md$outputObjects$events, function(e) "prepSpreadFitFireData" %in% e,
                                             logical(1))]
  expect_true(all(c("spreadFirePolys", "escapeFireSet") %in% made))
  expect_false(any(c("fireSense_annualSpreadFitCovariates", "fireSense_spreadFormula") %in% made))
})

test_that("prepare_SpreadFitFires runs the raster or the vector builder, as the full run did", {
  for (raster in c(TRUE, FALSE)) {
    sim <- suppressMessages(SpaDES.core::simInit())
    sim@params <- list(fireSense_dataPrepFit = list(useRasterizedFireForSpread = raster))
    sim@current <- data.table::data.table(moduleName = "fireSense_dataPrepFit")
    ran <- character()
    e <- new.env(parent = environment(prepare_SpreadFitFires))
    e$mod <- new.env()
    e$Par <- list(dataYears = 2000, fireYears = 2001:2002)
    e$prepare_SpreadFitFire_Raster <- function(sim) { ran <<- c(ran, "raster"); sim }
    e$prepare_SpreadFitFire_Vector <- function(sim) { ran <<- c(ran, "vector"); sim }
    f <- prepare_SpreadFitFires
    environment(f) <- e
    f(sim)
    expect_identical(ran, if (raster) "raster" else "vector")
    expect_false(is.null(e$mod$allYears))
  }
})
