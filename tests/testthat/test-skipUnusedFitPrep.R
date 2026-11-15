## `skipUnusedFitPrep`: a job that only predicts from the ledger's fits does not schedule the spread and escape
## fit-data events (`prepSpreadFitData`, `prepEscapeFitData`); `prepIgnitionFitData` stays (it makes `lightningMaps`
## and `ignitionFitRTM`, which the prediction modules read).

bothFits <- c("fireSense_ignitionFit", "fireSense_spreadFit")

## a simList with this module current, `doFit` as given, and optionally a spreadFit with its own parameters
prepSim <- function(doFit = FALSE, heldOutFold = NA, spreadFit = NULL, skip = TRUE) {
  sim <- suppressMessages(SpaDES.core::simInit())
  sim@params <- list(fireSense_dataPrepFit = list(doFit = doFit, heldOutFold = heldOutFold,
                                                  skipUnusedFitPrep = skip))
  if (!is.null(spreadFit)) sim@params$fireSense_spreadFit <- spreadFit
  sim@current <- data.table::data.table(moduleName = "fireSense_dataPrepFit")
  sim
}
scheduled <- function(sim) prepEventsToSchedule(bothFits, fitPrepNeeded(sim))

test_that("by default every prep event is scheduled, as before", {
  expect_identical(scheduled(prepSim(FALSE, skip = FALSE)),
                   c("prepSpreadFitData", "prepIgnitionFitData", "prepEscapeFitData"))
})

test_that("skipUnusedFitPrep: a predict-only job keeps only prepIgnitionFitData", {
  expect_identical(scheduled(prepSim(FALSE)), "prepIgnitionFitData")
  ## a spreadFit that only predicts (its default mode is "fit", which does nothing when the ledger has the fit)
  expect_identical(scheduled(prepSim(FALSE, spreadFit = list(mode = "fit", heldOutFold = NA))), "prepIgnitionFitData")
})

test_that("skipUnusedFitPrep: a job that may fit, validate, profile or hold out a fold prepares everything", {
  all3 <- c("prepSpreadFitData", "prepIgnitionFitData", "prepEscapeFitData")
  expect_identical(scheduled(prepSim(TRUE)), all3)
  expect_identical(scheduled(prepSim(NA)), all3)                 # the ledgers decide: not known here
  expect_identical(scheduled(prepSim(FALSE, heldOutFold = 2L)), all3)
  expect_identical(scheduled(prepSim(FALSE, spreadFit = list(mode = c("fit", "validate")))), all3)
  expect_identical(scheduled(prepSim(FALSE, spreadFit = list(mode = "profile"))), all3)
})

test_that("a project that prepares only the spread fit schedules nothing when it only predicts", {
  expect_length(prepEventsToSchedule("fireSense_spreadFit", FALSE), 0L)
})
