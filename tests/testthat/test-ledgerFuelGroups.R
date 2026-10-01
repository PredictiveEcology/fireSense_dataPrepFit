## With a SpreadFit in the ledger, the ignition and escape covariates of this run are still built here, from
## `sim$nonForestedLCCGroups` and `sim$missingLCCgroup` (dataPrepBuild). fireSense_dataPrepPredict builds the
## prediction covariates from the ledger's groups, so the fit must use them too. Before, they stayed the
## module default (`nf`), and fireSense_ignitionPredict stopped with "column not found: [nf]" (carbon run,
## ELF 4.2.2, 2026-09-29). This runs the module's own Init() with stand-ins for Drive, LandR and climate.

initBody <- function() {
  exprs <- parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE)
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                  identical(x[[2]], as.name("Init")), exprs)
  stopifnot(length(def) == 1L)
  as.list(body(eval(def[[1]][[3]])))[-1]
}

ledgerInitEnv <- function(ledger) {
  sim <- suppressMessages(SpaDES.core::simInit())
  sim@params <- list(fireSense_dataPrepFit = list(heldOutFold = NA, sppEquivCol = "LandR",
                                                  fuelClassCol = "FuelClass", fireYears = 1985:2024))
  sim@current <- data.table::data.table(moduleName = "fireSense_dataPrepFit")
  sim$studyArea <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1),
                                                                               c(0, 1), c(0, 0))))))
  spp <- data.table::data.table(LandR = c("Pice_mar", "Popu_tre"), FuelClass = c("Pice_mar", "Popu_tre"))
  sim$sppEquiv <- spp
  sim$nonForestedLCCGroups <- list(nf = c(40, 50, 60, 80, 100)) # the module default (.inputObjects)
  sim$missingLCCgroup <- "nf"
  sim$climateVariables <- list(historical_CMD = list())
  sim$spreadFitAdditionalColNames <- c("sppEquiv", "nonForestedLCCGroups", "missingLCCgroup", "params")
  e <- new.env()
  e$sim <- sim
  e$mod <- new.env()
  e$Par <- list(heldOutFold = NA, sppEquivCol = "LandR", spreadFitFilename = "latest",
                spreadFitGoogleDriveFolder = "folder", escapeSizeHa = 1)
  e$nonEscapedFireSizes <- function(...) NULL
  e$inputPath <- function(sim) tempdir()
  e$readSpreadFitLedger <- function(...) ledger
  e$sppHarmonize <- function(sppEquiv, sppNameVector, sppEquivCol, ...)
    list(sppEquiv = sppEquiv, sppNameVector = sppNameVector, sppEquivCol = sppEquivCol,
         sppColorVect = c(Pice_mar = "#000001", Popu_tre = "#000002", Mixed = "#000003"))
  e$addFitClimateVariables <- function(climateVariables, ...) climateVariables
  e$ledgerColumns <- function(ledger, colNames) ledger[intersect(colNames, names(ledger))]
  e$outputObjects <- function(sim) data.frame(objectName = c("sppEquiv", "sppNameVector", "sppColorVect",
                                                             "nonForestedLCCGroups", "missingLCCgroup"))
  e
}

oneELFLedger <- function() {
  groups <- list(nfLCC_100_60 = c(100, 60), nfLCC_40_50_80 = c(40, 50, 80))
  pars <- as.data.frame(matrix(0, nrow = 1, ncol = 0))
  for (nm in c(fireSenseUtils::logisticParamNames[[1]], "CMD", "youngAge", names(groups))) pars[[nm]] <- 0
  ledger <- data.frame(numIterations = 1L)
  ledger[[fireSenseUtils::polygonIDTxt]] <- "4.2.2"
  ledger$params <- list(pars)
  ledger$sppEquiv <- list(data.table::data.table(LandR = c("Pice_mar", "Popu_tre"),
                                                 FuelClass = c("Pice_mar", "Popu_tre")))
  ledger$nonForestedLCCGroups <- list(groups)
  ledger$missingLCCgroup <- "nfLCC_40_50_80"
  ledger
}

test_that("with one fitted ELF in the ledger, Init() uses that ELF's non-forest groups for this run's fit", {
  ledger <- oneELFLedger()
  e <- ledgerInitEnv(ledger)
  for (s in initBody()) suppressMessages(eval(s, e))
  expect_true(e$mod$haveSpreadFit)
  expect_identical(e$sim$nonForestedLCCGroups, ledger$nonForestedLCCGroups[[1]])
  expect_identical(e$sim$missingLCCgroup, "nfLCC_40_50_80")
  expect_false("nf" %in% names(e$sim$nonForestedLCCGroups))
})

test_that("without a ledger fit, Init() leaves the non-forest groups alone", {
  e <- ledgerInitEnv(NULL)
  for (s in initBody()) suppressMessages(eval(s, e))
  expect_false(e$mod$haveSpreadFit)
  expect_identical(e$sim$nonForestedLCCGroups, list(nf = c(40, 50, 60, 80, 100)))
  expect_identical(e$sim$missingLCCgroup, "nf")
})
