## With a SpreadFit in the ledger, the ignition and escape covariates of this run are still built here, from
## `sim$nonForestedLCCGroups` and `sim$missingLCCgroup` (dataPrepBuild). fireSense_dataPrepPredict builds the
## prediction covariates from the ledger's groups, so the fit must use them too. Before, they stayed the
## module default (`nf`), and fireSense_ignitionPredict stopped with "column not found: [nf]" (carbon run,
## ELF 4.2.2, 2026-09-29). This runs the module's own Init() with stand-ins for Drive, LandR and climate.

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

test_that("with own ELF and a neighbour in the ledger, Init() builds the single-ELF objects from the own row", {
  own <- oneELFLedger()
  nb <- oneELFLedger()
  nb[[fireSenseUtils::polygonIDTxt]] <- "4.2.1"
  nb$sppEquiv <- list(data.table::data.table(LandR = c("Abie_bal", "Pice_mar"),
                                             FuelClass = c("Abie_bal", "Pice_mar")))
  nb$nonForestedLCCGroups <- list(list(nfLCC_100 = 100, nfLCC_40 = 40))
  nb$missingLCCgroup <- "nfLCC_40"
  ledger <- rbind(nb, own)
  ledger$params <- list(own$params[[1]], own$params[[1]])
  ledger$sppEquiv <- list(nb$sppEquiv[[1]], own$sppEquiv[[1]])
  ledger$nonForestedLCCGroups <- list(nb$nonForestedLCCGroups[[1]], own$nonForestedLCCGroups[[1]])
  e <- ledgerInitEnv(ledger)
  e$sim$.ELFind <- "4.2.2"
  for (s in initBody()) suppressMessages(eval(s, e))
  expect_identical(e$sim$nonForestedLCCGroups, own$nonForestedLCCGroups[[1]])
  expect_identical(e$sim$missingLCCgroup, "nfLCC_40_50_80")
  expect_identical(e$sim$sppNameVector, c("Pice_mar", "Popu_tre"))
  expect_setequal(e$sim$sppEquiv$LandR, c("Pice_mar", "Popu_tre"))
  expect_length(e$sim$nonForestedLCCGroupsList, 2L)
  expect_length(e$sim$sppEquivs, 2L)
})
