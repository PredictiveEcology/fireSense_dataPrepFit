## ELF 6.2.2 (2026-10-21): its own ledger fit used CMD, its neighbours 6.3.1 and 6.3.2 CMDsm. Init() set
## climateVariablesForFire$spread to the union of every ledger row's variables, so `spread` was no longer
## "auto", no variable was chosen, and the spread formula had both CMD and CMDsm, two near-collinear terms.
## The union is still prepared (prediction blends neighbour ELFs); the formula takes this ELF's own.

ownAndNeighbourLedger <- function() {
  own <- oneELFLedger()                     # 4.2.2, fitted with CMD
  nb <- oneELFLedger()
  nb[[fireSenseUtils::polygonIDTxt]] <- "4.2.1"
  nbPars <- own$params[[1]]
  names(nbPars)[names(nbPars) == "CMD"] <- "CMDsm"
  ledger <- rbind(nb, own)
  ledger$params <- list(nbPars, own$params[[1]])
  ledger
}

initWithLedger <- function(doFit) {
  e <- ledgerInitEnv(ownAndNeighbourLedger())
  e$sim$.ELFind <- "4.2.2"
  e$sim$climateVariablesForFire <- list(ignition = c("CMD", "cumMDC", "CMDsm", "CMDsp"), spread = "auto")
  e$Par$doFit <- doFit
  for (s in initBody()) suppressMessages(eval(s, e))
  e
}

test_that("a reused ledger fit keeps its own climate variable; prediction still gets every ELF's", {
  e <- initWithLedger(doFit = NA)
  expect_identical(e$mod$spreadFormulaClimate, "CMD")
  expect_setequal(e$sim$climateVariablesForFire$spread, c("CMD", "CMDsm"))
})

test_that("a refit chooses its climate variable afresh, whatever the ledger rows hold", {
  e <- initWithLedger(doFit = TRUE)
  expect_identical(e$mod$spreadFormulaClimate, "auto")
  expect_setequal(e$sim$climateVariablesForFire$spread, c("CMD", "CMDsm"))
})

test_that("the spread formula has this ELF's climate variable only, not the ledger union", {
  stmts <- as.list(body(prepare_SpreadFit))[-1]
  txt <- vapply(stmts, function(s) paste(deparse(s), collapse = " "), character(1))
  patterns <- c("^spreadVars <-", "^sim\\$climateVariablesForFire\\$spread <- union", "^nonVegColnames <-",
                "^vegCols <- setdiff", "^dropCols <-", "^if \\(length\\(dropCols\\)", "^vegColsForRHS <-",
                "^RHS <-", "^if \\(is\\.null\\(sim\\$fireSense_spreadFormula\\)\\)")
  idx <- sort(unique(unlist(lapply(patterns, grep, x = txt))))
  expect_length(idx, length(patterns))

  env <- new.env(parent = environment(prepare_SpreadFit))
  env$sim <- new.env()
  env$sim$climateVariablesForFire <- list(spread = c("CMD", "CMDsm")) # Init()'s union
  env$sim$fireSense_spreadFormula <- NULL
  env$mod <- new.env()
  env$mod$spreadFormulaClimate <- "CMD"                               # this ELF's ledger fit
  env$Par <- list(spreadIntercept = FALSE)
  env$fireSenseVegData <- data.table::data.table(pixelID = 1:4, ids = 1L, burned = c(0L, 1L, 0L, 1L),
                                                 year = 2010L, fuelA = 1, youngAge = c(1, 0, 0, 0))
  for (i in idx) eval(stmts[[i]], env)

  expect_identical(env$sim$fireSense_spreadFormula, "~ 0 + CMD + youngAge + fuelA")
  expect_setequal(env$sim$climateVariablesForFire$spread, c("CMD", "CMDsm"))
})

## spreadClimateSplit: a reused ledger fit keeps its own split, whatever the parameter says
splitLedger <- function() {
  ledger <- ownAndNeighbourLedger()
  own <- ledger$params[[2]]
  names(own)[names(own) == "CMD"] <- "CMD_yr"
  own$CMD_dev <- 0
  ledger$params[[2]] <- own
  ledger
}
initSplit <- function(ledger, doFit, parSplit) {
  e <- ledgerInitEnv(ledger)
  e$sim$.ELFind <- "4.2.2"
  e$sim$climateVariablesForFire <- list(ignition = c("CMD", "cumMDC", "CMDsm", "CMDsp"), spread = "auto")
  e$Par$doFit <- doFit
  e$Par$spreadClimateSplit <- parSplit
  for (s in initBody()) suppressMessages(eval(s, e))
  e
}

test_that("a reused fit with split terms rebuilds the split formula even when spreadClimateSplit is FALSE", {
  e <- initSplit(splitLedger(), doFit = NA, parSplit = FALSE)
  expect_identical(e$mod$spreadFormulaClimate, "CMD")        # the variable, found from CMD_yr and CMD_dev
  expect_true(e$mod$spreadClimateSplit)
  expect_true(useClimateSplit(e$mod, e$Par))
})

test_that("a reused fit without split terms stays unsplit, also when spreadClimateSplit is TRUE", {
  e <- initSplit(ownAndNeighbourLedger(), doFit = NA, parSplit = TRUE)
  expect_false(e$mod$spreadClimateSplit)
  expect_false(useClimateSplit(e$mod, e$Par))
})

test_that("a refit follows the parameter", {
  e <- initSplit(splitLedger(), doFit = TRUE, parSplit = FALSE)
  expect_null(e$mod$spreadClimateSplit)
  expect_false(useClimateSplit(e$mod, e$Par))
  expect_true(useClimateSplit(e$mod, list(spreadClimateSplit = TRUE)))
})
