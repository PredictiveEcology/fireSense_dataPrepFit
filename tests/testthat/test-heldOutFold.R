## `heldOutFold` (the same parameter as fireSense_spreadFit's): a held-out fold is fitted without the
## SpreadFit ledger, so Init() must not read it. Init() needs a full simList (LandR, Google Drive), so this
## evaluates the statements of the module's own Init() up to `mod$haveSpreadFit <- ...` with stand-ins.

initLedgerStatements <- function() {
  exprs <- parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE)
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                  identical(x[[2]], as.name("Init")), exprs)
  stopifnot(length(def) == 1L)
  bod <- as.list(body(eval(def[[1]][[3]])))[-1]
  isEnd <- vapply(bod, function(x) any(grepl("mod$haveSpreadFit <-", deparse(x)[1], fixed = TRUE)), logical(1))
  bod[seq_len(which(isEnd)[1])]
}

## a simList with this module current and, optionally, another module carrying its own heldOutFold
initEnv <- function(ownFold, otherFold = NULL, ledgerCalls = NULL) {
  sim <- suppressMessages(SpaDES.core::simInit())
  sim@params <- list(fireSense_dataPrepFit = list(heldOutFold = ownFold))
  if (!is.null(otherFold)) sim@params$fireSense_spreadFit <- list(heldOutFold = otherFold)
  sim@current <- data.table::data.table(moduleName = "fireSense_dataPrepFit")
  e <- new.env()
  e$sim <- sim
  e$mod <- new.env()
  e$Par <- list(heldOutFold = ownFold, spreadFitFilename = "latest", spreadFitGoogleDriveFolder = "folder",
                escapeSizeHa = 1)
  e$nonEscapedFireSizes <- function(...) NULL
  e$inputPath <- function(sim) tempdir()
  e$ledgerCalls <- 0L
  e$readSpreadFitLedger <- function(...) {
    e$ledgerCalls <- e$ledgerCalls + 1L
    NULL
  }
  e
}

runInitStatements <- function(e) {
  for (s in initLedgerStatements()) eval(s, e)
  invisible(e)
}

test_that("with heldOutFold 1 or 2 the SpreadFit ledger is not read and the fold is treated as unfitted", {
  for (fold in 1:2) {
    e <- initEnv(fold, otherFold = fold)
    suppressMessages(runInitStatements(e))
    expect_identical(e$ledgerCalls, 0L)
    expect_null(e$sim$spreadFitPreRun)
    expect_false(e$mod$haveSpreadFit)
  }
})

test_that("with heldOutFold NA the SpreadFit ledger is read as before", {
  e <- initEnv(NA, otherFold = NA)
  suppressMessages(runInitStatements(e))
  expect_identical(e$ledgerCalls, 1L)
})

test_that("a heldOutFold that is not NA, 1 or 2 is an error", {
  e <- initEnv(3L, otherFold = 3L)
  expect_error(suppressMessages(runInitStatements(e)), "heldOutFold")
  expect_identical(e$ledgerCalls, 0L)
})

test_that("stops when fireSense_spreadFit has a different heldOutFold", {
  e <- initEnv(NA, otherFold = 1L)
  expect_error(suppressMessages(runInitStatements(e)), "multiple values for heldOutFold")
  expect_identical(e$ledgerCalls, 0L)
})
