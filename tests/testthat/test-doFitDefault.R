## doFit FALSE, the default since 2026-11-02: a user who adds fireSense only to predict must never start a
## fit by accident. Without a ledger fit for this ELF, Init() stops before any data preparation. This runs
## the module's own Init() with stand-ins for Drive and climate (helper-ledgerInit.R).

defaultDoFit <- function() {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  md$parameters$default[[match("doFit", md$parameters$paramName)]]
}

initDoFit <- function(ledger, elf, doFit) {
  e <- ledgerInitEnv(ledger)
  e$sim$.ELFind <- elf
  e$Par$doFit <- doFit
  e$prepared <- FALSE
  e$addFitClimateVariables <- function(climateVariables, ...) { e$prepared <- TRUE; climateVariables }
  out <- try(for (s in initBody()) suppressMessages(eval(s, e)), silent = TRUE)
  list(e = e, error = if (inherits(out, "try-error")) as.character(out) else "")
}

test_that("the default doFit is FALSE", {
  expect_identical(defaultDoFit(), FALSE)
})

test_that("default doFit, no ledger row for this ELF: Init() stops, naming the ELF and the ledger", {
  r <- initDoFit(oneELFLedger(), elf = "9.9.9", doFit = defaultDoFit())          # the ledger holds 4.2.2 only
  expect_match(r$error, "no usable stored fit for polygon 9\\.9\\.9 in the ledger.*the default.*doFit = NA")
  expect_false(r$e$prepared)                                                     # stopped before the ledger inputs
  ## an empty ledger is no fit either
  r <- initDoFit(NULL, elf = "4.2.2", doFit = defaultDoFit())
  expect_match(r$error, "no usable stored fit for polygon 4\\.2\\.2")
})

test_that("default doFit with this ELF's row, and doFit NA or TRUE without one, do not stop", {
  expect_identical(initDoFit(oneELFLedger(), elf = "4.2.2", doFit = defaultDoFit())$error, "")
  expect_identical(initDoFit(oneELFLedger(), elf = "9.9.9", doFit = NA)$error, "")
  expect_identical(initDoFit(oneELFLedger(), elf = "9.9.9", doFit = TRUE)$error, "")
})
