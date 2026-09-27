## fireSense_dataPrepFit.R ~774 (prepare_SpreadFit()): the spread formula's RHS was built as
## `climate + youngAgeTxt + vegCols`, but vegCols already includes a "youngAge" column whenever
## fireSenseUtils::fireSenseCovariatesCreate() finds young forest/non-forest pixels (real cache
## entries for ELFs 14.3 and 14.4, 2026-09-27). R's terms() silently drops the duplicate, so the
## resulting formula has one fewer distinct term than the raw covariate list that downstream code
## (the fitted parameter bounds) is built from.

test_that("the spread formula lists a veg-table youngAge column once, not twice", {
  body <- body(prepare_SpreadFit)
  stmts <- as.list(body)[-1]
  txt <- vapply(stmts, function(s) paste(deparse(s), collapse = " "), character(1))

  ## pull just the statements that build vegCols and the formula's RHS out of prepare_SpreadFit,
  ## in their original order, and eval them against a small fabricated fireSenseVegData -- the
  ## rest of that function (fire buffering, climate rasters, ...) is irrelevant to this defect
  patterns <- c("^nonVegColnames <-", "^vegCols <- setdiff", "^dropCols <-",
                "^if \\(length\\(dropCols\\)", "^vegColsForRHS <-", "^RHS <-",
                "^if \\(is\\.null\\(sim\\$fireSense_spreadFormula\\)\\)")
  idx <- sort(unique(unlist(lapply(patterns, grep, x = txt))))
  expect_true(length(idx) >= 5) # sanity: these statements are still present in prepare_SpreadFit

  env <- new.env(parent = environment(prepare_SpreadFit))
  env$sim <- new.env()
  env$sim$climateVariablesForFire <- list(spread = "CMD")
  env$sim$fireSense_spreadFormula <- NULL
  ## the shape fireSenseVegData has after joinFireBuffersToVeg() and setnames(..., "buffer",
  ## "burned"): a youngAge column that is not all-zero, as in ELFs 14.3/14.4
  env$fireSenseVegData <- data.table::data.table(
    pixelID = 1:10, ids = 1L, burned = rep(0:1, 5), year = 2010L,
    fuelA = c(rep(0, 5), rep(1, 5)), youngAge = c(rep(1, 3), rep(0, 7))
  )

  for (i in idx) eval(stmts[[i]], env)

  rhsTerms <- trimws(strsplit(env$RHS, "\\+")[[1]])
  expect_identical(anyDuplicated(rhsTerms), 0L)
  expect_identical(length(rhsTerms), length(env$vegCols) + 1L) # +1 for the climate variable
  expect_identical(env$sim$fireSense_spreadFormula, "~ 0 + CMD + youngAge + fuelA")
})
