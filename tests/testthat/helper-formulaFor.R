## prepare_SpreadFit()'s spread formula, from the statements that build it, for a toy fit table (fuelA,
## youngAge, and any `extraCols`, a named list of columns) and one spread climate variable, CMD
formulaFor <- function(spreadCovariates, intercept, spreadInteractions = NULL, extraCols = list()) {
  stmts <- as.list(body(prepare_SpreadFit))[-1]
  txt <- vapply(stmts, function(s) paste(deparse(s), collapse = " "), character(1))
  idx <- sort(unique(unlist(lapply(c("^spreadVars <-", "^nonVegColnames <-", "^vegCols <- setdiff", "^dropCols <-",
                                     "^if \\(length\\(dropCols\\)", "^vegColsForRHS <-", "^RHS <-",
                                     "^if \\(is\\.null\\(sim\\$fireSense_spreadFormula\\)\\)"), grep, x = txt))))
  env <- new.env(parent = environment(prepare_SpreadFit))
  env$sim <- new.env(); env$sim$climateVariablesForFire <- list(spread = "CMD")
  env$sim$fireSense_spreadFormula <- NULL
  env$mod <- new.env()
  env$Par <- list(spreadIntercept = intercept, spreadCovariates = spreadCovariates,
                  spreadInteractions = spreadInteractions)
  env$fireSenseVegData <- data.table::data.table(
    pixelID = 1:10, ids = 1L, burned = rep(0:1, 5), year = 2010L,
    fuelA = fireSenseUtils::logMinB(c(rep(0, 5), rep(500, 5))), youngAge = c(rep(1, 3), rep(0, 7)))
  for (nm in names(extraCols)) data.table::set(env$fireSenseVegData, NULL, nm, extraCols[[nm]])
  for (i in idx) eval(stmts[[i]], env)
  env$sim$fireSense_spreadFormula
}
