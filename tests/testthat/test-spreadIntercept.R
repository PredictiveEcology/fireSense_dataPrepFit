## `spreadIntercept` (default TRUE): the spread formula has an intercept, `~ 1 + RHS`, which
## fireSense_spreadFit fits with centred covariates. Off, the formula string is byte-identical to what it
## has always been, so every cache key that holds it is unchanged.

test_that("spreadFormula() gives '~ 0 + RHS' by default and '~ 1 + RHS' with an intercept", {
  rhs <- "CMD + youngAge + fuelA"
  expect_identical(spreadFormula(rhs), "~ 0 + CMD + youngAge + fuelA")
  expect_identical(spreadFormula(rhs, intercept = FALSE), "~ 0 + CMD + youngAge + fuelA")
  expect_identical(spreadFormula(rhs, intercept = TRUE), "~ 1 + CMD + youngAge + fuelA")
  ## the default string is the one the module built before the parameter existed
  expect_identical(spreadFormula(rhs), paste0("~ 0 + ", rhs))
  ## and the intercept is what terms() sees
  expect_equal(attr(terms(as.formula(spreadFormula(rhs, TRUE))), "intercept"), 1L)
  expect_equal(attr(terms(as.formula(spreadFormula(rhs))), "intercept"), 0L)
  expect_identical(fireSenseUtils::spreadDesignCols(spreadFormula(rhs, TRUE)),
                   c(fireSenseUtils::spreadInterceptTxt, "CMD", "youngAge", "fuelA"))
})

test_that("the module defines spreadIntercept, a logical that is TRUE by default", {
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("defineModule")),
                parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE))
  params <- eval(as.list(def[[1]][[3]])[["parameters"]], envir = asNamespace("SpaDES.core"))
  i <- match("spreadIntercept", params$paramName)
  expect_false(is.na(i))
  expect_identical(params$paramClass[[i]], "logical")
  expect_identical(params$default[[i]], TRUE)
})

test_that("prepare_SpreadFit() builds the formula with spreadFormula() from the spreadIntercept parameter", {
  src <- paste(deparse(body(prepare_SpreadFit)), collapse = "\n")
  expect_match(src, "spreadFormula\\(RHS, Par\\$spreadIntercept\\)")
  expect_no_match(src, '"~ 0 \\+ "', fixed = FALSE)
  ## the formula is only built when none was supplied
  expect_match(src, "is.null\\(sim\\$fireSense_spreadFormula\\)")
})

## `spreadCovariates = "none"`: the null spread model, the formula `~ 1`
test_that("an empty RHS is the formula '~ 1', and needs an intercept", {
  expect_identical(spreadFormula(character(), intercept = TRUE), "~ 1")
  expect_error(spreadFormula(character(), intercept = FALSE))
  expect_identical(fireSenseUtils::spreadDesignCols(spreadFormula(character(), TRUE)), fireSenseUtils::spreadInterceptTxt)
  ## "all" is what it was
  expect_identical(spreadFormula("CMD + youngAge"), "~ 0 + CMD + youngAge")
})

test_that("checkSpreadCovariates() allows 'all' and 'none', and 'none' only with an intercept", {
  expect_null(checkSpreadCovariates("all", FALSE))
  expect_null(checkSpreadCovariates("all", TRUE))
  expect_null(checkSpreadCovariates("none", TRUE))
  expect_error(checkSpreadCovariates("none", FALSE), "needs an intercept")
  expect_error(checkSpreadCovariates("some", TRUE), "must be")
  expect_error(checkSpreadCovariates(NULL, TRUE), "must be")
})

test_that("the module defines spreadCovariates, 'all' by default, and Init() checks it", {
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("defineModule")),
                parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE))
  params <- eval(as.list(def[[1]][[3]])[["parameters"]], envir = asNamespace("SpaDES.core"))
  i <- match("spreadCovariates", params$paramName)
  expect_false(is.na(i))
  expect_identical(params$paramClass[[i]], "character")
  expect_identical(params$default[[i]], "all")
  expect_match(paste(deparse(body(Init)), collapse = "\n"),
               "checkSpreadCovariates\\(Par\\$spreadCovariates, Par\\$spreadIntercept\\)")
})

test_that("prepare_SpreadFit() with spreadCovariates 'none' sets the formula '~ 1', and 'all' the one it always set", {
  stmts <- as.list(body(prepare_SpreadFit))[-1]
  txt <- vapply(stmts, function(s) paste(deparse(s), collapse = " "), character(1))
  idx <- sort(unique(unlist(lapply(c("^spreadVars <-", "^nonVegColnames <-", "^vegCols <- setdiff", "^dropCols <-",
                                     "^if \\(length\\(dropCols\\)", "^vegColsForRHS <-", "^RHS <-",
                                     "^if \\(is\\.null\\(sim\\$fireSense_spreadFormula\\)\\)"), grep, x = txt))))
  formulaFor <- function(spreadCovariates, intercept) {
    env <- new.env(parent = environment(prepare_SpreadFit))
    env$sim <- new.env(); env$sim$climateVariablesForFire <- list(spread = "CMD")
    env$sim$fireSense_spreadFormula <- NULL
    env$mod <- new.env()
    env$Par <- list(spreadIntercept = intercept, spreadCovariates = spreadCovariates)
    env$fireSenseVegData <- data.table::data.table(
      pixelID = 1:10, ids = 1L, burned = rep(0:1, 5), year = 2010L,
      fuelA = fireSenseUtils::logMinB(c(rep(0, 5), rep(500, 5))), youngAge = c(rep(1, 3), rep(0, 7)))
    for (i in idx) eval(stmts[[i]], env)
    env$sim$fireSense_spreadFormula
  }
  expect_identical(formulaFor("none", TRUE), "~ 1")
  expect_identical(formulaFor("all", FALSE), "~ 0 + CMD + youngAge + fuelA") # byte-identical, as before
  expect_identical(formulaFor("all", TRUE), "~ 1 + CMD + youngAge + fuelA")
})
