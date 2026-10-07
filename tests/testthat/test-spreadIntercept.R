## `spreadIntercept` (default FALSE): the spread formula has an intercept, `~ 1 + RHS`, which
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

test_that("the module defines spreadIntercept, a logical that is FALSE by default", {
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("defineModule")),
                parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE))
  params <- eval(as.list(def[[1]][[3]])[["parameters"]], envir = asNamespace("SpaDES.core"))
  i <- match("spreadIntercept", params$paramName)
  expect_false(is.na(i))
  expect_identical(params$paramClass[[i]], "logical")
  expect_identical(params$default[[i]], FALSE)
})

test_that("prepare_SpreadFit() builds the formula with spreadFormula() from the spreadIntercept parameter", {
  src <- paste(deparse(body(prepare_SpreadFit)), collapse = "\n")
  expect_match(src, "spreadFormula\\(RHS, Par\\$spreadIntercept\\)")
  expect_no_match(src, '"~ 0 \\+ "', fixed = FALSE)
  ## the formula is only built when none was supplied
  expect_match(src, "is.null\\(sim\\$fireSense_spreadFormula\\)")
})
