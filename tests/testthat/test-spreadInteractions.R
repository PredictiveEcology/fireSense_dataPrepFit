## `spreadInteractions` (default character(0)): "climate:fuel" adds one interaction term for each spread climate
## variable and each fuel covariate (not youngAge, nfLCC_* or treedWetland), after the main effects, for
## fireSense_spreadFit to fit with their own coefficients. Off, the formula string is byte-identical to what it was.

test_that("spreadRHS() without interactions is the string the module always built", {
  fuel <- c("dom_agb_Pice_mar", "dom_agb_Pinu_con", "nfLCC_1_2")
  old <- paste(paste0("CMD", collapse = " + "), "youngAge", paste0(fuel, collapse = " + "), sep = " + ")
  expect_identical(spreadRHS("all", "CMD", fuel), old)
  expect_identical(spreadRHS("all", "CMD", fuel, character(0)), old)
  expect_identical(spreadRHS("all", "CMD", fuel, NULL), old)
  expect_identical(spreadRHS("none", "CMD", fuel), character())
  expect_identical(spreadFormula(spreadRHS("all", "CMD", fuel), TRUE), paste0("~ 1 + ", old))
  expect_identical(spreadFormula(spreadRHS("all", "CMD", fuel)), paste0("~ 0 + ", old))
})

test_that("spreadRHS() with 'climate:fuel' adds climate:fuel for each fuel covariate, after the main effects", {
  fuel <- c("dom_agb_Pice_mar", "dom_agb_Pinu_con", "nfLCC_1_2", "treedWetland", "treedWetland_agb")
  expect_identical(
    spreadRHS("all", "CMD", fuel, "climate:fuel"),
    paste("CMD + youngAge", paste(fuel, collapse = " + "),
          "CMD:dom_agb_Pice_mar + CMD:dom_agb_Pinu_con + CMD:treedWetland_agb", sep = " + "))
  ## each climate variable
  expect_identical(
    spreadRHS("all", c("CMD", "Tmax"), "fuelA", "climate:fuel"),
    "CMD + Tmax + youngAge + fuelA + CMD:fuelA + Tmax:fuelA")
  ## the terms are design columns of the interaction, after the main effects
  f <- spreadFormula(spreadRHS("all", "CMD", c("fuelA", "fuelB"), "climate:fuel"), TRUE)
  expect_identical(f, "~ 1 + CMD + youngAge + fuelA + fuelB + CMD:fuelA + CMD:fuelB")
  expect_identical(fireSenseUtils::spreadDesignCols(f),
                   c(fireSenseUtils::spreadInterceptTxt, "CMD", "youngAge", "fuelA", "fuelB", "CMD:fuelA", "CMD:fuelB"))
  ## no fuel, no interaction
  expect_identical(spreadRHS("all", "CMD", "nfLCC_1_2", "climate:fuel"), "CMD + youngAge + nfLCC_1_2")
})

test_that("checkSpreadInteractions() allows character(0) and 'climate:fuel', and not with the null model", {
  expect_null(checkSpreadInteractions(character(0), "all"))
  expect_null(checkSpreadInteractions(NULL, "none"))
  expect_null(checkSpreadInteractions("climate:fuel", "all"))
  expect_error(checkSpreadInteractions("fuel:fuel", "all"), "must be")
  expect_error(checkSpreadInteractions(1, "all"), "must be")
  expect_error(checkSpreadInteractions("climate:fuel", "none"), "needs covariates")
})

test_that("spreadRHS() with a split climate variable interacts the between-year term (_yr) only, never _dev", {
  expect_identical(
    spreadRHS("all", "CMD", c("fuelA", "fuelB"), "climate:fuel", climateSplit = TRUE),
    "CMD_yr + CMD_dev + youngAge + fuelA + fuelB + CMD_yr:fuelA + CMD_yr:fuelB")
  expect_identical(
    spreadRHS("all", c("CMD", "Tmax"), "fuelA", "climate:fuel", climateSplit = TRUE),
    "CMD_yr + CMD_dev + Tmax_yr + Tmax_dev + youngAge + fuelA + CMD_yr:fuelA + Tmax_yr:fuelA")
  rhs <- spreadRHS("all", "CMD", c("fuelA", "fuelB"), "climate:fuel", climateSplit = TRUE)
  expect_false(grepl("_dev:|:[A-Za-z]*_dev", rhs))
  ## the design columns: the interaction's source columns include CMD_yr and the fuels, not CMD
  cols <- fireSenseUtils::spreadDesignCols(spreadFormula(rhs, TRUE))
  expect_identical(fireSenseUtils::spreadInteractionCols(cols), c("CMD_yr:fuelA", "CMD_yr:fuelB"))
  expect_false("CMD" %in% fireSenseUtils::spreadSourceCols(cols))
  ## a split without interactions, and interactions without a split, are as before
  expect_identical(spreadRHS("all", "CMD", "fuelA", character(0), climateSplit = TRUE), "CMD_yr + CMD_dev + youngAge + fuelA")
  expect_identical(spreadRHS("all", "CMD", "fuelA", "climate:fuel", climateSplit = FALSE), "CMD + youngAge + fuelA + CMD:fuelA")
})

test_that("prepare_SpreadFit() builds _yr:fuel interactions for a split fit, whatever spreadClimateSplit says for a reused fit", {
  stmts <- as.list(body(prepare_SpreadFit))[-1]
  txt <- vapply(stmts, function(s) paste(deparse(s), collapse = " "), character(1))
  idx <- sort(unique(unlist(lapply(c("^spreadVars <-", "^vegColsForRHS <-", "^checkSpreadInteractions", "^RHS <-"), grep, x = txt))))
  build <- function(modSplit, parSplit) {
    env <- new.env(parent = environment(prepare_SpreadFit))
    env$sim <- new.env(); env$sim$climateVariablesForFire <- list(spread = "CMD")
    env$mod <- new.env(); env$mod$spreadClimateSplit <- modSplit
    env$Par <- list(spreadIntercept = TRUE, spreadCovariates = "all", spreadInteractions = "climate:fuel",
                    spreadClimateSplit = parSplit)
    env$vegCols <- c("class1", "class2")
    for (i in idx) eval(stmts[[i]], env)
    env$RHS
  }
  want <- "CMD_yr + CMD_dev + youngAge + class1 + class2 + CMD_yr:class1 + CMD_yr:class2"
  expect_identical(build(modSplit = TRUE, parSplit = FALSE), want)
  expect_identical(build(modSplit = NULL, parSplit = TRUE), want)
  expect_match(build(modSplit = FALSE, parSplit = TRUE), "CMD:class1", fixed = TRUE)
})

test_that("the module defines spreadInteractions, character(0) by default, and Init() checks it", {
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("defineModule")),
                parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE))
  params <- eval(as.list(def[[1]][[3]])[["parameters"]], envir = asNamespace("SpaDES.core"))
  i <- match("spreadInteractions", params$paramName)
  expect_false(is.na(i))
  expect_identical(params$paramClass[[i]], "character")
  expect_identical(params$default[[i]], character(0))
  expect_match(paste(deparse(body(Init)), collapse = " "),
               "checkSpreadInteractions\\(Par\\$spreadInteractions, Par\\$spreadCovariates\\)")
})

test_that("prepare_SpreadFit() gives the formula it always gave, and with 'climate:fuel' the interactions of the fuels it keeps", {
  extra <- list(fuelB = fireSenseUtils::logMinB(c(rep(300, 5), rep(0, 5))), nfLCC_1_2 = rep(0:1, 5),
                fuelEmpty = fireSenseUtils::logMinB(rep(0, 10)))
  ## fuelEmpty is dropped (no biomass anywhere), so it has no interaction
  base <- formulaFor("all", TRUE, extraCols = extra)
  expect_identical(base, "~ 1 + CMD + youngAge + fuelA + fuelB + nfLCC_1_2")
  expect_identical(formulaFor("all", TRUE, character(0), extra), base)
  expect_identical(formulaFor("all", TRUE, "climate:fuel", extra),
                   "~ 1 + CMD + youngAge + fuelA + fuelB + nfLCC_1_2 + CMD:fuelA + CMD:fuelB")
  expect_identical(formulaFor("all", FALSE, "climate:fuel"), "~ 0 + CMD + youngAge + fuelA + CMD:fuelA")
})
