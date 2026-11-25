## With wholeLandscapeSpreadCovariates the spread-fit covariate tables have every burnable pixel, made by the
## same functions (annualSpreadCovariates(), nonAnnualSpreadCovariates()) called with all burnable pixels
## instead of the buffer pixels. They must have the buffer tables' names and columns, and the same values
## on the buffer pixels.

test_that("whole annual tables have the buffer tables' names and rows on the buffer pixels, and every burnable pixel", {
  skip_if_not_installed("fireSenseUtils")
  withr::local_options(reproducible.cachePath = withr::local_tempdir(), reproducible.useMemoise = FALSE)
  f <- wholeFixture()
  run <- function(pixels) annualSpreadCovariates(pixels, f$clim, f$allYears, f$burnablePixels,
                                                 tsds = f$tsds, firePixelsByYear = f$fires, cutoffForYoungAge = 15)
  buffer <- do.call(c, unname(run(f$buffers)))
  whole <- do.call(c, unname(run(lapply(f$burnablePixels, function(p) data.table::data.table(pixelID = p))))) 

  expect_identical(names(whole), names(buffer))
  expect_identical(names(whole), unname(unlist(f$allYears)))
  for (fy in names(buffer)) {
    expect_identical(names(whole[[fy]]), names(buffer[[fy]]))
    expect_setequal(whole[[fy]]$pixelID, f$burnablePixels[[fy]])
    w <- whole[[fy]][match(buffer[[fy]]$pixelID, whole[[fy]]$pixelID), ]
    expect_equal(w$CMD, buffer[[fy]]$CMD)
    expect_equal(w$youngAge, buffer[[fy]]$youngAge)
  }
  expect_identical(nrow(buffer$year2002), 0L)
  expect_false(17L %in% buffer$year2003$pixelID)   # a non-burnable buffer pixel is dropped from both
  expect_identical(whole$year2001$youngAge[whole$year2001$pixelID == 9L], 0L)
  expect_identical(whole$year2004$youngAge[whole$year2004$pixelID == 9L], 1L)  # burned in 2003
})

test_that("whole non-annual tables have the buffer tables' names and rows on the buffer pixels", {
  veg <- data.table::data.table(pixelID = rep(1:6, 2), year = rep(c(2000L, 2002L), each = 6),
                                fuel = c(1, 2, NA, 4, 5, 6, 7, 8, 9, 10, 11, 12), burned = 0L)
  cols <- c("pixelID", "fuel")
  inBuffer <- veg[pixelID %in% c(2L, 3L, 5L)]
  buffer <- Map(dy = c("2000", "2002"), function(dy) nonAnnualSpreadCovariates(inBuffer[year == as.integer(dy)], cols))
  whole <- Map(dy = c("2000", "2002"), function(dy) nonAnnualSpreadCovariates(veg[year == as.integer(dy)], cols))

  expect_identical(names(whole), names(buffer))
  expect_identical(whole[["2000"]]$pixelID, c(1L, 2L, 4L, 5L, 6L))   # the NA fuel row is dropped from both
  expect_identical(whole[["2002"]]$pixelID, 1:6)
  for (dy in names(buffer)) {
    w <- whole[[dy]][match(buffer[[dy]]$pixelID, whole[[dy]]$pixelID), ]
    expect_identical(w$fuel, buffer[[dy]]$fuel)
  }
})

test_that("prepare_SpreadFit picks the buffer or the burnable pixels by the parameter, into the existing outputs", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_false(md$parameters$default[md$parameters$paramName == "wholeLandscapeSpreadCovariates"][[1]])
  expect_false(any(grepl("Whole$", md$outputObjects$objectName)))
  body <- paste(deparse(body(prepare_SpreadFit)), collapse = " ")
  expect_match(body, "wholeLandscapeSpreadCovariates", fixed = TRUE)
  expect_match(body, "pixelsByYear <- if (wholeSpreadCovariates)", fixed = TRUE)
  expect_equal(lengths(regmatches(body, gregexpr("annualSpreadCovariates(", body, fixed = TRUE))), 1L) # one call, whichever pixels
})

test_that("spreadClimateSplit: the annual tables have CMD_yr (the year's burnable-pixel median) and CMD_dev (CMD minus it)", {
  skip_if_not_installed("fireSenseUtils")
  withr::local_options(reproducible.cachePath = withr::local_tempdir(), reproducible.useMemoise = FALSE)
  f <- wholeFixture()
  run <- function(split) do.call(c, unname(annualSpreadCovariates(f$buffers, f$clim, f$allYears, f$burnablePixels,
                                  tsds = f$tsds, firePixelsByYear = f$fires, cutoffForYoungAge = 15,
                                  climateSplit = split)))
  plain <- run(FALSE)
  split <- run(TRUE)
  expect_identical(names(split), names(plain))
  for (fy in names(plain)) {
    expect_identical(names(plain[[fy]]), c("pixelID", "CMD", "youngAge"))
    expect_true(all(c("pixelID", "CMD_yr", "CMD_dev", "youngAge") %in% names(split[[fy]])))
    expect_false("CMD" %in% names(split[[fy]]))
    expect_identical(split[[fy]]$pixelID, plain[[fy]]$pixelID)
  }
  ## pixels 1:16 burn: their values are y + 1:16, median y + 8.5 (not that of the fitted buffer pixels, which differ)
  for (fy in c("year2001", "year2003", "year2004")) {
    y <- as.integer(sub("year", "", fy))
    expect_equal(split[[fy]]$CMD_yr, rep(y + 8.5, nrow(split[[fy]])))
    expect_equal(split[[fy]]$CMD_dev, plain[[fy]]$CMD - (y + 8.5))
  }
  expect_false(isTRUE(all.equal(split$year2001$CMD_yr[1], median(plain$year2001$CMD))))
  expect_identical(nrow(split$year2002), 0L)
})

test_that("spreadClimateSplit is FALSE by default and the formula is then the same string; TRUE has the two terms", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_false(md$parameters$default[md$parameters$paramName == "spreadClimateSplit"][[1]])

  stmts <- as.list(body(prepare_SpreadFit))[-1]
  txt <- vapply(stmts, function(s) paste(deparse(s), collapse = " "), character(1))
  idx <- sort(unique(unlist(lapply(c("^spreadVars <-", "^vegColsForRHS <-", "^RHS <-"), grep, x = txt))))
  formulaFor <- function(Par, modSplit = NULL) {
    env <- new.env(parent = environment(prepare_SpreadFit))
    env$sim <- new.env(); env$sim$climateVariablesForFire <- list(spread = "CMD")
    env$mod <- new.env(); env$mod$spreadClimateSplit <- modSplit; env$Par <- Par
    env$vegCols <- c("class1", "class2")
    for (i in idx) eval(stmts[[i]], env)
    spreadFormula(env$RHS, Par$spreadIntercept)
  }
  expect_identical(formulaFor(list(spreadIntercept = TRUE)), "~ 1 + CMD + youngAge + class1 + class2")
  expect_identical(formulaFor(list(spreadIntercept = TRUE, spreadClimateSplit = FALSE)), "~ 1 + CMD + youngAge + class1 + class2")
  expect_identical(formulaFor(list(spreadIntercept = TRUE, spreadClimateSplit = TRUE)),
                   "~ 1 + CMD_yr + CMD_dev + youngAge + class1 + class2")
  ## a reused ledger fit's split (mod$spreadClimateSplit, from Init()) wins over the parameter
  expect_identical(formulaFor(list(spreadIntercept = TRUE, spreadClimateSplit = FALSE), modSplit = TRUE),
                   "~ 1 + CMD_yr + CMD_dev + youngAge + class1 + class2")
  expect_identical(formulaFor(list(spreadIntercept = TRUE, spreadClimateSplit = TRUE), modSplit = FALSE),
                   "~ 1 + CMD + youngAge + class1 + class2")
})
