## Reading the fitted-parameter ledger when rows predate a column.
##
## The ledger is shared cloud state that gains columns over time (covMinMax_spread, 2026-09-17). Init()
## selected the columns it reads with `df[colNames]`, which errors when any row set lacks one -- so adding a
## column to fireSenseUtils::spreadFitAdditionalColNamesTxt would have stopped every run until every row
## was rewritten.

test_that("ledgerColumns() keeps the requested columns that exist and ignores the missing ones", {
  old <- data.frame(numIterations = I(list(351L)), objFunVal = I(list(1:5)), polygonID = "6.2.1")
  wanted <- c("numIterations", "objFunVal", "covMinMax_spread")
  expect_identical(names(ledgerColumns(old, wanted)), c("numIterations", "objFunVal"))
  expect_no_error(ledgerColumns(old, wanted))
})

test_that("ledgerColumns() keeps an sf ledger's geometry, as the plain selection did", {
  skip_if_not_installed("sf")
  pt <- sf::st_sfc(sf::st_point(c(0, 0)), crs = 3978)
  old <- sf::st_sf(numIterations = 351L, polygonID = "6.2.1", geometry = pt)
  out <- ledgerColumns(old, c("numIterations", "covMinMax_spread"))
  expect_s3_class(out, "sf")
  expect_true("numIterations" %in% names(out))
  expect_false("covMinMax_spread" %in% names(out))
})

test_that("Init() reads the ledger through ledgerColumns()", {
  mainFiles <- list.files(moduleRoot, pattern = "\\.R$", full.names = TRUE)
  mainFile <- mainFiles[vapply(mainFiles, function(f)
    any(grepl("defineModule\\(", readLines(f, warn = FALSE))), logical(1))]
  src <- readLines(mainFile, warn = FALSE)
  expect_true(any(grepl("ledgerColumns(sim$spreadFitPreRun, colNames)", src, fixed = TRUE)))
})
