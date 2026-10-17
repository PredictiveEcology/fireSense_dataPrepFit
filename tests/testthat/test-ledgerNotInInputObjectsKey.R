## `spreadFitPreRun` was an expectsInput, so SpaDES.core put it in the `.inputObjects` cache key
## (simInit(), suppliedElsewhere()) whenever fireSense_ELFs had put the whole ledger in the sim: any
## ELF's new fit redid this ELF's data preparation. Init() always assigns it (this study area's
## rows from CacheGeo, or NULL for a held-out fold) before any cached event runs, so it is an output only
## and the cached events after init key on it through `.useCacheArgs`.

moduleParams <- function() {
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("defineModule")),
                parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE))
  eval(as.list(def[[1]][[3]])[["parameters"]], envir = asNamespace("SpaDES.core"))
}
fnDefs <- function() {
  files <- c(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"),
             list.files(testthat::test_path("..", "..", "R"), "\\.R$", full.names = TRUE))
  defs <- Filter(function(d) is.call(d) && as.character(d[[1]]) %in% c("<-", "=") && is.name(d[[2]]) &&
                   is.call(d[[3]]) && identical(d[[3]][[1]], as.name("function")),
                 unlist(lapply(files, function(f) as.list(parse(f, keep.source = FALSE))), recursive = FALSE))
  stats::setNames(lapply(defs, function(d) d[[3]]), vapply(defs, function(d) as.character(d[[2]]), ""))
}
cacheExtras <- function() {
  p <- moduleParams()
  p$default[[match(".useCacheArgs", p$paramName)]]
}

test_that("spreadFitPreRun is created by the module but not expected as an input", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_false("spreadFitPreRun" %in% md$inputObjects$objectName)
  expect_true("spreadFitPreRun" %in% md$outputObjects$objectName)
})

test_that(".inputObjects does not read spreadFitPreRun", {
  f <- fnDefs()
  expect_false("spreadFitPreRun" %in% all.names(f[[".inputObjects"]]))
})

test_that("Init() assigns spreadFitPreRun unconditionally, before the cached events", {
  txt <- deparse(fnDefs()[["Init"]])
  expect_true(any(grepl("^\\s*sim\\$spreadFitPreRun <- if \\(isTRUE\\(is.na\\(heldOutFold\\)\\)\\)", txt)))
})

test_that("every cached event after init keys on the ledger Init read, and only on it", {
  args <- cacheExtras()
  events <- c("dataPrepBuild", "prepIgnitionFitData", "prepEscapeFitData", "prepSpreadFitData")
  expect_true(all(events %in% names(args)))
  ## the same Init result gives the same key; this ELF's own new row gives another; NULL (held-out fold) another again
  keyOf <- function(extra, ledger) {
    sim <- new.env()
    sim$spreadFitPreRun <- ledger
    digest::digest(eval(extra, list2env(list(sim = sim), parent = asNamespace("fireSenseUtils"))))
  }
  own <- data.frame(polygonID = "5.1.1", objFunVal = 1)
  ownNew <- data.frame(polygonID = "5.1.1", objFunVal = 2)
  for (ev in events) {
    extra <- args[[ev]]$.cacheExtra
    expect_true(is.call(extra), info = ev)
    expect_identical(keyOf(extra, own), keyOf(extra, own), info = ev)
    expect_false(identical(keyOf(extra, own), keyOf(extra, ownNew)), info = ev)
    expect_false(identical(keyOf(extra, own), keyOf(extra, NULL)), info = ev)
  }
})

test_that("a ledger differing only in another ELF's rows gives the same Init result, hence the same key", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  sq <- function(x0) sf::st_sfc(sf::st_polygon(list(rbind(c(x0, 0), c(x0 + 1e4, 0), c(x0 + 1e4, 1e4),
                                                          c(x0, 1e4), c(x0, 0)))), crs = 3978)
  mk <- function(otherValue) {
    df <- data.frame(polygonID = c("4.1", "4.3"), objFunVal = c(otherValue, 2),
                     crs = I(rep(sf::st_crs(3978)$wkt, 2)))
    df$geometry <- c(sq(2e4), sq(4e4))
    df
  }
  value <- 1
  testthat::local_mocked_bindings(latestSpreadFits = function(...) mk(value), .package = "fireSenseUtils")
  domain <- sf::st_as_sf(sq(4e4))                            # 4.3
  first <- suppressMessages(readSpreadFitLedger("latest", "folder", domain, d))
  value <- 99                                                # ELF 4.1 refitted
  second <- suppressMessages(readSpreadFitLedger("latest", "folder", domain, d))
  expect_identical(digest::digest(first), digest::digest(second))
})
