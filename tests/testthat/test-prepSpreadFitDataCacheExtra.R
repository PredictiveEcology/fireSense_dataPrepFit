## A cached event's key covers this module's code, not the package functions it calls. The run
## scripts cache prepSpreadFitData, so a fix in fireSenseUtils (e.g. cleanUpSpreadFirePoints(),
## makeMutuallyExclusive()) came back as the old spread-fit data. The event's `.cacheExtra` must name
## every LandR and fireSenseUtils function the event reaches, through the module functions it calls.

## package functions reached from `start`, through the module's own functions
packageCallsFrom <- function(start, pkgs = c("LandR", "fireSenseUtils")) {
  defs <- moduleDefs()
  out <- unlist(lapply(moduleFunsFrom(start, defs), function(f) {
    txt <- paste(deparse(defs[[f]]), collapse = " ")
    qualified <- regmatches(txt, gregexpr("(LandR|fireSenseUtils):::?[A-Za-z0-9_.]+", txt))[[1]]
    nms <- all.names(defs[[f]])
    bare <- unlist(lapply(pkgs, function(p) {
      hits <- intersect(nms, getNamespaceExports(p))
      if (length(hits)) paste0(p, "::", hits)
    }))
    c(sub(":::", "::", qualified), bare)
  }))
  out <- unique(out)
  out[vapply(out, function(f) is.function(eval(str2lang(f))), logical(1))] # not data
}

test_that("prepSpreadFitData's cache key includes the package functions it reaches", {
  skip_if_not_installed("LandR")
  skip_if_not_installed("fireSenseUtils")
  ## the module's metadata, evaluated from source, so the checkout directory's name does not matter
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("defineModule")),
                parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE))
  params <- eval(as.list(def[[1]][[3]])[["parameters"]], envir = asNamespace("SpaDES.core"))
  extra <- params$default[[match(".useCacheArgs", params$paramName)]]$prepSpreadFitData$.cacheExtra
  expect_true(is.call(extra))
  txt <- paste(deparse(extra), collapse = " ")
  keyed <- regmatches(txt, gregexpr("(LandR|fireSenseUtils)::[A-Za-z0-9_.]+", txt))[[1]]

  called <- packageCallsFrom("prepare_SpreadFit")
  expect_true("fireSenseUtils::harmonizeFireData" %in% called) # the walk reaches the fire data
  missing <- setdiff(called, keyed)
  expect_identical(missing, character(0))
})

test_that("the cached harmonizeFireData() call is keyed on the functions it calls", {
  txt <- paste(deparse(moduleDefs()[["prepare_SpreadFitFire_Vector"]]), collapse = " ")
  expect_match(txt, ".cacheExtra = fireSenseUtils::harmonizeFireDataDeps()", fixed = TRUE)
})

## The escape response comes from fireSenseUtils::escapedIgnitions() (2026-10-25); a change to it must remake
## the cached escape data, so prepEscapeFitData's key names it and every other package function it reaches.
test_that("prepEscapeFitData's cache key includes the package functions it reaches", {
  skip_if_not_installed("LandR")
  skip_if_not_installed("fireSenseUtils")
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("defineModule")),
                parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE))
  params <- eval(as.list(def[[1]][[3]])[["parameters"]], envir = asNamespace("SpaDES.core"))
  extra <- params$default[[match(".useCacheArgs", params$paramName)]]$prepEscapeFitData$.cacheExtra
  txt <- paste(deparse(extra), collapse = " ")
  keyed <- regmatches(txt, gregexpr("(LandR|fireSenseUtils)::[A-Za-z0-9_.]+", txt))[[1]]
  called <- packageCallsFrom("prepare_EscapeFit")
  expect_true("fireSenseUtils::escapedIgnitions" %in% called)
  expect_identical(setdiff(called, keyed), character(0))
})
