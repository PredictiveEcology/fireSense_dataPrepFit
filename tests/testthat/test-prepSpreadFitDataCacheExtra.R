## A cached event's key covers this module's code, not the package functions it calls. The run
## scripts cache prepSpreadFitData, so a fix in fireSenseUtils (e.g. cleanUpSpreadFirePoints(),
## makeMutuallyExclusive()) came back as the old spread-fit data. The event's `.cacheExtra` must name
## every LandR and fireSenseUtils function the event reaches, through the module functions it calls.

moduleDefs <- function() {
  files <- c(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"),
             list.files(testthat::test_path("..", "..", "R"), "\\.R$", full.names = TRUE))
  isDef <- function(d) is.call(d) && as.character(d[[1]]) %in% c("<-", "=") && is.name(d[[2]]) &&
    is.call(d[[3]]) && identical(d[[3]][[1]], as.name("function"))
  defs <- Filter(isDef, unlist(lapply(files, function(f) as.list(parse(f, keep.source = FALSE))),
                               recursive = FALSE))
  stats::setNames(lapply(defs, function(d) d[[3]]), vapply(defs, function(d) as.character(d[[2]]), ""))
}

## the package functions a `.cacheExtra` expression names: `pkg::fun`, and the functions in the list a
## `pkg::...Deps()` call returns (matched to the package's exports)
keyedFunctions <- function(extra) {
  txt <- paste(deparse(extra), collapse = " ")
  named <- regmatches(txt, gregexpr("(LandR|fireSenseUtils)::[A-Za-z0-9_.]+", txt))[[1]]
  deps <- regmatches(txt, gregexpr("(LandR|fireSenseUtils)::[A-Za-z0-9_.]+Deps\\(\\)", txt))[[1]]
  fromDeps <- unlist(lapply(deps, function(d) {
    pkg <- sub("::.*", "", d)
    fns <- eval(str2lang(d))
    ex <- getNamespaceExports(pkg)
    unlist(lapply(fns, function(f) {
      hit <- ex[vapply(ex, function(e) identical(get(e, asNamespace(pkg)), f), logical(1))]
      if (length(hit)) paste0(pkg, "::", hit)
    }))
  }))
  unique(c(named, fromDeps))
}

## package functions reached from `start`, through the module's own functions
packageCallsFrom <- function(start, pkgs = c("LandR", "fireSenseUtils")) {
  defs <- moduleDefs()
  seen <- character()
  todo <- start
  out <- character()
  while (length(todo)) {
    f <- todo[1]
    todo <- todo[-1]
    if (f %in% seen) next
    seen <- c(seen, f)
    txt <- paste(deparse(defs[[f]]), collapse = " ")
    qualified <- regmatches(txt, gregexpr("(LandR|fireSenseUtils):::?[A-Za-z0-9_.]+", txt))[[1]]
    nms <- all.names(defs[[f]])
    bare <- unlist(lapply(pkgs, function(p) {
      hits <- intersect(nms, getNamespaceExports(p))
      if (length(hits)) paste0(p, "::", hits)
    }))
    out <- c(out, sub(":::", "::", qualified), bare)
    todo <- c(todo, setdiff(intersect(nms, names(defs)), seen))
  }
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
  keyed <- keyedFunctions(extra)

  called <- packageCallsFrom("prepare_SpreadFit")
  expect_true("fireSenseUtils::harmonizeFireData" %in% called) # the walk reaches the fire data
  missing <- setdiff(called, keyed)
  expect_identical(missing, character(0))
})

test_that("the cached harmonizeFireData() call is keyed on the functions it calls", {
  txt <- paste(deparse(moduleDefs()[["prepare_SpreadFitFire_Vector"]]), collapse = " ")
  expect_match(txt, ".cacheExtra = fireSenseUtils::harmonizeFireDataDeps()", fixed = TRUE)
})

## The escape response comes from fireSenseUtils::escapeIgnitionPoints() (2026-11-10); a change to it must remake
## the cached escape data, so prepEscapeFitData's key names it (with escapeFireSetDeps()) and every other package
## function it reaches. The same for prepIgnitionFitData, which makes those points.
test_that("prepEscapeFitData's cache key includes the package functions it reaches", {
  skip_if_not_installed("LandR")
  skip_if_not_installed("fireSenseUtils")
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("defineModule")),
                parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE))
  params <- eval(as.list(def[[1]][[3]])[["parameters"]], envir = asNamespace("SpaDES.core"))
  extra <- params$default[[match(".useCacheArgs", params$paramName)]]$prepEscapeFitData$.cacheExtra
  keyed <- keyedFunctions(extra)
  called <- packageCallsFrom("prepare_EscapeFit")
  expect_identical(setdiff(called, keyed), character(0))
})

test_that("prepIgnitionFitData's cache key includes the functions that make the core ignitions and escapes", {
  skip_if_not_installed("fireSenseUtils")
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("defineModule")),
                parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE))
  params <- eval(as.list(def[[1]][[3]])[["parameters"]], envir = asNamespace("SpaDES.core"))
  extras <- params$default[[match(".useCacheArgs", params$paramName)]]
  for (ev in c("prepIgnitionFitData", "prepEscapeFitData", "prepSpreadFitData"))
    expect_true(all(c("fireSenseUtils::escapeFireSet", "fireSenseUtils::escapeIgnitionPoints") %in%
                      keyedFunctions(extras[[ev]]$.cacheExtra)), info = ev)
})
