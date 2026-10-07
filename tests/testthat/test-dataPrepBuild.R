## init used to hold everything dataPrepBuild now does. init cannot be cached: it asks Google Drive
## whether a fit already exists, and caching would freeze that answer. So the land-cover, fuel-class
## and time-since-disturbance work re-ran on every run, 6.6 minutes per ELF on a warm cache. That work
## is now the cacheable dataPrepBuild event. These tests pin the split: a later merge must not fold
## the work back into init, and the cached event must not depend on anything its cache key does not
## see. test-dataPrepBuildEquivalence.R checks that the results are unchanged.

moduleExprs <- function() {
  parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE)
}

isAssign <- function(x) {
  is.call(x) && (identical(x[[1]], as.name("<-")) || identical(x[[1]], as.name("=")))
}

moduleFunctionBody <- function(name) {
  def <- Filter(function(x) isAssign(x) && identical(x[[2]], as.name(name)), moduleExprs())
  stopifnot(length(def) == 1L)
  body(eval(def[[1]][[3]]))
}

## the module's metadata, evaluated from source, so the checkout directory's name does not matter
moduleMetadataElement <- function(element) {
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("defineModule")), moduleExprs())
  stopifnot(length(def) == 1L)
  eval(as.list(def[[1]][[3]])[[element]], envir = asNamespace("SpaDES.core"))
}

## every call inside an expression
allCalls <- function(expr) {
  if (!is.call(expr)) return(list())
  args <- as.list(expr)[-1]
  args <- args[!vapply(args, function(a) missing(a), logical(1))] # the empty `i` in dt[, j]
  c(list(expr), unlist(lapply(args, allCalls), recursive = FALSE))
}

## one arm of doEvent's switch()
eventArm <- function(event) {
  sw <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("switch")),
               as.list(moduleFunctionBody("doEvent.fireSense_dataPrepFit")))
  stopifnot(length(sw) == 1L)
  as.list(sw[[1]])[[event]]
}

## names used as root$name or root[["name"]]
refNames <- function(expr, root) {
  refs <- Filter(function(x) {
    (identical(x[[1]], as.name("$")) && identical(x[[2]], as.name(root)) && is.name(x[[3]])) ||
      (identical(x[[1]], as.name("[[")) && identical(x[[2]], as.name(root)) && is.character(x[[3]]))
  }, allCalls(expr))
  unique(vapply(refs, function(x) as.character(x[[3]]), character(1)))
}

## names assigned as root$name <- value, including root$name[i] <- value and the like
assignedNames <- function(expr, root) {
  lhsName <- function(lhs) {
    while (is.call(lhs)) {
      if (identical(lhs[[1]], as.name("$")) && identical(lhs[[2]], as.name(root)))
        return(as.character(lhs[[3]]))
      lhs <- lhs[[2]]
    }
    NULL
  }
  unique(unlist(lapply(Filter(isAssign, allCalls(expr)), function(x) lhsName(x[[2]]))))
}

buildWork <- c("defineFlammable", "fuelClassPrep", "assessFuelClasses", "makeLandcoverDT",
               "correctMissingLCC", "makeTSD")

test_that("init does none of the land-cover and fuel-class work; dataPrepBuild does all of it", {
  initCode <- c(all.names(eventArm("init")), all.names(moduleFunctionBody("Init")))
  dpfExpectNone(intersect(buildWork, initCode), "build work still called from init")
  dpfExpectNone(setdiff(buildWork, all.names(moduleFunctionBody("dataPrepBuild"))),
                "build work missing from dataPrepBuild")
  expect_true("dataPrepBuild" %in% all.names(eventArm("dataPrepBuild")))
})

test_that("init schedules dataPrepBuild at the start, ahead of every other module's init", {
  scheduled <- Filter(function(x) identical(x[[1]], as.name("scheduleEvent")), allCalls(eventArm("init")))
  scheduled <- lapply(scheduled, function(x) match.call(SpaDES.core::scheduleEvent, x))
  build <- Filter(function(x) identical(x$eventType, "dataPrepBuild"), scheduled)
  expect_length(build, 1L)
  expect_identical(deparse(build[[1]]$eventTime), "start(sim)")
  ## init events are scheduled at .first(); anything lower runs before all of them, which is where
  ## this work ran while it was part of init
  expect_lt(eval(build[[1]]$eventPriority), SpaDES.core::.first())
})

test_that("dataPrepBuild reads only what its cache key covers, and writes only outputs", {
  inputs <- moduleMetadataElement("inputObjects")$objectName
  outputs <- moduleMetadataElement("outputObjects")$objectName
  build <- moduleFunctionBody("dataPrepBuild")

  ## A cached event's key covers this module's expected inputs, its functions and its mod$ contents.
  ## Outputs are read too, once written. An object that is neither is invisible to the key.
  printedOnly <- "nonForestedLCCGroupsList" # only print()ed, when a fit already exists
  dpfExpectNone(setdiff(refNames(build, "sim"), c(inputs, outputs, printedOnly)),
                "sim objects dataPrepBuild reads that its cache key does not cover")
  ## in particular, whether the user supplied the fuel objects is recorded by the uncached init
  expect_true(".userSuppliedObjNames" %in% refNames(moduleFunctionBody("Init"), "sim"))

  ## a cache hit restores declared outputs only, so anything else written here would be lost
  dpfExpectNone(setdiff(assignedNames(build, "sim"), outputs),
                "sim objects dataPrepBuild writes that are not declared outputs")

  ## mod$ values read here must already exist when the event's key is computed
  modSetEarlier <- c(assignedNames(moduleFunctionBody(".inputObjects"), "mod"),
                     assignedNames(moduleFunctionBody("Init"), "mod"))
  modRead <- setdiff(refNames(build, "mod"), assignedNames(build, "mod"))
  dpfExpectNone(setdiff(modRead, modSetEarlier), "mod$ values dataPrepBuild reads that are not set before it")
})

test_that("dataPrepBuild's cache key includes the LandR and fireSenseUtils functions it calls", {
  skip_if_not_installed("LandR")
  skip_if_not_installed("fireSenseUtils")
  params <- moduleMetadataElement("parameters")
  extra <- params$default[[match(".useCacheArgs", params$paramName)]]$dataPrepBuild$.cacheExtra
  expect_true(is.call(extra))

  isNamespaced <- function(x) is.call(x) && is.name(x[[1]]) && as.character(x[[1]]) %in% c("::", ":::")
  qualifiedName <- function(x) paste0(as.character(x[[2]]), "::", as.character(x[[3]]))
  keyed <- vapply(Filter(isNamespaced, as.list(extra)[-1]), qualifiedName, character(1))

  build <- moduleFunctionBody("dataPrepBuild")
  pkgs <- c("LandR", "fireSenseUtils")
  qualified <- vapply(Filter(function(x) isNamespaced(x) && as.character(x[[2]]) %in% pkgs, allCalls(build)),
                      qualifiedName, character(1))
  bare <- unlist(lapply(pkgs, function(p) {
    hits <- intersect(all.names(build), getNamespaceExports(p))
    if (length(hits)) paste0(p, "::", hits)
  }))
  called <- unique(c(qualified, bare))
  called <- called[vapply(called, function(f) is.function(eval(str2lang(f))), logical(1))] # not data
  dpfExpectNone(setdiff(called, keyed), "functions dataPrepBuild calls that are not in its .cacheExtra")
})
