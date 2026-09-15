## Helpers for test-dataPrepBuildEquivalence.R. They are plain functions, so a project can also
## source() this file to write a fixture with dpfFixtureFromSim().

dpfModule <- "fireSense_dataPrepFit"

## Expect nothing to be found, and name what was, so a failure says which object or function.
dpfExpectNone <- function(found, what) {
  testthat::expect(length(found) == 0L, paste0(what, ": ", paste(found, collapse = ", ")))
}

## terra objects are external pointers and do not survive serialisation, so they are wrapped with
## their values in memory. Plain lists are walked; every other object is kept as it is.
dpfPack <- function(x) {
  if (inherits(x, "SpatRaster"))
    return(structure(list(terra::wrap(x, proxy = FALSE)), class = "dpfPackedTerra"))
  if (inherits(x, "SpatVector"))
    return(structure(list(terra::wrap(x)), class = "dpfPackedTerra"))
  if (identical(class(x), "list")) x[] <- lapply(x, dpfPack)
  x
}

dpfUnpack <- function(x) {
  if (inherits(x, "dpfPackedTerra")) return(terra::unwrap(x[[1]]))
  if (identical(class(x), "list")) x[] <- lapply(x, dpfUnpack)
  x
}

## A simList's objects, without functions. SpaDES.core keeps its own state under dot-prefixed
## names, which ls() leaves out.
dpfSimObjects <- function(sim) {
  e <- SpaDES.core::envir(sim)
  nms <- ls(e)
  nms <- nms[!vapply(nms, function(n) is.function(get(n, envir = e)), logical(1))]
  mget(nms, envir = e)
}

## This module's mod$ contents. `mod` and `Par` are bindings back to the module, not contents.
dpfModObjects <- function(sim) {
  modEnv <- sim[[".modObjs"]][[dpfModule]]
  mget(setdiff(ls(modEnv, all.names = TRUE), c("mod", "Par")), envir = modEnv)
}

## Write a fixture from a simList stopped just before this module's init, e.g. from
## spades(sim, events = list(.stopBefore = list(fireSense_dataPrepFit = "init"))).
## By then this module's .inputObjects has run: its products are among the objects, and the mod$
## values it set are kept too.
dpfFixtureFromSim <- function(sim, file) {
  fixture <- list(
    objects = dpfPack(dpfSimObjects(sim)),
    mod = dpfPack(dpfModObjects(sim)),
    params = SpaDES.core::params(sim)[[dpfModule]],
    times = list(start = as.numeric(SpaDES.core::start(sim)), end = as.numeric(SpaDES.core::end(sim)),
                 timeunit = SpaDES.core::timeunit(sim)),
    userSuppliedObjNames = sim$.userSuppliedObjNames,
    inputPath = SpaDES.core::inputPath(sim),
    created = format(Sys.time()))
  qs2::qs_save(fixture, file)
  invisible(file)
}

## Run this module alone from a fixture. By default its events are cached as a project caches them,
## plus dataPrepBuild (the pre-split module has no such event, so the name has no effect there).
## Plots are off: they are side effects, not results.
dpfRunFromFixture <- function(fixture, modulePath, cachePath, outputPath, inputPath = fixture$inputPath,
                              useCache = c("dataPrepBuild", "prepIgnitionFitData", "prepEscapeFitData",
                                           "prepSpreadFitData")) {
  ## A fixture run uses the packages already installed. With Require on, simInit resolves the module's
  ## reqdPkgs and can install into the library that other running projects load from.
  op <- options(spades.useRequire = FALSE)
  on.exit(options(op), add = TRUE)
  ## prepSpreadFitData samples fire-buffer pixels inside forked workers (fireSenseUtils::bufferToArea.list
  ## uses parallel::mcMap). With R's default RNG each worker is seeded independently, so the same code gives
  ## different buffers on every run. With L'Ecuyer-CMRG each worker's stream derives from the session seed.
  rngOld <- RNGkind()
  on.exit(do.call(RNGkind, as.list(rngOld)), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  params <- fixture$params
  params$.useCache <- useCache
  params$.plots <- NA
  objects <- dpfUnpack(fixture$objects)
  ## .inputObjects already ran in the captured run. simInit skips it when every expected input is
  ## supplied, so the inputs that run did not have are passed as the NULL they were.
  md <- SpaDES.core::moduleMetadata(module = dpfModule, path = modulePath)
  objects[setdiff(md$inputObjects$objectName, names(objects))] <- list(NULL)
  sim <- SpaDES.core::simInit(
    times = fixture$times, modules = dpfModule,
    params = stats::setNames(list(params), dpfModule),
    objects = objects,
    paths = list(modulePath = modulePath, cachePath = cachePath,
                 inputPath = inputPath, outputPath = outputPath))
  list2env(dpfUnpack(fixture$mod), envir = sim[[".modObjs"]][[dpfModule]])
  ## simInit records every object passed in as supplied by the user; restore the captured run's
  ## record, which decides whether fuel classes are estimated
  sim$.userSuppliedObjNames <- fixture$userSuppliedObjNames
  messages <- character()
  set.seed(1)
  sim <- withCallingHandlers(
    SpaDES.core::spades(sim),
    message = function(m) {
      messages <<- c(messages, gsub("\033\\[[0-9;]*m", "", conditionMessage(m)))
      invokeRestart("muffleMessage")
    })
  list(sim = sim, messages = messages)
}

## Reduce an object to its result: raster values and geometry rather than file names or pointers,
## a data.table's contents rather than its internal attributes, and no Cache tags.
dpfNormalize <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.function(x) || is.environment(x)) return(class(x)[1])
  if (inherits(x, "SpatVector")) x <- sf::st_as_sf(x)
  if (inherits(x, "SpatRaster"))
    return(list(crs = terra::crs(x), extent = as.vector(terra::ext(x)), res = terra::res(x),
                names = names(x), cats = terra::cats(x),
                values = digest::digest(terra::values(x, mat = FALSE), algo = "xxhash64")))
  if (!is.language(x))
    for (a in c("tags", ".Cache", "call", "callInCache")) attr(x, a) <- NULL
  if (inherits(x, "data.table")) {
    x <- as.data.frame(x)
    attr(x, "sorted") <- NULL
    attr(x, "index") <- NULL
  }
  if (identical(class(x), "list")) x[] <- lapply(x, dpfNormalize)
  x
}

## A run's objects and this module's mod$ contents, normalised and in name order.
dpfResult <- function(sim) {
  objs <- dpfSimObjects(sim)
  modObjs <- dpfModObjects(sim)
  list(objects = lapply(objs[order(names(objs))], dpfNormalize),
       mod = lapply(modObjs[order(names(modObjs))], dpfNormalize))
}

## Names of whatever differs between two dpfResult()s.
dpfDifferences <- function(a, b) {
  diffIn <- function(x, y) {
    nms <- union(names(x), names(y))
    nms[!vapply(nms, function(n) isTRUE(all.equal(x[[n]], y[[n]])), logical(1))]
  }
  modDiffs <- diffIn(a$mod, b$mod)
  c(diffIn(a$objects, b$objects), if (length(modDiffs)) paste0("mod$", modDiffs))
}
