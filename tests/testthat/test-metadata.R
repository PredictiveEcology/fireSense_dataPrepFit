## The module's metadata is its public contract: a project using this module binds
## to these object names and classes. Renaming or retyping one breaks every caller,
## which is exactly the class of change the raster -> terra migration makes, so it is
## worth asserting here rather than discovering downstream.
##
## When a change is deliberate, update this file in the same commit and bump the
## module version to match: removed, renamed or retyped is a MAJOR bump.

test_that("module metadata parses", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_type(md, "list")
  expect_identical(md$name, moduleName)
})

test_that("inputs are the expected names and classes", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  inputs <- stats::setNames(md$inputObjects$objectClass, md$inputObjects$objectName)
  expect_identical(
    inputs[order(names(inputs))],
    c(climateVariables            = "list",
      climateVariablesForFire     = "list",
      cohortDatas                 = "list",
      firePolys                   = "list",
      firePolysForAge             = "list",
      historicalClimateRasters    = "list",
      historicalFireRaster        = "SpatRaster",
      ignitionFirePoints          = "SpatVector",
      missingLCCgroup             = "character",
      nonForestedLCCGroups        = "list",
      pixelGroupMaps              = "list",
      propFlammables              = "list",
      rasterToMatch               = "SpatRaster",
      rasterToMatch_biomassParam  = "SpatRaster",
      rasterToMatchLarge          = "SpatRaster",
      rstLCCs                     = "list",
      sppEquiv                    = "data.table",
      spreadFirePoints            = "list",
      spreadFirePolys             = "list",
      spreadFitAdditionalColNames = "character",
      standAgeMaps                = "list",
      studyArea                   = "SpatVector",
      studyArea_biomassParam      = "SpatVector")
  )
})

test_that("outputs are the expected names and classes", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  outputs <- stats::setNames(md$outputObjects$objectClass, md$outputObjects$objectName)
  expect_identical(
    outputs[order(names(outputs))],
    c(climateVariables                       = "list",
      climateVariablesForFire                = "list",
      fireBufferedListDT                     = "list",
      fireSense_annualSpreadFitCovariates    = "list",
      fireSense_escapeCovariates             = "data.table",
      fireSense_escapeFormula                = "character",
      fireSense_ignitionCovariates           = "data.table",
      fireSense_nonAnnualSpreadFitCovariates = "list",
      fireSense_spreadFormula                = "character",
      flammableRTM                           = "SpatRaster",
      flammableRTMs                          = "list",
      fuelClassTable                         = "data.table",
      ignitionFirePoints                     = "SpatVector",
      ignitionFitRTM                         = "SpatRaster",
      landcoverDT                            = "data.table",
      landcoverDTs                           = "list",
      lightningMaps                          = "SpatRaster",
      missingLCCgroup                        = "character",
      nonForest_timeSinceDisturbance         = "SpatRaster",
      nonForest_timeSinceDisturbances        = "list",
      nonForestedLCCGroups                   = "list",
      propFlammable                          = "SpatRaster",
      rstLCC                                 = "SpatRaster",
      rstLCC_RTM                             = "SpatRaster",
      rstLCCs                                = "list",
      sppColorVect                           = "character",
      sppEquiv                               = "data.table",
      sppNameVector                          = "character",
      spreadClimateSelection                 = "data.table",
      spreadFirePoints                       = "list",
      spreadFirePolys                        = "list",
      spreadFitPreRun                        = "data.frame",
      standAgeMap                            = "SpatRaster",
      studyAreaWithSpreadParams              = "sf")
  )
})

test_that("parameters are the expected names", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_identical(
    sort(md$parameters$paramName),
    sort(c(".studyAreaName",
           ".useCache", ".useCacheArgs", "areaMultiplier", "bufferForFireRaster", "cutoffForYoungAge",
           "dataYears", "estimateFuelClasses", "fireYears", "flammabilityThreshold",
           "forestedLCC", "fuelClassCol", "igAggFactor",
           "minBufferSize", "nonflammableLCC",
           "nonForestCanBeYoungAge", "sppEquivCol", "spreadFitFilename",
           "spreadFitGoogleDriveFolder", "targetFuelClasses",
           "useRasterizedFireForSpread", "whichModulesToPrepare"))
  )
})

test_that("every object .inputObjects() assigns is a declared input", {
  ## SpaDES restores only a module's declared inputs from a cached `.inputObjects`. An object assigned there
  ## but declared only as an output vanishes on a cache hit: that is how canClimateData got NULL
  ## `climateVariables` on the Mackenzie relaunch (2026-09-23), though the first run had worked.
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  src <- parse(file.path(modulePath, moduleName, paste0(moduleName, ".R")), keep.source = FALSE)
  io <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                 identical(as.character(e[[2]]), ".inputObjects"), as.list(src))[[1]][[3]]
  assigned <- character(0)
  walk <- function(e) {
    if (is.call(e)) {
      if (identical(e[[1]], as.name("<-")) && is.call(e[[2]]) && identical(e[[2]][[1]], as.name("$")) &&
          identical(e[[2]][[2]], as.name("sim")))
        assigned <<- c(assigned, as.character(e[[2]][[3]]))
      for (a in as.list(e)[-1]) if (!missing(a)) walk(a)
    }
  }
  walk(io)
  expect_true(length(assigned) > 0)
  expect_setequal(setdiff(unique(assigned), md$inputObjects$objectName), character(0))
})
