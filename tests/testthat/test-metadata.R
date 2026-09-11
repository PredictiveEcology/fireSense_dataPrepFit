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
    c(climateVariablesForFire     = "list",
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
      rstLCCs                     = "list",
      sppEquiv                    = "data.table",
      spreadFirePoints            = "list",
      spreadFirePolys             = "list",
      spreadFitAdditionalColNames = "character",
      standAgeMaps                = "list",
      studyArea                   = "SpatVector",
      studyArea_biomassParam      = "SpatVector",
      studyAreaReporting          = "sf")
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
      fireSense_ignitionFormula              = "character",
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
    sort(c(".plotInterval", ".saveInitialTime", ".saveInterval", ".studyAreaName",
           ".useCache", "areaMultiplier", "bufferForFireRaster", "cutoffForYoungAge",
           "dataYears", "estimateFuelClasses", "fireYears", "flammabilityThreshold",
           "forestedLCC", "fuelClassCol", "igAggFactor", "igFocalFactor",
           "minBufferSize", "modelAlgorithm", "nonflammableLCC",
           "nonForestCanBeYoungAge", "sppEquivCol", "spreadFitFilename",
           "spreadFitGoogleDriveFolder", "targetFuelClasses", "useCentroids",
           "useRasterizedFireForSpread", "whichModulesToPrepare"))
  )
})
