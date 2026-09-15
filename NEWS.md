# fireSense_dataPrepFit (development version)

- `init` is split in two. `init` keeps the Google Drive check for an existing fit (it cannot be cached,
  because caching would freeze that answer) and schedules a new `dataPrepBuild` event, which holds the
  land-cover, flammability, fuel-class, landcover-table and time-since-disturbance work. `dataPrepBuild`
  can be cached (add it to `.useCache`), so a warm cache no longer rebuilds all of it on every run
  (6.6 minutes per ELF observed). It runs at priority 0, straight after `init` and before other modules'
  `init` events, which is where the same work ran before. Results are unchanged with a cold or a warm
  cache. New parameter `.useCacheArgs` passes the LandR and fireSenseUtils functions `dataPrepBuild` calls
  as `.cacheExtra`, so a change to one of them re-runs the event. Version 1.2.0.9004.

- The cached `fuelClassPrep()` step passed `.omitArgs` (a typo for `omitArgs`) to `Cache()`. With caching on
  the two rasters it meant to omit were digested every time; with caching off (`spades.useCache = "eventsOnly"`)
  the stray argument made `Cache()`'s bypass call the result as a function (`could not find function "FUN"`).
  Version 1.2.0.9003.

- `fireRecordShapefile()` calls `reproducible::preProcess()` by its full name. caret, a
  fireSense_IgnitionFit dependency, defines its own `preProcess()` generic; attached after reproducible
  it masked the bare name, and the fire-record download stopped with 'argument "x" is missing, with no
  default'.
- Fire records are read with `fireregimetools` and downloaded with `reproducible::preProcess()`, so
  `reproducible::preProcessCheckURLs()` can recheck them against the server. NFDB ignition points come
  from `current_version/NFDB_point_shp.zip`; the URL used before (`NFDB_point.zip`) returns HTTP 404, so
  no NFDB release newer than the local copy could be fetched. NBAC perimeters still come from the newest
  release, and its URL is now part of the Cache key.

# fireSense_dataPrepFit 1.2.0

First release from `development` since `main` was last updated (2023-02-07). Full history: https://github.com/PredictiveEcology/fireSense_dataPrepFit/compare/4b4e921...v1.2.0

## Breaking changes

- Removed input `cohortData2001` (data.table).
- Removed input `cohortData2011` (data.table).
- Removed input `flammableRTM` (RasterLayer).
- Removed input `pixelGroupMap2001` (RasterLayer).
- Removed input `pixelGroupMap2011` (RasterLayer).
- Removed input `rstLCC` (RasterLayer).
- Removed input `standAgeMap2001` (RasterLayer).
- Removed input `standAgeMap2011` (RasterLayer).
- Input `ignitionFirePoints` is now `SpatVector` (was `list`).
- Input `rasterToMatch` is now `SpatRaster` (was `RasterLayer`).
- Input `studyArea` is now `SpatVector` (was `SpatialPolygonsDataFrame`).
- Removed output `firePolys` (list).
- Removed output `nonForest_timeSinceDisturbance2001` (RasterLayer).
- Removed output `nonForest_timeSinceDisturbance2011` (RasterLayer).
- Removed output `terrainDT` (data.table).
- Output `ignitionFitRTM` is now `SpatRaster` (was `RasterLayer`).
- Removed parameters: `.plotInitialTime`, `ignitionFuelClassCol`, `missingLCCgroup`, `spreadFuelClassCol`.

## New features

- New inputs: `climateVariablesForFire`, `cohortDatas`, `historicalFireRaster`, `missingLCCgroup`, `pixelGroupMaps`, `propFlammables`, `rasterToMatch_biomassParam`, `rstLCCs`, `spreadFirePolys`, `spreadFitAdditionalColNames`, `standAgeMaps`, `studyAreaReporting`, `studyArea_biomassParam`.
- New outputs: `climateVariables`, `climateVariablesForFire`, `flammableRTM`, `flammableRTMs`, `fuelClassTable`, `ignitionFirePoints`, `landcoverDTs`, `lightningMaps`, `missingLCCgroup`, `nonForest_timeSinceDisturbance`, `nonForest_timeSinceDisturbances`, `nonForestedLCCGroups`, `propFlammable`, `rstLCC`, `rstLCC_RTM`, `rstLCCs`, `sppColorVect`, `sppEquiv`, `sppNameVector`, `spreadFirePolys`, `spreadFitPreRun`, `standAgeMap`, `studyAreaWithSpreadParams`.
- New parameters: `bufferForFireRaster`, `dataYears`, `estimateFuelClasses`, `flammabilityThreshold`, `fuelClassCol`, `igFocalFactor`, `modelAlgorithm`, `spreadFitFilename`, `spreadFitGoogleDriveFolder`, `targetFuelClasses`, `useRasterizedFireForSpread`.

## Dependencies

- No longer depends on `spatialEco`.
- Now depends on `LandR`, `Require`, `SpaDES.project`, `climateData`, `reproducible`, `terra`.

## Testing

- testthat suite and CI (`testthat-module`), including a snapshot of the module's inputs, outputs and parameters in `tests/testthat/test-metadata.R`.
