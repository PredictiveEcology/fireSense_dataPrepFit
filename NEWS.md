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
