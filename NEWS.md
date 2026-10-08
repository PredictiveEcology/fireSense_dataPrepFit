# fireSense_dataPrepFit (development version)

# fireSense_dataPrepFit 1.3.0

This release reworks how the fire history, climate and fuel data are prepared for fitting fireSense. Fuels are now described by the dominant and secondary tree types in each place, plus treed wetland. Rare land-cover types are grouped with similar ones, and empty inputs are dropped instead of fitted. The climate variable that best separates bad fire years is chosen for each study area. Fire records come from the national fire databases through one shared tool, and the default fire years now run from 1985 to the latest year with climate data.

Several errors are fixed that could quietly feed a fit the wrong data: neighbouring regions' settings leaking into a run, old cached results surviving a refit or a change in a helper function, rock not counted as non-flammable, and fire points outside the study area. New options support validation (holding back part of the fire history) and fitting with an intercept. Settings that did nothing were removed, and the other fireSense modules are referred to by their new lower-case names.

- Fixed: the nested `Biomass_borealDataPrep` runs no longer receive the parent sim's `ecoregionLayer` or `ecoregionRst`. With `localEcozones` in the run, whether the fit used its local regions depended on module load order; the fit now always uses `Biomass_borealDataPrep`'s default ecoregions.

- Fixed: `Init()` builds the single-ELF objects (`sppEquiv`, `sppNameVector`, `sppColorVect`, `nonForestedLCCGroups`, `missingLCCgroup`) from the run's own ELF's ledger row (`sim$.ELFind`), not from every row that intersects the study area. Before, the neighbours' rows were merged into `sppEquiv` and `sppNameVector`, and `nonForestedLCCGroups`/`missingLCCgroup` stayed the module default (`nf`) whenever a neighbour was read, so `fireSense_ignitionFit` fitted `nf` and `fireSense_ignitionPredict` stopped with "column not found: [nf]". The per-ELF lists (`sppEquivs`, `sppNameVectors`, `nonForestedLCCGroupsList`, `missingLCCgroupList`) still hold every row. Without `.ELFind`, or with no row for it, the behaviour is unchanged.

- Fixed: `spreadFitPreRun` is declared as an input (it was only an output), so the ledger rows `Init()` reads are in the cache keys of the cached events (`dataPrepBuild`, `prepSpreadFitData`, ...). Before, a cache entry saved before an ELF was refitted restored its old ledger row over the new one on a hit, and prediction used stale parameters (ELF 14.3: the old per-class fit instead of the dom/sec fit, then `spreadPredict` "dom_agb_Pn_po.Ps_me not found"). `studyAreaWithSpreadParams` is derived from `spreadFitPreRun` in `Init()`, so keying the latter keys both.

# fireSense_dataPrepFit 1.2.0.9025

- `Init()` no longer stops with "sppColVect has unique colour values for a single species" when the SpreadFit ledger rows that intersect the study area (the own ELF and its neighbours) give one species different colours. It keeps one colour per species, from the own ELF's row (`sim$.ELFind`, now a declared optional input) or else the first row's, and messages the species that differed.

# fireSense_dataPrepFit 1.2.0.9024

- A spread covariate that is empty is dropped (`fireSenseUtils::emptySpreadCovariates()`): all zero, or, for a fuel class, all on the `logMinB()` floor, which a column sum never saw (a fuel class with no biomass in the buffers is 3.6 everywhere). Spread climate is no longer rounded to whole numbers by `fireSenseUtils::climateRasterToDataTable()`: the x1000 integer storage is the only rounding, as in `fireSense_dataPrepPredict`. Fits change and cached fits re-key (`emptySpreadCovariates` is in `prepSpreadFitData`'s `.cacheExtra`). Needs the `fireSenseUtils` change in PredictiveEcology/fireSenseUtils (floor to be set once it has a version).

# fireSense_dataPrepFit 1.2.0.9023

- New parameter `spreadIntercept` (default `FALSE`): with `TRUE` the spread formula is `~ 1 + ...`, so `fireSense_spreadFit` fits a free intercept with centred covariates (`fireSenseUtils::spreadInterceptTxt`) and the coefficients describe variation about the centre. Covariates are rescaled to 0-1, so without an intercept the level of the linear predictor is set by the coefficients alone and the climate coefficient trades off against the fuel and non-forest ones (ELF 13.1: the other coefficients explain a median 73% of CMD's variance in a final population). `FALSE` keeps the formula string `~ 0 + ...` byte for byte. Needs fireSenseUtils >= 0.2.3.9083 (PredictiveEcology/fireSenseUtils#131).
- New parameter `minCovariateProp` (default `0.05`): a land-cover class covering less than this share of the ELF's flammable pixels (most recent data year) gets no spread coefficient of its own. A rare non-forest class joins the non-forest group nearest its burn rate; with less treed wetland than this there is no `treedWetland_agb` and its tree AGB stays in `dom_agb_*`/`sec_agb_*`. Treed wetland is under 5% in 7 of the 10 fitted ELFs. Needs fireSenseUtils >= 0.2.3.9081 (PredictiveEcology/fireSenseUtils#129). Fits and cached `dataPrepBuild`/`prepSpreadFitData` results of ELFs whose covariates change are invalid.
- Fixed: with one fitted ELF in the SpreadFit ledger, `Init()` now sets `nonForestedLCCGroups` and `missingLCCgroup` to that
  ELF's groups, so this run's ignition and escape fits use the same non-forest columns (`nfLCC_*`) that
  `fireSense_dataPrepPredict` builds for prediction. Before, they kept the module default (`nf`), and
  `fireSense_ignitionPredict` stopped with "column not found: [nf]".
- `youngAge` is resolved for each fire year, as it was before February 2026, instead of once per data year. `prepare_SpreadFit()` builds the fuels with `fireSenseCovariatesCreate(youngAge = FALSE)` (no `youngAge`, nothing zeroed at the data year) and puts `youngAge` in each fire year's annual table from `fireSenseUtils::youngAgeAtYear()`: the data year's time since disturbance aged to that year, reset by every fire since (all of `firePolysForAge`, not only the fitted buffers); `NA` time since disturbance is not young. `prepare_IgnitionFit()` does the same for the ignition (and escape) covariates with `fireSenseUtils::prepare_FuelCovsCoarseByYear()`, aggregating the per-year values to the `igAggFactor` grid. `nonForestCanBeYoungAge = FALSE` now stops in both. Needs the `fireSenseUtils` change in PredictiveEcology/fireSenseUtils#119 (floor to be set once it has a version). Spread and ignition fits made with the data-year `youngAge`, and cached `prepSpreadFitData` and ignition covariates, are invalid.

# fireSense_dataPrepFit (development version)

- The pooled `other_agb` spread covariate is gone and `fuelCovariates = "domSecOther"` is renamed `"domSecWetland"` (now the default): spread fuels are `dom_agb_<class>`, `sec_agb_<class>` and `treedWetland_agb`. Needs the `fireSenseUtils` change that renames the value (PredictiveEcology/fireSenseUtils#116). Spread fits made with `other_agb` need refitting, and cached `prepSpreadFitData` results change.

- New parameter `heldOutFold` (`NA`, `1` or `2`; the same parameter as in `fireSense_spreadFit`). With `1` or `2`, `Init()` does not read the SpreadFit ledger: `sim$spreadFitPreRun` stays NULL and `mod$haveSpreadFit` is FALSE, so the fold derives its own species, fuel and climate objects. `paramCheckOtherMods()` stops if `fireSense_spreadFit` has a different value; set all three modules with `.globals = list(heldOutFold = ...)`.

- `snow` is no longer a `reqdPkgs`: nothing used it, and attaching it printed two "partial argument match of 'along'" warnings per run (from snow's `.onLoad()`) under `warnPartialMatchArgs = TRUE`.
- `fireSense_EscapeFit` no longer exists (`fireSense_ignitionFit` fits ignition and escape): it is removed from the `whichModulesToPrepare`
  default (now `fireSense_ignitionFit` and `fireSense_spreadFit`) and naming it stops with a message. Preparing `fireSense_ignitionFit`
  schedules `prepEscapeFitData` as well; the event and its caching are unchanged.

- `whichModulesToPrepare` default and comparisons use the renamed `fireSense_ignitionFit` and `fireSense_spreadFit` (formerly
  `fireSense_IgnitionFit`, `fireSense_SpreadFit`). A project setting `whichModulesToPrepare` must use the new names.


- `forestedLCC`, `cutoffForYoungAge`, `nonForestCanBeYoungAge`, `flammabilityThreshold`, `fuelClassCol`
  and `igAggFactor` now default to `fireSenseUtils::fireSenseForestedLCC`, `fireSenseYoungAgeCutoff`,
  `fireSenseNonForestCanBeYoungAge`, `fireSenseFlammabilityThreshold`, `fireSenseFuelClassCol` and
  `fireSenseIgAggFactor`, the same values as before but from the single source of truth also used by
  `fireSense_dataPrepPredict`. New parameter `scanfiVersion` (default `fireSenseUtils::fireSenseSCANFIVersion`,
  `"V3"`) is passed to `fireSenseUtils::makeFireSenseLCC()`. `fireSenseUtils:::fireSenseCovariatesCreate` is
  now called as `fireSenseUtils::fireSenseCovariatesCreate` (it is exported). Needs
  `fireSenseUtils@development (>= 0.2.3.9062)`. Version 1.2.0.9018.
- Fixed: `nonflammableLCC`'s default (`c(0, 20, 31, 32, 33)`) missed SCANFI's rock/exposed code
  (`30`), so rock entered fits as flammable non-forest. The default now comes from
  `fireSenseUtils::fireSenseNonflammableLCC`, the single source of truth `makeFireSenseLCC()`
  also uses. Needs `fireSenseUtils@development (>= 0.2.3.9060)`. Version 1.2.0.9017.
- New parameter `fuelCovariates` (default `"domSecOther"`): `prepare_SpreadFit()` now builds the
  spread covariates as `dom_agb_<class>`/`sec_agb_<class>` (the ELF's two fuel classes with the
  most total treed AGB), `other_agb` and `treedWetland_agb`, chosen once per ELF by the new
  `fireSenseUtils::chooseDomSecFuelClasses()` and recorded in `sim$fuelClassRoles`. `rstLCC` is
  now passed to `fireSenseCovariatesCreate()` (previously never passed, so `treedWetland` never
  appeared). `fuelCovariates = "species"` keeps the previous one-column-per-fuel-class behaviour.
  `chooseDomSecFuelClasses()` is added to the `prepSpreadFitData` cache key (`.useCacheArgs`).
  Needs `fireSenseUtils@development (>= 0.2.3.9057)`. Version 1.2.0.9016.
- A cached `prepSpreadFitData` event, or cached `harmonizeFireData()` call, now re-runs when a fireSenseUtils function it calls changes: they are keyed on those functions (`.useCacheArgs`, `fireSenseUtils::harmonizeFireDataDeps()`). Requires fireSenseUtils >= 0.2.3.9053. Version 1.2.0.9015.
- `spreadFitFilename` now defaults to `"latest"`: each polygon's fit comes from the most recent ledger file in
  `spreadFitGoogleDriveFolder` that has it (`fireSenseUtils::latestSpreadFits()`, which reads only the
  current model's files, `fireSenseParams_*<fireSenseUtils::spreadFitFileTag>.rds`). So "this ELF has a
  fit" means some such file has it. A named file is read as before. Needs reproducible >= 3.2.1.9042, whose `CacheGeo()` re-reads a local ledger file that has
  changed. Version 1.2.0.9014.
- `prepare_SpreadFit()` built the spread formula's RHS as `climate + youngAgeTxt + vegCols`, but `vegCols`
  (derived from `fireSenseVegData`) already includes a `youngAge` column whenever
  `fireSenseUtils::fireSenseCovariatesCreate()` finds young forest or non-forest pixels. The formula then
  listed `youngAge` twice; `terms()` silently drops the duplicate, leaving one fewer distinct term than the
  covariates actually used downstream. `youngAgeTxt` and the spread climate variable are now excluded from
  `vegCols` before building the RHS. Version 1.2.0.9013.
- `prepare_SpreadFitFire_Vector()` filtered `spreadFirePolys` by pixel size only, while `spreadFirePoints`
  was filtered with `escapedFires()` (pixel size and `escapeSizeHa`, PR #44). A fire between one pixel and
  `escapeSizeHa` then survived in the polygons but not the points, and `fireSenseUtils::harmonizeFireData()`
  stopped with "spread fire point and poly harmonization error in dataPrepFit". `spreadFirePolys` is now
  filtered with `escapedFires()` too. Version 1.2.0.9012.
- New parameter `escapeSizeHa` (default 50). A fire counts as escaped when it reached that size, in the
  escape model's response and in the fires the spread model is fitted to; before, any fire larger than one
  pixel (about 6 ha) counted. New output `nonEscapedFireSizesHa`: the study area's natural-cause fire sizes
  below it, for fireSense to size non-escaped ignitions. Version 1.2.0.9011.
- The default fire years are now 1985 to the latest year with historical climate
  (`climateData::latestHistoricalYear()`, 2024 today), and the default `dataYears` are 1985, 1990, 2000, 2010
  and 2020. The old default, 2002:2025, ran past the climate data, so every project had to set `fireYears`.
  fireregimetools is floored at 0.1.0.9008 (FOR-CAST main), which reads only the study area's part of the fire
  records. Version 1.2.0.9010.

- `climateVariables` and `rasterToMatchLarge` are now declared inputs. `.inputObjects` sets both, and a cached
  `.inputObjects` restores only declared inputs, so on a cache hit they vanished: canClimateData then failed with
  "`.l` must be a list, not NULL" (the Mackenzie relaunch, 2026-09-23), although the first run had worked. A new
  test fails if `.inputObjects` assigns any object that is not a declared input.

## Climate variables

- `climateVariablesForFire` now defaults to `ignition = c("CMD", "cumMDC", "CMD_sm", "CMD_sp")` (IgnitionFit's xgboost uses
  them all) and `spread = "auto"`. Unless supplied, `climateVariables` (what canClimateData prepares) is built from
  it, so the two agree and a user need not choose; canClimateData keeps its own default when used alone. Names
  are accepted with or without underscores.
- `spread = "auto"` picks, per study area, the candidate that best separates the bad fire years: the AUC of
  the top quarter of fire years (by area burned) against the rest, the Spearman correlation breaking ties. Below
  an AUC of 0.6 no variable separates them, and the first candidate is used, with a warning. The scores are
  the new output `spreadClimateSelection`.
- When predicting from existing fits of several ELFs, their spread variables may differ; every one is
  prepared (the union), instead of stopping.

- The non-xgboost ignition path is removed, as in fireSense_IgnitionFit. `prepare_IgnitionFit()` no longer
  builds `fireSense_ignitionFormula`; that branch used an object that was never defined, so it could only
  stop with an error. Removed with it: parameter `modelAlgorithm` and output `fireSense_ignitionFormula`.
- Removed, because the module declared them but never read them: parameters `igFocalFactor`,
  `useCentroids`, `.plotInterval`, `.saveInitialTime`, `.saveInterval`, and input `studyAreaReporting`.
  Stop setting them.
- Removed `dataPrepInit()`, which has had no caller since its event was removed. A length-one
  `climateVariablesForFire` is therefore not expanded to `ignition` and `spread`, and its names are not
  checked against `historicalClimateRasters`: supply both elements. Version 1.2.0.9007.

- Dead code removed (`Save()`, `rmMissingPixels()`, and the module copies of `putBackIntoRaster()` and
  `calcNonForestYoungAge()`, which live in `fireSenseUtils`), every function documented, metadata
  descriptions and the manual brought up to date. No change in behaviour.

- `dataPrepBuild` now always clips `ignitionFirePoints` to the study-area polygon, via the new
  `clipPointsToStudyArea()`. The clip used to run only inside an `if (!same.crs(points, rasterToMatch))`
  branch, and clipped to the `rasterToMatch` rectangle rather than the polygon, so points that arrived
  already in the raster's CRS were never clipped at all. `prepare_IgnitionFit()` asserts that every
  ignition point is within the study area, and ELF 5.1.2 failed it three times on a single NFDB point
  618 m outside the polygon but inside the raster extent. Version 1.2.0.9006.

- `.inputObjects` no longer stops with `object 'LCC' not found` when `rstLCCs` is supplied. Its per-year
  step returned the land cover and flammable proportion it builds whether or not it had built them; a
  supplied `rstLCCs` (with its `propFlammables`) is now left as it is. Version 1.2.0.9005.

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
