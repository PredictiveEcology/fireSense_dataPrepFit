defineModule(sim, list(
  name = "fireSense_dataPrepFit",
  description = "Prepare data required by `fireSense_IginitionFit`, `fireSense_EscapeFit`, and `fireSense_SpreadFit`.",
  keywords = "fireSense",
  authors = c(
    person("Ian", "Eddy", role = c("aut", "cre"), email = "ian.eddy@nrcan-rncan.gc.ca"),
    person("Eliot", "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = "aut"),
    person(c("Alex", "M"), "Chubaty", role = "ctb", email = "achubaty@for-cast.ca")
  ),
  childModules = character(0),
  version = list(fireSense_dataPrepFit = "1.2.0.9010"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = deparse(list("README.md", "fireSense_dataPrepFit.Rmd")),
  loadOrder = list(before = c("Biomass_speciesData", "Biomass_borealDataPrep", "Biomass_speciesParameters")),
  reqdPkgs = list("data.table", "fastDummies", "reproducible", "Require",
                  "PredictiveEcology/climateData@development (>= 2.2.3.9006)",
                  "PredictiveEcology/fireSenseUtils@development (>= 0.2.3.9024)",
                  "FOR-CAST/fireregimetools@main (>= 0.1.0.9008)",
                  "ggplot2", "parallel", "purrr", "raster", "sf", "sp",
                  "PredictiveEcology/LandR@development (>= 1.2.0.9015)",
                  "PredictiveEcology/SpaDES.core@development (>= 2.0.2.9006)",
                  "PredictiveEcology/SpaDES.project@development",
                  "PredictiveEcology/SpaDES.tools@development (>= 2.1.1.9000)",
                  "snow", "terra"),
  parameters = bindrows(
    defineParameter("areaMultiplier", c("numeric", "name"), quote(fireSenseUtils::multiplier), NA, NA,
                    paste("Size of the unburned buffer sampled around each fire: a scalar (buffer area is",
                          "`areaMultiplier * fireSize`) or a quoted function of `fireSize`.",
                          "See `?fireSenseUtils::bufferToArea`.")),
    defineParameter("bufferForFireRaster", "numeric", 1000, 0, NA,
                    paste("Buffer distance within which separate patches of burned pixels count as one fire.",
                          "Only used when `useRasterizedFireForSpread = TRUE`.")),
    defineParameter("cutoffForYoungAge", "numeric", 15, NA, NA,
                    "Age at and below which pixels are considered 'young' (`young <- age <= cutoffForYoungAge`)"),
    defineParameter("dataYears", "integer", c(1985L, 1990L, 2000L, 2010L, 2020L), NA_integer_, NA_integer_,
                    paste("Two or more increasing years for which vegetation, land cover and stand age are built",
                          "(`cohortDatas`, `rstLCCs`, `standAgeMaps`, ...).",
                          "Each fire year uses the data year at or before it, so no `fireYears` may precede the first,",
                          "and every data year needs at least one fire year before the next data year.")),
    defineParameter("estimateFuelClasses", "logical", TRUE, NA, NA,
                    paste("Estimate fuel classes with `fireSenseUtils::assessFuelClasses`? Skipped if the user supplies",
                          "`nonForestedLCCGroups`, `fuelClassTable` or a `FuelClass` column that differs from LandR's,",
                          "or if a previous SpreadFit exists for the study area.")),
    defineParameter("fireYears", "integer", 1985L:climateData::latestHistoricalYear(), NA, NA,
                    paste("Years of fire records to use for fitting. None may precede the first of `dataYears`,",
                          "and `historicalClimateRasters` must cover all of them. The default runs from 1985, the",
                          "first SCANFI V2 year, to the latest year with historical climate for every tile",
                          "(`climateData::latestHistoricalYear()`); climate is the last of the inputs to reach a year.")),
    defineParameter("flammabilityThreshold", "numeric", 0.1, 0, 1,
                    paste("Minimum proportion of flammable fine-resolution land cover for a `rasterToMatch`",
                          "pixel to be flammable. Only used when `rstLCCs` is built here.")),
    defineParameter("forestedLCC", "numeric", c(81, 210, 220, 230, 240), NA, NA,
                    paste("Forested land cover classes - these differ from non-forest because the biomass",
                          "and composition of fuels are taken into account by fireSense, while non-forest",
                          "classes are treated categorically")),
    defineParameter("igAggFactor", "numeric", 4, 1, NA,
                    "Aggregation factor (number of `rasterToMatch` cells per side) for the ignition and escape covariates."),
    defineParameter("fuelClassCol", "character", "FuelClass", NA, NA,
                    "the column in `sppEquiv` that defines unique fuel classes. A column ",
                    "named `FuelClass` exists in the `LandR::sppEquivalencies_CA` and will be used ",
                    "by default. To change the `FuelClass` classifications, add a column to that table, ",
                    "or to `sim$sppEquiv` and then modify this `fuelClassCol` parameter"),
    defineParameter("minBufferSize", "numeric", 5000, NA, NA,
                    paste("Minimum number of cells in each fire's burned-plus-buffer sample, applied after `areaMultiplier`.")),
    defineParameter("nonflammableLCC", "numeric", c(0, 20, 31, 32, 33), NA, NA,
                    "Non-flammable classes in `rstLCCs`; the default is water, snow/ice, rock and barren land."),
    defineParameter("nonForestCanBeYoungAge", "logical", TRUE, NA, NA,
                    paste("if TRUE, burned non-forest will be treated as `youngAge`. Recommended to be TRUE",
                          "as burned forest is often classified as non-forest")),
    defineParameter("sppEquivCol", "character", "LandR", NA, NA,
                    "column name in `sppEquiv` object that defines unique species in `cohortData`"),
    defineParameter("spreadFitGoogleDriveFolder", "character",
                    "https://drive.google.com/drive/folders/1X9-mRjyLMNpgkP_cfqhbr_AQEPOsVCHf",
                    NA, NA, paste("URL of the Google Drive folder holding the ledger of previous SpreadFit results",
                                  "(`spreadFitFilename`), read with `reproducible::CacheGeo`.")),
    defineParameter("spreadFitFilename", "character", "fireSenseParams.rds",
                    NA, NA, paste("Name of the ledger file in `spreadFitGoogleDriveFolder`: study area polygons with",
                                  "their fitted SpreadFit parameters.")),
    defineParameter("targetFuelClasses", "numeric", 5, 1, 7,
                    "the target number of unique fuel classes when using semi-automated approach"),
    defineParameter("useRasterizedFireForSpread", "logical", FALSE, NA, NA,
                    paste("Use `historicalFireRaster` in place of fire polygons for spread?",
                          "Not currently supported: `TRUE` stops with an error when preparing SpreadFit.")),
    defineParameter("whichModulesToPrepare", "character",
                    c("fireSense_IgnitionFit", "fireSense_SpreadFit", "fireSense_EscapeFit"),
                    NA, NA, "Which fireSense fit modules to prep? defaults to all 3"),
    defineParameter(".studyAreaName", "character", NULL, NA, NA,
                    "`studyArea` name used in file names and cache tags; `NULL` derives it from `sim$studyArea`."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    paste("Should this entire module be run with caching activated? This is intended",
                          "for data-type modules, where stochasticity and time are not relevant")),
    defineParameter(".useCacheArgs", "list",
                    list(dataPrepBuild = list(.cacheExtra = quote(list(
                      LandR::.compareRas, LandR::asInt, LandR::defineFlammable, LandR::isInt,
                      fireSenseUtils::assessFuelClasses, fireSenseUtils::fuelClassPrep,
                      fireSenseUtils::makeLandcoverDT, fireSenseUtils::makeTSD)))),
                    NA, NA,
                    paste("Extra `reproducible::Cache()` arguments, by event. A cached event's digest covers",
                          "this module's code but not the package functions it calls, so `dataPrepBuild`",
                          "passes those in `.cacheExtra`: a changed function then re-runs the event."))
  ),
  inputObjects = bindrows(
    expectsInput("climateVariables", "list", sourceURL = NA,
                 paste("Climate variable definitions for `climateData::prepClimateLayers` (canClimateData). Unless",
                       "supplied, `.inputObjects` builds them from `climateVariablesForFire`. Declared as an input",
                       "because `.inputObjects` sets it: a cached `.inputObjects` restores only declared inputs,",
                       "so without this a cache hit returned no `climateVariables` and canClimateData failed.")),
    expectsInput("climateVariablesForFire", "list", sourceURL = NA,
                 paste("List with elements `ignition` and `spread`, each a character vector of climate variable",
                       "names, with or without underscores (e.g., `CMD_sm` or `CMDsm`). IgnitionFit uses all of",
                       "`ignition`; SpreadFit uses `spread`. Default: `ignition = c('CMD', 'cumMDC', 'CMD_sm', 'CMD_sp')`,",
                       "`spread = 'auto'`: the `ignition` variable that best separates the study area's worst fire years",
                       "(see `spreadClimateSelection`). Unless supplied, `climateVariables` is built from these.")),
    expectsInput("cohortDatas", "list", sourceURL = NA,
                 paste("List of `cohortData` data.tables, one per `dataYears`, named `year<year>`.",
                       "If not supplied, built by running Biomass_borealDataPrep for each data year.")),
    expectsInput("spreadFirePoints", "list", sourceURL = NA,
                 paste("List of spatial points, one per fire year, named `year<year>`; each point is the ignition",
                       "location of one fire in `firePolys`. The default is the polygon centroids.")),
    expectsInput("spreadFitAdditionalColNames", "character",
                 desc = paste("Names of the ledger (`spreadFitPreRun`) columns to read.",
                              "The default is `fireSenseUtils::spreadFitAdditionalColNamesTxt`.")),
    expectsInput("firePolys", "list", sourceURL = NA,
                 paste("List of fire polygons, one per year in `fireYears`, named `year<year>`.",
                       "The default is the latest NBAC release.")),
    expectsInput("firePolysForAge", "list", sourceURL = NA,
                 paste("List of annual fire polygons used for time since disturbance; as `firePolys`, but",
                       "starting `cutoffForYoungAge` years before the first of `fireYears`.")),
    expectsInput("historicalFireRaster", "SpatRaster",
                 sourceURL = "https://opendata.nfis.org/downloads/forest_change/CA_Forest_Fire_1985-2020.zip",
                 paste("Optional raster of fire year, 1985-2020. If supplied it replaces `firePolysForAge` for time since",
                       "disturbance. Only downloaded when `useRasterizedFireForSpread = TRUE`.")),
    expectsInput("historicalClimateRasters", "list", sourceURL = NA,
                 paste("List of SpatRasters of historical climate, named by climate variable,",
                       "with layers named `year<year>`. Must be supplied, and must cover `fireYears`.")),
    expectsInput("ignitionFirePoints", "SpatVector", sourceURL = NA,
                 paste("Points of annual ignitions, of every fire size, with columns `YEAR` and `SIZE_HA`.",
                       "The default is the lightning- and natural-caused fires of the current NFDB release.")),
    expectsInput("missingLCCgroup", "character", sourceURL = NA,
                 paste("The `nonForestedLCCGroups` name given to forested pixels that are absent from `cohortData`.",
                       "The default is the first name; it is replaced if fuel classes are estimated.")),
    expectsInput("nonForestedLCCGroups", "list",
                 paste("Named list of non-forest land cover classes, e.g. `list(wetland = c(19, 23, 32))`.",
                       "Each group becomes a fuel covariate. The default is one group, `nf`, of every class that is",
                       "neither forested nor non-flammable; it is replaced if fuel classes are estimated.")),
    expectsInput("pixelGroupMaps", "list", sourceURL = NA,
                 "List of `pixelGroupMap` SpatRasters matching `cohortDatas`, named `year<year>`."),
    expectsInput("propFlammables", "list", sourceURL = NA,
                 paste("List of SpatRasters of the proportion of flammable land cover in a pixel, one per `dataYears`.",
                       "Built with `rstLCCs` when that is not supplied.")),
    expectsInput("rasterToMatch", "SpatRaster", sourceURL = NA,
                 "Template raster for `studyArea`. The default is 240 m, from SCANFI land cover."),
    expectsInput("rasterToMatchLarge", "SpatRaster", sourceURL = NA,
                 paste("Optional larger template. If supplied it defines `studyArea_biomassParam`; if not,",
                       "`.inputObjects` sets it to `rasterToMatch` for Biomass_speciesData. Declared as an input",
                       "because `.inputObjects` sets it (a cached `.inputObjects` restores only declared inputs).")),
    expectsInput("rasterToMatch_biomassParam", "SpatRaster", sourceURL = NA,
                 "Template raster for `studyArea_biomassParam`, passed to Biomass_borealDataPrep. Expected to ",
                 "cover at least `rasterToMatch` (formerly `rasterToMatchLarge`)."),
    expectsInput("rstLCCs", "list", sourceURL = NA,
                 paste("List of land cover SpatRasters, one per `dataYears`, named `year<year>`, on",
                       "`rasterToMatch_biomassParam`. The default is from `fireSenseUtils::makeFireSenseLCC`.")),
    expectsInput("sppEquiv", "data.table", sourceURL = NA,
                 "Table of LandR species equivalencies. The default is from `LandR::speciesInStudyArea`."),
    expectsInput("standAgeMaps", "list", sourceURL = NA,
                 "List of stand age SpatRasters, one per `dataYears`, named `year<year>`;",
                 " used to create `cohortDatas` and time since disturbance."),
    expectsInput("spreadFirePolys", "list", sourceURL = NA,
                  "Not needed from the user: `firePolys` in the CRS of `rasterToMatch`, declared as an input",
                  " because a later event modifies it."),
    expectsInput("studyArea", "SpatVector", sourceURL = NA,
                 "Study area for all data. Should be buffered to limit edge effects on fire spread."),
    expectsInput("studyArea_biomassParam", "SpatVector", sourceURL = NA,
                 "study area passed to Biomass_borealDataPrep for vegetation calibration")
  ),
  outputObjects = bindrows(
    createsOutput("climateVariablesForFire", "list",
                  paste("As the input. If a previous SpreadFit exists, `spread` becomes the climate variables",
                        "of that fit, and they are added to `ignition`.")),
    createsOutput("spreadFitPreRun", "data.frame",
                  desc = paste("Ledger rows of previous SpreadFit results that overlap `studyArea`, from `CacheGeo`:",
                               "a geometry column (convert with `sf::st_as_sf`), plus `polygonID` and the columns in",
                               "`spreadFitAdditionalColNames`. `NULL` if there is no previous fit.")),
    createsOutput("studyAreaWithSpreadParams", "sf",
                  desc = "Same as `spreadFitPreRun`; not created if there is no previous fit."),
    createsOutput("sppColorVect", "character",
                  desc = "Named vector of hex colours, one per species. ",
                  "Only created if a previous SpreadFit exists, from the `sppEquiv` stored with it."),
    createsOutput("sppNameVector", "character",
                 desc = paste("Sorted species names (`sppEquivCol`) from the `sppEquiv` stored with a previous SpreadFit;",
                              "only created if one exists.")),
    createsOutput("spreadClimateSelection", "data.table",
                  paste("With `climateVariablesForFire$spread = 'auto'`: each candidate's AUC for separating the",
                        "worst quarter of fire years (by area burned), its Spearman correlation, and which was chosen.")),
    createsOutput("climateVariables", "list",
                  paste("Climate variable definitions, as used by `climateData::prepClimateLayers` (canClimateData).",
                        "Unless supplied, built from `climateVariablesForFire` for `fireYears` (and projected years",
                        "unless canClimateData's `climateGCM` is 'NRV'). If a previous SpreadFit exists, the",
                        "variables of that fit are added.")),
    createsOutput("fireBufferedListDT", "list",
                  "list of data.tables with fire id, `pixelID`, and buffer status"),
    createsOutput("fuelClassTable", "data.table",
                  "Fuel class assigned to each tree species; only created if fuel classes are estimated."),
    createsOutput("spreadFirePolys", "list",
                  "List of annual fire polygons used for SpreadFit: larger than one pixel and matched to `spreadFirePoints`."),
    createsOutput("fireSense_annualSpreadFitCovariates", "list",
                  "List of data.tables, one per fire year, of `pixelID` and the spread climate covariates in the fire buffers."),
    createsOutput("fireSense_escapeCovariates", "data.table",
                  "ignition covariates with added column of escapes"),
    createsOutput("fireSense_escapeFormula", "character",
                  "formula for escape, using fuel classes and landcover, as character"),
    createsOutput("fireSense_ignitionCovariates", "data.table",
                  "table of aggregated ignition covariates with annual ignitions"),
    createsOutput("fireSense_nonAnnualSpreadFitCovariates", "list",
                  "List of data.tables, one per `dataYears`, of `pixelID` and the fuel covariates in the fire buffers."),
    createsOutput("fireSense_spreadFormula", "character",
                  "formula for spread, using climate and vegetation covariates, as character"),
    createsOutput("ignitionFirePoints", "SpatVector",
                  paste("The input, in the CRS of `rasterToMatch` and clipped to `studyArea`.")),
    createsOutput("ignitionFitRTM", "SpatRaster",
                  paste("Template raster with the resolution and extent of `fireSense_ignitionCovariates`.",
                        "Attributes: `nonNAs`, the number of rows in that table, and `meanForestB`,",
                        "the mean forest biomass per pixel.")),
    createsOutput("landcoverDTs", "list",
                  "List of data.tables, one per `dataYears`, of `pixelID` and a 0/1 column per non-forest group, for flammable pixels."),
    createsOutput("lightningMaps", "SpatRaster",
                  paste("A 4-layer SpatRaster of lightning: lightningDays, lightningDensity, positiveCG, positiveCGdensity")),
    createsOutput("missingLCCgroup", "character",
                  "As the input, or the estimated group if fuel classes are estimated."),
    createsOutput("nonForestedLCCGroups", "list",
                  "As the input, or the estimated groups if fuel classes are estimated."),
    createsOutput("nonForest_timeSinceDisturbances", "list",
                  paste("List of SpatRasters, one per `dataYears`, of years since disturbance in flammable pixels,",
                        "forested or not.")),
    createsOutput("rstLCCs", "list",
                  "The input, on `rasterToMatch`."),
    createsOutput("flammableRTMs", "list", "List of binary SpatRasters of flammable land cover on `rasterToMatch`, one per `dataYears`."),
    createsOutput("sppEquiv", "data.table", "sppEquiv table potentially modified with new or overwritten fuel class"),
    createsOutput("spreadFirePoints", "list",
                  paste("List of ignition points, one per fire year, for fires larger than one pixel,",
                        "harmonized with `spreadFirePolys`.")),
    
    # For fireSense_**Predict modules
    createsOutput("propFlammable", "SpatRaster",
                 "Last element of `propFlammables`."),
    createsOutput("standAgeMap", "SpatRaster", "Last element of `standAgeMaps`."),
    createsOutput("rstLCC_RTM", "SpatRaster",
                  paste0("Last element of the output `rstLCCs`, i.e., on `rasterToMatch`.")),
    createsOutput("rstLCC", "SpatRaster",
                  paste0("Last element of the input `rstLCCs`, i.e., on `rasterToMatch_biomassParam`.")),
    createsOutput("flammableRTM", "SpatRaster", "Last element of `flammableRTMs`."),
    createsOutput("landcoverDT", "data.table",
                  "Last element of `landcoverDTs`."),
    createsOutput("nonForest_timeSinceDisturbance", "SpatRaster", 
                  "Last element of `nonForest_timeSinceDisturbances`.")
    
  )
))

#' Event dispatcher
#'
#' @param sim a `simList`.
#' @param eventTime current simulation time.
#' @param eventType character, the event to run.
#' @return the `simList`, invisibly.
doEvent.fireSense_dataPrepFit = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      if (!all(P(sim)$whichModulesToPrepare %in%
               c("fireSense_SpreadFit", "fireSense_IgnitionFit", "fireSense_EscapeFit"))) {
        stop("unrecognized module to prepare - review parameter whichModulesToPrepare")
      }

      ## schedule future event(s)
      if ("fireSense_IgnitionFit" %in% P(sim)$whichModulesToPrepare)
        sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "prepIgnitionFitData", eventPriority = 1)
      if ("fireSense_EscapeFit" %in% P(sim)$whichModulesToPrepare)
        sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "prepEscapeFitData", eventPriority = 1)
      if ("fireSense_SpreadFit" %in% P(sim)$whichModulesToPrepare) {
        sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "prepSpreadFitData", eventPriority = 1)
      }

      sim <- Init(sim)
      ## Priority 0 runs dataPrepBuild straight after this event and ahead of every other module's
      ## init (those are at .first(), i.e. 1).
      sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "dataPrepBuild", eventPriority = 0)

      sim <- scheduleEvent(sim, end(sim), "fireSense_dataPrepFit", "plotAndMessage", eventPriority = 9)
      sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "cleanUp", eventPriority = 10)
    },
    dataPrepBuild = {
      sim <- dataPrepBuild(sim)
    },
    prepIgnitionFitData = {
      sim <- prepare_IgnitionFit(sim)
    },
    prepEscapeFitData = {
      sim <- prepare_EscapeFit(sim)
    },
    prepSpreadFitData = {
      sim <- prepare_SpreadFit(sim)
    },
    plotAndMessage = {
      sim <- plotAndMessage(sim)
    },
    cleanUp = {
      sim <- cleanUpMod(sim)
    },
    warning(paste("Undefined event type: \"", current(sim)[1, "eventType", with = FALSE],
                  "\' in module \'", current(sim)[1, "moduleName", with = FALSE], "\'", sep = ""))
  )
  return(invisible(sim))
}

#' Look for a previous SpreadFit for this study area
#'
#' Reads the SpreadFit ledger with `CacheGeo`. If it has rows for `sim$studyArea`, sets the species,
#' fuel and climate variable objects to those of that fit. Not cacheable: it asks Google Drive.
#'
#' @param sim a `simList`.
#' @return the `simList`, invisibly, with `mod$haveSpreadFit` and `mod$userSuppliedFuelObjs` set.
Init <- function(sim) {

  sa <- sim$studyArea
  if (inherits(sa, "SpatVector")) sa <- st_as_sf(sa)
  
  prepInputsFSURL <- SpaDES.core::paramCheckOtherMods(sim, "spreadFitGoogleDriveFolder")
  sim$spreadFitPreRun <- CacheGeo(cloudFolderID = Par$spreadFitGoogleDriveFolder,
                              targetFile = Par$spreadFitFilename, purge = 7,
                              domain = sa, action = "nothing", useCache = FALSE,
                              destinationPath = inputPath(sim), bufferOK = TRUE)
  mod$haveSpreadFit <- is(sim$spreadFitPreRun, "sf") || is(sim$spreadFitPreRun, "data.frame")
  if (mod$haveSpreadFit) {
    sim$studyAreaWithSpreadParams <- sim$spreadFitPreRun

    # remove the column called "params" ... this just allows for partial matching, with or without "s"
    #   in case somebody uses `parameters`, `param`, or `params`
    colNames <- setdiff(sim$spreadFitAdditionalColNames,
                        grep(value = TRUE, "param", sim$spreadFitAdditionalColNames))
    df <- ledgerColumns(sim$spreadFitPreRun, colNames)
    
    assignToSim <- outputObjects(sim)$objectName
    outs2 <- Map(nam = sim$spreadFitPreRun[[fireSenseUtils::polygonIDTxt]], ind = seq_len(NROW(df)), function(ind, nam) {
      dfList <- lapply(df, function(x) x[[ind]]) # Take only first entry, which if within 1 ELF, will be correct
      list2env(dfList, environment()) # nolint: vars numIterations objFunVal sppEquiv nonForestedLCCGroups missingLCCgroup geometry
      sim$sppNameVector <- unique(sppEquiv[[Par$sppEquivCol]])
      sppOuts <- sppHarmonize(sppEquiv, sim$sppNameVector, P(sim)$sppEquivCol, sppColorVect = NULL,
                              vegLeadingProportion = NULL, studyArea = sim$studyArea)
      list2env(sppOuts, envir = environment())  # nolint: vars sppEquiv sppNameVector sppEquivCol sppColorVect

      pars <- sim$spreadFitPreRun$params[[ind]]
      lpn <- fireSenseUtils::logisticParamNames
      allMatched <- sapply(lpn, function(lpn) all(lpn %in% colnames(pars)))
      numMatches <- sapply(lpn, function(lpn) sum(colnames(pars) %in% lpn))
      whLogistic <- which(allMatched & numMatches == max(numMatches))  
      
      FuelNames <- c(sim$sppEquiv[[P(sim)$fuelClassCol]], names(nonForestedLCCGroups))
      hasYoungAge <- youngAgeTxt %in% colnames(pars)
      if (isTRUE(hasYoungAge))
        FuelNames <- c(youngAgeTxt, FuelNames)
      
      ClimateNames <- setdiff(setdiff(colnames(pars), lpn[[whLogistic]]), FuelNames)
      
      allVars <- getFromNamespace(".allowedClimateVars", ns = asNamespace("climateData"))
      allVarsNoUnderscore <- gsub("_", "", allVars)
      whClimateVar <- which(allVarsNoUnderscore %in% ClimateNames)
      theseClimVars <- allVars[whClimateVar]
      theseClimVarsNoUnderscore <- allVarsNoUnderscore[whClimateVar]
      
      # Append the ones needed in FireSense_spreadFit to whatever was supplied by user
      if (is.null(sim$climateVariables)) {
        stop("You need to provide sim$climateVariables; perhaps try adding the PredictiveEcology/canClimateData module")
      }
      keepObjs <- intersect(ls(), assignToSim)
      append(mget(keepObjs), 
             list(theseClimVars = theseClimVars, theseClimVarsNoUnderscore = theseClimVarsNoUnderscore))
    })

    # Ok. With the multi-ELF reality, many objects must become lists --> the names of objects will become plural if we need them    
    keepCols <- setdiff(colnames(outs2[[1]]$sppEquiv), P(sim)$fuelClassCol)
    
    # sppEquiv
    sim$sppEquiv <- Map(o = outs2, function(o) o[["sppEquiv"]]) |>
      rbindlist(use.names = TRUE) |> 
      unique(by = keepCols) |> 
      setorderv(Par$sppEquivCol)
    sim$sppEquivs <- Map(o = outs2, function(o) o[["sppEquiv"]])
    sim$sppColorVects <- Map(o = outs2, function(o) o[["sppColorVect"]])
    sim$sppNameVectors <- Map(o = outs2, function(o) o[["sppNameVector"]]) 
    for (ind in seq_len(NROW(sim$spreadFitPreRun$sppEquiv))) {
      if (!isTRUE(all.equal(sim$spreadFitPreRun$sppEquiv[[ind]], sim$sppEquivs[[ind]], check.attributes = FALSE)  )) {
        sim$spreadFitPreRun$sppEquiv[[ind]] <- sim$sppEquivs[[ind]]
      }
    }
    
    # sppNameVector # not used in FS modules
    sim$sppNameVector <- Map(o = outs2, function(o) o[["sppNameVector"]]) |> 
      unlist() |> 
      unique() |> 
      sort()
    
    # sppColourVector
    scv <- sort(lapply(seq_along(outs2), function(x) outs2[[x]][["sppColorVect"]]) |> unlist())
    scvUnique <- data.frame(scv = scv, scn = names(scv)) |> unique() 
    scvUnique <- scvUnique[order(scvUnique$scn),]
    whMixed <- which(scvUnique$scn %in% "Mixed")
    scvUnique1 <- scvUnique[-whMixed,]
    scvUnique <- rbind(scvUnique1, scvUnique[whMixed,])
    dups <- duplicated(scvUnique$scn)
    if (any(dups)) {
      stop("sppColVect has unique colour values for a single species; this needs to be reworked and rethought")
    }
    sim$sppColorVect <- scvUnique$scv |> setNames(scvUnique$scn)
    
    sim$nonForestedLCCGroupsList <- Map(x = outs2, function(x) x[["nonForestedLCCGroups"]])
    sim$missingLCCgroupList <- Map(x = outs2, function(x) x[["missingLCCgroup"]])
    
    ## Each ELF's fit names its own spread climate variable(s): with spread = "auto" two ELFs can
    ## differ. Prepare the union; each ELF's own formula picks its variables from it. Variables already
    ## being prepared keep their definition (and years); only missing ones are added, like the rest.
    fitClimVars <- unique(unlist(Map(x = outs2, function(x) x[["theseClimVars"]]), use.names = FALSE))
    fitClimVarsNoUnderscore <- gsub("_", "", fitClimVars)
    gcm <- tryCatch(P(sim, module = "canClimateData")$climateGCM, error = function(e) NULL)
    projYears <- tryCatch(P(sim, module = "canClimateData")$projectedClimateYears, error = function(e) NULL)
    sim$climateVariables <- addFitClimateVariables(sim$climateVariables, fitClimVars,
                                                   historicalYears = P(sim)$fireYears,
                                                   projected = !identical(gcm, "NRV"),
                                                   projectedYears = if (is.null(projYears)) 2011:2100 else projYears)
    sim$climateVariablesForFire[["spread"]] <- fitClimVarsNoUnderscore
    ## ignition (xgboost) can use every one of them
    sim$climateVariablesForFire[["ignition"]] <-
      sort(unique(c(sim$climateVariablesForFire[["ignition"]], fitClimVarsNoUnderscore)))
  
  }
  
  ## dataPrepBuild is cached. Its digest covers expected inputs, this module's functions and
  ## mod$ contents, but not sim$.userSuppliedObjNames, so record what it needs from that here.
  fuelObjs <- c("nonForestedLCCGroups", "fuelClassTable")
  mod$userSuppliedFuelObjs <- stats::setNames(fuelObjs %in% sim$.userSuppliedObjNames, fuelObjs)
  return(invisible(sim))
}

#' Build land cover, flammability, fuel class and time-since-disturbance objects
#'
#' The expensive, cacheable part of initialization, as its own event. Runs whether or not a
#' previous SpreadFit exists; `mod$haveSpreadFit`, set by `Init()`, decides whether fuel classes
#' are assessed, and is in the event's digest.
#'
#' @param sim a `simList`.
#' @return the `simList`, invisibly.
dataPrepBuild <- function(sim) {
  sim$sppEquiv <- copy(sim$sppEquiv) # `:= NULL` below must not alter the user's table

  #because BBDP wants objects potentially larger than studyArea,
  #crop rstLCC and standAgeMap to create smaller objects before their derived objects
  # (landcoverDT/flammableMap and nonForest_timeSinceDisturbance, respectively).
  # This is at rasterToMatch_biomassParam which is needed by Biomass_borealDataPrep
  sim$rstLCC <- tail(sim$rstLCCs, 1)[[1]]
  
  objs <- sim$rstLCCs
  if (!LandR::.compareRas(sim$rasterToMatch, sim$rasterToMatch_biomassParam, stopOnError = FALSE)) {
    objs2 <- lapply(sim$rstLCCs, function(x) {
      postProcess(x, to = sim$rasterToMatch, method = "near")})
    objs <- objs2
  }
  if (!all(vapply(objs, isInt, logical(1)))) {
    objs <- Map(obj = objs, function(obj) LandR::asInt(obj))
  }
  # This makes rstLCCs same as sim$rasterToMatch instead of rasterToMatch_biomassParam
  sim$rstLCCs <- objs
  
  sim$flammableRTMs <- Map(dy = mod$dyChars, function(dy) {
    defineFlammable(sim$rstLCCs[[dy]],
                    nonFlammClasses = P(sim)$nonflammableLCC,
                    to = sim$rasterToMatch)
  })
  # recover the factors; the defineFlammable needed it to be integer; so this is after that
  sim$rstLCCs <- Map(r = sim$rstLCCs, function(r) {
    levels(r) <- terra::cats(sim$rstLCC)
    r
  })
  
  digFlammableRTMs <- .robustDigest(sim$flammableRTMs)

  # Create the "objects for prediction cases
  #this object is still at biomassParam size
  sim$standAgeMap <- tail(sim$standAgeMaps, 1)[[1]]
  sim$propFlammable <- tail(sim$propFlammables, 1)[[1]]
  
  # This is now RTM
  sim$rstLCC_RTM <- tail(sim$rstLCCs, 1)[[1]]
  sim$flammableRTM <- tail(sim$flammableRTMs, 1)[[1]]
  sim$landcoverDT <- tail(sim$landcoverDTs, 1)[[1]]
  
  ## ensure studyArea consists of a single polygon
  mod$studyAreaUnion <- if (is(sim$studyArea, "sf")) {
    sf::st_union(sim$studyArea)
  } else if (is(sim$studyArea, "SpatVector")) {
    terra::aggregate(sim$studyArea)
  } else {
    stop("studyArea must be either an sf or SpatVector object")
  }

  if (!terra::same.crs(mod$studyAreaUnion, sim$rasterToMatch)) {
    mod$studyAreaUnion <- projectTo(mod$studyAreaUnion, terra::crs(sim$rasterToMatch))
  }

  sim$ignitionFirePoints <- clipPointsToStudyArea(sim$ignitionFirePoints, mod$studyAreaUnion)

  ## possible, if user-supplied
  if (!terra::same.crs(sim$firePolys[[1]], sim$rasterToMatch)) {
    sim$spreadFirePolys <- Map(fp = sim$firePolys,
                               function(fp)
                                 projectTo(fp, st_crs(sim$rasterToMatch))) |>
      Cache(.functionName = "projectTo_for_firePolys")
  } else {
    sim$spreadFirePolys <- sim$firePolys
  }

  fires <- Reduce(rbind, sim$spreadFirePolys)
  #this must ensure landcover overrides species - it does not currently

  fuelObjs <- names(mod$userSuppliedFuelObjs)
  userSupplied <- unname(mod$userSuppliedFuelObjs)
  sppFCSupplied <- LandR::sppEquivalencies_CA[sim$sppEquiv, on = "LandR"]
  userSuppliedFC <- sppFCSupplied[, FuelClass == i.FuelClass]

  # Determine whether user has supplied their own FuelClass col in sppEquiv; their own nonForestedLCCGroups,
  #  their own fuelClassTable. If so, then don't estimate them here.
  needToEstimateFuelClasses <- P(sim)$estimateFuelClasses && all(userSupplied %in% FALSE) && all(userSuppliedFC %in% TRUE)

  ageGroups <- yearGroups(mod$dys, Par$fireYears)

  digRstLCC <- .robustDigest(sim$rstLCCs)
  digRTMs <- .robustDigest(sim$pixelGroupMaps)
  if (!mod$haveSpreadFit) {
    landscape <- Map(f = fuelClassPrep,
                     pixelGroupMap = sim$pixelGroupMaps,
                     cohortData = sim$cohortDatas,
                     rstLCC = sim$rstLCCs,
                     yearRange = ageGroups,
                     MoreArgs = list(nonflammableLCC = P(sim)$nonflammableLCC,
                                     fires = fires,
                                     nonforestLCC = sim$nonForestedLCCGroups)) |>
      Cache(.functionName = "fuelClassPrep", userTags = c("fireSenseDataPrepFit", "fuelClassPrep"),
            omitArgs = c("pixelGroupMap", "rstLCC"), .cacheExtra = list(digRTMs, digRstLCC))

    # Combine landscapes and finalize data
    landscape <- rbindlist(landscape)
    landscape <- landscape[, .(cell, speciesCode, lcc, B, totalBiomass, burned, year)]

    nonforestLCC <- unlist(sim$nonForestedLCCGroups) #unlist for backwards compatibility

    #set up plots
    speciesGG_DF <- landscape[!is.na(B), .(percentBurn = sum(burned)/.N * 100), .(speciesCode)]
    lccGG_DF <- landscape[is.na(B)]
    lccGG_DF[!lcc %in% nonforestLCC, newLCC := paste("not LandR forest", lcc)]
    lccGG_DF[is.na(newLCC), newLCC := lcc]
    lccGG_DF <- lccGG_DF[, .(percentBurn = sum(burned)/.N * 100), .(newLCC)]

    #TODO: add this to a Plots call?
    #It is a bar plot of the percent burned of each species and lcc

    if (needToEstimateFuelClasses) {

      fuelClassObjects <- assessFuelClasses(landscape = landscape,
                                            fuelCol = P(sim)$fuelClassCol,
                                            sppEquiv = sim$sppEquiv,
                                            sppEquivCol = P(sim)$sppEquivCol,
                                            targetFuelClasses = P(sim)$targetFuelClasses,
                                            nonforestLCC = nonforestLCC) |>
        Cache(userTags = c("assessFuelClasses", P(sim)$fuelClassCol))

      # sppEquiv
      temp <- fuelClassObjects$modSppEquiv[, .(species, assignedFuelClass)]
      setnames(temp, c(P(sim)$sppEquivCol, P(sim)$fuelClassCol))
      sim$sppEquiv[, P(sim)$fuelClassCol := NULL]
      sim$missingLCCgroup <- fuelClassObjects$missingLCCgroup

      sim$fuelClassTable <- fuelClassObjects$modSppEquiv
      sim$nonForestedLCCGroups <- fuelClassObjects$nonForestedLCCGroups
      message("Estimated fuel classes for this study area:")
      messageDF(sim$fuelClassTable)
      sim$sppEquiv <- sim$sppEquiv[temp, on = P(sim)$sppEquivCol]
    } else {
      us <- if (any(userSupplied %in% TRUE)) {
        paste(fuelObjs[userSupplied], collapse = ", ")
      } else {
        ""
      }
      usFC <- if (any(userSuppliedFC %in% TRUE)) {
        paste0("custom ", Par$fuelClassCol, " column in `sppEquiv`")
      } else {
        ""
      }
      message("User has supplied ", if (nzchar(us)) us, if (nzchar(usFC)) usFC, "; ",
              "Not running `assessFuelClasses` to determine nonforest fuel classes. Using:")
      print(sim$nonForestedLCCGroups)
    }
  } else {
    message("User is using pre-estimated SpreadFit parameters",
            "Not running `assessFuelClasses` to determine nonforest fuel classes. Using:")
    print(sim$nonForestedLCCGroupsList)
  }
  
  #make this object small (used by the landcoverDT and time-since-disturbance steps)
  standAgeMaps <- lapply(sim$standAgeMaps, reproducible::postProcess, to = sim$rasterToMatch)

  # One landcover table per year for THIS study area, keyed by pixelID only. A fit is
  # always single-ELF: it uses the fuel objects assessed or supplied above
  # (nonForestedLCCGroups, missingLCCgroup), never ledger geometry.
  sim$landcoverDTs <- Map(dy = names(sim$rstLCCs), function(dy) {
    ll <- makeLandcoverDT(rstLCC = sim$rstLCCs[[dy]],
                          flammableRTM = sim$flammableRTMs[[dy]],
                          forestedLCC = P(sim)$forestedLCC,
                          nonForestedLCCGroups = sim$nonForestedLCCGroups)
    correctMissingLCC(ll, sim[["pixelGroupMaps"]][[dy]], sim$missingLCCgroup)
  }) |>
    Cache(.functionName = "makeLandcoverDT",
          .cacheExtra = list(rstLCC = digRstLCC, flammableRTM = digFlammableRTMs,
                             rasterToMatchs = digRTMs,
                             sim$missingLCCgroup, P(sim)$forestedLCC, sim$nonForestedLCCGroups))

  ## TODO: this object is used to track annual youngAge of all pixels, forested or not
  ## so "nonForest" is a poor choice of name. It should not have values for non-flammable pixels.
  
  sim$nonForest_timeSinceDisturbances <- Map(dy = mod$dyChars, dyNum = mod$dys, function(dy, dyNum) {
    tsd <- makeTSD(
      year = dyNum,
      fireRaster = sim$historicalFireRaster, ## can be NULL
      firePolys = sim$firePolysForAge,
      standAgeMap = standAgeMaps[[dy]],
      lcc = sim$landcoverDTs[[dy]],
      cutoffForYoungAge = P(sim)$cutoffForYoungAge
    )
    tsd[sim$flammableRTMs[[dy]][] == 0] <- NA
    tsd
  })

  ## Until youngAge treatment is identical between spread and ignition, no point in prepping veg here
  ## Currently youngAge is resolved annually in spread, but only once in ignition
  ## e.g. if a pixel ignited in 2008, its youngAge status in ignition is still determined by whether it was 15 in 2010,
  ## but its youngAge status for spread is deterimined by whether standAge < 15 in 2008

  # Create the "objects for prediction cases
  sim$flammableRTM <- tail(sim$flammableRTMs, 1)[[1]]
  sim$landcoverDT <- tail(sim$landcoverDTs, 1)[[1]]
  sim$nonForest_timeSinceDisturbance <- tail(sim$nonForest_timeSinceDisturbances, 1)[[1]]
  return(invisible(sim))
}


#' Prepare the covariates and formula for fireSense_SpreadFit
#'
#' Buffers each fire, joins the buffers to the fuel covariates of their data year and to annual
#' climate, and splits the result into annual and non-annual tables.
#'
#' @param sim a `simList`.
#' @return the `simList`, invisibly, with `fireSense_annualSpreadFitCovariates`,
#'   `fireSense_nonAnnualSpreadFitCovariates`, `fireSense_spreadFormula`, `fireBufferedListDT`,
#'   `spreadFirePoints` and `spreadFirePolys`.
prepare_SpreadFit <- function(sim) {

  ## prep veg data ---------------------------------------------------------------------------------
  doAssertion <- getOption("fireSenseUtils.assertions", TRUE)

  ## sanity check the inputs
  lapply(sim$historicalClimateRasters, compareGeom, x = sim$rasterToMatch)

  ## when landcoverDT is included, as is the case here, non-forest pixels in cohortData are masked out
  ## this is necessary when LandR and fireSense have differing concepts of non-forest

  dig1 <- .robustDigest(list(sim$landcoverDTs, sim$flammableRTMs))
  dig1a <- .robustDigest(list(sim$cohortDatas, sim$pixelGroupMaps, sim$nonForest_timeSinceDisturbances))
  dig2 <- append(dig1, dig1a)
  
  # This adds youngAge
  vegData <- Map(f = fireSenseUtils:::fireSenseCovariatesCreate,
                        cohortData = sim$cohortDatas,
                        pixelGroupMap = sim$pixelGroupMaps,
                        flammableRTM = sim$flammableRTMs,
                        landcoverDT = sim$landcoverDTs,
                        nonForest_timeSinceDisturbance = sim$nonForest_timeSinceDisturbances,
                        MoreArgs = list(sppEquiv = sim$sppEquiv,
                                        sppEquivCol = P(sim)$sppEquivCol,
                                        fuelClassCol = P(sim)$fuelClassCol,
                                        cutoffForYoungAge = P(sim)$cutoffForYoungAge,
                                        missingLCCgroup = sim$missingLCCgroup,
                                        nonForestedLCCGroups = sim$nonForestedLCCGroups,
                                        nonForestCanBeYoungAge = P(sim)$nonForestCanBeYoungAge,
                                        studyAreaName = P(sim)$.studyAreaName
                        ) 
  ) |>
    Cache(.cacheExtra = dig2, 
          omitArgs = c("landcoverDT", "flammableRTM", "cohortData", "pixelGroupMap", "nonForest_timeSinceDisturbance"),
          .functionName = "spreadCovariatesCreate")
  # Add "year" column
  vegData <- Map(v = vegData, n = names(vegData), function(v, n) {
    v[, year := gsub("[[:alpha:]+]", "", n)]
    v
  })
  vegData <- rbindlist(vegData)
  
  # prep the fire data ####
  # sim$fireBufferedListDT is made in these functions
  if (P(sim)$useRasterizedFireForSpread) {
    sim <- prepare_SpreadFitFire_Raster(sim)
  } else {
    sim <- prepare_SpreadFitFire_Vector(sim)
  }

  ## join fire and veg data
  fireSenseVegData <- joinFireBuffersToVeg(sim$fireBufferedListDT, vegData, mod$allYears)

  rm(vegData)
  gc()

  ## TODO: lines from creation of vegData onwards should be reviewed. Seems redundant.
  ## TODO: should column ids be in vegData? currently pixels appear > 1 time, within a fire period
  #because they can be burned or unburned, or in >1 fire years that are < 10 years apart
  setnames(fireSenseVegData, "buffer", "burned")

  nonVegColnames <- c("pixelID", "burned", "ids", grep(fireSenseUtils::yearTxt, ignore.case = TRUE, colnames(fireSenseVegData), value = TRUE))
  vegCols <- setdiff(names(fireSenseVegData),
                     nonVegColnames)
  dropCols <- names(which(colSums(fireSenseVegData[, ..vegCols], na.rm = TRUE) == 0))

  ## spreadFit will fail if there are empty (all zero) columns
  if (length(dropCols) > 0) {
    message("Dropping column(s) from spreadFit covariate table: ",
            paste(dropCols, collapse = ", "))
    vegCols <- vegCols[!vegCols %in% dropCols]
    set(fireSenseVegData, NULL, dropCols, NULL)
  }

  if (isTRUE(doAssertion)) {
    ttt <- table(fireSenseVegData$burned)
    ratioZeroToOne <- ttt[1] / ttt[2]
    if (ratioZeroToOne < 5) {
      stop("The number of pixels in the fire buffers should be at least 5x the number of burned pixels\n",
           "Please create larger buffers around fires in fireBufferedListDT, e.g., via ",
           "fireSenseUtils::bufferToArea(..., areaMultiplier = multiplier)")
    }
  }

  ## spread = "auto": the climate variable that best separates this study area's bad fire years
  ## (R/fireClimateVariables.R)
  if (identical(sim$climateVariablesForFire$spread, "auto")) {
    fy <- paste0(fireSenseUtils::yearTxt, P(sim)$fireYears)
    burned <- setNames(numeric(length(fy)), fy)
    b <- vapply(sim$fireBufferedListDT, function(d) as.numeric(sum(d$buffer == 1)), numeric(1))
    b <- b[names(b) %in% fy]
    burned[names(b)] <- b
    sim$spreadClimateSelection <- selectSpreadClimateVariable(sim$historicalClimateRasters, burned,
                                                              candidates = sim$climateVariablesForFire$ignition)
    sim$climateVariablesForFire$spread <- sim$spreadClimateSelection[chosen == TRUE, var]
    message("Spread climate variable (auto): ", sim$climateVariablesForFire$spread,
            "; bad-fire-year AUC: ", paste0(sim$spreadClimateSelection$var, " ",
                                            round(sim$spreadClimateSelection$auc, 2), collapse = ", "))
  }

  RHS <- paste(paste0(sim$climateVariablesForFire$spread, collapse = " + "), youngAgeTxt,
               paste0(vegCols, collapse = " + "), sep =  " + ")

  ## this is a funny way to get years but avoids years with 0 fires
  allYears <- unname(unlist(mod$allYears))
  yearsWithFire <- allYears[allYears %in% names(sim$fireBufferedListDT)]

  #### prep climate data ####

  spreadClimate <- sim$historicalClimateRasters[sim$climateVariablesForFire$spread]
  #don't need climate data for years outside fire years
  spreadClimate <- lapply(spreadClimate, FUN = function(x){
    x <- terra::subset(x, subset = names(x) %in% paste0(fireSenseUtils::yearTxt, P(sim)$fireYears))
  })

  climateDT <- climateRasterToDataTable(historicalClimateRasters = spreadClimate) |>
    Cache(userTags = c("climateRasterToDataTable", names(spreadClimate)))

  fbl <- rbindlist(sim$fireBufferedListDT, idcol = fireSenseUtils::yearTxt)
  rmCols <- setdiff(colnames(fbl), c("pixelID", fireSenseUtils::yearTxt))
  set(fbl, NULL, rmCols, NULL)
  fbl <- climateDT[fbl, on = c(fireSenseUtils::yearTxt, "pixelID"), nomatch = NULL]

  fireSense_annualSpreadFitCovariates <- split(fbl, by = fireSenseUtils::yearTxt, keep.by = FALSE)

  ## prepare non-annual spread fit covariates by getting the youngAge
  # have to remove years that have no fires and so no climate data needed
  missingYears <- setdiff(unlist(mod$allYears), names(fireSense_annualSpreadFitCovariates))
  for (my in missingYears) # keep the colnames even though NROW is 0
    fireSense_annualSpreadFitCovariates[[my]] <-  fireSense_annualSpreadFitCovariates[[1]][0]
  fireSense_annualSpreadFitCovariates <- fireSense_annualSpreadFitCovariates[order(names(fireSense_annualSpreadFitCovariates))]

  colsToExtract <- c("pixelID", vegCols)

  ## group the annual covariates by data year
  annualCovariates <- Map(yrsChar = mod$allYears, function(yrsChar) {
    fireSense_annualSpreadFitCovariates[yrsChar]
  })

  ## youngAge is deliberately NOT added to annualCovariates: it is a non-annual covariate

  ## get rid of nonflammable pixels (here because the calcYoungAge function assigns ages to NA values,
  ## due to inconsistent treatment of non-forest age pixels  in kNN years and other products (0 vs NA)
  annualCovariates <- Map(ac = annualCovariates, yrNum = names(annualCovariates),
      function(ac,
               yrNum) {
        yrChar <- paste0(fireSenseUtils::yearTxt, yrNum)
    lapply(ac, function(x) {
      x[pixelID %in% sim$landcoverDTs[[yrChar]]$pixelID, ]
    })
  })

  # Confirm that YoungAge is the right amount. If there are "normal" amount of fires,
  #  then it should be < 0.2
  propYoungAge <- unlist(unlist(unname(
    lapply(annualCovariates, function(ac) 
      lapply(ac, function(x) if (!is.null(x$youngAge)) mean(x$youngAge)))), recursive = FALSE))
  lotsOfYA <- propYoungAge > 0.2
  if (any(lotsOfYA, na.rm = TRUE))
    warning("There are individual years with >20% of the pixels in YoungAge; confirm this is ",
            "expected... ", paste(names(propYoungAge)[lotsOfYA %in% TRUE], collapse = ", ")
            )


  if (!P(sim)$nonForestCanBeYoungAge) {
    ## TODO: test this inversion of makeMutuallyExclusive's regular use
    args <- as.list(rep(youngAgeTxt, length = length(sim$nonForestedLCCGroups)))
    names(args) <- names(sim$nonForestedLCCGroups)
  } else {
    ## this is done later in spreadFit - but done here for accuracy of outputs
    args <- list(names(sim$nonForestedLCCGroups)) |> setNames(youngAgeTxt)
  }

  annualCovariates <- lapply(annualCovariates, makeMutuallyExclusive, mutuallyExclusiveCols = args)

  sim$fireSense_annualSpreadFitCovariates <- do.call(c, unname(annualCovariates))

  # nonAnnuals are all the fuels; fire occurrence; youngAge and climate are annual (as of Jan 23, 2026)
  nonAnnuals <- Map(yrsChar = mod$allYears, yrGroup = names(mod$allYears), function(yrsChar, yrGroup) {
    yrsNum <- gsub("[^0-9]", "", yrsChar) |> as.integer()
    fireSenseVegData[year < max(yrsNum) & year >= as.integer(yrGroup), .SD, .SDcols = colsToExtract] %>%
      na.omit(.) %>%
      as.data.table(.) %>%
      .[!duplicated(pixelID), ]
  })
  sim$fireSense_nonAnnualSpreadFitCovariates <- nonAnnuals

  if (is.null(sim$fireSense_spreadFormula)) {
    sim$fireSense_spreadFormula <- paste0("~ 0 + ", RHS)
  }

  return(invisible(sim))
}

#' Fire buffers and ignition points from `historicalFireRaster`
#'
#' Not currently supported: stops with an error. The code after the `stop()` predates the
#' per-data-year `flammableRTMs` and `landcoverDTs` and must be revised before use.
#'
#' @param sim a `simList`.
#' @return the `simList`, invisibly, with `fireBufferedListDT` and `spreadFirePoints`.
prepare_SpreadFitFire_Raster <- function(sim) {
  stop("these methods need to be revised for the two flammable RTMs, two landcoverDTs")
  historicalFireRaster <- sim$historicalFireRaster

  ## build initial burn IDs by buffering  - then using clump(raster) or patches(terra)
  ## historical fire Raster is currently not in outputs - if assigned to sim here, it should be added
  ## as we modify it by removing non-flammable fires.

  historicalFireRaster <- mask(historicalFireRaster, sim$flammableRTM,
                               maskvalues = 0, updatevalue = NA)

  nCores <- ifelse(grepl("Windows", Sys.info()[["sysname"]]), 1L, length(Par$fireYears))

  ## this is analogous to buffer to area but for raster datasets as opposed to polygon
  ## the inner looping function is very similar - one difference is that non-flammable
  ## pixels do not count toward the buffer size, unlike the polygonal version.
  sim$fireBufferedListDT <- rasterFireBufferDT(years =  P(sim)$fireYears,
                                               fireRaster = historicalFireRaster,
                                               flammableRTM = sim$flammableRTM,
                                               bufferForFireRaster = P(sim)$bufferForFireRaster, verb = 1,
                                               areaMultiplier = eval(P(sim)$areaMultiplier),
                                               minSize = P(sim)$minBufferSize,
                                               cores = nCores) |>
    Cache(userTags = c(currentModule(sim), "rasterFireBufferDT"))
  ## TODO: test that this is the correct method for missing years
  missingYears <- unlist(lapply(sim$fireBufferedListDT, is.null))

  if (any(missingYears)) {
    actualFireYears <- P(sim)$fireYears[!missingYears]
    sim$fireBufferedListDT <- sim$fireBufferedListDT[!missingYears]
  }

  ## next up: generate spread fire points. no harmonization is needed with this approach :)
  sim$spreadFirePoints <- lapply(sim$fireBufferedListDT,
                                 rasterFireSpreadPoints,
                                 flammableRTM = sim$flammableRTM)

  ## TODO: this is temporary while we migrate out of spatial/raster constructs
  ## the fireSenseUtils::yearTxt prefix is added by fireSenseUtils::makeLociList - discuss what to do
  tempFun <- function(pts, year){
    pts$YEAR <- year
    return(pts)
  }

  sim$spreadFirePoints <- Map(pts = sim$spreadFirePoints,
                              year = P(sim)$fireYears[!missingYears], f = tempFun)
  names(sim$spreadFirePoints) <- names(sim$fireBufferedListDT)

  return(invisible(sim))
}

#' Fire buffers and ignition points from fire polygons
#'
#' Drops fires of one pixel or less, then, per data year, buffers the fires and harmonizes polygons
#' and points with `fireSenseUtils::harmonizeFireData`. Data years with no fires are dropped from
#' `mod$allYears`.
#'
#' @param sim a `simList`.
#' @return the `simList`, invisibly, with `fireBufferedListDT`, `spreadFirePoints` and
#'   `spreadFirePolys`.
prepare_SpreadFitFire_Vector <- function(sim) {

  # prep fire data -----------------------------------------------------------------
  IDvar <- grep("ID", names(sim$spreadFirePolys[[1]]), value = TRUE) |> setdiff("GID")
  if (is.null(sim$spreadFirePolys[[1]][[IDvar]])) {
    stop("firePolys needs a numeric an ID column, such as FIRE_ID, NFIREID")
  }

  if (!is.numeric(sim$spreadFirePolys[[1]][[IDvar]]) | !is.numeric(sim$spreadFirePoints[[1]][[IDvar]])) {

    message("need numeric FIRE_ID column in fire polygons and points. Coercing to numeric...")
    #this is true of the current NFBB
    origNames <- names(sim$spreadFirePolys)
    PointsAndPolys <- Map(year = names(sim$spreadFirePolys),
                          function(year, polys = sim$spreadFirePolys, points = sim$spreadFirePoints) {
                            polys <- polys[[year]]
                            points <- points[[year]]
                            ## ensure matching IDs
                            points <- points[points[[IDvar]] %in% polys[[IDvar]],]
                            polys <- polys[polys[[IDvar]] %in% points[[IDvar]],]
                            points[[IDvar]] <- as.numeric(as.factor(points[[IDvar]][[IDvar]]))
                            polys[[IDvar]] <- as.numeric(as.factor(polys[[IDvar]][[IDvar]]))
                            return(list(polys = polys, points = points))
                          })
    sim$spreadFirePoints <- Map(x = PointsAndPolys, function(x) x[["points"]])
    sim$spreadFirePolys <- Map(x = PointsAndPolys, function(x) x[["polys"]])
    rm(PointsAndPolys)
  }

  ## drop fires less than 1 px in size
  pixSizeHa <- prod(res(sim$flammableRTMs[[1]])) / 1e4
  ## using x[x$SIZE_HA] will work with terra or sf, while subset will not (I believe...)
  haColname <- grep("_HA$", names(sim$spreadFirePoints[[1]]), value = TRUE)[1] # has been SIZE_HA, POLY_HA

  sim$spreadFirePoints <- lapply(sim$spreadFirePoints, function(x, minSize = pixSizeHa) {
    x <- x[x[[haColname]] > minSize,]
    if (NROW(x) > 0) x else NULL
    x
  })

  sim$spreadFirePoints[sapply(sim$spreadFirePoints, is.null)] <- NULL ## silly R

  sim$spreadFirePolys <- lapply(sim$spreadFirePolys, function(x) {
    x <- x[x[[haColname]] > pixSizeHa,]
    if (nrow(x) > 0) x else NULL
  })
  sim$spreadFirePolys[sapply(sim$spreadFirePolys, is.null)] <- NULL
  sim$spreadFirePoints <- sim$spreadFirePoints[names(sim$spreadFirePoints) %in% names(sim$spreadFirePolys)]
  ## this covers when years are NA, which are caused by fire years with no available data

  ## years run separately because flammableRTM is different
  allYearsVect <- unlist(mod$allYears)
  nCores <- ifelse(grepl("Windows", Sys.info()[["sysname"]]), 1L,
                   sum(names(sim$spreadFirePolys) %in% allYearsVect))
  pointsIDcolumn <- grep("ID$", names(sim$spreadFirePolys[[1]]), value = TRUE)[1]

  # There may be entire decades with no fires
  yearsWeHave <- names(sim$spreadFirePolys)
  ays <- mod$allYears
  aysAll <- Map(ay = ays, function(ay) isTRUE(any(yearsWeHave %in% ay)))
  mod$allYears <- mod$allYears[unlist(aysAll)]
  
  harmonizeds <- Map(yr = mod$allYears, nam = names(mod$allYears),
                     function(yr, nam) {
    yrNam <- grep(nam, names(sim$flammableRTMs), value= TRUE)
    harmonizeFireData(
      firePolys = sim$spreadFirePolys[names(sim$spreadFirePolys) %in% yr], ## protects from missing years
      flammableRTM = sim$flammableRTMs[[yrNam]],
      spreadFirePoints = sim$spreadFirePoints[names(sim$spreadFirePoints) %in% yr], ## protects from missing years
      areaMultiplier = eval(P(sim)$areaMultiplier),
      minSize = P(sim)$minBufferSize,
      pointsIDcolumn = pointsIDcolumn,
      cores = nCores
    ) |>
      Cache(.functionName = paste0("harmonizedFireData_", yrNam),
      userTags = c("harmonizeFireData", P(sim)$.studyAreaName))
  }
  )

  sim$fireBufferedListDT <- Map(pp = harmonizeds, function(pp) pp[["fireBufferedListDT"]]) |>
    unname() |> unlist(recursive = FALSE)
  sim$spreadFirePoints <- Map(pp = harmonizeds, function(pp) pp[["spreadFirePoints"]]) |>
    unname() |> unlist(recursive = FALSE)
  sim$spreadFirePolys <- Map(pp = harmonizeds, function(pp) pp[["firePolys"]]) |>
    unname() |> unlist(recursive = FALSE)

  omitYears <- sapply(sim$spreadFirePoints, is.null)
  if (any(omitYears)) {
    sim$fireBufferedListDT[omitYears] <- NULL
    sim$spreadFirePolys[omitYears] <- NULL
    sim$spreadFirePoints[omitYears] <- NULL
  }

  ## drop fire years from these lists that don't have any buffer points post-harmonization
  sim$spreadFirePolys <- Map(year = names(sim$spreadFirePolys), function(year) {
    poly <- sim$spreadFirePolys[[year]]
    point <- sim$spreadFirePoints[[year]]
    poly <- poly[poly[[IDvar]] %in% point[[IDvar]],]
    return(poly)
  })
  return(invisible(sim))

}

#' Prepare the covariates for fireSense_IgnitionFit
#'
#' Aggregates fuel, climate and lightning covariates by `igAggFactor` and counts the ignitions in
#' each coarse pixel and year. Sets `mod$allYears`, which the spread preparation uses. No formula is
#' built: fireSense_IgnitionFit fits with xgboost, which does not use one.
#'
#' @param sim a `simList`.
#' @return the `simList`, invisibly, with `fireSense_ignitionCovariates`, `ignitionFitRTM` and
#'   `lightningMaps`.
prepare_IgnitionFit <- function(sim) {

  stopifnot(
    "all ignitionFirePoints are not within studyArea" = identical(
      nrow(st_as_sf(sim$ignitionFirePoints)),
      nrow(st_intersection(st_as_sf(sim$ignitionFirePoints), st_as_sf(mod$studyAreaUnion)))
    )
  )

  dig1 <- .robustDigest(list(sim$landcoverDTs, sim$flammableRTMs))
  dig1a <- .robustDigest(list(sim$cohortDatas, sim$pixelGroupMaps, sim$nonForest_timeSinceDisturbances))
  dig2 <- append(dig1, dig1a)

  #  Makes youngAge, amongst other things
  fuelCovsCoarse <- Map(
    f = function(..., fact = P(sim)$igAggFactor, rasTemplate = sim$flammableRTM) {
      prepare_FuelCovsCoarse(..., rasTemplate = rasTemplate, fact = fact)
    }, 
    cohortData = sim$cohortDatas,
    pixelGroupMap = sim$pixelGroupMaps,
    flammableRTM = sim$flammableRTMs,
    landcoverDT = sim$landcoverDTs,
    nonForest_timeSinceDisturbance = sim$nonForest_timeSinceDisturbances,
    MoreArgs = list(sppEquiv = sim$sppEquiv,
                    sppEquivCol = P(sim)$sppEquivCol,
                    fuelClassCol = P(sim)$fuelClassCol,
                    cutoffForYoungAge = P(sim)$cutoffForYoungAge,
                    missingLCCgroup = sim$missingLCCgroup,
                    nonForestedLCCGroups = sim$nonForestedLCCGroups,
                    nonForestCanBeYoungAge = P(sim)$nonForestCanBeYoungAge,
                    studyAreaName = P(sim)$.studyAreaName
    ))  |>
    Cache(.cacheExtra = list(prepare_FuelCovsCoarse = prepare_FuelCovsCoarse, dig2, fireSenseCovariatesCreate = fireSenseCovariatesCreate), # add the inner function
          omitArgs = c("landcoverDT", "flammableRTM", "cohortData", "pixelGroupMap", "nonForest_timeSinceDisturbance"),
          .functionName = "ignitionCovariatesCreate")


  # Climate data
  ignitionClimateCoarse <- prepare_ignitionClimate(
    ignitionClimateList = sim$historicalClimateRasters[sim$climateVariablesForFire$ignition],
    fact = P(sim)$igAggFactor,
    digest = dig2)

  sim$lightningMaps <- prepare_LightningData(sim$rasterToMatch, P(sim)$igAggFactor,
                                             dPath = inputPath(sim))

  do.call(compareGeom, unname(Reduce(append, list(sim$lightningMaps, ignitionClimateCoarse, fuelCovsCoarse))))

  ## ignition won't have same years as spread so we do not use names of init objects
  ## The reason is some years may have ignitions but no fires, e.g. 2010 in RIA
  years <- yearGroups(Par$dataYears, Par$fireYears, FALSE)
  years <- Map(y = years, function(y) paste0(fireSenseUtils::yearTxt, y))
  mod$allYears <- years
  allYearsVect <- unlist(mod$allYears)
  #assume that if multiple climate variables are present, they are of equal length
  #else bigger problems exist
  whAvailable <- allYearsVect %in% names(ignitionClimateCoarse[[1]])
  yearsNotAvailable <- allYearsVect[!whAvailable]

  ## Explicitly requested fire years must never be dropped quietly
  checkClimateYears(allYearsVect, whAvailable)

  ## join fuel class, LCC, and climate, subsetting to flamIndex, calculating n of ignitions

  sim$fireSense_ignitionCovariates <-
    mergePreparedCovs(years, fuelCovsCoarse, sim$ignitionFirePoints, sim$nonForestedLCCGroups,
                      ignitionClimateCoarse, sim$lightningMaps["lightningDays"], 
                      digest = append(dig2, list(P(sim)$igAggFactor)))

  ## make new ignition object, ignitionFitRTM
  sim$ignitionFitRTM <- rast(fuelCovsCoarse[[1]][[1]])
  sim$ignitionFitRTM <- setValues(sim$ignitionFitRTM, 1) ## avoids a warning
  attributes(sim$ignitionFitRTM)$nonNAs <- nrow(sim$fireSense_ignitionCovariates)

  ## assign mean forest biomass- for use in plotting in ignitionFit
  tempCD <- LandR::addPixels2CohortData(tail(sim$cohortDatas, 1)[[1]], tail(sim$pixelGroupMaps, 1)[[1]])
  bPerPixel <- tempCD[age > 0, .(bPerPixel = sum(B)), .(pixelIndex)]
  meanForestB <- mean(bPerPixel$bPerPixel)
  attr(sim$ignitionFitRTM, "meanForestB") <- meanForestB
  rm(tempCD, bPerPixel)

  return(invisible(sim))
}

#' Prepare the covariates and formula for fireSense_EscapeFit
#'
#' Adds to the ignition covariates the number of escapes: ignitions that grew beyond one
#' `flammableRTMs` pixel. Needs `prepare_IgnitionFit()` to have run.
#'
#' @param sim a `simList`.
#' @return the `simList`, invisibly, with `fireSense_escapeCovariates` and `fireSense_escapeFormula`.
prepare_EscapeFit <- function(sim) {

  if (is.null(sim$fireSense_ignitionCovariates)) {
    ## the datasets are essentially the same, with one column difference
    stop("Please include ignitionFit in parameter 'whichModulesToPrepare' if running EscapeFit")
  }

  escapeThreshHa <- prod(res(sim$flammableRTMs[[1]])) / 10000
  escapes <- sim$ignitionFirePoints[sim$ignitionFirePoints$SIZE_HA > escapeThreshHa, ]

  ## make a template aggregated raster - values are irrelevant, only need pixelID
  aggregatedRas <- terra::aggregate(sim$historicalClimateRasters[[1]][[1]],
                                    fact = P(sim)$igAggFactor, fun = mean) |>
    Cache(.functionName = "aggregate_historicalClimateRasters_forTemplate")

  if (is(escapes, "SpatVector")) {
    escapesOrig <- escapes
    coords <- pointCoords(escapesOrig)
  } else {
    escapes <- sf::st_as_sf(escapes)
    coords <- pointCoords(escapes)
  }
  escapeCells <- cellFromXY(aggregatedRas, coords)
  escapeDT <- as.data.table(escapes)
  setnames(escapeDT, "YEAR", fireSenseUtils::yearTxt)
  escapeDT[, pixelID := escapeCells]
  escapeDT <- escapeDT[, .(year, pixelID)]
  escapeDT <- escapeDT[, .(escapes = .N), .(year, pixelID)]
  escapeDT[, year := as.numeric(year)]
  escapeDT <- escapeDT[sim$fireSense_ignitionCovariates, on = c("pixelID", fireSenseUtils::yearTxt)]
  escapeDT[is.na(escapes), escapes := 0]

  escapeVars <- names(escapeDT)[!names(escapeDT) %in% c(fireSenseUtils::yearTxt, "pixelID", "escapes",
                                                        sim$climateVariablesForFire$ignition,
                                                        "ignitions", ranEffsLabel)]

  interactionsDF <- as.data.table(expand.grid(escapeVars, sim$climateVariablesForFire$ignition))
  interactionsDF[, interaction := do.call(paste, c(.SD, sep = ":")), .SDcols = names(interactionsDF)]
  interactions <- interactionsDF$interaction

  ## sanity check for base::abbreviate
  if (!length(unique(interactions)) == length(escapeVars) * length(sim$climateVariablesForFire$ignition)) {
    warning("automated escape formula construction needs review")
  }
  if (is.null(sim$fireSense_escapeFormula)) {
    sim$fireSense_escapeFormula <- paste0("cbind(escapes, ignitions - escapes) ~ ",
                                          paste0("(1|", ranEffsLabel, ")"), " + ",
                                          paste0(interactions, collapse = " + "))
  }

  if (any(escapeDT$escapes > escapeDT$ignitions)) {
    stop("issue with escapes outnumbering ignitions in a pixel - contact module creators")
  }
  sim$fireSense_escapeCovariates <- escapeDT

  return(invisible(sim))
}

#' Free memory held in `mod`
#'
#' @param sim a `simList`.
#' @return the `simList`, invisibly.
cleanUpMod <- function(sim) {
  mod$firePolysForAge <- NULL
  mod$fireSenseVegData <- NULL

  return(invisible(sim))
}

#' Placeholder for the `plotAndMessage` event; does nothing
#'
#' @param sim a `simList`.
#' @return the `simList`, invisibly.
plotAndMessage <- function(sim) {
  ## TODO: this could plot the ignition/spread covariates
  return(invisible(sim))
}

#' Build `cohortData`, `pixelGroupMap` and `standAgeMap` for each data year
#'
#' Runs a nested `simInitAndSpades()` of Biomass_borealDataPrep (and Biomass_speciesData, if it is
#' in the project) once per `dataYears`, cached on the inputs and on those modules' code. Missing
#' modules are downloaded to this module's `submodules` folder.
#'
#' @param sim a `simList`.
#' @return the `simList`, with `cohortData<year>`, `pixelGroupMap<year>` and `standAgeMap<year>`
#'   for each data year.
runBorealDP_forCohortData <- function(sim) {
  ## Biomass_species should be only run if it is already used in project
  ## TODO: this should really be specified by the user by a param, which would allow additional mods to be run
  modules <- modules(sim)
  modulesInProject <- list.dirs(modulePath(sim), full.names = TRUE, recursive = FALSE) |> as.list()
  names(modulesInProject) <- modulesInProject
  modulesInProject <- lapply(modulesInProject, basename)
  modules <- modifyList(modules, modulesInProject)

  neededModule <- "Biomass_borealDataPrep"
  if ("Biomass_speciesData" %in% modules) {
    neededModule <- c("Biomass_borealDataPrep", "Biomass_speciesData")
  }

  pathsLocal <- paths(sim)
  if (any(!neededModule %in% modules)) {
    ## NOTE: don't install pkgs mid-stream; use module metadata to declare pkgs for installation
    modulePathLocal <- file.path(modulePath(sim), currentModule(sim), "submodules")
    getModule(paste0("PredictiveEcology/", neededModule, "@development"),
              modulePath = modulePathLocal, overwrite = FALSE)
    pathsLocal$modulePath <- modulePathLocal
  }
  neededYears <- Par$dataYears

  ecoFile <- ifelse(is.null(sim$ecoregionRst), "ecoregionLayer", "ecoregionRst")
  objsNeeded <- c(ecoFile,
                  "firePerimeters",
                  "rasterToMatch", "studyArea",
                  "rstLCCs",
                  "standAgeMaps",
                  "studyArea_biomassParam", "rasterToMatch_biomassParam", #needed by BBDP
                  "species", "speciesTable", "sppEquiv")
  objsNeeded <- intersect(ls(sim), objsNeeded)
  objsNeeded <- mget(objsNeeded, envir = envir(sim))
  ## simInit applies `objects` after the modules' .inputObjects, so a NULL passed here would
  ## replace what those modules build (e.g. Biomass_speciesData's sppEquiv) with NULL.
  objsNeeded <- objsNeeded[!vapply(objsNeeded, is.null, logical(1))]
  cds <- Map(ny = neededYears, function(ny, objs = objsNeeded) {
    messageColoured(colour = "yellow", "Running Biomass_borealDataPrep for year ", ny)
    messageColoured(colour = "yellow", "  inside fireSense_dataPrepFit to estimate cohortData", ny)
    yrChar <- names(sim$rstLCCs)
    yrChar <- grep(ny, yrChar, value = TRUE)
    rstLCC <- objs[["rstLCCs"]][[yrChar]]
    standAgeMap <- objs[["standAgeMaps"]][[yrChar]]
    objs[["standAgeMaps"]] <- NULL
    objs <- c(objs, 
              "rstLCC" = rstLCC, 
              "standAgeMap" = standAgeMap)
    #now duplicated
    objs[[yrChar]] <- NULL
    objs[["rstLCCs"]] <- NULL
    parms <- list()
    for (nm in neededModule) {
      parms[[nm]] <- P(sim, module = nm)
      parms[[nm]][["dataYear"]] <- ny
      parms[[nm]][["forestedLCCClasses"]] <- P(sim)$forestedLCC
    }
    parms$Biomass_borealDataPrep$exportModels <- "none"

    # Digest the source code of modules; in case they change
    sourceCodeDig <- moduleCodeDigest(pathsLocal$modulePath, neededModule)

    outNY <- SpaDES.core::simInitAndSpades(paths = pathsLocal,
                                           params = parms,
                                           times = list(start = ny, end = ny),
                                           modules = neededModule,
                                           objects = objs) |>
      Cache(.functionName = paste0("simInitAndSpades_insideFireSenseDataPrepFit", ny),
            omitArgs = "paths",  # paths includes temp files that are always different
            .cacheExtra = sourceCodeDig)
    cohDatObj <- paste0(cohDat, ny)
    pixGrpMap <- paste0(pixGM, ny)
    saObj <- paste0(saMap, ny)
    outNY[[cohDatObj]] <- outNY[[cohDat]]
    outNY[[pixGrpMap]] <- outNY[[pixGM]]
    outNY[[saObj]] <- outNY[[saMap]]
    mget(c(cohDatObj, pixGrpMap, saObj), envir = envir(outNY))
  })
  lapply(cds, function(cd) list2env(cd, envir = envir(sim))) # nolint: vars cohortDatas pixelGroupMaps standAgeMaps
  sim
}

#' Default inputs
#'
#' Builds whatever is not supplied: study areas and templates, `sppEquiv`, land cover and stand age
#' per data year, `cohortDatas` and `pixelGroupMaps` (a nested simulation, see
#' `runBorealDP_forCohortData()`), NBAC fire polygons, NFDB ignition points and the non-forest
#' land cover groups. `historicalClimateRasters` must be supplied.
#'
#' @param sim a `simList`.
#' @return the `simList`, invisibly.
.inputObjects <- function(sim) {
  if (!suppliedElsewhere("studyArea", sim)) {
    sim$studyArea <- LandR::randomStudyArea(size = 10000 * 6.25 * 20000)
  }

  mod$dys <- P(sim)$dataYears
  mod$dataYears <- P(sim)$dataYears

  if (!suppliedElsewhere("studyArea_biomassParam", sim)) {
    if (is.null(sim$studyAreaLarge)) {
      if (!is.null(sim$rasterToMatchLarge)) {
        message("Creating studyArea_biomassParam from sim$rasterToMatchLarge")
        sim$studyArea_biomassParam <- terra::as.polygons(sim$rasterToMatchLarge) |> terra::aggregate()
      } else {
        warning("studyArea_biomassParam is not supplied; using studyArea, which is likely wrong as ",
                "several modules in the Biomass_** family expect studyArea_biomassParam to be from ",
                "the ***Large family of rasterToMatchLarge/studyAreaLarge")
        sim$studyArea_biomassParam <- sim$studyArea
      }
    } else {
      warning("please replace studyAreaLarge with studyArea_biomassParam")
      sim$studyArea_biomassParam <- sim$studyAreaLarge
    }
  }

  ## suppliedElsewhere() is TRUE whenever another module declares sppEquiv as an output, even
  ## when that module has not run yet and sim$sppEquiv is still NULL, so also build it when there
  ## is NO table. A table with zero rows is not "no table": fireSense_ELFs supplies one for an ELF
  ## with no tree species, and that must not be rebuilt here.
  if (!suppliedElsewhere("sppEquiv", sim, where = c("user", "initEvent")) || is.null(sim$sppEquiv)) {
    ## the same table fireSense_ELFs uses: no _Spp genus entries, only species with LANDIS
    ## traits, Engelmann spruce merged into Pice_eng
    sim$sppEquiv <- LandR::speciesInStudyArea(studyArea = sim$studyArea, sppEquivCol = Par$sppEquivCol,
                                              dPath = inputPath(sim))$sppEquiv
  }

  SpaDES.core::paramCheckOtherMods(sim, paramToCheck = "sppEquivCol")

  if (is.null(P(sim)$.studyAreaName)) {
    P(sim)$.studyAreaName <- studyAreaName(sim$studyArea)
  }
  cacheTags <- c(currentModule(sim), P(sim)$.studyAreaName)
  dPath <- asPath(inputPath(sim), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  if (!suppliedElsewhere("rasterToMatch", sim)) {
    sim$rasterToMatch <-
      {
        LandR::prepInputs_SCANFI_LCC_FAO(year = 2020,
                                         destinationPath = dPath,
                                         to = sim$studyArea) |>
          terra::aggregate(fact = 240/30)
      } |>
      Cache(.functionName = "prepInputs_rasterToMatch")
  }

  if (!suppliedElsewhere("rasterToMatch_biomassParam", sim)) {
    if (!is.null(sim$rasterToMatchLarge)) {
      warning("please use rasterToMatch_biomassParam in place of rasterToMatchLarge")
      sim$rasterToMatch_biomassParam <- sim$rasterToMatchLarge
    } else {
      sim$rasterToMatch_biomassParam <- sim$rasterToMatch
      sim$rasterToMatchLarge <- sim$rasterToMatch # needed for Biomass_speciesData
    }
  }

  ## Climate variables for the fire models (R/fireClimateVariables.R): the default here, unless supplied.
  ## `climateVariables` follows from them, so canClimateData prepares exactly these layers.
  if (!suppliedElsewhere("climateVariablesForFire", sim, where = c("sim", "user"))) {
    sim$climateVariablesForFire <- defaultClimateVariablesForFire
  }
  if (!suppliedElsewhere("climateVariables", sim, where = c("sim", "user"))) {
    gcm <- tryCatch(P(sim, module = "canClimateData")$climateGCM, error = function(e) NULL)
    projYears <- tryCatch(P(sim, module = "canClimateData")$projectedClimateYears, error = function(e) NULL)
    sim$climateVariables <- fireClimateLayers(sim$climateVariablesForFire, historicalYears = P(sim)$fireYears,
                                              projected = !identical(gcm, "NRV"),
                                              projectedYears = if (is.null(projYears)) 2011:2100 else projYears)
  }
  ## the rest of the module names climate layers without underscores ("CMDsm")
  sim$climateVariablesForFire <- lapply(sim$climateVariablesForFire, function(v) gsub("_", "", v))

  doRstLCCs <- !suppliedElsewhere("rstLCCs", sim)
  doStandAgeMaps <- !suppliedElsewhere("standAgeMaps", sim)

  dyChars <- paste0(fireSenseUtils::yearTxt, mod$dys)
  outs <- Map(dyChar = dyChars, dy = mod$dys, function(dyChar, dy) {
    if (doRstLCCs) {
      #use a threshold to to assign non-flammable cover (e.g. if < 10% flammable cover)
      opts11 <- options(reproducible.prepInputsUrlTiles = FALSE)
      on.exit(options(opts11))
      LCC <- fireSenseUtils::makeFireSenseLCC(
        neededYear = dy,
        writeTo = .suffix("rstLCC.tif",
                          paste0(dy, "_", P(sim)$.studyAreaName)),
        destinationPath = inputPath(sim),
        maskTo = sim$studyArea_biomassParam,
        to = sim$rasterToMatch_biomassParam,
        nonflammableLCC = P(sim)$nonflammableLCC,
        flammabilityThreshold = P(sim)$flammabilityThreshold) |>
        Cache(userTags = c("makeFireSenseLCC", dy),
              .functionName = paste0("makeFireSenseLCC", dy),
              ## Cache digests only makeFireSenseLCC's own code, so the functions it CALLS must be named
              ## here. Which ones depends on lccSource, a run-time option, so ask rather than list.
              .cacheExtra = fireSenseUtils::makeFireSenseLCCDeps())
    }
    
    if (doStandAgeMaps) {
      standAgeMap <- prepInputsStandAgeMap(rasterToMatch = sim$rasterToMatch_biomassParam,
                                           destinationPath = dPath,
                                           dataYear = dy) |>
        Cache(.functionName = paste0("prepInputsStandAgeMap", dy),
              userTags = c(cacheTags, "prepInputsStandAgeMap"),
              ## Cache digests only prepInputsStandAgeMap's own code; it adjusts ages in fires with this function
              .cacheExtra = list(LandR::replaceAgeInFires))
    } else {
      standAgeMap <- sim$standAgeMaps[[dyChar]]
    }
    
    ## land cover only when it was built here: a supplied rstLCCs (and its propFlammables) stays as it is
    out <- list(standAgeMaps = standAgeMap)
    if (doRstLCCs) {
      out$rstLCCs <- LCC$lcc
      out$propFlammables <- LCC$flammableProp
    }
    return(out)
  })
  outsRev <- Require::invertList(outs)
  list2env(outsRev, envir = envir(sim)) # nolint: vars standAgeMaps propFlammables rstLCCs
  
  mod$dyChars <- sapply(mod$dys, function(dy) grep(dy, names(sim$rstLCCs), value = TRUE))
  
  if (!all(suppliedElsewhere("cohortDatas", sim),
           suppliedElsewhere("pixelGroupMaps", sim)
  )) {
    ## This runs simInitAndSpades if needed
    sim <- runBorealDP_forCohortData(sim)
    objsHere <- c(cohDat, pixGM)
    for (nam in objsHere) {
      namPlural <- paste0(nam, "s")
      # This next line puts them in the ascending order
      cdnames <- sapply(mod$dys, function(dy) grep(paste0(nam, dy), names(sim), value = TRUE))
      sim[[namPlural]] <- mget(cdnames, envir(sim)) # nolint: vars cohortDatas pixelGroupMaps # nolint: unresolved_accessor
      names(sim[[namPlural]]) <- mod$dyChars # nolint: unresolved_accessor
      rm(list = cdnames, envir = envir(sim))
    }
  }

  if (!P(sim)$useRasterizedFireForSpread) {
    if (!suppliedElsewhere("firePolys", sim) | !suppliedElsewhere("firePolysForAge", sim)) {
      ## don't want to needlessly postProcess the same firePolys objects

      fireYears <- c(min(P(sim)$fireYears - P(sim)$cutoffForYoungAge):max(P(sim)$fireYears))
      ## the newest NBAC release; the shapefile's name carries the release, so it keys the Cache
      nbacShp <- fireRecordShapefile(fireSenseUtils::latestNBACUrl(), destinationPath = dPath)
      allFirePolys <- firePolysByYear(
        shp = nbacShp,
        years = fireYears,
        studyArea = postProcessTo(sim$studyArea, projectTo = sim$rasterToMatch)) |>
        Cache(omitArgs = "shp", .cacheExtra = basename(nbacShp),
              userTags = c(cacheTags, "firePolys", paste0(fireYears, collapse = ":")))
    }
    if (anyPlotting(Par$.plots)) {
      fp <- allFirePolys[!sapply(allFirePolys, is.null)]
      yrsDone <- names(fp)
      r1 <- {
        Map(nam = yrsDone, function(nam) {
          r <- sim$rasterToMatch
          wh <- r[] > 0
          # force the scale on the individual layers to be 1:2; must do the rtm colors larger so that they always plot
          r[wh] <- 2
          p <- rasterize(allFirePolys[[nam]], sim$rasterToMatch)
          r[p[] > 0] <- p[p[] > 0 ]
          r
        }) |>
          terra::rast() |>
          Plots(
            .plotInitialTime = NULL, # needs to be NULL or else won't do it
            maxnl = length(fp),
            legend = FALSE,
            maxcell = 1e6, zlim = c(1, 2), col = c("red", "grey"),
            deviceArgs = list(width = 11, height = 8, units = "in", res = 300),
            filename = "Historical Fire Maps") } |>
        Cache(.cacheExtra = attr(allFirePolys, "tags"),
              omitArgs = "data",
              .functionName = "Plots_fireMaps") # uses the cacheId of the firePolysByYear; only plot if changed
    }

    if (!suppliedElsewhere("firePolys", sim)) {
      sim$firePolys <- allFirePolys[names(allFirePolys) %in% paste0(fireSenseUtils::yearTxt, P(sim)$fireYears)]
      if (sum(lengths(sim$firePolys)) == 0) {
        stop("There are no fires in this study area during these years:\n",
                paste(paste0(fireSenseUtils::yearTxt, P(sim)$fireYears), collapse = ", "))
      }
    }

    if (!suppliedElsewhere("firePolysForAge", sim)) {
      sim$firePolysForAge <- allFirePolys
    }

    if (!suppliedElsewhere("spreadFirePoints", sim)) {
      message("... preparing polyCentroids")
      centerFun <- function(x) {
        if (is.null(x)) {
          return(NULL)
        } else {
          cent <- terra::centroids(x)
          return(cent)
        }
      }

      sim$spreadFirePoints <- suppressWarnings(lapply(sim$firePolys, centerFun))

      names(sim$spreadFirePoints) <- names(sim$firePolys)
    }

    if (all(!is.null(sim$spreadFirePoints), !is.null(sim$firePolys))) {
      ## may be NULL if passed by objects - add to Init?
      ## this is necessary because centroids may be fewer than fires if fire polys were small
      min1Fire <- lapply(sim$spreadFirePoints, length) > 0
      sim$spreadFirePoints <- sim$spreadFirePoints[min1Fire]
      sim$firePolys <- sim$firePolys[min1Fire]
    }

    if (length(sim$firePolys) != length(sim$spreadFirePoints)) {
      stop("mismatched years between firePolys and firePoints")
      ## TODO: need to implement a better approach that matches each year's IDS
      ## these are mostly edge cases if a user passes only one of spreadFirePoints/firePolys
    }
  }

  if (!suppliedElsewhere("ignitionFirePoints", sim)) {
    ## the URL is the same for every NFDB release; the shapefile's name carries the release, so it keys the Cache
    nfdbShp <- fireRecordShapefile(
      "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_pnt/current_version/NFDB_point_shp.zip",
      destinationPath = dPath)
    ignitionFirePoints <- nfdbFirePoints(
      shp = nfdbShp,
      years = P(sim)$fireYears,
      studyArea = sim$studyArea) |>
      Cache(omitArgs = "shp", .cacheExtra = basename(nfdbShp),
            userTags = c("ignitionFirePoints", P(sim)$.studyAreaName)) |>
      postProcessTo(projectTo = sim$rasterToMatch)
    sim$ignitionFirePoints <- ignitionFirePoints[ignitionFirePoints$CAUSE %in% c("L", "N"),]
    if (nrow(sim$ignitionFirePoints) == 0) {
      stop("no lightning- or natural-caused (CAUSE L or N) NFDB fire points in the study area during fireYears")
    }
  }

  if (!suppliedElsewhere("historicalClimateRasters", sim)) {
    stop("please supply sim$historicalClimateRasters")
  }

  if (P(sim)$useRasterizedFireForSpread) {
    if (!suppliedElsewhere("historicalFireRaster", sim)) {
      sim$historicalFireRaster <- prepInputs(url = extractURL("historicalFireRaster", sim),
                                             rasterToMatch = sim$rasterToMatch,
                                             destinationPath = dPath,
                                             studyArea = sim$studyArea,
                                             method = "near", ## only use near or ngb; bilinear is wrong!
                                             filename2 = paste0("wildfire_", P(sim)$.studyAreaName, ".tif")) |>
        Cache(userTags = c("historicalFireRaster", P(sim)$.studyAreaName))
    }
  }

  if (!suppliedElsewhere("nonForestedLCCGroups", sim)) {
    # Compare in CODE space: `nonflammableLCC` and `forestedLCC` are land-cover
    # codes, but freq() on a categorical raster (makeFireSenseLCC() attaches the
    # SCANFI levels) returns the LABELS, so drop the levels first.
    rstLCC <- terra::rast(sim$rstLCCs)
    if (any(terra::is.factor(rstLCC))) {
      rstLCC <- terra::deepcopy(rstLCC)
      levels(rstLCC) <- NULL
    }
    vals <- freq(rstLCC)
    forestOrNonFlamm <- sort(unique(c(Par$nonflammableLCC, Par$forestedLCC)))
    sim$nonForestedLCCGroups <- list(nf = sort(setdiff(unique(vals$value), forestOrNonFlamm)))
    ## TODO: consider moving this to init - and checking if unsupplied
  }

  if (!suppliedElsewhere("missingLCCgroup", sim)) {
    sim$missingLCCgroup <- names(sim$nonForestedLCCGroups)[1]
  }

  if (!suppliedElsewhere("spreadFitAdditionalColNames")) {
    sim$spreadFitAdditionalColNames <- fireSenseUtils::spreadFitAdditionalColNamesTxt
  }

  return(invisible(sim))
}

youngAgeTxt <- fireSenseUtils::youngAgeTxt
ranEffsLabel <- fireSenseUtils::yearTxt
cohDat <- "cohortData"
pixGM <- "pixelGroupMap"
saMap <- "standAgeMap"

#' Group fire years by the data year whose vegetation they use
#'
#' Each fire year belongs to the latest data year at or before it. Stops if a fire year precedes
#' the first data year, or if a data year gets no fire years.
#'
#' @param dataYears increasing integer years with vegetation data.
#' @param fireYears integer years of fire records.
#' @param minmaxOnly if `TRUE`, return only the first and last fire year of each group.
#' @return list with one integer vector per data year, named by it.
yearGroups <- function(dataYears, fireYears, minmaxOnly = TRUE) {
  if (any(fireYears < min(dataYears))) {
    stop("fireYears before the first dataYear (", min(dataYears), ") have no vegetation data: ",
         paste(fireYears[fireYears < min(dataYears)], collapse = ", "))
  }
  ageGroups <- split(fireYears, factor(dataYears[findInterval(fireYears, dataYears)], levels = dataYears))
  empty <- lengths(ageGroups) == 0
  if (any(empty)) {
    stop("dataYears with no fireYears before the next dataYear: ", paste(dataYears[empty], collapse = ", "))
  }
  if (isTRUE(minmaxOnly))
    ageGroups <- Map(ag = ageGroups, function(ag) c(min(ag), max(ag)))
  ageGroups
}

#' Stop when the climate rasters do not cover every requested fire year
#'
#' Deliberately not a warning plus truncation: that would fit a different window than the one
#' requested, with nothing in the results to say so.
#'
#' @param allYearsVect character, the requested fire years as `year<year>`.
#' @param whAvailable logical, same length: is there climate for that year?
#' @return `allYearsVect`, invisibly, if all are available; otherwise an error naming the missing years.
checkClimateYears <- function(allYearsVect, whAvailable) {
  if (all(whAvailable))
    return(invisible(allYearsVect))
  missingYears <- gsub(fireSenseUtils::yearTxt, "", allYearsVect[!whAvailable])
  stop("sim$historicalClimateRasters has no climate for ", sum(!whAvailable), " of the ",
       length(allYearsVect), " requested fireYears: ", paste(missingYears, collapse = ", "),
       "\nEither supply climate for those years or set P(sim)$fireYears to a window the ",
       "climate data covers. It is NOT truncated automatically: that would silently fit a ",
       "different window than the one requested.")
}

#' Join fire buffers to the vegetation of their data year
#'
#' @param fireBufferedListDT list of data.tables of `pixelID`, `ids` and `buffer`, named `year<year>`.
#' @param vegData data.table of fuel covariates by `pixelID` and `year` (the data year, as
#'   character). Modified by reference: `year` becomes integer.
#' @param allYears list of `year<year>` fire years, named by data year, as from `yearGroups()`.
#' @return data.table with one row per buffered pixel and fire year (`fireYear`), with the
#'   covariates of that fire year's data year only.
joinFireBuffersToVeg <- function(fireBufferedListDT, vegData, allYears) {
  vegData[, year := gsub(fireSenseUtils::yearTxt, "", year) |> as.integer()]
  indices <- Map(yrs = allYears, dataYear = as.integer(names(allYears)), function(yrs, dataYear) {
    dt <- rbindlist(fireBufferedListDT[yrs], idcol = "fireYear")
    dt2 <- dt[vegData[year == dataYear], on = c("pixelID")]
    dt2[!is.na(buffer)]
  })
  rbindlist(indices)
}

#' Coordinates of points, for `terra::cellFromXY()`
#'
#' @param points `SpatVector` or `sf` points.
#' @return two-column matrix of x and y, also for a single point (hence `drop = FALSE`:
#'   `cellFromXY()` rejects a plain vector).
pointCoords <- function(points) {
  if (is(points, "SpatVector")) {
    terra::geom(points)[, c("x", "y"), drop = FALSE]
  } else {
    sf::st_coordinates(points)
  }
}

#' Project points to the study area's CRS, if needed, and drop those outside it
#'
#' `prepare_IgnitionFit()` asserts that every ignition point is within the study area, so the
#' clip is unconditional and against the polygon, not the `rasterToMatch` rectangle.
#'
#' @param points `SpatVector` or `sf` points.
#' @param studyAreaUnion single study area polygon (`SpatVector` or `sf`).
#' @return the points inside `studyAreaUnion`, in its CRS.
clipPointsToStudyArea <- function(points, studyAreaUnion) {
  if (!terra::same.crs(points, studyAreaUnion)) {
    points <- projectTo(points, terra::crs(studyAreaUnion))
  }
  maskTo(points, studyAreaUnion)
}

#' Digest of the code of modules, for a Cache key that changes when that code does
#'
#' @param modulePath directory containing the modules.
#' @param modules character, module names.
#' @return digest of each module's `.R` file and the files in its `R/` folder. The paths must be
#'   full (`full.names = TRUE`): the digest of a bare file name does not read the file.
moduleCodeDigest <- function(modulePath, modules) {
  outerDirs <- file.path(modulePath, modules)
  files <- dir(c(outerDirs, file.path(outerDirs, "R")), pattern = "\\.R$", full.names = TRUE)
  .robustDigest(asPath(files))
}
