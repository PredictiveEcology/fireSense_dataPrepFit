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
  version = list(fireSense_dataPrepFit = "1.2.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = deparse(list("README.md", "fireSense_dataPrepFit.Rmd")),
  loadOrder = list(before = c("Biomass_speciesData", "Biomass_borealDataPrep", "Biomass_speciesParameters")),
  reqdPkgs = list("data.table", "fastDummies", "reproducible", "Require",
                  "PredictiveEcology/climateData@development (>= 2.2.3)",
                  "PredictiveEcology/fireSenseUtils@development (>= 0.2.3.9005)",
                  "FOR-CAST/fireregimetools@main (>= 0.1.0.9006)",
                  "ggplot2", "parallel", "purrr", "raster", "sf", "sp",
                  "PredictiveEcology/LandR@development (>= 1.2.0.9012)",
                  "PredictiveEcology/SpaDES.core@development (>= 2.0.2.9006)",
                  "PredictiveEcology/SpaDES.project@development",
                  "PredictiveEcology/SpaDES.tools@development (>= 2.1.1.9000)",
                  "snow", "terra"),
  parameters = bindrows(
    defineParameter("areaMultiplier", c("numeric", "name"), quote(fireSenseUtils::multiplier), NA, NA,
                    paste("Either a scalar that will buffer `areaMultiplier * fireSize` or ",
                          "a quoted function of `fireSize`. See `?fireSenseUtils::bufferToArea`.")),
    defineParameter("bufferForFireRaster", "numeric", 1000, 0, NA,
                    paste("The distance that determine whether separate patches of burned pixels originated",
                          "from the same fire. Only relevant when `useRasterizedFireForSpread = TRUE`.",
                          "This param is separate from `minBufferSize`, which is used to determine the",
                          "minimum sample of burned and unburned pixels to include in each fire.")),
    defineParameter("cutoffForYoungAge", "numeric", 15, NA, NA,
                    "Age at and below which pixels are considered 'young' (`young <- age <= cutoffForYoungAge`)"),
    defineParameter("dataYears", "integer", c(2000L, 2010L, 2020L), NA_integer_, NA_integer_,
                    paste("A numeric vector of length 2 or more (only tested with 2 and 3) indicating",
                          "which years should be used for standAgeMaps, rstLCCs etc.")),
    defineParameter("estimateFuelClasses", "logical", TRUE, NA, NA,
                    paste("estimate fuel classes from combination of data and P(sim)$fuelClassCol?")),
    defineParameter("fireYears", "integer", 2002:2025, NA, NA,
                    paste("A numeric vector indicating which years should be extracted",
                          "from the fire databases to use for fitting.",
                          "Should *not* include years prior to 2002, to ensure correct intialization from data.")),
    defineParameter("flammabilityThreshold", "numeric", 0.1, 0, 1,
                    paste("Minimum proportion of flammable old pixel needed to define a new pixel
                          as flammable when upscaling the default flammable maps`.")),
    defineParameter("forestedLCC", "numeric", c(81, 210, 220, 230, 240), NA, NA,
                    paste("Forested land cover classes - these differ from non-forest because the biomass",
                          "and composition of fuels are taken into account by fireSense, while non-forest",
                          "classes are treated categorically")),
    defineParameter("igAggFactor", "numeric", 4, 1, NA, # was 4 before xgboost, Jun 19, 2025
                    "aggregation factor for rasters during ignition prep."),
    defineParameter("igFocalFactor", "integer", 1L, NA, NA, # was prior to xgboost Jun 19, 2025 as.integer(3), as.integer(3), NA,
                    paste("Use focal statistics at the base resolution as an alternative to aggregating ignition covariates",
                          "This will occur if `P(sim)$igAggFactor is <= 1, and igFocalFactor is > 1.",
                          "The parameter is the `w` in the `terra::focal` function, i.e. the number of cells. It must be odd,",
                          "thus 3 is the minimum")),
    defineParameter("fuelClassCol", "character", "FuelClass", NA, NA,
                    "the column in `sppEquiv` that defines unique fuel classes. A column ",
                    "named `FuelClass` exists in the `LandR::sppEquivalencies_CA` and will be used ",
                    "by default. To change the `FuelClass` classifications, add a column to that table, ",
                    "or to `sim$sppEquiv` and then modify this `fuelClassCol` parameter"),
    defineParameter("modelAlgorithm", "character", "xgboost", NA, NA,
                    "Can be `xgboost`, `glmmtmb`, `glm.nb`, `glmmadaptive`, `glm`; only `xgboost` is supported currently",
                    "Should be the same as modelAlgorithm used in fireSense_IgnitionFit"),
    defineParameter("minBufferSize", "numeric", 5000, NA, NA,
                    paste("Minimum number of cells in buffer and nonbuffer. This is imposed after the",
                          "multiplier on the `bufferToArea` fn")),
    defineParameter("nonflammableLCC", "numeric", c(0, 20, 31, 32, 33), NA, NA,
                    "non-flammable LCC in rstLCC layers - defaulting to water, snow/ice, rock, barren land."),
    defineParameter("nonForestCanBeYoungAge", "logical", TRUE, NA, NA,
                    paste("if TRUE, burned non-forest will be treated as `youngAge`. Recommended to be TRUE",
                          "as burned forest is often classified as non-forest")),
    defineParameter("sppEquivCol", "character", "LandR", NA, NA,
                    "column name in `sppEquiv` object that defines unique species in `cohortData`"),
    defineParameter("spreadFitGoogleDriveFolder", "character",
                    "https://drive.google.com/drive/folders/1X9-mRjyLMNpgkP_cfqhbr_AQEPOsVCHf",
                    # KNN pre Oct 2025: "https://drive.google.com/drive/u/0/folders/1spxq7CnL4kNcJoUQlRek2CmBJ1InAmbP",
                    NA, NA, "A Googledrive folder url where a file with fireSense studyArea exists as an 'sf' class object"),
    defineParameter("spreadFitFilename", "character", "fireSenseParams.rds",
                    NA, NA, "A Googledrive folder url where a file with fireSense studyArea exists as an 'sf' class object"),
    defineParameter("targetFuelClasses", "numeric", 5, 1, 7,
                    "the target number of unique fuel classes when using semi-automated approach"),
    defineParameter("useCentroids", "logical", TRUE, NA, NA,
                    paste("Should fire ignitions start at the `sim$firePolygons` centroids",
                          "or at the ignition points in `sim$firePoints`?")),
    defineParameter("useRasterizedFireForSpread", "logical", FALSE, NA, NA,
                    paste("Should rasterized fire be used in place of a vectorized fire dataset?",
                          "This method attributes burned pixels to specific fires,",
                          "only examines the latest fire in a pixel, and may be subject to temporal error.",
                          "is therefore more appropriate in areas with low rates of fire,",
                          "or where the NFDB dataset may be incomplete (e.g., northern Ontario).")),
    defineParameter("whichModulesToPrepare", "character",
                    c("fireSense_IgnitionFit", "fireSense_SpreadFit", "fireSense_EscapeFit"),
                    NA, NA, "Which fireSense fit modules to prep? defaults to all 3"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "Describes the simulation time interval between plot events."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first save event should occur."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events."),
    defineParameter(".studyAreaName", "character", NULL, NA, NA,
                    "`studyArea` name that will be appended to file-backed rasters"),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    paste("Should this entire module be run with caching activated? This is intended",
                          "for data-type modules, where stochasticity and time are not relevant"))
  ),
  inputObjects = bindrows(
    expectsInput("climateVariablesForFire", "list", sourceURL = NA,
                 paste("A list detailing which climate variables in `sim$historicalClimateRasters`",
                       "to use for which fire processes (ignition and spread). If the list is length one,",
                       "both processes will use the same variables. The default is to use 'MDC'.")),
    expectsInput("cohortDatas", "list", sourceURL = NA,
                 paste0("List of (2) data.tables with cohortData that defines the ",
                        "cohorts by pixelGroup in the years represented by the names of the list")),
    expectsInput("spreadFirePoints", "list", sourceURL = NA,
                 paste("named list of spatial points for each fire year",
                       "with each point denoting an ignition location.")),
    expectsInput("spreadFitAdditionalColNames", "character",
                 desc = paste0("The column names used to attach the spreadFit object and several ancilliary objects")),
    expectsInput("firePolys", "list", sourceURL = NA,
                 paste0("List of sf polygon objects representing annual fire polygons.",
                        "List must be named with followign convention: `year<numeric year>`")),
    expectsInput("firePolysForAge", "list", sourceURL = NA,
                 "list of fire polygons used to classify `timeSinceDisturbance` in nonforest LCC"),
    expectsInput("historicalFireRaster", "SpatRaster",
                 sourceURL = "https://opendata.nfis.org/downloads/forest_change/CA_Forest_Fire_1985-2020.zip",
                 "a raster with values representing fire year 1985-2020"),
    expectsInput("historicalClimateRasters", "list", sourceURL = NA,
                 paste("length-one list of containing a raster stack of historical climate",
                       "list named after the variable and raster layers named as `year<numeric year>`")),
    expectsInput("ignitionFirePoints", "SpatVector", sourceURL = NA,
                 paste("SpatVector of points representing annual ignition locations.",
                       "This includes all fires regardless of size. It should have the same CRS",
                       "as `sim$rasterToMatch`")),
    expectsInput("missingLCCgroup", "character", NA,
                 paste("if a pixel is forested but is absent from `cohortData`, it will be grouped in this class.",
                       "It can be estimated if `P(sim)$estimateFuelClasses` is TRUE.",
                       "If supplied, it must be one of the names in `sim$nonForestedLCCGroups`")),
    expectsInput("nonForestedLCCGroups", "list",
                 paste("a named list of non-forested landcover groups, e.g. list('wetland' = c(19, 23, 32))",
                       "These will become fuel covariates, and the groups will be estimated if",
                       "`P(sim)$estimateFuelClasses` is TRUE")),
    expectsInput("pixelGroupMaps", "list", sourceURL = NA,
                 paste0("List of (2) SpatRaster that goes with cohortDatas and that defines the ",
                        "cohorts by pixelGroup in the years represented by the names of the list")),
    expectsInput("propFlammables", "list", sourceURL = NA,
                 "Lists of (2) SpatRasters with proportion of flammable landcover in a pixel - for post-hoc analysis"),
    expectsInput("rasterToMatch", "SpatRaster", sourceURL = NA,
                 "template raster for study area. Assumes some buffering of core area to limit edge effect of fire."),
    expectsInput("rasterToMatch_biomassParam", "SpatRaster", sourceURL = NA,
                 "template raster for studyArea_biomassParam. Passed to Biomass_borealDataPrep. This is expected to be ",
                 "*large* as in previously it was rasterToMatchLarge. Several downstream steps assume that it is ",
                 "larger than studyArea or rasterToMatch"),
    expectsInput("rstLCCs", "list", sourceURL = NA,
                 paste0("List of (2) SpatRasters of land cover - updated so that pixels above `P(sim)$flammabilityThreshold",
                        "have an assigned flammable landcover")),
    expectsInput("sppEquiv", "data.table", sourceURL = NA,
                 "table of LandR species equivalencies"),
    expectsInput("standAgeMaps", "list", sourceURL = NA,
                 "list of length 2 of maps of stand age in dataYear[[1]] and dataYear[[2]]",
                 " used to create `cohortDatas`"),
    expectsInput("spreadFirePolys", "list", sourceURL = NA,
                  "list of sf polygon objects representing annual fires; this is an 'input' because ",
                  " the object is modified in a subsequent event; this is not required at 'init'"),
    expectsInput("studyArea", "SpatVector", sourceURL = NA,
                 "study area that determines spatial boundaries of all data. Should be buffered to accomodate edge effects"),
    expectsInput("studyArea_biomassParam", "SpatVector", sourceURL = NA,
                 "study area passed to Biomass_borealDataPrep for vegetation calibration"),
    expectsInput("studyAreaReporting", "sf", sourceURL = NA,
                 desc = paste("(optional) study area used for reporting purposes, specifically whether fires inside",
                              "the studyAreaReporting polygon are being removed for also falling partially outside",
                              "the studyArea polygon, indicating the buffered studyArea shoudl be expanded."))
  ),
  outputObjects = bindrows(
    createsOutput("climateVariablesForFire", "list", # sourceURL = NA,
                  paste("A list detailing which climate variables in `sim$historicalClimateRasters`",
                        "to use for which fire processes (ignition and spread). If the list is length one,",
                        "both processes will use the same variables. This will be an output ",
                        "if there is an existing studyAreaWithParams")),
    createsOutput("spreadFitPreRun", "data.frame",
                  desc = paste("This is a data.frame that has a geometry list column, so it can be ",
                               "converted to a sf or SpatVector (e.g., `terra::vect(sf::st_as_sf(sim$spreadFitPreRun))` ",
                               " , plus other mostly list columns:",
                               "numIterations, objFunVal (not list), params, sppEquiv, ", 
                               "nonForestedLCCGroups, missingLCCgroup, and polygonID. These are from ",
                               "previously fitted SpreadFit. If no pre-existing object exists from ",
                               "CacheGeo, this will be NULL")),
    createsOutput("studyAreaWithSpreadParams", "sf",
                  desc = paste("This is the studyArea, but with parameters from a previously ",
                               "fitted SpreadFit. If no pre-existing object exists from ",
                               "CacheGeo, this will be NULL")),
    createsOutput("sppColorVect", "character",
                  desc = "named character vector of hex colour codes corresponding to each species. ",
                  "This may be updated if CacheGeo has different values than the inputted version of this object"),
    createsOutput("sppNameVector", "character",
                 desc = paste("a vector of species names pulled from `sppEquiv`.",
                              "Is sets it to sim$sppEquiv[[Par$sppEquivCol]]")),
    createsOutput("climateVariables", "list",
                  paste("a list, named by climate variable using 'projected_' or 'historical_'",
                        "prefixes, with each list element containing a list of three arguments:",
                        "vars - the raw variables used to derive the target variable,",
                        "fun - the quoted function used to derive the target variable, where",
                        "'quote(calcAsIs)' denotes target variables that ARE the raw variable,",
                        "and dots - additional arguments passed to 'fun'. See the .inputObjects",
                        "for examples of how to build this object and ?climateData::prepClimateLayers",
                        "for how it is used . The GCM, SSP, and selected years, whether projected or",
                        "historical, are set by module parameters and must be identical for all variables")),
    createsOutput("fireBufferedListDT", "list",
                  "list of data.tables with fire id, `pixelID`, and buffer status"),
    createsOutput("fuelClassTable", "data.table",
                  "table with assigned fuel class of each tree species, after running assessFuelClasses"),
    createsOutput("spreadFirePolys", "list",
                  "list of sf polygon objects representing annual fires"),
    createsOutput("fireSense_annualSpreadFitCovariates", "list",
                  "list of tables with climate covariates, `youngAge`, burn status, `polyID`, and `pixelID`"),
    createsOutput("fireSense_escapeCovariates", "data.table",
                  "ignition covariates with added column of escapes"),
    createsOutput("fireSense_escapeFormula", "character",
                  "formula for escape, using fuel classes and landcover, as character"),
    createsOutput("fireSense_ignitionCovariates", "data.table",
                  "table of aggregated ignition covariates with annual ignitions"),
    createsOutput("fireSense_ignitionFormula", "character",
                  "formula for ignition, using climate and vegetation covariates, as character"),
    createsOutput("fireSense_nonAnnualSpreadFitCovariates", "list",
                  "list of two tables with vegetation covariates, burn status, polyID, and `pixelID`"),
    createsOutput("fireSense_spreadFormula", "character",
                  "formula for spread, using climate and vegetation covariates, as character"),
    createsOutput("ignitionFirePoints", "SpatVector",
                  paste("Same as object that is an input, but possibly changed CRS")),
    createsOutput("ignitionFitRTM", "SpatRaster",
                  paste("A (template) raster with information with regards to the spatial",
                        "resolution and geographical extent of `fireSense_ignitionCovariates`.",
                        "Used to pass this information onto `fireSense_ignitionFitted`",
                        "Needs to have number of non-NA cells as attribute (`attributes(ignitionFitRTM)$nonNAs`).")),
    createsOutput("landcoverDTs", "list",
                  "List of (2) `data.table`s with `pixelID` and relevant landcover classes for flammable pixels in each layer"),
    createsOutput("lightningMaps", "SpatRaster",
                  paste("A 4-layer SpatRaster of lightning: lightningDays, lightningDensity, positiveCG, positiveCGdensity")),
    createsOutput("missingLCCgroup", "character",
                  "if estimating fuel classes, the nonforest class to assign forested pixels absent from `sim$cohortData`"),
    createsOutput("nonForestedLCCGroups", "list",
                  paste("a named list of non-forested landcover groups forming distinct fuel classes",
                        "e.g. list('wetland' = c(19, 23, 32))")),
    createsOutput("nonForest_timeSinceDisturbances", "list",
                  "List of (2) SpatRaster with time since burn for non-forested pixels in each layer"),
    createsOutput("rstLCCs", "list", # sourceURL = NA,
                  paste0("At rasterToMatch geoms. List of (3) SpatRasters of land cover - updated so that pixels above `P(sim)$flammabilityThreshold",
                         "have an assigned flammable landcover")),
    createsOutput("flammableRTMs", "list", "List of (3 binary SpatRaster of flammable landcover for years given by the list names"),
    createsOutput("sppEquiv", "data.table", "sppEquiv table potentially modified with new or overwritten fuel class"),
    createsOutput("spreadFirePoints", "list",
                  paste("Named list of `sf` polygon objects representing annual fire centroids.",
                        "This only includes fires that escaped (e.g. `size > res(flammableRTM)`.")),
    
    # For fireSense_**Predict modules
    createsOutput("propFlammable", "SpatRaster", # sourceURL = NA,
                 "SpatRaster with proportion of flammable landcover in a pixel - for post-hoc analysis"),
    createsOutput("standAgeMap", "SpatRaster", "Single layer, which will be taken from the last of standAgeMaps"),
    createsOutput("rstLCC_RTM", "SpatRaster", # sourceURL = NA,
                  paste0("Same as rstLCC, but at rasterToMatch geoms")),
    createsOutput("rstLCC", "SpatRaster", # sourceURL = NA,
                  paste0("At rasterToMatch_biomassParams geoms. Final SpatRaster of land cover, i.e, conditions at start(sim). Taken from last layer of rstLCCs")),
    createsOutput("flammableRTM", "SpatRaster", "Flammable landcover i.e, conditions at start(sim). Taken from last layer of rstLCCs"),
    createsOutput("landcoverDT", "data.table",
                  "`pixelID` and relevant landcover classes for flammable pixels in each layer, ",
                  "i.e, conditions at start(sim). Taken from last layer of landcoverDTs"),
    createsOutput("nonForest_timeSinceDisturbance", "SpatRaster", 
                  paste("raster tracking the time since disturbance (e.g. burn) in each pixel, forested or not.", 
                        "Forested pixels are tracked based on the age in cohortData"))
    
  )
))

doEvent.fireSense_dataPrepFit = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      if (!all(P(sim)$whichModulesToPrepare %in%
               c("fireSense_SpreadFit", "fireSense_IgnitionFit", "fireSense_EscapeFit"))) {
        stop("unrecognized module to prepare - review parameter whichModulesToPrepare")
        ## NOTE: the camelcase is still different with FS from LandR Biomass
      }

      ## schedule future event(s)
      if ("fireSense_IgnitionFit" %in% P(sim)$whichModulesToPrepare)
        sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "prepIgnitionFitData", eventPriority = 1)
      if ("fireSense_EscapeFit" %in% P(sim)$whichModulesToPrepare)
        sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "prepEscapeFitData", eventPriority = 1)
      if ("fireSense_SpreadFit" %in% P(sim)$whichModulesToPrepare) {
        sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "prepSpreadFitData", eventPriority = 1)
      }

      ## do stuff for this event
      sim <- Init(sim)

      sim <- scheduleEvent(sim, end(sim), "fireSense_dataPrepFit", "plotAndMessage", eventPriority = 9)
      sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "cleanUp", eventPriority = 10)
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

Init <- function(sim) {

  sa <- sim$studyArea
  if (inherits(sa, "SpatVector")) sa <- st_as_sf(sa)
  
  prepInputsFSURL <- SpaDES.core::paramCheckOtherMods(sim, "spreadFitGoogleDriveFolder")
  sim$spreadFitPreRun <- CacheGeo(cloudFolderID = Par$spreadFitGoogleDriveFolder,
                              targetFile = Par$spreadFitFilename, purge = 7,
                              domain = sa, action = "nothing", useCache = FALSE,
                              destinationPath = inputPath(sim), bufferOK = TRUE) # |> Cache()
  mod$haveSpreadFit <- is(sim$spreadFitPreRun, "sf") || is(sim$spreadFitPreRun, "data.frame")
  if (mod$haveSpreadFit) {
    sim$studyAreaWithSpreadParams <- sim$spreadFitPreRun
    
    
    
    # remove the column called "params" ... this just allows for partial matching, with or without "s"
    #   in case somebody uses `parameters`, `param`, or `params`
    colNames <- setdiff(sim$spreadFitAdditionalColNames,
                        grep(value = TRUE, "param", sim$spreadFitAdditionalColNames))
    df <- sim$spreadFitPreRun
    df <- df[colNames]
    
    #affected objects sppNameVector
    assignToSim <- outputObjects(sim)$objectName
    outs2 <- Map(nam = sim$spreadFitPreRun[[fireSenseUtils::polygonIDTxt]], ind = seq_len(NROW(df)), function(ind, nam) {
      dfList <- lapply(df, function(x) x[[ind]]) # Take only first entry, which if within 1 ELF, will be correct
      list2env(dfList, environment()) # nolint: vars numIterations objFunVal sppEquiv nonForestedLCCGroups missingLCCgroup geometry
      # list2env(dfList, envir(sim)) # nolint: vars numIterations objFunVal sppEquiv nonForestedLCCGroups missingLCCgroup geometry
      sim$sppNameVector <- unique(sppEquiv[[Par$sppEquivCol]])
      sppOuts <- sppHarmonize(sppEquiv, sim$sppNameVector, P(sim)$sppEquivCol, sppColorVect = NULL,
                              vegLeadingProportion = NULL, studyArea = sim$studyArea)
      list2env(sppOuts, envir = environment())  # nolint: vars sppEquiv sppNameVector sppEquivCol sppColorVect
      # list2env(sppOuts, envir = envir(sim))  # nolint: vars sppEquiv sppNameVector sppEquivCol sppColorVect
      
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
    
    climateVariablesList <- Map(x = outs2, function(x) 
      modifyList(sim$climateVariables,
                 climateLayers(.climVars = x[["theseClimVars"]])))
    theseClimVarsNoUnderscore <- Map(x = outs2, function(x) 
      x[["theseClimVarsNoUnderscore"]])
    theseClimVarsNoUnderscore <- unique(theseClimVarsNoUnderscore)[[1]]
    
    sim$climateVariables <- unique(climateVariablesList)
    if (length(sim$climateVariables) > 1)
      stop("Currently, can't have different climate variables by ELF; must all be the same")
    sim$climateVariables <- sim$climateVariables[[1]]

    # Take exactly the ones in the existing object
    sim$climateVariablesForFire[["spread"]] <- theseClimVarsNoUnderscore

    # Append the ones in the object as a decent guess. There can be more variables for ignitionfit
    sim$climateVariablesForFire[["ignition"]] <-
      sort(unique(c(sim$climateVariablesForFire[["ignition"]], theseClimVarsNoUnderscore)))
  
  }
  
  sim$sppEquiv <- copy(sim$sppEquiv) #debugging error where FuelClass disappears

  #because BBDP wants objects potentially larger than studyArea,
  #crop rstLCC and standAgeMap to create smaller objects before their derived objects
  # (landcoverDT/flammableMap and nonForest_timeSinceDisturbance, respectively).
  ## sanity checks
  # This is at rasterToMatch_biomassParam which is needed by Biomass_borealDataPrep
  sim$rstLCC <- tail(sim$rstLCCs, 1)[[1]]
  
  objs <- sim$rstLCCs
  if (!LandR::.compareRas(sim$rasterToMatch, sim$rasterToMatch_biomassParam, stopOnError = FALSE)) {
    objs2 <- lapply(sim$rstLCCs, function(x) {
      postProcess(x, to = sim$rasterToMatch, method = "near")})
    objs <- objs2
  }
  if (!isInt(objs$year2000)) {
    objs <- Map(obj = objs, function(obj) LandR::asInt(obj))
  }
  # This makes rstLCCs same as sim$rasterToMatch instead of rasterToMatch_biomassParam
  sim$rstLCCs <- objs
  
  ## TODO: make this table multidimensional or a list
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

  if (!terra::same.crs(sim$ignitionFirePoints, sim$rasterToMatch)) {
    ## project it first, faster than the postProcessTo sequence pre-crop, project, mask, crop
    sim$ignitionFirePoints <- projectTo(sim$ignitionFirePoints, sim$rasterToMatch) |>
      cropTo(sim$rasterToMatch) |>
      maskTo(sim$rasterToMatch)
  }

  ## possible, if user-supplied
  if (!terra::same.crs(sim$firePolys[[1]], sim$rasterToMatch)) {
    sim$spreadFirePolys <- Map(fp = sim$firePolys,
                               function(fp)
                                 projectTo(fp, st_crs(sim$rasterToMatch))) |>
      Cache(.functionName = "projectTo_for_firePolys")
  } else {
    sim$spreadFirePolys <- sim$firePolys
  }

  #glm of burned ~ biomass (of each species) link = logit
  #which species are similar enough to lump

  #make two firemaps 2002-2020, and 2012-2020 -
  # use fun = min to ignore repeated burns (since these would become youngAge anyway)


  # Apply the function for each time period
  fires <- Reduce(rbind, sim$spreadFirePolys)
  #this must ensure landcover overrides species - it does not currently

  fuelObjs <- c("nonForestedLCCGroups", "fuelClassTable")
  userSupplied <- fuelObjs %in% sim$.userSuppliedObjNames
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
            .omitArgs = c("pixelGroupMap", "rstLCC"), .cacheExtra = list(digRTMs, digRstLCC))

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
    # print(sim$nonForestedLCCGroups)

  }
  
  #make this object small (used by the landcoverDT and time-since-disturbance steps)
  standAgeMaps <- lapply(sim$standAgeMaps, reproducible::postProcess, to = sim$rasterToMatch)

  # One landcover table per year for THIS study area, keyed by pixelID only. A fit is
  # always single-ELF: it uses the fuel objects assessed or supplied above
  # (nonForestedLCCGroups, missingLCCgroup), never ledger geometry. The per-polygon
  # build introduced for multi-ELF prediction (2026-06) depended on
  # studyAreaWithSpreadParamsNoOverlap, which exists only when the ledger already
  # covers the study area, so every fresh fit failed here with
  # "object 'studyAreaWithSpreadParamsNoOverlap' not found".
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
  
  
  

  ## cannot merge because before subsetting due to column differences over time

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


##TODO NAME OF THIS FUNCTION MUST CHANGE
dataPrepInit <- function(sim) {

  #TODO: test that this is mistake-proof
  if (length(sim$climateVariablesForFire) == 1) {
    sim$climateVariablesForFire <- list(
      ignition = sim$climateVariablesForFire[[1]],
      spread = sim$climateVariablesForFire[[1]]
    )
  }
  
  if (!all(unlist(sim$climateVariablesForFire) %in% names(sim$historicalClimateRasters))) {
    stop("mismatch between sim$climateVariablesForFire and sim$historicalClimateRasters")
  }
  return(invisible(sim)) 
}

prepare_SpreadFit <- function(sim) {

  ## Put in format for DEOptim that distinguishes annual and nonannual covariates
  ## Prepare annual spread fit covariates
  ## prep veg data ---------------------------------------------------------------------------------
  doAssertion <- getOption("fireSenseUtils.assertions", TRUE)

  ## sanity check the inputs
  lapply(sim$historicalClimateRasters, compareGeom, x = sim$rasterToMatch)

  ## when landcoverDT is included, as is the case here, non-forest pixels in cohortData are masked out
  ## this is necessary when LandR and fireSense have differing concepts of non-forest

  dig1 <- .robustDigest(list(sim$landcoverDTs, sim$flammableRTMs))
  dig1a <- .robustDigest(list(sim$cohortDatas, sim$pixelGroupMaps, sim$nonForest_timeSinceDisturbances))
  dig2 <- append(dig1, dig1a)
  
  # This should be an exact replicate of ignitionCovariatesCreate
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
    #TODO: fix this approach to work with two flammableRTMs
    sim <- prepare_SpreadFitFire_Raster(sim)
  } else {
    sim <- prepare_SpreadFitFire_Vector(sim)
  }

  ## join fire and veg data
  fireSenseVegData <- joinFireBuffersToVeg(sim$fireBufferedListDT, vegData, mod$allYears)

  ## TODO: discuss if this is expected (as far as I can tell, it is)

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

  RHS <- paste(paste0(sim$climateVariablesForFire$spread, collapse = " + "), youngAgeTxt,
               paste0(vegCols, collapse = " + "), sep =  " + ")

  ## this is a funny way to get years but avoids years with 0 fires
  allYears <- unname(unlist(mod$allYears))
  yearsWithFire <- allYears[allYears %in% names(sim$fireBufferedListDT)]

  #### prep climate data ####

  ## index removed as argument - as flammable pixels change between 2010 and 2011 (mainly water)
  ## the wide layout of this object is incompatible (or we have NAs in rows)
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

  ## join the climate variables with the other annual covariate - youngAge (a.k.a. time since fire)
  ## pmap allows for internal debugging when there are large lists that are passed in; Map does not
  annualCovariates <- Map(yrsChar = mod$allYears, function(yrsChar) {
    fireSense_annualSpreadFitCovariates[yrsChar]
  })

  # This should put climate and youngAge in every list element
  # This adds youngAge to annualCovariates: this is wrong; it should not as youngAge is non-annual
  #  KEEP THIS OUT (Eliot Jan 2026)

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

prepare_SpreadFitFire_Raster <- function(sim) {
  stop("these methods need to be revised for the two flammable RTMs, two landcoverDTs")
  ## TODO: do this, obviously
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
  ## ultimately this function should combine the climate data to avoid needless iteration,
  ## and even this duplicated step should be a function of "fire period" for >2 periods
  ## however, the rasterized fire prep is significantly different, and needs review first
  allYearsVect <- unlist(mod$allYears)
  nCores <- ifelse(grepl("Windows", Sys.info()[["sysname"]]), 1L,
                   sum(names(sim$spreadFirePolys) %in% allYearsVect))
  if (FALSE) {
    # This chunk visualizes the largest fire in each year, along with the buffers
    sizes <- lapply(harmonized2010$firePolys, function(x) max(x$SIZE_HA))
    biggestFires <- Map(size = sizes, fires = harmonized2010$firePolys,
                        function(size, fires) fires[fires$SIZE_HA == size, c("FIRE_ID", "SIZE_HA")])
    par(mfrow = c(3,4))
    dd <- Map(bf = biggestFires, buff = harmonized2010$fireBufferedListDT,
              function(bf, buff) {
                a <- sf::st_buffer(bf, dist = 15000)
                b <- sim$flammableRTM2010
                b[] <- NA
                b[buff$pixelID] <- 1
                b <- terra::crop(b, a)
              })
    # terra::plot(b, add = TRUE)
    Map(p = dd, r = biggestFires, function(p, r) {
      terra::plot(p$flammable)
      terra::plot(r, add = TRUE)
    })
  }

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

prepare_IgnitionFit <- function(sim) {

  stopifnot(
    "all ignitionFirePoints are not within studyArea" = identical(
      nrow(st_as_sf(sim$ignitionFirePoints)),
      nrow(st_intersection(st_as_sf(sim$ignitionFirePoints), st_as_sf(mod$studyAreaUnion)))
    )
  )

  ## account for forested pixels that aren't in cohortData
  ## TODO: make this elegant
  ## first put landcover into raster stack
  ## non-flammable pixels require zero values for non-forest landcover, not NA
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

  #instead of aggregating, take the focal --> Eliot -- this doesn't make sense, since we need to aggregate to coarser res

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
  yearsInClimateRast <- allYearsVect[whAvailable]
  yearsNotAvailable <- allYearsVect[!whAvailable]

  if (length(yearsNotAvailable)) {
    warning("P(sim)$fireYears includes more years than are available in ",
            "sim$historicalClimateRasters; \nmissing: ", paste(yearsNotAvailable, collapse = ", "),
            "\ntruncating P(sim)$fireYears to: ",
            paste0(min(yearsInClimateRast), ":", max(yearsInClimateRast)))
    years <- Map(y = years, function(y) intersect(y, yearsInClimateRast))
    P(sim)$fireYears <- intersect(as.numeric(gsub(fireSenseUtils::yearTxt, "", yearsInClimateRast)),
                                  P(sim)$fireYears)
  }

  ## join fuel class, LCC, and climate, subsetting to flamIndex, calculating n of ignitions

  sim$fireSense_ignitionCovariates <-
    mergePreparedCovs(years, fuelCovsCoarse, sim$ignitionFirePoints, sim$nonForestedLCCGroups,
                      ignitionClimateCoarse, sim$lightningMaps["lightningDays"], 
                      digest = append(dig2, list(P(sim)$igAggFactor)))
  # END OF CREATING ignitionCovariates
  
  
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

  
  if (grepl("xgb", Par$modelAlgorithm) %in% FALSE) {
    ## build formula
    igCovariates <- names(sim$fireSense_ignitionCovariates)
    igCovariates <- igCovariates[!igCovariates %in%
                                   c(names(ignitionClimate),
                                     fireSenseUtils::yearTxt, "ignitions", "ignitionsNoGT1", "pixelID")]
    ## this is safer for multiple climate variables
    interactionsDF <- as.data.table(expand.grid(igCovariates, sim$climateVariablesForFire$ignition))
    interactionsDF[, interaction := do.call(paste, c(.SD, sep = ":")), .SDcols = names(interactionsDF)]
    interactions <- interactionsDF$interaction

    ## sanity check for base::abbreviate
    if (!length(unique(interactions)) == length(igCovariates) * length(sim$climateVariablesForFire$ignition)) {
      warning("automated ignition formula construction needs review")
    }
    if (is.null(sim$fireSense_ignitionFormula)) {
      sim$fireSense_ignitionFormula <- paste0("ignitions ~ ",
                                              paste0("(1|", ranEffsLabel, ")"), " + ",
                                              # this longer formula has had more unrealistic results 12/12/2024
                                              paste0(interactions, collapse = " + "))
    }
  } else {
    # Won't have sim$fireSense_ignitionFormula for xgboost
  }

  return(invisible(sim))
}

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

cleanUpMod <- function(sim) {
  mod$firePolysForAge <- NULL
  mod$fireSenseVegData <- NULL

  return(invisible(sim))
}

### template for save events
Save <- function(sim) {
  sim <- saveFiles(sim)
  return(invisible(sim))
}

### template for plot events
plotAndMessage <- function(sim) {
  ## TODO: this could plot the ignition/spread covariates
  return(invisible(sim))
}

rmMissingPixels <- function(fbldt, pixelIDsAllowed)  {
  fbldt <- rbindlist(fbldt, idcol = fireSenseUtils::yearTxt)
  fbldt <- fbldt[pixelID %in% unique(pixelIDsAllowed)]
  fireBufferedListDT <- split(fbldt, by = fireSenseUtils::yearTxt, keep.by = FALSE)
}

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
  # rstLCC * see below
  neededYears <- Par$dataYears # c(2010, 2020)
  
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
    outerDirs <- file.path(pathsLocal$modulePath, neededModule)
    innerRDirs <- file.path(outerDirs, "R")
    allModuleFiles <- dir(c(outerDirs, innerRDirs ), pattern = ".R$")
    sourceCodeDig <- .robustDigest(asPath(allModuleFiles))

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
  ## when that module has already run and left it NULL (fireSense_ELFs does, for an ELF without
  ## Engelmann spruce), so also build it when there is no table.
  if (!suppliedElsewhere("sppEquiv", sim, where = c("user", "initEvent")) || NROW(sim$sppEquiv) == 0) {
    sp <- LandR::speciesInStudyArea(studyArea = sim$studyArea, dPath = inputPath(sim))
    sp <- LandR::equivalentName(sp$speciesList, df = sppEquivalencies_CA, column = Par$sppEquivCol)
    sp <- sp[nzchar(sp)]
    sim$sppEquiv <- sppEquivalencies_CA[get(Par$sppEquivCol) %in% sp]
    sim$sppEquiv <- sim$sppEquiv[LANDIS_traits != "",] # ONLY USE THE SPECIES THAT HAVE TRAITS
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

  if (!suppliedElsewhere("climateVariablesForFire", sim)) {
    sim$climateVariablesForFire <- defaultClimateVariablesForFire
  }

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
              .functionName = paste0("makeFireSenseLCC", dy))
    }
    
    if (doStandAgeMaps) {
      standAgeMap <- prepInputsStandAgeMap(rasterToMatch = sim$rasterToMatch_biomassParam,
                                           # studyArea = sim$studyArea_biomassParam,
                                           destinationPath = dPath,
                                           dataYear = dy) |>
        Cache(.functionName = paste0("prepInputsStandAgeMap", dy),
              userTags = c(cacheTags, "prepInputsStandAgeMap"))
    } else {
      standAgeMap <- sim$standAgeMaps[[dyChar]]
    }
    
    return(list(standAgeMaps = standAgeMap, rstLCCs = LCC$lcc, propFlammables = LCC$flammableProp))
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
          #TODO: switch to the above when terra conversion is complete
          return(cent)
        }
      }

      sim$spreadFirePoints <- suppressWarnings(lapply(sim$firePolys, centerFun))
      # st_centroid assumes attributes are constant over geometries

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
      #this was happening with data update - the module will still run with no fire
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
    # SCANFI levels) returns the LABELS, so nothing matched and every class present
    # -- coniferous, broadleaf, mixedwood included -- became "non-forest".
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

defaultClimateVariablesForFire <- list("spread" = "MDC",
                                       "ignition" = "MDC")


yearGroups <- function(dataYears, fireYears, minmaxOnly = TRUE) {
  mm <- match(dataYears, fireYears)
  mm[is.na(mm)] <- 1
  mm1 <- rep(mm[1:2], diff(mm))
  mm2 <- c(mm1, rep(length(mm1) + 1, length(fireYears) - length(mm1) ))
  ageGroups <- split(fireYears, mm2)
  names(ageGroups) <- dataYears
  if (isTRUE(minmaxOnly))
    ageGroups <- Map(ag = ageGroups, function(ag) c(min(ag), max(ag)))
  ageGroups
}

## Join each group of fire years (`allYears`, named by data year, from `yearGroups()`) to the
## vegetation of that data year only. `vegData` is modified by reference (`year` becomes integer).
joinFireBuffersToVeg <- function(fireBufferedListDT, vegData, allYears) {
  vegData[, year := gsub(fireSenseUtils::yearTxt, "", year) |> as.integer()]
  indices <- Map(yrs = allYears, dataYear = as.integer(names(allYears)), function(yrs, dataYear) {
    dt <- rbindlist(fireBufferedListDT[yrs], idcol = "fireYear")
    dt2 <- dt[vegData[year == dataYear], on = c("pixelID")]
    dt2[!is.na(buffer)]
  })
  rbindlist(indices)
}

## x/y of each point as a two-column matrix, for terra::cellFromXY(). `terra::geom(points)[, c("x", "y")]`
## without drop = FALSE turned a single point into a plain vector, which cellFromXY() rejects
## ("unable to find an inherited method ... xy = numeric"; ELF 12.1, one escape).
pointCoords <- function(points) {
  if (is(points, "SpatVector")) {
    terra::geom(points)[, c("x", "y"), drop = FALSE]
  } else {
    sf::st_coordinates(points)
  }
}

# polygonIDTxt <- "polygonID"
