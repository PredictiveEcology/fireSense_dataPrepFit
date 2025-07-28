defineModule(sim, list(
  name = "fireSense_dataPrepFit",
  description = "Prepare data required by `fireSense_IginitionFit`, `fireSense_EscapeFit`, and `fireSense_SpreadFit`.",
  keywords = "fireSense",
  authors = c(
    person("Ian", "Eddy", role = c("aut", "cre"), email = "ian.eddy@nrcan-rncan.gc.ca"),
    person(c("Alex", "M"), "Chubaty", role = c("ctb"), email = "achubaty@for-cast.ca")
  ),
  childModules = character(0),
  version = list(fireSense_dataPrepFit = "1.1.1"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = deparse(list("README.md", "fireSense_dataPrepFit.Rmd")),
  loadOrder = list(before = c("Biomass_speciesData", "Biomass_borealDataPrep", "Biomass_speciesParameters"),
                   after = c("canClimateData")),
  reqdPkgs = list("data.table", "fastDummies", "reproducible",
                  "PredictiveEcology/fireSenseUtils@development (>= 0.0.5.9091)",
                  "ggplot2", "parallel", "purrr", "raster", "sf", "sp",
                  "PredictiveEcology/LandR@development (>= 1.1.5.9029)",
                  "PredictiveEcology/SpaDES.core@development (>= 2.0.2.9006)",
                  "PredictiveEcology/SpaDES.project@development",
                  "PredictiveEcology/SpaDES.tools (>= 2.0.4.9002)",
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
    defineParameter("estimateFuelClasses", "logical", TRUE, NA, NA,
                    paste("estimate fuel classes from combination of data and P(sim)$fuelClassCol?")),
    defineParameter("fireYears", "integer", 2002:2021, NA, NA,
                    paste("A numeric vector indicating which years should be extracted",
                          "from the fire databases to use for fitting.",
                          "Should *not* include years prior to 2002, to ensure correct intialization from data.")),
    defineParameter("flammabilityThreshold", "numeric", 0.1, 0, 1,
                    paste("Minimum proportion of flammable old pixel needed to define a new pixel
                          as flammable when upscaling the default flammable maps`.")),
    defineParameter("forestedLCC", "numeric", c(0, 81, 210, 220, 230, 240), NA, NA,
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
    defineParameter("spreadFitGoogleDriveFolder", "character", "https://drive.google.com/drive/u/0/folders/1spxq7CnL4kNcJoUQlRek2CmBJ1InAmbP",
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
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first plot event should occur."),
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
    expectsInput("cohortData2001", "data.table", sourceURL = NA,
                 paste0("Table that defines the cohorts by pixelGroup in 2001")),
    expectsInput("cohortData2011", "data.table", sourceURL = NA,
                 paste0("Table that defines the cohorts by pixelGroup in 2011")),
    expectsInput("spreadFirePoints", "list", sourceURL = NA,
                 paste("named list of spatial points for each fire year",
                       "with each point denoting an ignition location.")),
    # expectsInput("fireSense_spreadLogisticTermNames", "character",
    #              desc = paste0("The term names for the logistic terms in the spread fit")),
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
    expectsInput("ignitionFirePoints", "sf", sourceURL = NA,
                 paste("list of sf polygon objects representing annual ignition locations.",
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
    expectsInput("pixelGroupMap2001", "SpatRaster", sourceURL = NA,
                 "defines the `pixelGroups` for cohortData table in 2001"),
    expectsInput("pixelGroupMap2011", "SpatRaster",
                 "defines the `pixelGroups` for cohortData table in 2011"),
    expectsInput("propFlammable2001", "SpatRaster", sourceURL = NA,
                 "proportion of flammable landcover in a pixel - for post-hoc analysis"),
    expectsInput("propFlammable2011", "SpatRaster", sourceURL = NA,
                 "proportion of flammable landcover in a pixel - for post-hoc analysis"),
    expectsInput("rasterToMatch", "SpatRaster", sourceURL = NA,
                 "template raster for study area. Assumes some buffering of core area to limit edge effect of fire."),
    expectsInput("rasterToMatch_biomassParam", "SpatRaster", sourceURL = NA,
                 "template raster for studyArea_biomassParam. Passed to Biomass_borealDataPrep."),
    expectsInput("rstLCC2001", "SpatRaster", sourceURL = NA,
                 paste0("Raster of 2001 land cover - updated so that pixels above `P(sim)$flammabilityThreshold",
                        "have an assigned flammable landcover")),
    expectsInput("rstLCC2011", "SpatRaster", sourceURL = NA,
                 paste0("Raster of 2011 land cover - updated so that pixels above `P(sim)$flammabilityThreshold",
                        "have an assigned flammable landcover")),
    expectsInput("sppEquiv", "data.table", sourceURL = NA,
                 "table of LandR species equivalencies"),
    expectsInput("standAgeMap2001", "SpatRaster", sourceURL = NA,
                 "map of stand age in 2001 used to create `cohortData2001`"),
    expectsInput("standAgeMap2011", "SpatRaster", sourceURL = NA,
                 "map of stand age in 2011 used to create `cohortData2011`"),
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
    createsOutput("ignitionFirePoints", "sf",
                 paste("Same as object that is an input, but possibly changed CRS")),
    createsOutput("ignitionFitRTM", "SpatRaster",
                  paste("A (template) raster with information with regards to the spatial",
                        "resolution and geographical extent of `fireSense_ignitionCovariates`.",
                        "Used to pass this information onto `fireSense_ignitionFitted`",
                        "Needs to have number of non-NA cells as attribute (`attributes(ignitionFitRTM)$nonNAs`).")),
    createsOutput("landcoverDT2001", "data.table",
                  "data.table with `pixelID` and relevant landcover classes for flammable pixels in 2001 "),
    createsOutput("landcoverDT2011", "data.table",
                  "data.table with `pixelID` and relevant landcover classes for flammable pixels in 2011"),
    createsOutput("lightningMaps", "SpatRaster",
                  paste("A 4-layer SpatRaster of lightning: lightningDays, lightningDensity, positiveCG, positiveCGdensity")),
    createsOutput("missingLCCgroup", "character",
                  "if estimating fuel classes, the nonforest class to assign forested pixels absent from `sim$cohortData`"),
    createsOutput("nonForestedLCCGroups", "list",
                  paste("a named list of non-forested landcover groups forming distinct fuel classes",
                        "e.g. list('wetland' = c(19, 23, 32))")),
    createsOutput("nonForest_timeSinceDisturbance2001", "SpatRaster",
                  "time since burn for non-forested pixels in 2001"),
    createsOutput("nonForest_timeSinceDisturbance2011", "SpatRaster",
                  "time since burn for non-forested pixels in 2011"),
    createsOutput("rstLCC2001", "SpatRaster",
                 paste0("Raster of 2001 land cover - updated so that pixels above `P(sim)$flammabilityThreshold",
                        "have an assigned flammable landcover")),
    createsOutput("rstLCC2011", "SpatRaster",
                 paste0("Raster of 2011 land cover - updated so that pixels above `P(sim)$flammabilityThreshold",
                        "have an assigned flammable landcover")),
    createsOutput("flammableRTM2001", "SpatRaster", "binary raster of flammable landcover for 2001"),
    createsOutput("flammableRTM2011", "SpatRaster", "binary raster of flammable landcover for 2011"),
    createsOutput("sppEquiv", "data.table", "sppEquiv table potentially modified with new or overwritten fuel class"),
    createsOutput("spreadFirePoints", "list",
                  paste("Named list of `sf` polygon objects representing annual fire centroids.",
                        "This only includes fires that escaped (e.g. `size > res(flammableRTM)`."))
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

  sim$sppEquiv <- copy(sim$sppEquiv) #debugging error where FuelClass disappears

  #because BBDP wants objects potentially larger than studyArea,
  #crop rstLCC and standAgeMap to create smaller objects before their derived objects
  # (landcoverDT/flammableMap and nonForest_timeSinceDisturbance, respectively).
  ## sanity checks
  objs <- c(sim$standAgeMap2001, sim$standAgeMap2011, sim$rstLCC2001, sim$rstLCC2011)
  if (!LandR::.compareRas(sim$rasterToMatch, sim$rasterToMatch_biomassParam, stopOnError = FALSE)) {
    objs <- lapply(objs, FUN = postProcess, to = sim$rasterToMatch)
  }
  if (!isInt(objs[[1]]) | !isInt(objs[[3]])) {
    objs <- lapply(objs, LandR::asInt)
  }
  standAgeMap2001 <- objs[[1]]
  standAgeMap2011 <- objs[[2]]
  rstLCC2001 <- objs[[3]]
  rstLCC2011 <- objs[[4]]

  sim$flammableRTM2001 <- defineFlammable(rstLCC2001,
                                          nonFlammClasses = P(sim)$nonflammableLCC,
                                          to = sim$rasterToMatch)

  sim$flammableRTM2011 <- defineFlammable(rstLCC2011,
                                          nonFlammClasses = P(sim)$nonflammableLCC,
                                          to = sim$rasterToMatch)
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
    # projectTo(fp, crs(sim$rasterToMatch))) |> #terrarize
    sim$spreadFirePolys <- Map(fp = sim$firePolys, function(fp)
      projectTo(fp, st_crs(sim$rasterToMatch))) |>
      Cache(.functionName = "projectTo_for_firePolys")
  } else {
    sim$spreadFirePolys <- sim$firePolys
  }

  #glm of burned ~ biomass (of each species) link = logit
  #which species are similar enough to lump

  #make two firemaps 2002-2011, and 2012-2020 -
  # use fun = min to ignore repeated burns (since these would become youngAge anyway)


  # Apply the function for each time period
  fires <- do.call(rbind, sim$spreadFirePolys)
  #this must ensure landcover overrides species - it does not currently

  # Determine whether user has supplied their own FuelClass col in sppEquiv; their own nonForestedLCCGroups,
  #  their own fuelClassTable. If so, then don't estimate them here.
  fuelObjs <- c("nonForestedLCCGroups", "fuelClassTable")
  userSupplied <- fuelObjs %in% sim$.userSuppliedObjNames
  sppFCSupplied <- LandR::sppEquivalencies_CA[sim$sppEquiv, on = "LandR"]
  userSuppliedFC <- sppFCSupplied[, FuelClass == i.FuelClass]

  needToEstimateFuelClasses <- P(sim)$estimateFuelClasses && all(userSupplied %in% FALSE) && all(userSuppliedFC %in% TRUE)

  if (needToEstimateFuelClasses) {

    # THIS COMMENTED SECTION REBUILT THE OBJECT; NEEDED BECAUSE OF A BUG THAT MAYBE IS NO LONGER
    # df <- prepInputs(url = "https://drive.google.com/file/d/1KSexouuUvxEFuCdu173i5R7xyRDrKDqI/view?usp=drive_link")
    #
    # df1 <- lapply(df, function(x) {
    #   y <- parse(text = x)
    #   z <- try(eval(y), silent = TRUE)
    #   if (is(z, "try-error"))
    #     x
    #   else
    #     z
    #   })
    #
    # # I(list(df1$focal_modal)),
    #
    # df <- data.frame(I(list(as.data.table(df1$params))),
    #                  I(list(as.data.table(df1$sppEquiv))),
    #                  I(list(df1$nonForestedLCCGroups)),
    #                  I(list(df1$missingLCCgroup))
    #                  ) |> setNames(sim$spreadFitAdditionalColNames)
    #
    # sim$studyAreaWithSpreadParams <- sf::st_as_sf(sim$studyArea) |>
    #   dplyr::mutate(df)
    #
    # le <- function(x) {x}
    # spreadFitPreRun <- CacheGeo(cloudFolderID = Par$spreadFitGoogleDriveFolder,
    #                             targetFile = Par$spreadFitFilename,
    #                             domain = sim$studyArea,
    #                             destinationPath = inputPath(sim),
    #                             FUN = le(studyAreaFireSense),
    #                             le = le,
    #                             studyAreaFireSense = sim$studyAreaWithSpreadParams,
    #                             action = "update")
    sa <- sim$studyArea
    if (inherits(sa, "SpatVector")) sa <- st_as_sf(sa)
    spreadFitPreRun <- CacheGeo(cloudFolderID = Par$spreadFitGoogleDriveFolder,
                                targetFile = Par$spreadFitFilename,
                                domain = sa, action = "nothing",
                                destinationPath = getPaths()$inputPath, bufferOK = TRUE) # |> Cache()
    haveSpreadFit <- is(spreadFitPreRun, "sf") || is(spreadFitPreRun, "data.frame")

    if (haveSpreadFit) {

      if (TRUE) {
        # remove the column called "params" ... this just allows for partial matching, with or without "s"
        #   in case somebody uses `parameters`, `param`, or `params`
        colNames <- setdiff(sim$spreadFitAdditionalColNames,
                            grep(value = TRUE, "param", sim$spreadFitAdditionalColNames))
        df <- as.data.frame(spreadFitPreRun)
        df <- df[colNames]
        dfList <- lapply(df, function(x) x[[1]])
        list2env(dfList, envir(sim)) # sim$nonForestedLCCGroups, sim$sppEquiv, sim$missingLCCgroup

      } else {
        # This is to rebuild the objects if they are not contained within the spreadFitPreRun; defunct

        # speciesInCohortData <- lapply(list(sim$cohortData2001$speciesCode,
        #                                    sim$cohortData2011$speciesCode), levels) |>
        #   unlist() |> unique()
        #
        #
        # cols <- colnames(spreadFitPreRun$params[[1]])
        # allTerms <- setdiff(cols, sim$fireSense_spreadLogisticTermNames)
        # climateTerms <- sim$climateVariablesForFire$spread
        # fuelTerms <- setdiff(allTerms, climateTerms)
        # matureTerms <- setdiff(fuelTerms, "youngAge")
        # forestFuelClasses <- grep("^nf", matureTerms, invert = TRUE, value = TRUE)
        # dt <- data.table(species = speciesInCohortData)
        # # dt <- data.table(species = temp[[Par$sppEquivCol]])
        # dt2 <- abbreviateSpNames(dt)
        # forestFuelClassesSplit <- strsplit(forestFuelClasses, split = "\\.")
        # forestFuelClasses2 <- lapply(forestFuelClassesSplit, function(fc) {
        #   b <- lapply(fc, function(fcInner) {
        #     colHere <- c("species", "assignedFuelClass")
        #     a <- lapply(colHere, function(ch) {
        #       if (!is.null(dt2[[ch]]))
        #         if (any(grepl(fcInner, dt2[[ch]]))) {
        #           a <- dt2[data.table(fcInner) |> setNames(ch), on = ch]
        #           # set(a, NULL, "fuelType", paste(a[[ch]], sep = "."))
        #         }
        #     }) |> rbindlist(fill = TRUE)
        #     return(a)
        #   })
        #   d <- rbindlist(b)
        #   if (NROW(d) > 1) {
        #     ft <- paste(d[[2]], collapse = ".")
        #   } else {
        #     ft <- d[[1]]
        #   }
        #   set(d, NULL, "FuelClass", ft)
        #   d
        # }) |> rbindlist() |> set(NULL, "assignedFuelClass", NULL)
        #
        # sim$sppEquiv[, P(sim)$fuelClassCol := NULL]
        # setnames(forestFuelClasses2, old = "species", new = Par$sppEquivCol)
        # sim$sppEquiv <- sim$sppEquiv[forestFuelClasses2, on = Par$sppEquivCol]
        # # temp[sim$sppEquiv, on = Par$sppEquivCol]
        #
        # # now non-forest
        # nonForestFuelClasses <- grep("^nf", matureTerms, value = TRUE)
        # nonForestFuelClasses2 <- strsplit(nonForestFuelClasses, split = "_")
        # nonForestFuelClasses3 <- lapply(nonForestFuelClasses2, function(x) as.integer(x[-1]))
        # names(nonForestFuelClasses3) <- nonForestFuelClasses
        # sim$nonForestedLCCGroups <- nonForestFuelClasses3
      }
    }
  }


  if (!haveSpreadFit) {
    landscape <- Cache(
      Map,
      f = fuelClassPrep,
      pixelGroupMap = list(sim$pixelGroupMap2001, sim$pixelGroupMap2011),
      cohortData = list(sim$cohortData2001, sim$cohortData2011),
      rstLCC = list(rstLCC2001, rstLCC2001),
      yearRange = list(c(2002, 2011), c(2012, 2020)),
      MoreArgs = list(nonflammableLCC = P(sim)$nonflammableLCC,
                      fires = fires,
                      nonforestLCC = sim$nonForestedLCCGroups),
      userTags = c("fireSenseDataPrepFit", "fuelClassPrep")
    )

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
    # lccGG <- ggplot(lccGG_DF, aes(x = newLCC, y = percentBurn)) +
    #   geom_bar(stat = "identity") +
    #   labs(y = "% in study area burned", x = "nonForest LCC")
    #
    # sppGG <- ggplot(data = speciesGG_DF, aes(x  = speciesCode, y = percentBurn)) +
    #   geom_bar(stat = "identity") +
    #   labs(x = "fuel covariate", y = "% burned")

    # # Determine whether user has supplied their own FuelClass col in sppEquiv; their own nonForestedLCCGroups,
    # #  their own fuelClassTable. If so, then don't estimate them here.
    # fuelObjs <- c("nonForestedLCCGroups", "fuelClassTable")
    # userSupplied <- fuelObjs %in% sim$.userSuppliedObjNames
    # sppFCSupplied <- LandR::sppEquivalencies_CA[sim$sppEquiv, on = "LandR"]
    # userSuppliedFC <- sppFCSupplied[, FuelClass == i.FuelClass]

    if (needToEstimateFuelClasses) {

      # spreadFitPreRun <- CacheGeo(cloudFolderID = Par$spreadFitGoogleDriveFolder,
      #                             targetFile = Par$spreadFitFilename,
      #                             domain = sim$studyArea, action = "nothing",
      #                             destinationPath = getPaths()$inputPath, bufferOK = TRUE)

      fuelClassObjects <- Cache(
        assessFuelClasses(
          landscape = landscape,
          fuelCol = P(sim)$fuelClassCol,
          sppEquiv = sim$sppEquiv,
          sppEquivCol = P(sim)$sppEquivCol,
          targetFuelClasses = P(sim)$targetFuelClasses,
          nonforestLCC = nonforestLCC),
        userTags = c("assessFuelClasses", P(sim)$fuelClassCol))

      # sppEquiv
      temp <- fuelClassObjects$modSppEquiv[, .(species, assignedFuelClass)]
      setnames(temp, c(P(sim)$sppEquivCol, P(sim)$fuelClassCol))
      sim$sppEquiv[, P(sim)$fuelClassCol := NULL]
      sim$missingLCCgroup <- fuelClassObjects$missingLCCgroup

      # if (!(is(spreadFitPreRun, "sf") || is(spreadFitPreRun, "data.frame"))) {
      # set(temp, NULL, "FuelClass", NULL)
      # temp <- forestFuelClasses2[temp, on = c("species" = Par$sppEquivCol)]

      sim$fuelClassTable <- fuelClassObjects$modSppEquiv
      sim$nonForestedLCCGroups <- fuelClassObjects$nonForestedLCCGroups
      message("Estimated fuel classes for this study area:")
      messageDF(sim$fuelClassTable)

      # } else {
      # cols <- colnames(spreadFitPreRun$params[[1]])
      # allTerms <- setdiff(cols, sim$fireSense_spreadLogisticTermNames)
      # climateTerms <- sim$climateVariablesForFire$spread
      # fuelTerms <- setdiff(allTerms, climateTerms)
      # matureTerms <- setdiff(fuelTerms, "youngAge")
      # forestFuelClasses <- grep("^nf", matureTerms, invert = TRUE, value = TRUE)
      # dt <- data.table(species = speciesInCohortData)
      # # dt <- data.table(species = temp[[Par$sppEquivCol]])
      # dt2 <- abbreviateSpNames(dt)
      # forestFuelClassesSplit <- strsplit(forestFuelClasses, split = "\\.")
      # forestFuelClasses2 <- lapply(forestFuelClassesSplit, function(fc) {
      #   b <- lapply(fc, function(fcInner) {
      #     colHere <- c("species", "assignedFuelClass")
      #     a <- lapply(colHere, function(ch) {
      #       if (!is.null(dt2[[ch]]))
      #         if (any(grepl(fcInner, dt2[[ch]]))) {
      #           a <- dt2[data.table(fcInner) |> setNames(ch), on = ch]
      #           # set(a, NULL, "fuelType", paste(a[[ch]], sep = "."))
      #         }
      #     }) |> rbindlist(fill = TRUE)
      #     return(a)
      #   })
      #   d <- rbindlist(b)
      #   if (NROW(d) > 1) {
      #     ft <- paste(d[[2]], collapse = ".")
      #   } else {
      #     ft <- d[[1]]
      #   }
      #   set(d, NULL, "FuelClass", ft)
      #   d
      # }) |> rbindlist() |> set(NULL, "assignedFuelClass", NULL)
      #
      # forestFuelClasses2[sim$sppEquiv, on = c("species" = Par$sppEquivCol)]
      # setnames(temp, old = "species", new = Par$sppEquivCol)
      # # temp[sim$sppEquiv, on = Par$sppEquivCol]
      #
      # # now non-forest
      # nonForestFuelClasses <- grep("^nf", matureTerms, value = TRUE)
      # nonForestFuelClasses2 <- strsplit(nonForestFuelClasses, split = "_")
      # nonForestFuelClasses3 <- lapply(nonForestFuelClasses2, function(x) as.integer(x[-1]))
      # names(nonForestFuelClasses3) <- nonForestFuelClasses
      # sim$nonForestedLCCGroups <- nonForestFuelClasses3
      # }
      #overwrite sim$sppEquiv
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
    print(sim$nonForestedLCCGroups)

  }


  ## TODO: make this table multidimensional or a list
  sim$landcoverDT2001 <- makeLandcoverDT(rstLCC = rstLCC2001, flammableRTM = sim$flammableRTM2001,
                                         forestedLCC = P(sim)$forestedLCC, sim$nonForestedLCCGroups)
  sim$landcoverDT2011 <- makeLandcoverDT(rstLCC = rstLCC2001, flammableRTM = sim$flammableRTM2011,
                                         forestedLCC = P(sim)$forestedLCC, sim$nonForestedLCCGroups)
  sim$landcoverDT2001 <- correctMissingLCC(sim$landcoverDT2001, sim$pixelGroupMap2001, sim$missingLCCgroup)
  sim$landcoverDT2011 <- correctMissingLCC(sim$landcoverDT2011, sim$pixelGroupMap2011, sim$missingLCCgroup)

  ## cannot merge because before subsetting due to column differences over time

  ## TODO: this object is used to track annual youngAge of all pixels, forested or not
  ## so "nonForest" is a poor choice of name. It should not have values for non-flammable pixels.
  sim$nonForest_timeSinceDisturbance2001 <- makeTSD(
    year = 2001,
    fireRaster = sim$historicalFireRaster, ## can be NULL
    firePolys = sim$firePolysForAge,
    standAgeMap = standAgeMap2001,
    lcc = sim$landcoverDT2001,
    cutoffForYoungAge = P(sim)$cutoffForYoungAge
  )
  sim$nonForest_timeSinceDisturbance2001[sim$flammableRTM2001[] == 0] <- NA

  sim$nonForest_timeSinceDisturbance2011 <- makeTSD(
    year = 2011,
    fireRaster = sim$historicalFireRaster, ## can be NULL
    firePolys = sim$firePolysForAge,
    standAgeMap = standAgeMap2011,
    lcc = sim$landcoverDT2011,
    cutoffForYoungAge = P(sim)$cutoffForYoungAge
  )
  sim$nonForest_timeSinceDisturbance2011[sim$flammableRTM2011[] == 0] <- NA

  ## Until youngAge treatment is identical between spread and ignition, no point in prepping veg here
  ## Currently youngAge is resolved annually in spread, but only once in ignition
  ## e.g. if a pixel ignited in 2008, its youngAge status in ignition is still determined by whether it was 15 in 2001,
  ## but its youngAge status for spread is deterimined by whether standAge < 15 in 2008

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
  vegData <- Map(
    f = cohortsToFuelClasses,
    cohortData = list(sim$cohortData2001, sim$cohortData2011),
    pixelGroupMap = list(sim$pixelGroupMap2001, sim$pixelGroupMap2011),
    landcoverDT = list(sim$landcoverDT2001, sim$landcoverDT2011),
    flammableRTM = list(sim$flammableRTM2001, sim$flammableRTM2011),
    MoreArgs = list(sppEquiv = sim$sppEquiv,
                    sppEquivCol = P(sim)$sppEquivCol,
                    fuelClassCol = P(sim)$fuelClassCol,
                    cutoffForYoungAge = -1)
  ) |>
    Cache(.functionName = "cohortsToFuelClasses") #youngAge will be resolved annually downstream

  ## log transform the biomass values, setting zeroes to zero
  ## ideally we move this into cohortsToFuelClasses...
  ## there is no young age here unlike ignition
  vegData <- lapply(vegData, FUN = function(x) {
    #to lessen the leverage of zeroes where there is no biomass
    #TODO: discuss implications for spread
    #change the zeroes to one log below the minimum in the data (in this case 100 g/m2)
    if (FALSE) {
      minimumB <- exp(log(100) - 1)
      x[x < minimumB] <- minimumB
      # x[x == 0] <- NA # for NAs
      x <- log(x)

    } else {
      x[] <- logMinB(x[]) # same as in dataPrepPredict
    }
    dt <- as.data.table(values(x))
    dt[, pixelID := 1:ncell(x)]
    return(dt)
  })
  gc()
  vegData[[1]][, year := 2002]
  vegData[[2]][, year := 2012]

  vegData[[1]] <- vegData[[1]][sim$landcoverDT2001, on = c("pixelID")]
  vegData[[2]] <- vegData[[2]][sim$landcoverDT2011, on = c("pixelID")]
  vegData <- rbindlist(vegData)

  #### prep the fire data ####
  if (P(sim)$useRasterizedFire) {
    #TODO: fix this approach to work with two flammableRTMs
    sim <- prepare_SpreadFitFire_Raster(sim)
  } else {
    sim <- prepare_SpreadFitFire_Vector(sim)
  }

  ## join fire and veg data
  pre2012 <- paste0("year", min(P(sim)$fireYears):2011)
  post2012 <- paste0("year", 2012:max(P(sim)$fireYears))

  pre2012Indices <- sim$fireBufferedListDT[pre2012] %>%
    rbindlist(.) %>%
    .[vegData[year < 2012], on = c("pixelID")] %>%
    .[!is.na(buffer)] ## some buffered pixels are non-flammable...
  ## TODO: discuss if this is expected (as far as I can tell, it is)

  post2012Indices <- sim$fireBufferedListDT[post2012] %>%
    rbindlist(.) %>%
    .[vegData[year >= 2012], on = c("pixelID")] %>%
    .[!is.na(buffer)]
  rm(vegData)
  gc()

  ## TODO: lines from creation of vegData onwards should be reviewed. Seems redundant.
  ## TODO: should column ids be in vegData? currently pixels appear > 1 time, within a fire period
  #because they can be burned or unburned, or in >1 fire years that are < 10 years apart
  fireSenseVegData <- rbind(pre2012Indices, post2012Indices)
  setnames(fireSenseVegData, "buffer", "burned")

  vegCols <- setdiff(names(fireSenseVegData), c("pixelID", "burned", "ids", "year"))
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

  RHS <- paste(paste0(sim$climateVariablesForFire$spread, collapse = " + "), "youngAge",
               paste0(vegCols, collapse = " + "), sep =  " + ")

  ## this is a funny way to get years but avoids years with 0 fires
  years <- paste0("year", P(sim)$fireYears)
  yearsWithFire <- years[years %in% names(sim$fireBufferedListDT)]
  pre2012int <- as.integer(min(P(sim)$fireYears):2011)
  post2012int <- as.integer(2012:max(P(sim)$fireYears))
  pre2012 <- yearsWithFire[yearsWithFire %in% paste0("year", pre2012int)]
  post2012 <- yearsWithFire[yearsWithFire %in% paste0("year", post2012int)]

  #### prep climate data ####

  ## index removed as argument - as flammable pixels change between 2001 and 2011 (mainly water)
  ## the wide layout of this object is incompatible (or we have NAs in rows)
  spreadClimate <- sim$historicalClimateRasters[sim$climateVariablesForFire$spread]
  #don't need climate data for years outside fire years
  spreadClimate <- lapply(spreadClimate, FUN = function(x){
    x <- terra::subset(x, subset = names(x) %in% paste0("year", P(sim)$fireYears))
  })

  climateDT <- Cache(climateRasterToDataTable,
                     historicalClimateRasters = spreadClimate,
                     userTags = c("climateRasterToDataTable", names(spreadClimate)))

  fbl <- rbindlist(sim$fireBufferedListDT, idcol = "year")
  rmCols <- setdiff(colnames(fbl), c("pixelID", "year"))
  set(fbl, NULL, rmCols, NULL)
  fbl <- climateDT[fbl, on = c("year", "pixelID"), nomatch = NULL]

  fireSense_annualSpreadFitCovariates <- split(fbl, by = "year", keep.by = FALSE)

  ## prepare non-annual spread fit covariates by getting the youngAge
  pre2012Indices <- sim$fireBufferedListDT[names(sim$fireBufferedListDT) %in% pre2012]
  post2012Indices <- sim$fireBufferedListDT[!names(sim$fireBufferedListDT) %in% pre2012]
  colsToExtract <- c("pixelID", vegCols)

  nonAnnualpre2012 <- fireSenseVegData[year < 2012, .SD, .SDcols = colsToExtract] %>%
    na.omit(.) %>%
    as.data.table(.) %>%
    .[!duplicated(pixelID), ]

  nonAnnualpost2012 <- fireSenseVegData[year >= 2012, .SD, .SDcols = colsToExtract] %>%
    na.omit(.) %>%
    as.data.table(.) %>%
    .[!duplicated(pixelID)] ## remove duplicates from same pixel diff year

  ## join the climate variables with the other annual covariate - youngAge (a.k.a. time since fire)
  ## pmap allows for internal debugging when there are large lists that are passed in; Map does not
  annualCovariates <- list(fireSense_annualSpreadFitCovariates[pre2012],
                           fireSense_annualSpreadFitCovariates[post2012])

  annualCovariates <- Cache(
    purrr::pmap,
    .l = list(
      #years = list(c(2001:2010), c(2011:max(P(sim)$fireYears))),
      years = list(pre2012int, post2012int),
      annualCovariates = annualCovariates,
      standAgeMap = list(sim$nonForest_timeSinceDisturbance2001,
                         sim$nonForest_timeSinceDisturbance2011)
    ),
    .f = calcYoungAge,
    fireBufferedListDT = sim$fireBufferedListDT,
    cutoffForYoungAge = P(sim)$cutoffForYoungAge
  )

  ## get rid of nonflammable pixels (here because the calcYoungAge function assigns ages to NA values,
  ## due to inconsistent treatment of non-forest age pixels  in kNN years and other products (0 vs NA)
  annualCovariates[[1]] <- lapply(annualCovariates[[1]], function(x) {
    x[pixelID %in% sim$landcoverDT2001$pixelID, ]
  })
  annualCovariates[[2]] <- lapply(annualCovariates[[2]], function(x) {
    x[pixelID %in% sim$landcoverDT2011$pixelID, ]
  })

  if (!P(sim)$nonForestCanBeYoungAge) {
    ## TODO: test this inversion of makeMutuallyExclusive's regular use
    args <- as.list(rep("youngAge", length = length(sim$nonForestedLCCGroups)))
    names(args) <- names(sim$nonForestedLCCGroups)
  } else {
    ## this is done later in spreadFit - but done here for accuracy of outputs
    args <- list("youngAge" = names(sim$nonForestedLCCGroups))
  }

  annualCovariates <- lapply(annualCovariates, makeMutuallyExclusive, mutuallyExclusiveCols = args)

  sim$fireSense_annualSpreadFitCovariates <- do.call(c, annualCovariates)

  sim$fireSense_nonAnnualSpreadFitCovariates <- list(nonAnnualpre2012, nonAnnualpost2012)
  names(sim$fireSense_nonAnnualSpreadFitCovariates) <- c(paste(names(pre2012Indices), collapse = "_"),
                                                         paste(names(post2012Indices), collapse = "_"))

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
  sim$fireBufferedListDT <- Cache(rasterFireBufferDT, years =  P(sim)$fireYears,
                                  fireRaster = historicalFireRaster,
                                  flammableRTM = sim$flammableRTM,
                                  bufferForFireRaster = P(sim)$bufferForFireRaster, verb = 1,
                                  areaMultiplier = eval(P(sim)$areaMultiplier),
                                  minSize = P(sim)$minBufferSize,
                                  cores = nCores, userTags = c(currentModule(sim), "rasterFireBufferDT"))
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
  ## the "year" prefix is added by fireSenseUtils::makeLociList - discuss what to do
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
  pre2012 <- paste0("year", min(P(sim)$fireYears):2011)
  post2012 <- paste0("year", 2012:max(P(sim)$fireYears))

  ## prep fire data -----------------------------------------------------------------
  if (is.null(sim$spreadFirePolys[[1]]$FIRE_ID)) {
    stop("firePolys needs a numeric FIRE_ID column")
  }

  if (!is.numeric(sim$spreadFirePolys[[1]]$FIRE_ID) | !is.numeric(sim$spreadFirePoints[[1]]$FIRE_ID)) {

    message("need numeric FIRE_ID column in fire polygons and points. Coercing to numeric...")
    #this is true of the current NFBB
    origNames <- names(sim$spreadFirePolys)
    PointsAndPolys <- lapply(names(sim$spreadFirePolys),
                             function(year, polys = sim$spreadFirePolys, points = sim$spreadFirePoints) {
                               polys <- polys[[year]]
                               points <- points[[year]]
                               ## ensure matching IDs
                               points <- points[points$FIRE_ID %in% polys$FIRE_ID,]
                               polys <- polys[polys$FIRE_ID %in% points$FIRE_ID,]
                               points$FIRE_ID <- as.numeric(as.factor(points$FIRE_ID))
                               polys$FIRE_ID <- as.numeric(as.factor(polys$FIRE_ID))
                               return(list(polys = polys, points = points))
                             })
    sim$spreadFirePoints <- lapply(PointsAndPolys, FUN = function(x) x[["points"]])
    sim$spreadFirePolys <- lapply(PointsAndPolys, FUN = function(x) x[["polys"]])
    rm(PointsAndPolys)
    names(sim$spreadFirePolys) <- origNames
    names(sim$spreadFirePoints) <- origNames
  }

  ## drop fires less than 1 px in size
  pixSizeHa <- prod(res(sim$flammableRTM2011)) / 1e4
  ## using x[x$SIZE_HA] will work with terra or sf, while subset will not (I believe...)
  sim$spreadFirePoints <- lapply(sim$spreadFirePoints, function(x, minSize = pixSizeHa) {
    x <- x[x$SIZE_HA > minSize,]
    if (nrow(x) > 0) x else NULL
  })

  sim$spreadFirePoints[sapply(sim$spreadFirePoints, is.null)] <- NULL ## silly R

  sim$spreadFirePolys <- lapply(sim$spreadFirePolys, function(x) {
    x <- x[x$SIZE_HA > pixSizeHa,]
    if (nrow(x) > 0) x else NULL
  })
  sim$spreadFirePolys[sapply(sim$spreadFirePolys, is.null)] <- NULL
  sim$spreadFirePoints <- sim$spreadFirePoints[names(sim$spreadFirePoints) %in% names(sim$spreadFirePolys)]
  ## this covers when years are NA, which are caused by fire years with no available data

  ## years run separately because flammableRTM is different
  ## ultimately this function should combine the climate data to avoid needless iteration,
  ## and even this duplicated step should be a function of "fire period" for >2 periods
  ## however, the rasterized fire prep is significantly different, and needs review first
  nCores <- ifelse(grepl("Windows", Sys.info()[["sysname"]]), 1L,
                   sum(names(sim$spreadFirePolys) %in% pre2012))
  if (FALSE) {
    # This chunk visualizes the largest fire in each year, along with the buffers
    sizes <- lapply(harmonized2001$firePolys, function(x) max(x$SIZE_HA))
    biggestFires <- Map(size = sizes, fires = harmonized2001$firePolys,
                        function(size, fires) fires[fires$SIZE_HA == size, c("FIRE_ID", "SIZE_HA")])
    par(mfrow = c(3,4))
    dd <- Map(bf = biggestFires, buff = harmonized2001$fireBufferedListDT,
              function(bf, buff) {
      a <- sf::st_buffer(bf, dist = 15000)
      b <- sim$flammableRTM2001
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

  harmonized2001 <- harmonizeFireData(
    firePolys = sim$spreadFirePolys[names(sim$spreadFirePolys) %in% pre2012], ## protects from missing years
    flammableRTM = sim$flammableRTM2001,
    spreadFirePoints = sim$spreadFirePoints[names(sim$spreadFirePoints) %in% pre2012], ## protects from missing years
    areaMultiplier = eval(P(sim)$areaMultiplier),
    minSize = P(sim)$minBufferSize,
    pointsIDcolumn = "FIRE_ID",
    cores = nCores
  ) |>
    Cache(userTags = c("harmonizeFireData", P(sim)$.studyAreaName, "2001"))
  harmonized2011 <- harmonizeFireData(
    firePolys = sim$spreadFirePolys[names(sim$spreadFirePolys) %in% post2012],
    flammableRTM = sim$flammableRTM2011,
    spreadFirePoints = sim$spreadFirePoints[names(sim$spreadFirePoints) %in% post2012],
    areaMultiplier = eval(P(sim)$areaMultiplier), minSize = P(sim)$minBufferSize,
    pointsIDcolumn = "FIRE_ID",
    cores = nCores
  ) |>
    Cache(userTags = c("harmonizeFireData", P(sim)$.studyAreaName, "2011"))

  sim$fireBufferedListDT <- append(harmonized2001$fireBufferedListDT,
                                   harmonized2011$fireBufferedListDT)
  sim$spreadFirePoints <- append(harmonized2001$spreadFirePoints, harmonized2011$spreadFirePoints)
  sim$spreadFirePolys <- append(harmonized2001$firePolys, harmonized2011$firePolys)

  omitYears <- sapply(sim$spreadFirePoints, is.null)
  sim$fireBufferedListDT[omitYears] <- NULL
  sim$spreadFirePolys[omitYears] <- NULL
  sim$spreadFirePoints[omitYears] <- NULL

  ## drop fire years from these lists that don't have any buffer points post-harmonization
  sim$spreadFirePolys <- lapply(names(sim$spreadFirePolys), FUN = function(year) {
    poly <- sim$spreadFirePolys[[year]]
    point <- sim$spreadFirePoints[[year]]
    poly <- poly[poly$FIRE_ID %in% point$FIRE_ID,]
    return(poly)
  })
  names(sim$spreadFirePolys) <- names(sim$spreadFirePoints)
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
  LCCras <- Map(
    f = putBackIntoRaster,
    landcoverDT = list(sim$landcoverDT2001, sim$landcoverDT2011),
    flammableMap = list(sim$flammableRTM2001, sim$flammableRTM2011),
    MoreArgs = list(lcc = names(sim$nonForestedLCCGroups))
  ) |>
    Cache(.functionName = "putBackIntoRaster",
          userTags = c("putBackIntoRaster", P(sim)$.studyAreaName))

  fuelClasses <- Map(
    f = cohortsToFuelClasses,
    cohortData = list(sim$cohortData2001, sim$cohortData2011),
    flammableRTM = list(sim$flammableRTM2001, sim$flammableRTM2011),
    landcoverDT = list(sim$landcoverDT2001, sim$landcoverDT2011),
    pixelGroupMap = list(sim$pixelGroupMap2001, sim$pixelGroupMap2011),
    MoreArgs = list(sppEquiv = sim$sppEquiv,
                    sppEquivCol = P(sim)$sppEquivCol,
                    fuelClassCol = P(sim)$fuelClassCol,
                    cutoffForYoungAge = P(sim)$cutoffForYoungAge)
  ) |>
    Cache(.functionName = "cohortsToFuelClasses")

  fuelClasses <- lapply(fuelClasses, FUN = function(x){
    bCols <- unique(sim$sppEquiv[[P(sim)$fuelClassCol]])
    xYA <- terra::subset(x, !names(x) %in% bCols)
    xBiomass <- terra::subset(x, names(x) %in% bCols)
    #to lessen the leverage of zeroes where there is no biomass
    #change the zeroes to one log below the minimum in the data (in this case 100 g/m2)
    # (this assumes every class has one pixel with minimum B - probably a safe assumption)
    minimumB <- exp(log(100) - 1)
    xBiomass[xBiomass < minimumB] <- minimumB
    xBiomass <- log(xBiomass)
    x <- c(xBiomass, xYA)
    return(x)
  })

  if (P(sim)$nonForestCanBeYoungAge) {
    ## this modifies the NF landcover by converting some NF to a new YA layer
    ## it must be done before aggregating

    LCCras <- Map(
      f = calcNonForestYoungAge,
      landcoverDT = list(sim$landcoverDT2001, sim$landcoverDT2011),
      NFTSD = list(sim$nonForest_timeSinceDisturbance2001,
                   sim$nonForest_timeSinceDisturbance2011),
      LCCras = list(LCCras[[1]], LCCras[[2]]),
      MoreArgs = list(cutoffForYoungAge = P(sim)$cutoffForYoungAge)
    ) |>
      Cache(.functionName = "calcNonForestYoungAge")

    for (i in c(1:2)) {
      if ("youngAge" %in% names(fuelClasses[[i]])) {

        YA1 <- fuelClasses[[i]]$youngAge
        YA2 <- LCCras[[i]]$youngAge
        bothYA <- YA1 + YA2
        fuelClasses[[i]]$youngAge <- bothYA
      }  else {
        fuelClasses[[i]]$youngAge <- LCCras[[i]]$youngAge
      }
      toKeep <- setdiff(names(LCCras[[i]]), "youngAge")
      LCCras[[i]] <- terra::subset(LCCras[[i]], toKeep) ## to avoid double-counting
    }
  }

  ignitionClimate <- sim$historicalClimateRasters[sim$climateVariablesForFire$ignition]
  #instead of aggregating, take the focal
  if (P(sim)$igAggFactor > 1) {

  LCCras <- lapply(LCCras, aggregate, fact = P(sim)$igAggFactor, fun = mean) |>
    Cache(.functionName = "aggregate_LCCras_to_coarse")

  ## must specify terra::aggregate to avoid conflict with stats::aggregate
  fuelClasses <- lapply(fuelClasses, FUN = terra::aggregate, fact = P(sim)$igAggFactor, fun = mean) |>
    Cache(.functionName = "aggregate_fuelClasses_to_coarse")

  ignitionClimate <- lapply(X = ignitionClimate, FUN = terra::aggregate,
                            fact = P(sim)$igAggFactor, fun = mean) |>
    Cache(.functionName = "aggregate_historicalClimateRasters_to_coarse")
  } else if (P(sim)$igFocalFactor > 2) {
    #two will trigger, 1 does nothing.
    igSpatial <- lapply(X = list(ignitionClimate, fuelClasses, LCCras), FUN = function(x, size = P(sim)$igFocalFactor) {
        #these are all lists due to time
        x <- lapply(x, FUN = terra::focal, w =  size, fun = mean, na.rm = TRUE)
      })
    ignitionClimate <- igSpatial[[1]]
    fuelClasses <- igSpatial[[2]]
    LCCras <- igSpatial[[3]]
  }
  names(LCCras) <- c("year2001", "year2011")
  names(fuelClasses) <- c("year2001", "year2011")

  lightningUrls <- list(lightningDays = "1jeKJquhVJsesoNk2EPP1QZkttX3Zwp5c",
                        lightningDensity = "12fnhfKtER-JXkl06M4_yZ3GvpZWtQlIr" ,
                        positiveCG = "1bn6cQ23tvPicFLHn1tz4Z3AqDzJI4r60",
                        positiveCGdensity = "1GNixhXj1Ex1jT0tWXfhmxef-dX3ze1a4")
  digRTM <- sim$rasterToMatch
  sim$lightningMaps <- Map(url = lightningUrls, nam = names(lightningUrls),
                       function(url, nam) {
                         {
                           prepInputs(url = url,
                                      fun = readLightningData(targetFile, to = sim$rasterToMatch),
                                      destinationPath = inputPath(sim))  |>
                             terra::aggregate(fact = P(sim)$igAggFactor)} |>
                           Cache(.functionName = paste0("prepInputs_lightning_", nam), omitArgs = "...",
                                 .cacheExtra = list(url = url, igAggFactor = P(sim)$igAggFactor,
                                                    rtm = digRTM))
                       })



  compareGeom(ignitionClimate[[1]], fuelClasses[[1]], fuelClasses[[2]]) ## safety check

  ## ignition won't have same years as spread so we do not use names of init objects
  ## The reason is some years may have ignitions but no fires, e.g. 2001 in RIA
  pre2012 <- paste0("year", min(P(sim)$fireYears):2011)
  post2012 <- paste0("year", 2012:max(P(sim)$fireYears))
  allYears <- c(pre2012, post2012)

  #assume that if multiple climate variables are present, they are of equal length
  #else bigger problems exist
  whAvailable <- allYears %in% names(ignitionClimate[[1]])
  yearsInClimateRast <- allYears[whAvailable]
  yearsNotAvailable <- allYears[!whAvailable]
  if (length(yearsNotAvailable)) {
    warning("P(sim)$fireYears includes more years than are available in ",
            "sim$historicalClimateRasters; \nmissing: ", paste(yearsNotAvailable, collapse = ", "),
            "\ntruncating P(sim)$fireYears to: ",
            paste0(min(yearsInClimateRast), ":", max(yearsInClimateRast)))
    pre2012 <- intersect(pre2012, yearsInClimateRast)
    post2012 <- intersect(post2012, yearsInClimateRast)
    P(sim)$fireYears <- intersect(as.numeric(gsub("year", "", yearsInClimateRast)),
                                  P(sim)$fireYears)
  }

  ## join fuel class, LCC, and climate, subsetting to flamIndex, calculating n of ignitions
  fireSense_ignitionCovariates <- Map(
    f = fireSenseUtils::stackAndExtract,
    years = list(pre2012, post2012),
    fuel = list(fuelClasses$year2001, fuelClasses$year2011),
    LCC = list(LCCras$year2001, LCCras$year2011),
    MoreArgs = list(climate = ignitionClimate,
                    fires = sim$ignitionFirePoints)
  ) |>
    Cache(
      .functionName = "stackAndExtract",
      userTags = names(ignitionClimate)
    )

  fireSense_ignitionCovariates <- rbindlist(fireSense_ignitionCovariates)

  ## remove any pixels that are 0 for all classes
  fireSense_ignitionCovariates[, coverSums := rowSums(.SD),
                               .SD = setdiff(names(fireSense_ignitionCovariates),
                                             c(names(ignitionClimate), "cell", "ignitions", "year"))]
  fireSense_ignitionCovariates <- fireSense_ignitionCovariates[coverSums > 0]
  set(fireSense_ignitionCovariates, NULL, "coverSums", NULL)

  ## rename cells to pixelID - though aggregated raster is not saved
  setnames(fireSense_ignitionCovariates, old = "cell", new = "pixelID")
  fireSense_ignitionCovariates[, year := as.numeric(year)]

  # add lightning
  # https://www.tandfonline.com/doi/full/10.1080/07055900.2020.1845117
  set(fireSense_ignitionCovariates, NULL, "lightning",
      sim$lightningMaps[[2]][fireSense_ignitionCovariates[["pixelID"]]])

  ## for random effect
  if (grepl("xgb", Par$modelAlgorithm) %in% FALSE) {
    # ranEffs <- "yearChar"
    set(fireSense_ignitionCovariates, NULL, ranEffsLabel, as.character(fireSense_ignitionCovariates$year))
  }
  firstCols <- c("pixelID", "ignitions", names(ignitionClimate), "youngAge")
  firstCols <- firstCols[firstCols %in% names(fireSense_ignitionCovariates)]
  setcolorder(fireSense_ignitionCovariates, neworder = firstCols)

  sim$fireSense_ignitionCovariates <- fireSense_ignitionCovariates

  ## make new ignition object, ignitionFitRTM
  sim$ignitionFitRTM <- rast(fuelClasses$year2001[[1]])
  sim$ignitionFitRTM <- setValues(sim$ignitionFitRTM, 1) ## avoids a warning
  attributes(sim$ignitionFitRTM)$nonNAs <- nrow(sim$fireSense_ignitionCovariates)

  ## assign mean forest biomass- for use in plotting in ignitionFit
  tempCD <- LandR::addPixels2CohortData(sim$cohortData2011, sim$pixelGroupMap2011)
  bPerPixel <- tempCD[age > 0, .(bPerPixel = sum(B)), .(pixelIndex)]
  meanForestB <- mean(bPerPixel$bPerPixel)
  attributes(sim$ignitionFitRTM)$meanForestB <- meanForestB
  rm(tempCD, bPerPixel)

  ## build formula
  igCovariates <- names(sim$fireSense_ignitionCovariates)
  igCovariates <- igCovariates[!igCovariates %in%
                                 c(names(ignitionClimate),
                                   "year", "yearChar", "ignitions", "ignitionsNoGT1", "pixelID")]

  if (grepl("xgb", Par$modelAlgorithm) %in% FALSE) {
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
                                              # paste0(sim$climateVariablesForFire$ignition, collapse = " + "), " + ",
                                              # paste0(igCovariates, collapse = " + "), " + ",
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

  escapeThreshHa <- prod(res(sim$flammableRTM2001)) / 10000
  escapes <- sim$ignitionFirePoints[sim$ignitionFirePoints$SIZE_HA > escapeThreshHa, ]

  ## make a template aggregated raster - values are irrelevant, only need pixelID
  aggregatedRas <- terra::aggregate(sim$historicalClimateRasters[[1]][[1]],
                                    fact = P(sim)$igAggFactor, fun = mean) |>
    Cache(.functionName = "aggregate_historicalClimateRasters_forTemplate")

  coords <- st_coordinates(escapes)
  escapeCells <- cellFromXY(aggregatedRas, coords)
  escapeDT <- as.data.table(escapes)
  setnames(escapeDT, "YEAR", "year")
  escapeDT[, pixelID := escapeCells]
  escapeDT <- escapeDT[, .(year, pixelID)]
  escapeDT <- escapeDT[, .(escapes = .N), .(year, pixelID)]
  escapeDT[, year := as.numeric(year)]
  escapeDT <- escapeDT[sim$fireSense_ignitionCovariates, on = c("pixelID", "year")]
  escapeDT[is.na(escapes), escapes := 0]

  # ranEffs <- "yearChar"
  escapeVars <- names(escapeDT)[!names(escapeDT) %in% c("year", "pixelID", "escapes",
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
  fbldt <- rbindlist(fbldt, idcol = "year")
  fbldt <- fbldt[pixelID %in% unique(pixelIDsAllowed)]
  fireBufferedListDT <- split(fbldt, by = "year", keep.by = FALSE)
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
  cohDat <- "cohortData"
  pixGM <- "pixelGroupMap"
  saMap <- "standAgeMap"
  # rstLCC * see below
  neededYears <- c(2001, 2011)

  ecoFile <- ifelse(is.null(sim$ecoregionRst), "ecoregionLayer", "ecoregionRst")
  objsNeeded <- c(ecoFile,
                  "firePerimeters",
                  "rasterToMatch", "studyArea",
                  "rstLCC2011", "rstLCC2001",
                  "studyArea_biomassParam", "rasterToMatch_biomassParam", #needed by BBDP
                  "species", "speciesTable", "sppEquiv")
  objsNeeded <- intersect(ls(sim), objsNeeded)
  objsNeeded <- mget(objsNeeded, envir = envir(sim))
  cds <- lapply(neededYears, function(ny, objs = objsNeeded) {
    messageColoured(colour = "yellow", "Running Biomass_borealDataPrep for year ", ny)
    messageColoured(colour = "yellow", "  inside fireSense_dataPrepFit to estimate cohortData", ny)
    rstLCC <- objs[[paste0("rstLCC", ny)]]
    objs <- c(objs, "rstLCC" = rstLCC)
    #now duplicated
    objs[[paste0("rstLCC", ny)]] <- NULL
    otherYr <- ifelse(ny == 2001, 2011, 2001)
    objs[[paste0("rstLCC", otherYr)]] <- NULL
    parms <- list()
    for (nm in neededModule) {
      parms[[nm]] <- P(sim, module = nm)
      parms[[nm]][["dataYear"]] <- ny
      parms[[nm]][["forestedLCCClasses"]] <- P(sim)$forestedLCC
    }
    parms$Biomass_borealDataPrep$exportModels <- "none"
    outNY <- do.call(SpaDES.core::simInitAndSpades, list(paths = pathsLocal,
                                                         params = parms,
                                                         times = list(start = ny, end = ny),
                                                         modules = neededModule,
                                                         objects = objs)) |>
      Cache(.functionName = paste0("simInitAndSpades_insideFireSenseDataPrepFit", ny))
    cohDatObj <- paste0(cohDat, ny)
    pixGrpMap <- paste0(pixGM, ny)
    saObj <- paste0(saMap, ny)
    outNY[[cohDatObj]] <- outNY[[cohDat]]
    outNY[[pixGrpMap]] <- outNY[[pixGM]]
    outNY[[saObj]] <- outNY[[saMap]]
    mget(c(cohDatObj, pixGrpMap, saObj), envir = envir(outNY))
  })
  lapply(cds, function(cd) list2env(cd, envir = envir(sim)))
  sim
}

.inputObjects <- function(sim) {
  if (!suppliedElsewhere("studyArea", sim)) {
    sim$studyArea <- LandR::randomStudyArea(size = 10000 * 6.25 * 20000)
  }


  if (!suppliedElsewhere("studyArea_biomassParam", sim)) {
    if (is.null(sim$studyAreaLarge)) {
      sim$studyArea_biomassParam <- sim$studyArea
    } else {
      warning("please replace studyAreaLarge with studyArea_biomassParam")
      sim$studyArea_biomassParam <- sim$studyAreaLarge
    }
  }

  if (!suppliedElsewhere("sppEquiv", sim)) {
    sp <- LandR::speciesInStudyArea(studyArea = sim$studyArea)
    sp <- LandR::equivalentName(sp$speciesList, df = sppEquivalencies_CA, column = Par$sppEquivCol)
    sim$sppEquiv <- sppEquivalencies_CA[get(Par$sppEquivCol) %in% sp]
  }

  SpaDES.core::paramCheckOtherMods(sim, paramToCheck = "sppEquivCol")

  if (is.null(P(sim)$.studyAreaName)) {
    P(sim)$.studyAreaName <- studyAreaName(sim$studyArea)
  }
  cacheTags <- c(currentModule(sim), P(sim)$.studyAreaName)
  dPath <- asPath(inputPath(sim), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  if (!suppliedElsewhere("rasterToMatch", sim)) {
    sim$rasterToMatch <- Cache(LandR::prepInputsLCC(year = 2005, ## TODO: use 2010
                                                    destinationPath = dPath,
                                                    studyArea = sim$studyArea,
                                                    useCache = TRUE))
  }

  if (!suppliedElsewhere("rasterToMatch_biomassParam", sim)) {
    if (!is.null(sim$rasterToMatchLarge)) {
      warning("please use rasterToMatch_biomassParam in place of rasterToMatchLarge")
      sim$rasterToMatch_biomassParam <- sim$rasterToMatchLarge
    } else {
      sim$rasterToMatch_biomassParam <- sim$rasterToMatch
    }
  }

  if (!suppliedElsewhere("climateVariablesForFire", sim)) {
    sim$climateVariablesForFire <- list("spread" = "MDC",
                                        "ignition" = "MDC")
  }

  if (!suppliedElsewhere("rstLCC2001", sim)) {
    #use a threshold to to assign non-flammable cover (e.g. if < 10% flammable cover)
    LCC2001 <- Cache(makeFireSenseLCC,
                            neededYear = 2001,
                            writeTo = .suffix("rstLCC.tif",
                                              paste0(2001, "_", P(sim)$.studyAreaName)),
                            destinationPath = inputPath(sim),
                            studyArea = sim$studyArea_biomassParam,
                            rasterToMatch = sim$rasterToMatch_biomassParam,
                            nonflammableLCC = P(sim)$nonflammableLCC,
                            flammabilityThreshold = P(sim)$flammabilityThreshold,
                            userTags = c("makeFireSenseLCC", 2001))
    sim$rstLCC2001 <- LCC2001$lcc
    sim$propFlammable2001 <- LCC2001$flammableProp
  }

  if (!suppliedElsewhere("rstLCC2011", sim)) {
    #use a threshold to to assign non-flammable cover (e.g. if < 10% flammable cover)
    LCC2011 <- Cache(makeFireSenseLCC,
                     neededYear = 2011,
                     writeTo = .suffix("rstLCC.tif",
                                       paste0(2011, "_", P(sim)$.studyAreaName)),
                     destinationPath = inputPath(sim),
                     studyArea = sim$studyArea_biomassParam,
                     rasterToMatch = sim$rasterToMatch_biomassParam,
                     nonflammableLCC = P(sim)$nonflammableLCC,
                     flammabilityThreshold = P(sim)$flammabilityThreshold,
                     userTags = c("makeFireSenseLCC", 2011))
    sim$rstLCC2011 <- LCC2011$lcc
    sim$propFlammable2011 <- LCC2011$flammableProp
  }

  if (!suppliedElsewhere("standAgeMap2001", sim)) {
    sim$standAgeMap2001 <- Cache(prepInputsStandAgeMap,
                                 rasterToMatch = sim$rasterToMatch_biomassParam,
                                 studyArea = sim$studyArea_biomassParam,
                                 destinationPath = dPath,
                                 startTime = 2001,
                                 userTags = c(cacheTags, "prepInputsStandAgeMap2001"))
  }

  if (!suppliedElsewhere("standAgeMap2011", sim)) {
    sim$standAgeMap2011 <- Cache(prepInputsStandAgeMap,
                                 rasterToMatch = sim$rasterToMatch_biomassParam,
                                 studyArea = sim$studyArea_biomassParam,
                                 destinationPath = dPath,
                                 startTime = 2011,
                                 userTags = c(cacheTags, "prepInputsStandAgeMap2011"))
  }

  if (!all(suppliedElsewhere("cohortData2001", sim),
           suppliedElsewhere("cohortData2011", sim),
           suppliedElsewhere("pixelGroupMap2011", sim),
           suppliedElsewhere("pixelGroupMap2001", sim))) {
    ## This runs simInitAndSpades if needed
    sim <- runBorealDP_forCohortData(sim)
  }

  if (!P(sim)$useRasterizedFireForSpread) {
    if (!suppliedElsewhere("firePolys", sim) | !suppliedElsewhere("firePolysForAge", sim)) {
      ## don't want to needlessly postProcess the same firePolys objects

      saNotLatLong <- if (isTRUE(sf::st_is_longlat(sim$studyArea))) {
        terra::project(sim$studyArea, terra::crs(sim$rasterToMatch))
      } else {
        sim$studyArea
      }

      fireYears <- c(min(P(sim)$fireYears - P(sim)$cutoffForYoungAge):max(P(sim)$fireYears))
      ## TODO: check why this isn't resulting in identical crs between firePolys, studyArea
      allFirePolys <- Cache(fireSenseUtils::getFirePolygons,
                            fun = "sf::st_read",
                            years = fireYears,
                            useInnerCache = TRUE,
                            destinationPath = dPath,
                            cropTo = sim$rasterToMatch,
                            maskTo = saNotLatLong,
                            projectTo = sim$rasterToMatch,
                            userTags = c(cacheTags, "firePolys", paste0(fireYears, collapse = ":")))
    }

    if (!suppliedElsewhere("firePolys", sim)) {
      sim$firePolys <- allFirePolys[names(allFirePolys) %in% paste0("year", P(sim)$fireYears)]
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
          # cent <- terra::centroids(x)
          #TODO: switch to the above when terra conversion is complete
          cent <- sf::st_centroid(x)
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
    ignitionFirePoints <- Cache(
      getFirePoints_NFDB_V2,
      studyArea = sim$studyArea,
      years = P(sim)$fireYears,
      NFDB_pointPath = dPath,
      userTags = c("ignitionFirePoints", P(sim)$.studyAreaName),
      plot = !is.na(P(sim)$.plotInitialTime)
    ) ## default redownload means it will update annually - I think this is fine?
    ignitionFirePoints <- postProcess(ignitionFirePoints, projectTo = sim$rasterToMatch)
    sim$ignitionFirePoints <- ignitionFirePoints[ignitionFirePoints$CAUSE %in% c("L", "N"),]
    if (nrow(sim$ignitionFirePoints) == 0) {
      stop("no ignitions present - review getFirePoints-NFDB_V2")
      #this was happening with data update - the module will still run with no fire
    }
  }

  if (!suppliedElsewhere("historicalClimateRasters", sim)) {
    stop("please supply sim$historicalClimateRasters")
  }

  if (P(sim)$useRasterizedFire) {
    if (!suppliedElsewhere("historicalFireRaster", sim)) {
      sim$historicalFireRaster <- Cache(prepInputs,
                                        url = extractURL("historicalFireRaster", sim),
                                        rasterToMatch = sim$rasterToMatch,
                                        destinationPath = dPath,
                                        studyArea = sim$studyArea,
                                        method = "near", ## only use near or ngb; bilinear is wrong!
                                        filename2 = paste0("wildfire_", P(sim)$.studyAreaName, ".tif"),
                                        userTags = c("historicalFireRaster", P(sim)$.studyAreaName))
    }
  }

  if (!suppliedElsewhere("nonForestedLCCGroups", sim)) {
    vals <- freq(sim$rstLCC2011)
    forestOrNonFlamm <- sort(unique(c(Par$nonflammableLCC, Par$forestedLCC)))
    sim$nonForestedLCCGroups <- list(nf = setdiff(vals$value, forestOrNonFlamm))
    ## TODO: consider moving this to init - and checking if unsupplied
    # sim$nonForestedLCCGroups <- list(
    #   #"nf_dryland" = c(50, 100, 40), # shrub, herbaceous, bryoid
    #   #"nf_wetland" = c(80)), #non-treed wetland.
    #   "nf_highFlam" = c(50, 100, 40), # shrub, herbaceous
    #   "nf_lowFlam" = c(80)) # bryoids + non-treed wetland.
  }


  if (!suppliedElsewhere("missingLCCgroup", sim)) {
    sim$missingLCCgroup <- names(sim$nonForestedLCCGroups)[1]
  }

  if (!suppliedElsewhere("spreadFitAdditionalColNames")) {
    sim$spreadFitAdditionalColNames <- fireSenseUtils::spreadFitAdditionalColNames
  }

  # if (!suppliedElsewhere("fireSense_spreadLogisticTermNames")) {
  #   sim$fireSense_spreadLogisticTermNames <- c("inflectionPoint1", "inflectionPoint2",
  #                                              "maxAsymptote", "hillSlope1", "hillSlope2")
  #
  # }

  return(invisible(sim))
}

ranEffsLabel <- "yearChar"
