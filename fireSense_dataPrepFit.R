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
  loadOrder = list(after = c("Biomass_borealDataPrep", "Biomass_speciesParameters")),
  reqdPkgs = list("data.table", "fastDummies",
                  "PredictiveEcology/fireSenseUtils@lccFix (>= 0.0.5.9062)",
                  "ggplot2", "parallel", "purrr", "raster", "sf", "sp",
                  "PredictiveEcology/LandR@lccFix (>= 1.1.0.9081)",
                  "PredictiveEcology/SpaDES.core@development (>= 2.0.2.9006)",
                  "PredictiveEcology/SpaDES.project@transition",
                  "PredictiveEcology/SpaDES.tools (>= 2.0.4.9002)",
                  "snow", "terra"),
  parameters = bindrows(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter("areaMultiplier", c("numeric", "function"), fireSenseUtils::multiplier, NA, NA,
                    paste("Either a scalar that will buffer `areaMultiplier * fireSize` or a function",
                          "of `fireSize`. See `?fireSenseUtils::bufferToArea`.")),
    defineParameter("bufferForFireRaster", "numeric", 1000, 0, NA,
                    paste("The distance that determine whether separate patches of burned pixels originated",
                          "from the same fire. Only relevant when `useRasterizedFireForSpread = TRUE`.",
                          "This param is separate from `minBufferSize`, which is used to determine the",
                          "minimum sample of burned and unburned pixels to include in each fire.")),
    defineParameter("cutoffForYoungAge", "numeric", 15, NA, NA,
                    "Age at and below which pixels are considered 'young' (`young <- age <= cutoffForYoungAge`)"),
    defineParameter("fireYears", "integer", 2002:2021, NA, NA,
                    paste("A numeric vector indicating which years should be extracted",
                          "from the fire databases to use for fitting")),
    defineParameter("forestedLCC", "numeric", c(210, 220, 230, 240), NA, NA,
                    paste("Forested land cover classes - these differ from non-forest because the biomass", 
                          "and composition of fuels are taken into account by fireSense, while non-forest",
                          "classes are treated categorically")),
    defineParameter("igAggFactor", "numeric", 40, 1, NA,
                    "aggregation factor for rasters during ignition prep."),
    defineParameter("ignitionFuelClassCol", "character", "FuelClass", NA, NA,
                    "the column in `sppEquiv` that defines unique fuel classes for ignition. A column ",
                    "named 'FuelClass' exists in the LandR::sppEquivalencies_CA and will be used ",
                    "by default. To change the FuelClass classifications, add a column to that table, ",
                    "or to `sim$sppEquiv` and then modify this `ignitionFuelClassCol` parameter"),
    defineParameter("minBufferSize", "numeric", 5000, NA, NA,
                    paste("Minimum number of cells in buffer and nonbuffer. This is imposed after the",
                          "multiplier on the `bufferToArea` fn")),
    defineParameter("missingLCCgroup", "character", "nf_highFlam", NA, NA,
                    paste("if a pixel is forested but is absent from `cohortData`, it will be grouped in this class.",
                          "Must be one of the names in `sim$nonForestedLCCGroups`")),
    defineParameter("nonflammableLCC", "numeric", c(20, 31, 32, 33), NA, NA,
                    "non-flammable LCC in rstLCC layers - defaulting to water, snow/ice, rock, barren land."),
    defineParameter("nonForestCanBeYoungAge", "logical", TRUE, NA, NA,
                    paste("if TRUE, burned non-forest will be treated as `youngAge`. Recommended to be TRUE", 
                          "as burned forest is often classified as non-forest")),
    defineParameter("sppEquivCol", "character", "LandR", NA, NA,
                    "column name in `sppEquiv` object that defines unique species in `cohortData`"),
    defineParameter("spreadFuelClassCol", "character", "FuelClass", NA, NA,
                    "if using fuel classes for spread, the column in `sppEquiv` that defines unique ",
                    "fuel classes. The column 'FuelClass' in LandR::sppEquivalencies_CA and will be used ",
                    "by default. To change the FuelClass classifications, add a column to that table, ",
                    "or to `sim$sppEquiv` and then modify this `spreadFuelClassCol` parameter"),
    defineParameter("useCentroids", "logical", TRUE, NA, NA,
                    paste("Should fire ignitions start at the `sim$firePolygons` centroids",
                          "or at the ignition points in `sim$firePoints`?")),
    defineParameter("usePiecewiseRegression", "logical", FALSE, NA, NA,
                    paste("Should fire ignitions fitting use the hockey stick 'piecewise regression' ",
                          "approach or the newer glmmTMB with zero inflated poisson mixed effect.")),
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
                       "This includes all fires regardless of size")),
    expectsInput("nonForestedLCCGroups", "list",
                 paste("a named list of non-forested landcover groups",
                       "e.g. list('wetland' = c(19, 23, 32))",
                       "These will become covariates in `fireSense_IgnitionFit`")),
    expectsInput("pixelGroupMap2001", "SpatRaster", sourceURL = NA,
                 "defines the `pixelGroups` for cohortData table in 2001"),
    expectsInput("pixelGroupMap2011", "SpatRaster",
                 "defines the `pixelGroups` for cohortData table in 2011"),
    expectsInput("rasterToMatch", "SpatRaster", sourceURL = NA,
                 "template raster for study area. Assumes some buffering of core area to limit edge effect of fire."),
    expectsInput("rstLCC2001", "SpatRaster", sourceURL = NA,
                 "Raster of land cover - will use Biomass_borealDataPrep to generate if missing."),
    expectsInput("rstLCC2011", "SpatRaster", sourceURL = NA,
                 "Raster of land cover - will use Biomass_borealDataPrep to generate if missing."),
    expectsInput("sppEquiv", "data.table", sourceURL = NA,
                 "table of LandR species equivalencies"),
    expectsInput("standAgeMap2001", "SpatRaster", sourceURL = NA,
                 "map of stand age in 2001 used to create `cohortData2001`"),
    expectsInput("standAgeMap2011", "SpatRaster", sourceURL = NA,
                 "map of stand age in 2011 used to create `cohortData2011`"),
    expectsInput("studyArea", "sf", sourceURL = NA,
                 "studyArea that determines spatial boundaries of all data")
  ),
  outputObjects = bindrows(
    createsOutput("fireBufferedListDT", "list",
                  "list of data.tables with fire id, `pixelID`, and buffer status"),
    createsOutput("firePolys", "list",
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
                  "data.table with `pixelID` and relevant landcover classes in 2001"),
    createsOutput("landcoverDT2011", "data.table",
                  "data.table with `pixelID` and relevant landcover classes in 2011"),
    createsOutput("nonForest_timeSinceDisturbance2001", "SpatRaster",
                  "time since burn for non-forested pixels in 2001"),
    createsOutput("nonForest_timeSinceDisturbance2011", "SpatRaster",
                  "time since burn for non-forested pixels in 2011"),
    createsOutput("flammableRTM2001", "SpatRaster", "binary raster of flammable landcover for 2001"),
    createsOutput("flammableRTM2011", "SpatRaster", "binary raster of flammable landcover for 2011"),
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
        #the camelcase is still different with FS from LandR Biomass
      }

      # schedule future event(s)
      if ("fireSense_IgnitionFit" %in% P(sim)$whichModulesToPrepare)
        sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "prepIgnitionFitData", eventPriority = 1)
      if ("fireSense_EscapeFit" %in% P(sim)$whichModulesToPrepare)
        sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "prepEscapeFitData", eventPriority = 1)
      if ("fireSense_SpreadFit" %in% P(sim)$whichModulesToPrepare) {
        sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "prepSpreadFitData", eventPriority = 1)
      }

      # do stuff for this event
      sim <- Init(sim)

      sim <- scheduleEvent(sim, end(sim), "fireSense_dataPrepFit", "plotAndMessage", eventPriority = 9)
      sim <- scheduleEvent(sim, start(sim), "fireSense_dataPrepFit", "cleanUp", eventPriority = 10) #cleans up Mod objects
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

## event functions
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
Init <- function(sim) {
 
  #TODO: standardize sim$climateVariablesForFire if user provided
  #this approach will be wrong if they pass a list length one... 
  if (length(sim$climateVariablesForFire) == 1) {
    sim$climateVariablesForFire <- list("ignition" = sim$climateVariablesForFire,
                                        "spread" = sim$climateVariablesForFire)
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
    # project it first, faster than the postProcessTo sequence pre-crop, project, mask, crop
    sim$ignitionFirePoints <- projectTo(sim$ignitionFirePoints, sim$rasterToMatch) |>
      cropTo(sim$rasterToMatch) |>
      maskTo(sim$rasterToMatch)
  }
  
  #possible if user-supplied
  if (!terra::same.crs(sim$firePolys[[1]], sim$rasterToMatch)) {
    # projectTo(fp, crs(sim$rasterToMatch))) |> #terrarize
    sim$firePolys <- Map(fp = sim$firePolys, function(fp)
      projectTo(fp, st_crs(sim$rasterToMatch))) |>
      Cache(.functionName = "projectTo_for_firePolys")
  }

  
  ## sanity checks
  if (!LandR::.compareRas(sim$standAgeMap2001, sim$standAgeMap2011, sim$rasterToMatch,
                          stopOnError = FALSE)) {
    sim$standAgeMap2001 <- postProcess(sim$standAgeMap2001, cropTo = sim$rasterToMatch,
                                       projectTo = sim$rasterToMatch, maskTo = mod$studyAreaUnion)
    sim$standAgeMap2011 <- postProcess(sim$standAgeMap2011, cropTo = sim$rasterToMatch,
                                        projectTo = sim$rasterToMatch, maskTo = mod$studyAreaUnion)
  }

  if (!P(sim)$missingLCCgroup %in% names(sim$nonForestedLCCGroups)) {
    stop("please review names of sim$nonForestedLCCGroups and param missingLCCGroup")
  }
  
  igFuels <- sim$sppEquiv[[P(sim)$ignitionFuelClassCol]]
  spreadFuels <- sim$sppEquiv[[P(sim)$spreadFuelClassCol]]

  if (any(c(is.null(spreadFuels), is.na(spreadFuels),
            is.null(igFuels), is.na(igFuels)))) {
    stop("All species must have spread and ignition fuelClasses defined")
  }

  #TODO: make this table multidimensional or a list
  sim$landcoverDT2001 <- makeLandcoverDT(rstLCC = sim$rstLCC2001, flammableRTM = sim$flammableRTM2001,
                                         forestedLCC = P(sim)$forestedLCC, sim$nonForestedLCCGroups)
  sim$landcoverDT2011 <- makeLandcoverDT(rstLCC = sim$rstLCC2011, flammableRTM = sim$flammableRTM2011,
                                         forestedLCC = P(sim)$forestedLCC, sim$nonForestedLCCGroups)

  sim$landcoverDT2001 <- correctMissingLCC(sim$landcoverDT2001, sim$pixelGroupMap2001, P(sim)$missingLCCgroup)
  sim$landcoverDT2011 <- correctMissingLCC(sim$landcoverDT2011, sim$pixelGroupMap2011, P(sim)$missingLCCgroup)
  
  # cannot merge because before subsetting due to column differences over time

  ## TODO: this object  be used to track annual youngAge of all pixels, forested or not
  #so "nonForest" is a poor choice of name 
  sim$nonForest_timeSinceDisturbance2001 <- makeTSD(year = 2001,
                                                    fireRaster = sim$historicalFireRaster, #can be NULL
                                                    firePolys = sim$firePolysForAge,
                                                    standAgeMap = sim$standAgeMap2001, lcc = sim$landcoverDT2001,
                                                    cutoffForYoungAge = P(sim)$cutoffForYoungAge)
  sim$nonForest_timeSinceDisturbance2011 <- makeTSD(year = 2011, 
                                                    fireRaster = sim$historicalFireRaster, #can be NULL 
                                                    firePolys = sim$firePolysForAge,
                                                    standAgeMap = sim$standAgeMap2011,
                                                    lcc = sim$landcoverDT2011,
                                                    cutoffForYoungAge = P(sim)$cutoffForYoungAge)

  #until youngAge treatment is identical between spread and ignition, no point in prepping veg here
  #Currently youngAge is resolved annually in spread, but only once in ignition
  #e.g. if a pixel ignited in 2008, its youngAge status in ignition is still determined by whether it was 15 in 2001,
  #but its youngAge status for spread is deterimined by whether standAge < 15 in 2008

  return(invisible(sim))
}

prepare_SpreadFit <- function(sim) {

  ## Put in format for DEOptim that distinguishes annual and nonannual covariates
  ## Prepare annual spread fit covariates
  ####prep veg data####
  doAssertion <- getOption("fireSenseUtils.assertions", TRUE)

  ## sanity check the inputs
  compareGeom(sim$rasterToMatch, sim$flammableRTM2011, sim$flammableRTM2001)
  compareGeom(sim$rasterToMatch, sim$standAgeMap2001, sim$standAgeMap2011)
  lapply(sim$historicalClimateRasters, compareGeom, x = sim$rasterToMatch)
  
  #when landcoverDT is included, as is the case here, non-forest pixels in cohortData are masked out
  #this is necessary when LandR and fireSense have differing concepts of non-forest
  vegData <- Map(f = cohortsToFuelClasses,
                 cohortData = list(sim$cohortData2001, sim$cohortData2011),
                 pixelGroupMap = list(sim$pixelGroupMap2001, sim$pixelGroupMap2011),
                 landcoverDT = list(sim$landcoverDT2001, sim$landcoverDT2011), 
                 flammableRTM = list(sim$flammableRTM2001, sim$flammableRTM2011),
                 MoreArgs = list(sppEquiv = sim$sppEquiv,
                                 sppEquivCol = P(sim)$sppEquivCol,
                                 fuelClassCol = P(sim)$spreadFuelClassCol,
                                 cutoffForYoungAge = -1)) |>
    Cache(.functionName = "cohortsToFuelClasses") #youngAge will be resolved annually downstream
    
  vegData <- lapply(vegData, FUN = function(x) {
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

  lccNames <- setdiff(names(vegData), c("pixelID", "year"))
  

  #prep the fire data
  if (P(sim)$useRasterizedFire) {
    sim <- prepare_SpreadFitFire_Raster(sim)
  } else {
    sim <- prepare_SpreadFitFire_Vector(sim)
  }

  ####join fire and veg data ####
  pre2012 <- paste0("year", min(P(sim)$fireYears):2011)
  post2012 <- paste0("year", 2012:max(P(sim)$fireYears))
  #these fire year objects get a bit tedious because they need to be in integer, string, 
  #with and without missing fire years
  
  pre2012Indices <- sim$fireBufferedListDT[pre2012] %>%
    rbindlist(.) %>%
    .[vegData[year < 2012], on = c("pixelID")] %>%
    .[!is.na(buffer)] #some buffered pixels are non-flammable...
  #TODO: discuss if this is expected (as far as I can tell, it is)
 
  post2012Indices <- sim$fireBufferedListDT[post2012] %>%
    rbindlist(.) %>%
    .[vegData[year >= 2012], on = c("pixelID")] %>%
    .[!is.na(buffer)] 
  rm(vegData)
  gc()

  #TODO: lines from creation of vegData onwards should be reviewed. Seems redundant.. 
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
    ratioZeroToOne <- ttt[1]/ttt[2]
    if (ratioZeroToOne < 5)
      stop("The number of pixels in the fire buffers should be at least 5x the number of burned pixels\n",
           "Please create larger buffers around fires in fireBufferedListDT, e.g., via ",
           "fireSenseUtils::bufferToArea(..., areaMultiplier = multiplier)")
  }

  RHS <- paste(paste0(names(sim$historicalClimateRasters), collapse = " + "), "youngAge",
               paste0(vegCols, collapse = " + "), sep =  " + ")

  ## this is a funny way to get years but avoids years with 0 fires
  years <- paste0("year", P(sim)$fireYears)
  yearsWithFire <- years[years %in% names(sim$fireBufferedListDT)]
  pre2012int <- as.integer(min(P(sim)$fireYears):2010)
  post2012int <- as.integer(2011:max(P(sim)$fireYears))
  pre2012 <- yearsWithFire[yearsWithFire %in% paste0("year", pre2012int)]
  post2012 <- yearsWithFire[yearsWithFire %in% paste0("year", post2012int)]

  ##climate ###
  
  #TODO: index removed as argument - so if this breaks - this is the issue
  #we would need to modify this object to account for two sets of flammable indices
  # as flammable pixels change between 2001 and 2011 (mainly water) 
  spreadClimate <- sim$historicalClimateRasters[sim$climateVariablesForFire$spread]

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
    .[!duplicated(pixelID),]
 

  nonAnnualpost2012 <- fireSenseVegData[year >= 2012, .SD, .SDcols = colsToExtract] %>%
    na.omit(.) %>%
    as.data.table(.) %>%
    .[!duplicated(pixelID)] ## remove duplicates from same pixel diff year

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
 
  if (!P(sim)$nonForestCanBeYoungAge) {
    #TODO: test this inversion of makeMutuallyExclusive's regular use
    args <- as.list(rep("youngAge", length = length(sim$nonForestedLCCGroups)))
    names(args) <- names(sim$nonForestedLCCGroups)
  } else {
    #this is done later in spreadFit - but done here for accuracy of outputs
    args <- list("youngAge" = names(sim$nonForestedLCCGroups))
  }  
  annualCovariates <- lapply(annualCovariates, makeMutuallyExclusive, 
                             mutuallyExclusiveCols = args)

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
  #TODO: do this, obviously
  historicalFireRaster <- sim$historicalFireRaster

  #build initial burn IDs by buffering  - then using clump(raster) or patches(terra)
  #historical fire Raster is currently not in outputs - if assigned to sim here, it should be added
  # as we modify it by removing non-flammable fires.

  historicalFireRaster <- mask(historicalFireRaster, sim$flammableRTM,
                               maskvalues = 0, updatevalue = NA)

  nCores <- ifelse(grepl("Windows", Sys.info()[["sysname"]]), 1L, length(sim$fireYears))

  #this is analogous to buffer to area but for raster datasets as opposed to polygon
  #the inner looping function is very similar - one difference is that non-flammable
  #pixels do not count toward the buffer size, unlike the polygonal version.
  sim$fireBufferedListDT <- Cache(rasterFireBufferDT, years =  P(sim)$fireYears,
                                  fireRaster = historicalFireRaster, flammableRTM = sim$flammableRTM,
                                  bufferForFireRaster = P(sim)$bufferForFireRaster, verb = 1,
                                  areaMultiplier = P(sim)$areaMultiplier, minSize = P(sim)$minBufferSize,
                                  cores = nCores, userTags = c(currentModule(sim), "rasterFireBufferDT"))
  #TODO: test that this is the correct method for missing years
  missingYears <- unlist(lapply(sim$fireBufferedListDT, is.null))

  if (any(missingYears)) {
    actualFireYears <- P(sim)$fireYears[!missingYears]
    sim$fireBufferedListDT <- sim$fireBufferedListDT[!missingYears]
  }

  #next up: generate spread fire points. no harmonization is needed with this approach :)
  sim$spreadFirePoints <- lapply(sim$fireBufferedListDT,
                                 rasterFireSpreadPoints,
                                 flammableRTM = sim$flammableRTM)

  #TODO: this is temporary while we migrate out of spatial/raster constructs
  #the "year" prefix is added by fireSenseUtils::makeLociList - discuss what to do
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
  ## sanity check
  ## TODO: is there a terra version of st_contains?
  #this takes a long time and shouldn't be necessary if postProcess is functional
  # stopifnot(
  #   "all annual firePolys are not within studyArea" = all(unlist(lapply(sim$firePolys, function(x) {
  #     SA <- st_as_sf(mod$studyAreaUnion)
  #     x <- st_as_sf(x)
  #     length(sf::st_contains(SA, x)) == 1
  #   })))
  # )

  ####prep fire data ####
  if (is.null(sim$firePolys[[1]]$FIRE_ID)) {
    stop("firePolys needs a numeric FIRE_ID column")
  }

  if (!is.numeric(sim$firePolys[[1]]$FIRE_ID) | !is.numeric(sim$spreadFirePoints[[1]]$FIRE_ID)) {

    message("need numeric FIRE_ID column in fire polygons and points. Coercing to numeric...")
    #this is true of the current NFBB
    origNames <- names(sim$firePolys)
    PointsAndPolys <- lapply(names(sim$firePolys),
                             function(year, polys = sim$firePolys, points = sim$spreadFirePoints) {
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
    sim$firePolys <- lapply(PointsAndPolys, FUN = function(x) x[["polys"]])
    rm(PointsAndPolys)
    names(sim$firePolys) <- origNames
    names(sim$spreadFirePoints) <- origNames
  }

  ## drop fires less than 1 px in size
  pixSizeHa <- prod(res(sim$flammableRTM2011)) / 1e4
  #using x[x$SIZE_HA] will work with terra or sf, while subset will not (I believe...)
  sim$spreadFirePoints <- lapply(sim$spreadFirePoints, function(x, minSize = pixSizeHa) {
    x <- x[x$SIZE_HA > minSize,]
    if (nrow(x) > 0) x else NULL
  })
  
  #silly R  
  sim$spreadFirePoints[sapply(sim$spreadFirePoints, is.null)] <- NULL

  sim$firePolys <- lapply(sim$firePolys, function(x) {
    x <- x[x$SIZE_HA > pixSizeHa,]
    if (nrow(x) > 0) x else NULL
  })
  sim$firePolys[sapply(sim$firePolys, is.null)] <- NULL

  #years run separately because flammableRTM is different
  #ultimately this function should combine the climate data to avoid needless iteration, 
  #and even this duplicated step should be a function of "fire period" for >2 periods 
  #however, the rasterized fire prep is significantly different, and needs review first 
  harmonized2001 <- Cache(harmonizeFireData,
                          firePolys = sim$firePolys[names(sim$firePolys) %in% pre2012], #protects from missing years
                          flammableRTM = sim$flammableRTM2001,
                          spreadFirePoints = sim$spreadFirePoints[pre2012], 
                          areaMultiplier = P(sim)$areaMultiplier, minSize = P(sim)$minBufferSize,
                          pointsIDcolumn = "FIRE_ID",
                          userTags = c("harmonizeFireData", P(sim)$.studyAreaName, "2001"))
  harmonized2011 <- Cache(harmonizeFireData,
                          sim$firePolys[names(sim$firePolys) %in% post2012], 
                          sim$flammableRTM2011, 
                          spreadFirePoints = sim$spreadFirePoints[names(sim$spreadFirePoints) %in% post2012],
                          areaMultiplier = P(sim)$areaMultiplier, minSize = P(sim)$minBufferSize,
                          pointsIDcolumn = "FIRE_ID",
                          userTags = c("harmonizeFireData", P(sim)$.studyAreaName, "2011"))
  
  sim$fireBufferedListDT <- append(harmonized2001$fireBufferedListDT, 
                                   harmonized2011$fireBufferedListDT)
  sim$spreadFirePoints <- append(harmonized2001$spreadFirePoints, harmonized2011$spreadFirePoints)
  sim$firePolys <- append(harmonized2001$firePolys, harmonized2011$firePolys)
  
  ## drop fire years from these lists that don't have any buffer points post-harmonization
  #TOOD: is this necessary?
  omitYears <- sapply(sim$spreadFirePoints, is.null)
  sim$fireBufferedListDT[omitYears] <- NULL
  sim$firePolys[omitYears] <- NULL
  sim$spreadFirePoints[omitYears] <- NULL

  return(invisible(sim))
}

prepare_IgnitionFit <- function(sim) {
  
  stopifnot(
    "all ignitionFirePoints are not within studyArea" = identical(
      nrow(st_as_sf(sim$ignitionFirePoints)),
      nrow(st_intersection(st_as_sf(sim$ignitionFirePoints), st_as_sf(mod$studyAreaUnion)))
    )
  )

  # account for forested pixels that aren't in cohortData
  #TODO: make this elegant 
  ## first put landcover into raster stack
  ## non-flammable pixels require zero values for non-forest landcover, not NA
  LCCras <- Map(
    f = putBackIntoRaster,
    landcoverDT = list(sim$landcoverDT2001, sim$landcoverDT2011),
    flammableMap = list(sim$flammableRTM2001, sim$flammableRTM2011),
    MoreArgs = list(lcc = names(sim$nonForestedLCCGroups))) |>
    Cache(.functionName = "putBackIntoRaster",
          userTags = c("putBackIntoRaster", P(sim)$.studyAreaName))
  
  fuelClasses <- Map(f = cohortsToFuelClasses,
                     cohortData = list(sim$cohortData2001, sim$cohortData2011),
                     flammableRTM = list(sim$flammableRTM2001, sim$flammableRTM2011),
                     landcoverDT = list(sim$landcoverDT2001, sim$landcoverDT2011),
                     pixelGroupMap = list(sim$pixelGroupMap2001, sim$pixelGroupMap2011),
                     MoreArgs = list(sppEquiv = sim$sppEquiv,
                                     sppEquivCol = P(sim)$sppEquivCol,
                                     fuelClassCol = P(sim)$ignitionFuelClassCol,
                                     cutoffForYoungAge = P(sim)$cutoffForYoungAge)) |>
    Cache(.functionName = "cohortsToFuelClasses")
  
  if (P(sim)$nonForestCanBeYoungAge) {
    ## this modifies the NF landcover by converting some NF to a new YA layer
    ## it must be done before aggregating

    LCCras <- Map(f = calcNonForestYoungAge,
                  landcoverDT = list(sim$landcoverDT2001, sim$landcoverDT2011),
                  NFTSD = list(sim$nonForest_timeSinceDisturbance2001,
                               sim$nonForest_timeSinceDisturbance2011),
                  LCCras = list(LCCras[[1]], LCCras[[2]]),
                  MoreArgs = list(cutoffForYoungAge = P(sim)$cutoffForYoungAge)) |>
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
  
  LCCras <- lapply(LCCras, aggregate, fact = P(sim)$igAggFactor, fun = mean) |>
    Cache(.functionName = "aggregate_LCCras_to_coarse")
  names(LCCras) <- c("year2001", "year2011")
  fuelClasses <- lapply(fuelClasses, FUN = aggregate, fact = P(sim)$igAggFactor, fun = mean) |>
    Cache(.functionName = "aggregate_fuelClasses_to_coarse")
  names(fuelClasses) <- c("year2001", "year2011")
  
  ignitionClimate <- sim$historicalClimateRasters[sim$climateVariablesForFire$ignition]
  ignitionClimate <- lapply(X = ignitionClimate, FUN = terra::aggregate, 
                            fact = P(sim)$igAggFactor, fun = mean) |>
    Cache(.functionName = "aggregate_historicalClimateRasters_to_coarse")
  #safety
  compareGeom(ignitionClimate[[1]], fuelClasses[[1]], fuelClasses[[2]])
  
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
  
  #TODO: the cache bevavior is too permissive
  #this is joining fuel class, LCC, and climate, subsetting to flamIndex, calculating n of ignitions
  fireSense_ignitionCovariates <- Map(f = fireSenseUtils::stackAndExtract,
                                      years = list(pre2012, post2012),
                                      fuel = list(fuelClasses$year2001, fuelClasses$year2011),
                                      LCC = list(LCCras$year2001, LCCras$year2011),
                                      MoreArgs = list(climate = ignitionClimate,
                                                      fires = sim$ignitionFirePoints
                                      )) |> Cache(.functionName = "stackAndExtract", 
                                                  userTags = names(ignitionClimate))
  
  fireSense_ignitionCovariates <- rbindlist(fireSense_ignitionCovariates)
  
  #remove any pixels that are 0 for all classes
  fireSense_ignitionCovariates[, coverSums := rowSums(.SD),
                               .SD = setdiff(names(fireSense_ignitionCovariates),
                                             c(names(ignitionClimate), "cell", "ignitions", "year"))]
  fireSense_ignitionCovariates <- fireSense_ignitionCovariates[coverSums > 0]
  set(fireSense_ignitionCovariates, NULL, "coverSums", NULL)
  
  
  #rename cells to pixelID - though aggregated raster is not saved
  setnames(fireSense_ignitionCovariates, old = "cell", new = "pixelID")
  fireSense_ignitionCovariates[, year := as.numeric(year)]
  
  # for random effect
  ranEffs <- "yearChar"
  set(fireSense_ignitionCovariates, NULL, ranEffs, as.character(fireSense_ignitionCovariates$year))
  firstCols <- c("pixelID", "ignitions", names(ignitionClimate), "youngAge")
  firstCols <- firstCols[firstCols %in% names(fireSense_ignitionCovariates)]
  setcolorder(fireSense_ignitionCovariates, neworder = firstCols)
  
  if (isTRUE(P(sim)$usePiecewiseRegression)) {
    response <- "ignitions"
  } else {
    response <- "ignitionsNoGT1"
    set(fireSense_ignitionCovariates, NULL, response, pmin(fireSense_ignitionCovariates$ignitions, 1))
    # fireSense_ignitionCovariates[, ignitionsNoGT1 := ifelse(ignitions > 1, 1, ignitions)]
  }
  
  sim$fireSense_ignitionCovariates <- fireSense_ignitionCovariates
  
  #make new ignition object, ignitionFitRTM
  sim$ignitionFitRTM <- rast(fuelClasses$year2001[[1]])
  sim$ignitionFitRTM <- setValues(sim$ignitionFitRTM, 1) #avoids a warning
  attributes(sim$ignitionFitRTM)$nonNAs <- nrow(sim$fireSense_ignitionCovariates)
  
  #assign mean forest biomass- for use in plotting in ignitionFit# 
  tempCD <- LandR::addPixels2CohortData(sim$cohortData2011, sim$pixelGroupMap2011)
  bPerPixel <- tempCD[age > 0, .(bPerPixel = sum(B)), .(pixelIndex)]
  meanForestB <- mean(bPerPixel$bPerPixel)
  attributes(sim$ignitionFitRTM)$meanForestB <- meanForestB
  rm(tempCD, bPerPixel)
  
  #build formula
  igCovariates <- names(sim$fireSense_ignitionCovariates)
  igCovariates <- igCovariates[!igCovariates %in%
                                 c(names(ignitionClimate), 
                                   "year", "yearChar", "ignitions", "ignitionsNoGT1", "pixelID")]
  
  #this is safer for multiple climate variables 
  interactionsDF <- as.data.table(expand.grid(igCovariates, sim$climateVariablesForFire$ignition))
  interactionsDF[, interaction := do.call(paste, c(.SD, sep = ":")), .SDcols = names(interactionsDF)]
  interactions <- interactionsDF$interaction
  
  #sanity check for base::abbreviate
  if (isTRUE(P(sim)$usePiecewiseRegression)) {
    if (length(sim$climateVariablesForFire$ignition) > 1) {
      stop("cannot use multiple climate variables with piecewise ignition fit formula")
    }
    
    pwNames <- abbreviate(igCovariates, minlength = 3, use.classes = TRUE, strict = FALSE)
    
    pw <- paste0(igCovariates, ":", "pw(", sim$climateVariablesForFire$ignition, ", k_", pwNames, ")")
    if (!all(c(length(unique(pw)), length(unique(interactions))) == length(igCovariates))) {
      warning("automated ignition formula construction needs review")
    }
    #don't overwrite if it exists
    if (is.null(sim$fireSense_ignitionFormula)) {
      sim$fireSense_ignitionFormula <- paste0(response, " ~ ", paste0(interactions, collapse = " + "), " + ",
                                              paste0(pw, collapse  = " + "), "- 1")
    }
    
  } else {
    if (!length(unique(interactions)) == length(igCovariates) * length(sim$climateVariablesForFire$ignition)) {
      warning("automated ignition formula construction needs review")
    }
    if (is.null(sim$fireSense_ignitionFormula)) {
      sim$fireSense_ignitionFormula <- paste0(response, " ~ ",
                                              paste0("(1|", ranEffs, ")"), " + ",
                                              paste0(sim$climateVariablesForFire$ignition, collapse = " + "), " + ",
                                              paste0(igCovariates, collapse = " + "), " + ",
                                              paste0(interactions, collapse = " + "))
    }
  }
  return(invisible(sim))
}

prepare_EscapeFit <- function(sim) {
  if (is.null(sim$fireSense_ignitionCovariates)) {
    #the datasets are essentially the same, with one column difference
    stop("Please include ignitionFit in parameter 'whichModulesToPrepare' if running EscapeFit")
  }

  escapeThreshHa <- prod(res(sim$flammableRTM2001))/10000
  escapes <- sim$ignitionFirePoints[sim$ignitionFirePoints$SIZE_HA > escapeThreshHa,]

  #make a template aggregated raster - values are irrelevant, only need pixelID
  aggregatedRas <- aggregate(sim$historicalClimateRasters[[1]][[1]],
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

  sim$fireSense_escapeCovariates <- escapeDT

  escapeVars <- names(escapeDT)[!names(escapeDT) %in% c("year", "pixelID", "escapes", "ignitions")]
  LHS <- paste0("cbind(escapes, ignitions - escapes) ~ ")
  RHS <- paste0(escapeVars, collapse = " + ")
  
  if (is.null(sim$fireSense_escapeFormula)) {
    sim$fireSense_escapeFormula <- paste0(LHS, RHS, " - 1")
  }

  if (any(sim$fireSense_escapeCovariates$escapes > sim$fireSense_escapeCovariates$ignitions)) {
    stop("issue with escapes outnumbering ignitions in a pixel - contact module creators")
  }

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
}

### template for plot events
plotAndMessage <- function(sim) {
  #TODO: this could plot the ignition/spread covariates
  return(invisible(sim))
}

rmMissingPixels <- function(fbldt, pixelIDsAllowed)  {
  fbldt <- rbindlist(fbldt, idcol = "year")
  fbldt <- fbldt[pixelID %in% unique(pixelIDsAllowed)]
  fireBufferedListDT <- split(fbldt, by = "year", keep.by = FALSE)
}

runBorealDP_forCohortData <- function(sim) {
  ## Biomass_species should be only run if it is already in the simList
  neededModule <- "Biomass_borealDataPrep"
  if ("Biomass_speciesData" %in% modules(sim)) {
    neededModule <- c("Biomass_borealDataPrep", "Biomass_speciesData")
  }
  
  # neededModule <- "Biomass_borealDataPrep"
  pathsLocal <- paths(sim)
  if (any(!neededModule %in% modules(sim))) {
    ## don't install pkgs mid-stream; already use module metadata to declare pkgs for installation
    # Require::Install("PredictiveEcology/SpaDES.project@transition")
    modulePathLocal <- file.path(modulePath(sim), currentModule(sim), "submodules")
    getModule(file.path("PredictiveEcology", paste0(neededModule, "@development")),
              modulePath = modulePathLocal, overwrite = FALSE)
    pathsLocal$modulePath <- modulePathLocal
  }
  cohDat <- "cohortData"
  pixGM <- "pixelGroupMap"
  saMap <- "standAgeMap"
  rstLCC <- "rstLCC"
  neededYears <- c(2001, 2011)
  
  if (!is.null(sim$cohortData)) {
    alreadyDone <- P(sim, "dataYear", "Biomass_borealDataPrep")
    cohDatObj <- paste0(cohDat, alreadyDone)
    pixGrpMap <- paste0(pixGM, alreadyDone)
    saObj <- paste0(saMap, alreadyDone)
    rstLCCobj <- paste0(rstLCC, alreadyDone)

    sim[[cohDatObj]] <- sim[[cohDat]]
    sim[[pixGrpMap]] <- sim[[pixGM]]
    sim[[saObj]] <- sim[[saMap]]
    sim[[rstLCCobj]] <- sim[[rstLCC]]
    
    messageColoured(colour = "yellow", "fireSense_dataPrepFit will use estimates of ",
                    paste0("cohortData", alreadyDone, collapse = ", "), " from modules already run")
    
    neededYears <- setdiff(neededYears, alreadyDone)
  }
  
  if (is.null(sim$studyAreaLarge)) {
    sim$studyAreaLarge <- sim$studyArea
  }
  
  ecoFile <- ifelse(is.null(sim$ecoregionRst), "ecoregionLayer", "ecoregionRst")
  objsNeeded <- c(ecoFile, 
                  "firePerimeters", 
                  #TODO: fire perimeters is an optional object but it should be passed
                        # because the function call is written incorrectly (isn't unqique to 2001 2011 sims)
                  # and I believe this is triggering some kind of bug somewhere
                  "rasterToMatchLarge", "rasterToMatch",
                  "studyAreaLarge", "studyArea",
                  "species", "speciesTable", "sppEquiv")
  objsNeeded <- intersect(ls(sim), objsNeeded)
  objsNeeded <- mget(objsNeeded, envir = envir(sim))
  
  cds <- lapply(neededYears, function(ny, objs = objsNeeded) {
    messageColoured(colour = "yellow", "Running Biomass_borealDataPrep for year ", ny)
    messageColoured(colour = "yellow", "  inside fireSense_dataPrepFit to estimate cohortData", ny)
    
    parms <- list()
    ## if needModule is vectorized - we will have to rethink
    for (nm in neededModule) {
      parms[[nm]] <- P(sim, module = nm)
      parms[[nm]][["dataYear"]] <- ny
      parms[[nm]][["exportModels"]] <- "none"
      parms[[nm]] <- parms[[nm]][!names(parms[[nm]]) %in% SpaDES.core:::paramsDontCacheOn]
    }
    
    if (".globals" %in% names(params(sim))) {
      parms[".globals"] <- params(sim)[".globals"]
    }

    outNY <- Cache(do.call(SpaDES.core::simInitAndSpades, list(paths = pathsLocal,
                                                             params = parms,
                                                             times = list(start = ny, end = ny),
                                                             modules = neededModule,
                                                             objects = objs)),
                 .functionName = "simInitAndSpades")
    cohDatObj <- paste0(cohDat, ny)
    pixGrpMap <- paste0(pixGM, ny)
    saObj <- paste0(saMap, ny)
    rstLCCobj <- paste0(rstLCC, ny)
    outNY[[cohDatObj]] <- outNY[[cohDat]]
    outNY[[pixGrpMap]] <- outNY[[pixGM]]
    outNY[[saObj]] <- outNY[[saMap]]
    outNY[[rstLCCobj]] <- outNY[[rstLCC]]
    mget(c(cohDatObj, pixGrpMap, saObj, rstLCCobj), envir = envir(outNY))
  })
  lapply(cds, function(cd) list2env(cd, envir = envir(sim)))
  sim
}

.inputObjects <- function(sim) {
  if (!suppliedElsewhere("studyArea", sim)) {
    stop("Please supply study area - this object is key")
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
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")
  
  if (!suppliedElsewhere("rasterToMatch", sim)) {
    sim$rasterToMatch <- Cache(LandR::prepInputsLCC(year = 2005, ## TODO: use 2010
                                                    destinationPath = dPath,
                                                    studyArea = sim$studyArea,
                                                    useCache = TRUE))
  }
  
  if (!suppliedElsewhere("climateVariablesForFire", sim)) {
    sim$climateVariablesForFire <- list("spread" = "MDC", 
                                        "ignition" = "MDC")
  }
  
  if (!all(suppliedElsewhere("cohortData2011", sim),
           suppliedElsewhere("pixelGroupMap2011", sim),
           suppliedElsewhere("pixelGroupMap2001", sim),
           suppliedElsewhere("cohortData2001", sim),
           suppliedElsewhere("rstLCC2001", sim),
           suppliedElsewhere("rstLCC2011", sim))) {
    # This runs simInitAndSpades if needed
    sim <- runBorealDP_forCohortData(sim)
  }
  
  if (!P(sim)$useRasterizedFireForSpread) {
    if (!suppliedElsewhere("firePolys", sim) | !suppliedElsewhere("firePolysForAge", sim)) {
      # don't want to needlessly postProcess the same firePolys objects
      
      saNotLatLong <- if (isTRUE(sf::st_is_longlat(sim$studyArea))) {
        terra::project(sim$studyArea, terra::crs(sim$rasterToMatch))
      } else {
        sim$studyArea
      }
      
      fireYears <- c(min(P(sim)$fireYears - P(sim)$cutoffForYoungAge):max(P(sim)$fireYears))
      #check why this isn't resulting in identical crs between firePolys, studyArea
      allFirePolys <- Cache(fireSenseUtils::getFirePolygons, fun = "st_read",
                            years = fireYears, useInnerCache = TRUE, destinationPath = dPath,
                            cropTo = sim$rasterToMatch, maskTo = saNotLatLong, projectTo = sim$rasterToMatch,
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
  
  if (!suppliedElsewhere("standAgeMap2001", sim)) {
    sim$standAgeMap2001 <- Cache(prepInputsStandAgeMap,
                                 rasterToMatch = sim$rasterToMatch,
                                 studyArea = sim$studyArea,
                                 destinationPath = dPath,
                                 filename2 = "standAgeMap2001.tif",
                                 startTime = 2001,
                                 userTags = c(cacheTags, 'prepInputsStandAgeMap2001'))
  }
  
  if (!suppliedElsewhere("standAgeMap2011", sim)) {
    sim$standAgeMap2011 <- Cache(prepInputsStandAgeMap,
                                 rasterToMatch = sim$rasterToMatch,
                                 studyArea = sim$studyArea,
                                 destinationPath = dPath,
                                 filename2 = 'standAgeMap2011.tif',
                                 startTime = 2011,
                                 userTags = c(cacheTags, 'prepInputsStandAgeMap2011'))
  }
  
  if (!suppliedElsewhere("ignitionFirePoints", sim)) {
    ignitionFirePoints <- Cache(
      getFirePoints_NFDB_V2,
      studyArea = sim$studyArea,
      years = P(sim)$fireYears,
      NFDB_pointPath = dPath,
      userTags = c("ignitionFirePoints", P(sim)$.studyAreaName),
      plot = !is.na(P(sim)$.plotInitialTime)
    ) #default redownload means it will update annually - I think this is fine?
    sim$ignitionFirePoints <- ignitionFirePoints[ignitionFirePoints$CAUSE == "L",]
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
                                        method = "near",  ## make sure this is near or ngb; bilinear is wrong
                                        filename2 = paste0("wildfire_", P(sim)$.studyAreaName, ".tif"),
                                        userTags = c("historicalFireRaster", P(sim)$.studyAreaName))
    }
  }
  if (!suppliedElsewhere("nonForestedLCCGroups", sim)) {
    
    #TODO: consider moving this to init - and checking if unsupplied
    sim$nonForestedLCCGroups <- list(
      "nf_highFlam" = c(50, 100), #shrub, herbaceous, non-treed wetland
      "nf_lowFlam" = c(40, 81)) #bryoids + non-treed wetland. 
  }
  
  if (!suppliedElsewhere("rstLCC2011", sim) & !suppliedElsewhere("rstLCC2001", sim)) {
    #these are optional improvements to lcc that are not necessary for BBDP
    sim$rstLCC2011 <- Cache(LandR::prepInputs_NTEMS_Nonforest, rstLCC = sim$rstLCC2011, 
                            destinationPath = dPath, 
                            userTags = c("rstLCC2011", P(sim)$.studyAreaName))
    sim$rstLCC2001 <- Cache(LandR::prepInputs_NTEMS_Nonforest, rstLCC = sim$rstLCC2001, 
                            destinationPath = dPath,
                            userTags = c("rstLCC2001", P(sim)$.studyAreaName))
  }
  
  if (!suppliedElsewhere("flammableRTM2001", sim) & !suppliedElsewhere("flammableRTM2011", sim)) {
    
    if (!isInt(sim$rstLCC2001)) sim$rstLCC2001 <- asInt(sim$rstLCC2001)
    if (!isInt(sim$rstLCC2011)) sim$rstLCC2011 <- asInt(sim$rstLCC2011)
    
    sim$flammableRTM2001 <- defineFlammable(sim$rstLCC2001,
                                            nonFlammClasses = P(sim)$nonflammableLCC,
                                            to = sim$rasterToMatch)
    sim$flammableRTM2011 <- defineFlammable(sim$rstLCC2011,
                                            nonFlammClasses = P(sim)$nonflammableLCC,
                                            to = sim$rasterToMatch)
  }
  
  return(invisible(sim))
}

