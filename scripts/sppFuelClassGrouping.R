library(reproducible)
library(terra)
library(sf)
library(data.table)
library(magrittr)
library(raster)
library(LandR)

#hard-codey stuff

cacheDir <- "C:/Ian/testing/cache"
options("reproducible.cachePath" = cacheDir)
knnPath <- 'C:/Ian/data/kNN'
inputPath <-  "C:/Ian/testing/"
Ecozones <- terra::vect("C:/Ian/Data/Ecoregions/ecozone_shp/Ecozones/ecozones.shp")
#the zones are only used to get the WBI study area
#the EcoProvinces do not have a legend, but one exists in the downloadable FDGB -
# https://agriculture.canada.ca/atlas/data_donnees/nationalEcologicalFramework/data_donnees/fgdb/ep/nef_ca_ter_ecoprovince_v2_2.gdb.zip
legend <- fread("C:/Ian/Data/Ecoregions/Ecoprovince_info.txt")

###

minCov <- 10
getOption("reproducible.useTerra") #must be true
options("reproducible.rasterRead" = "terra::rast")

####data prep #####
sppEquiv <- LandR::sppEquivalencies_CA
sppEquiv <- sppEquiv[KNN %in% c("Pice_Mar", "Pice_Gla", "Pice_Eng", "Abie_Bal", "Abie_Las",
                                  "Pinu_Con", "Pinu_Ban", "Popu_Tre", "Pseu_Men", "Popu_Bal")]
#TODO: KNN is not finding Pseu_Men_Gla and Pseu_Men_Men
#also maybe use .studyAreaName


Canada <- prepInputs(url = paste0("https://www12.statcan.gc.ca/census-recensement/2011/",
                                  "geo/bound-limit/files-fichiers/2016/lpr_000b16a_e.zip"),
                     destinationPath = inputPath,
                     fun = "terra::vect")

SA <- Canada[Canada$PRUID %in% c(46,47,48,49,59,60,61),]
Ecozones <- postProcess(Ecozones, studyArea = SA)
Ecozones$ECOZONE
EcozoneSA <- Ecozones[Ecozones$ECOZONE %in% c(9, 6, 12, 14, 4, 5, 11),] %>%
  project(., SA)
SA <- terra::aggregate(SA) %>%
  intersect(., EcozoneSA)
SA <- terra::aggregate(SA)

Ecoprovinces <- prepInputs(url = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/province/ecoprovince_shp.zip",
                           destinationPath = inputPath,
                           useCache = "overwrite",
                           fun = "terra::vect")
#using SA as studyArea was losing the province attribute data - bug?
Ecoprovinces <- project(Ecoprovinces, SA)
Ecoprovinces <- intersect(Ecoprovinces, SA)
Ecoprovinces$ECOPROVINC_num <- as.numeric(Ecoprovinces$ECOPROVINC) #for rasterize


#Add a legend
legend[, ECOPROVINC := as.character(as.factor(round(ECOPROVINCE_ID, digits = 1)))]
legend[, c("ECOPROVINCE_NAME_FR", "SHAPE_Length", "SHAPE_Area", "OBJECTID") := NULL] #sorry Celine no francais
legend <- unique(legend) #remove duplicate entries
temp <- as.data.table(Ecoprovinces)[, .(PROVINCE_, ECOPROVINC, ECOPROVINC_num)] #the PROVINCE_ field is not in the original table
legend <- legend[temp, on = c("ECOPROVINC")]
legend <- unique(legend)
rm(temp)

options("reproducible.rasterFun" = "raster::raster") #ugh I don't have terra LandR yet
rasterToMatch <- Cache(prepInputsStandAgeMap, studyArea = Ecoprovinces, destinationPath = inputPath,
                       userTags = c("standAge", "chickenscratch"))
Ecoprovinces <- st_as_sf(Ecoprovinces)
Ecoprovinces <- st_transform(Ecoprovinces, crs(rasterToMatch))
SAsf <- st_as_sf(SA)

kNN2001 <- Cache(prepSpeciesLayers_KNN, destinationPath = knnPath, studyArea = SAsf,
                                 rasterToMatch = rasterToMatch_rast,
                                 outputPath = inputPath, sppEquiv = sppEquiv,
                                 year = 2001, sppEquivCol = "LandR",
                 userTags = c("knn2001", "chickenscratch"))

kNN2011 <- Cache(prepSpeciesLayers_KNN, destinationPath = knnPath, studyArea = SAsf,
                                 rasterToMatch = rasterToMatch_rast,
                                 outputPath = inputPath, sppEquiv = sppEquiv,
                                 year = 2011, sppEquivCol = "LandR",
                 useCache = "overwrite",
                 userTags = c("kNN2011", "chickenscratch"))

fireRas2001filename <- file.path(inputPath, "fireRas2001.tif")
fireRas2011filename <- file.path(inputPath, "fireRas2011.tif")
if (!file.exists(fireRas2001filename)){
  fireUrl <- "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip"

  #do not change this chunk to  avoid z error
  ###
  firePolys <- prepInputs(url = fireUrl, fun = "terra::vect",
                             destinationPath = inputPath, useCache = TRUE,
                             userTags = c("firePoly", "chickenscratch")) %>%
    project(., SA)
  firePolys <- terra::intersect(firePolys, SA)
  firePolys <- sf::st_as_sf(firePolys)
  firePolys <- st_transform(firePolys, crs(rasterToMatch))
  firePolys$foo <- 1 #for fasterize
  firePolys <- st_cast(firePolys, "MULTIPOLYGON") #for fasterize
  ####

  #make two groups - the fires after kNN 2001 and after kNN 2011
  #this follows fireSense_dataPrepFit and could be done in the Init.
  #In fact the only difference is we are joining fire data with ALL pixels here,
  # while in fireSense_dataPrepFit we take only a fire buffer.
  firePolys2000_2010 <- firePolys[firePolys$YEAR > 2001 & firePolys$YEAR < 2011,]
  firePolys2011_Present <- firePolys[firePolys$YEAR > 2011,]
  fireRas2001 <- fasterize::fasterize(firePolys2000_2010, raster = rasterToMatch,
                                      field = "foo", fun = "sum") #collapse multiple years
  fireRas2011 <- fasterize::fasterize(firePolys2011_Present, raster = rasterToMatch,
                                      field = "foo", fun = "sum")
  writeRaster(fireRas2001, fireRas2001filename, overwrite = TRUE)
  writeRaster(fireRas2011, fireRas2011filename, overwrite = TRUE)
  rm(firePolys2000_2010)
  rm(firePolys2011_Present)
} else {
  fireRas2001 <- raster(fireRas2001filename)
  fireRas2011 <- raster(fireRas2011filename)
}
#to save time

#some ecoprovinces are dropped because they are intersection errors
ecoprovinceRas <- terra::rasterize(Ecoprovinces, rasterToMatch, field = "ECOPROVINC_num")

####build table - there are two measurements of species + fire occurence ####
makeDT <- function(kNNRas, ecoprovinceRas, fireRas, minCover) {
  burnsBySpp <- as.data.table(values(kNNRas))
  spp <- names(burnsBySpp)
  burnsBySpp[, totalCover := rowSums(.SD, na.rm = TRUE), .SDcol = spp]
  burnsBySpp[, pixelID := 1:ncell(ecoprovinceRas)]
  burnsBySpp <- burnsBySpp[!is.na(totalCover) & totalCover > minCover] #must have at least some relevant cover
  gc()
  #some covers exceed 100 - likely rounding?
  burnsBySpp[, burn := fireRas[pixelID]]
  burnsBySpp[, ecoprovince := ecoprovinceRas[pixelID]]
  return(burnsBySpp)
}

burnsBySpp2001 <- makeDT(kNNRas = kNN2001, ecoprovinceRas = ecoprovinceRas, fireRas = fireRas2001, minCover = minCov)
burnsBySpp2001[is.na(burn), burn := 0]
burnsBySpp2011 <- makeDT(kNNRas = kNN2011, ecoprovinceRas, fireRas = fireRas2011, minCover = minCov)
burnsBySpp2011[is.na(burn), burn := 0]
burnsBySpp <- rbind(burnsBySpp2001, burnsBySpp2011)

### for now
fwrite(burnsBySpp, "C:/Ian/testing/burnsBySpp.csv")



rm(burnsBySpp2001, burnsBySpp2011, fireRas2001, fireRas2011)
gc()
burnsBySpp[, pixelID := NULL]
burnsBySpp[, ECOPROVINC_num := as.character(as.factor(round(ecoprovince, digits = 1)))]

spp <- names(burnsBySpp)[!names(burnsBySpp) %in% c("totalCover", "burn", "ecoprovince")]

#write summarizing functions
propBurn <- function(x, Burn, minCover){
  n <- x > minCover
  denom <- sum(n, na.rm = TRUE)
  num <- sum(Burn[n], na.rm = TRUE)
  propBurned <- round(num/denom, digits = 5)
  return(propBurned)  #as percent
}
count <- function(x, minCover){
  n <- x > minCover
  N <- sum(n, na.rm = TRUE)
  return(N)
}
#calculate percent of each vegetation that burned in each ecoprovince. Multiple burns allowed
burnSummary <- burnsBySpp[, lapply(.SD, propBurn, Burn = burn, minCover = minCov),
                          .SDcol = spp, .(ecoprovince)]
legend$ECOPROVINC_num <- as.character(legend$ECOPROVINC_num)
legend <- legend[, .(ECOPROVINCE_NAME_EN, ECOZONE_ID,ECOPROVINC_num)]
burnSummary <- legend[burnSummary, on = c("ECOPROVINC_num" = "ECOPROVINC_num")]
burnLong <- melt(burnSummary, measure.vars = spp, value.name = "propBurn", variable.name = "spp")
setkey(burnLong, spp)
burnLong <- burnLong[!duplicated(burnLong)]

#calculate sample sizes -
Samps <- burnsBySpp[, lapply(.SD, count, minCover = minCov), .SDcol = spp, .(ecoprovince)]

Samps[, sums := sum(.SD, na.rm = TRUE), .SDcol = spp, .(ecoprovince)]
Samps <- Samps[, lapply(.SD, FUN = function(x){round(x/sums, digits = 4)} * 100), .SDcol = spp, .(ecoprovince)]
Samps <- melt(Samps, id.vars = "ecoprovince", variable.name = "spp", value.name = "PctCover")

lenend$
Samps[, ECOPROVINC_num := NULL]
Samps1 <- legend[Samps, on = c("ecoprovince")]

#TODO fix this join
burnLong <- Samps[burnLong, on = c("ECOPROVINCE_NAME_EN", "ECOZONE_ID", "ECOPROVINC_num", "spp")]

# function for classifying percent burn into clusters


clusterFun <- function(province, dt, clusterGroup, minSamp){
  dt <- dt[ECOPROVINC_num == province,]
  dt <- dt[!is.nan(propBurn)]
  dt <- dt[PctCover > minSamp]
  if (length(unique(dt$propBurn)) <= clusterGroup) {
    return(NULL)
  } else {
    kMeans <- kmeans(dt$propBurn, centers = clusterGroup)
    dt[, numGroups := clusterGroup]
    dt[, cluster := kMeans$cluster]
    dt[, clusterSS := kMeans$betweenss/kMeans$totss]
    return(dt)
  }
}

#run for 2 and 3 group clusters
kmeans2 <- rbindlist(lapply(unique(burnLong$ECOPROVINC_num), FUN = clusterFun,
                            dt = burnLong, minSamp = 2.5, clusterGroup = 2))
kmeans3 <- rbindlist(lapply(unique(burnLong$ECOPROVINC_num), FUN = clusterFun,
                            dt = burnLong, minSamp =2.5, clusterGroup = 3))
means <- rbind(kmeans2, kmeans3)
means[, numGroups := as.factor(numGroups)]
means[, cluster := as.factor(cluster)]
fwrite(means, file.path(inputPath, "kmeans_burnability_by_spp_Ecoprovince.csv"))


sppCols <- LandR::sppColors(sppEquiv, "LandR", palette ="Accent")




c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
library(ggplot2)
ggplot(data = means[numGroups == 2,], aes(y = PctCover, x = cluster, fill = spp)) +
  geom_bar(stat = "identity", position = 'stack') +
  # geom_bar(data = means[numGroups == 2], position = "dodge", stat = "identity") +
  theme_bw() +
  labs(x = "group", y = "cover", fill = "spp") +
  theme(strip.text.x = element_text(size = 8)) +
  theme(axis.text = element_text(angle = 90)) +
  facet_wrap(~ ECOPROVINCE_NAME_EN) +
  scale_fill_discrete(type = c25[1:10])

# EliotsMeans <- dcast(data = means, formula = ECOPROVINC + numGroups ~ spp, value.var = "cluster")
