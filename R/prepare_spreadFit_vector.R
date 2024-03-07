prepare_spreadPeriod_vector <- function(flammableRTM, firePolys, minBufferSize, 
                                        period, areaMultiplier, minBufferSize, nCores, 
                                        studyAreaName, spreadFirePoints) {
  
  firePolys <- firePolys[period]
  spreadFirePoints <- spreadFirePoints[period]
  
  fireBufferedListDT <- Cache(bufferToArea,
                              poly = firePolys,
                              polyName = names(firePolys),
                              rasterToMatch = flammableRTM, 
                              verb = TRUE,
                              areaMultiplier = areaMultiplier,
                              field = "FIRE_ID",
                              cores = nCores,
                              minSize = minBufferSize,
                              userTags = c("bufferToArea"),
                              omitArgs = "cores")
  
  if (!any(sapply(fireBufferedListDT, is.data.table)))
    fireBufferedListDT <- lapply(fireBufferedListDT, as.data.table)
  
  if (!any(sapply(listDT, is.data.table)))
    listDT <- lapply(listDT, as.data.table)
  ## drop fire years from these lists that don't have any buffer points pre-harmonization
  omitYears <- sapply(listDT, function(x) nrow(x) == 0)
  listDT[omitYears] <- NULL
  firePolys[omitYears] <- NULL
  spreadFirePoints[omitYears] <- NULL
  
  ## Post buffering, new issues --> must make sure points and buffers match
  pointsIDColumn <- "FIRE_ID"
  
  spreadFirePoints <- Cache(harmonizeBufferAndPoints,
                            cent = spreadFirePoints,
                            buff = fireBufferedListDT,
                            ras = flammableRTM,
                            idCol = pointsIDColumn,
                            userTags = c("harmonizeBufferAndPoints", studyAreaName))
  
  ## drop fire years from these lists that don't have any buffer points post-harmonization
  omitYears <- sapply(spreadFirePoints, is.null)
  fireBufferedListDT[omitYears] <- NULL
  firePolys[omitYears] <- NULL
  spreadFirePoints[omitYears] <- NULL
  
  ## Also 2 other problems:
  ## 1. Big fire, but ignition is in non-flammable pixels e.g., lake -- bad;
  ##    solution -- pick nearest pixel in burned polygon
  ## 2. Small fire, ignition in non-flammable pixel, but NO pixel in burned polygon
  ##    is actually flammable -- remove this from data
  #FIRE_ID is hardcoded here but since we enforce its presence, it should be permissible..
  out22 <- Map(f = cleanUpSpreadFirePoints,
               firePoints = spreadFirePoints,
               bufferDT = fireBufferedListDT,
               MoreArgs = list(flammableRTM = flammableRTM)) |>
    Cache(userTags = studyAreaName, .functionName = "cleanUpSpreadFirePoints")
  
  
  out22 <- purrr::transpose(out22)
  
  
  return(list(spreadFirePoints = out22$SpatialPoints,
              fireBufferedListDT = out22$FireBuffered))
}
