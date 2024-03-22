harmonizeFireData <- function(firePolys, flammableRTM, spreadFirePoints, 
                              areaMultiplier, minSize, pointsIDcolumn = "FIRE_ID") {
  fireYears <- names(firePolys)
  fireBufferedListDT <- Cache(bufferToArea, 
                              poly = firePolys,
                              polyName = fireYears,
                              rasterToMatch = flammableRTM, 
                              verb = TRUE,
                              areaMultiplier = areaMultiplier,
                              field = pointsIDcolumn,
                              minSize = minSize,
                              cacheTags = c("bufferToArea", fireYears[1]))
  
  if (!any(sapply(fireBufferedListDT, is.data.table)))
    fireBufferedListDT <- lapply(fireBufferedListDT, as.data.table)
  
  harmonized <- harmonizeBufferAndPoints(
    cent = spreadFirePoints,
    buff = fireBufferedListDT,
    ras = flammableRTM,
    idCol = pointsIDcolumn
  ) %>%
    Map(f = cleanUpSpreadFirePoints, firePoints = ., bufferDT = fireBufferedListDT, 
        MoreArgs = list(flammableRTM = flammableRTM)) %>%
    purrr::transpose(.)
  
  return(list(fireBufferedListDT = harmonized$FireBuffered, 
              firePolys = firePolys, #should be returned because some years may have been converted to NULL
              spreadFirePoints = harmonized$SpatialPoints))
}
