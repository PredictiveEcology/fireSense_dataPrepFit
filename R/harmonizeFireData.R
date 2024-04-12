## TODO: move to fireSenseUtils to get documentation and tests!!
harmonizeFireData <- function(firePolys, flammableRTM, spreadFirePoints,
                              areaMultiplier, minSize, pointsIDcolumn = "FIRE_ID") {
  ## safety to ensure missing years actually removed
  spreadFirePoints[is.na(names(spreadFirePoints))] <- NULL

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

  if (!any(sapply(fireBufferedListDT, is.data.table))) {
    fireBufferedListDT <- lapply(fireBufferedListDT, as.data.table)
  }

  harmonized <- harmonizeBufferAndPoints(
    cent = spreadFirePoints,
    buff = fireBufferedListDT,
    ras = flammableRTM,
    idCol = pointsIDcolumn
  )

  ## ensure mismatched (e.g. points w/ no polys) and now missing years actually removed.
  ## TODO: move this check code into the corresponding fireSenseUtils functions
  emptyYearsPoints <- which(vapply(harmonized, is.null, logical(1))) |> names()
  emptyYearsPolys <- which(vapply(fireBufferedListDT, function(x) nrow(x) == 0, logical(1))) |> names()
  stopifnot(emptyYearsPoints == emptyYearsPolys)
  emptyYears <- unique(c(emptyYearsPoints, emptyYearsPolys))
  if (length(emptyYears > 0)) {
    harmonized[[emptyYears]] <- NULL
    fireBufferedListDT[[emptyYears]] <- NULL
  }

  harmonized <- Map(f = cleanUpSpreadFirePoints,
                    firePoints = harmonized,
                    bufferDT = fireBufferedListDT,
                    MoreArgs = list(flammableRTM = flammableRTM)) |>
    purrr::transpose()

  return(list(fireBufferedListDT = harmonized$FireBuffered,
              firePolys = firePolys, #should be returned because some years may have been converted to NULL
              spreadFirePoints = harmonized$SpatialPoints))
}
