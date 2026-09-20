#' Assign flammable pixels with no fuel to the `missingLCC` group
#'
#' These are pixels whose land cover is forest (no non-forest group) but that have no `pixelGroup`.
#'
#' @param landcoverDT data.table of `pixelID` and one 0/1 column per non-forest group, from
#'   `fireSenseUtils::makeLandcoverDT`. Modified by reference.
#' @param pixelGroupMap `SpatRaster` of pixel groups for the same year.
#' @param missingLCC character, the `landcoverDT` column to set to 1 for those pixels.
#' @return `landcoverDT`, invisibly.
correctMissingLCC <- function(landcoverDT, pixelGroupMap, missingLCC) {
  landcoverDT[, rowSums := rowSums(.SD), .SDcols = setdiff(names(landcoverDT), "pixelID")]
  forestPix <- landcoverDT[rowSums == 0, ]$pixelID
  problemPix <- forestPix[is.na(pixelGroupMap[forestPix])]
  ## The non-forests aren't the same between years, due to cohortData being different
  landcoverDT[pixelID %in% problemPix & rowSums == 0, eval(missingLCC) := 1]

  set(landcoverDT, NULL, "rowSums", NULL)
}
