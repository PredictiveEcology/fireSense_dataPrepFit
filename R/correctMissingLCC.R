fixMissingLCC <- function(landcoverDT, pixelGroupMap, missingLCC) {
  
  landcoverDT[, rowSums := rowSums(.SD), .SD = setdiff(names(landcoverDT), "pixelID")]
  problemPix <- forestPix[is.na(pixelGroupMap[forestPix])]
  set(landcoverDT, NULL, 'rowSums', NULL)
  
  ## The non-forests aren't the same between years, due to cohortData being different
  landcoverDT[pixelID %in% problemPix, eval(missingLCC) := 1]
  
}