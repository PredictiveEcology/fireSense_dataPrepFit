correctMissingLCC <- function(landcoverDT, pixelGroupMap, missingLCC) {
  
  landcoverDT[, rowSums := rowSums(.SD), .SD = setdiff(names(landcoverDT), "pixelID")]
  forestPix <- landcoverDT[rowSums == 0,]$pixelID
  problemPix <- forestPix[is.na(pixelGroupMap[forestPix])]
  ## The non-forests aren't the same between years, due to cohortData being different
  landcoverDT[pixelID %in% problemPix & rowSums == 0, eval(missingLCC) := 1]
  
  set(landcoverDT, NULL, 'rowSums', NULL)
}
