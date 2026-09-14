## Each fire-year group must be joined to the vegetation of its own data year only. The join
## used to keep every data year at or before the group's first fire year, so the 2010 group got
## the 2000 and 2010 vegetation and the 2020 group got all three. Every pixel then appeared once
## per data year, and for ELF 14.3 the join overflowed:
##   Error in vecseq(f__, len__, limit): Join results in 5708751 rows; more than 5704496
test_that("fire buffers join only to their own data year's vegetation", {
  dataYears <- c(2000L, 2010L, 2020L)
  fireYears <- 2002:2025
  pixels <- 1:10

  vegData <- data.table::rbindlist(lapply(dataYears, function(dy) {
    data.table::data.table(pixelID = pixels, fuel = dy, year = as.character(dy))
  }))
  ## every pixel is in the buffer of every fire year: enough repeats to overflow the old join
  fireBufferedListDT <- lapply(fireYears, function(fy) {
    data.table::data.table(pixelID = pixels, buffer = rep(0:1, length.out = length(pixels)))
  })
  names(fireBufferedListDT) <- paste0(fireSenseUtils::yearTxt, fireYears)

  allYears <- yearGroups(dataYears, fireYears, FALSE)
  allYears <- Map(y = allYears, function(y) paste0(fireSenseUtils::yearTxt, y))

  out <- joinFireBuffersToVeg(fireBufferedListDT, data.table::copy(vegData), allYears)

  expect_identical(nrow(out), length(fireYears) * length(pixels))
  expect_false(anyDuplicated(out, by = c("fireYear", "pixelID")) > 0)
  fireYearNum <- as.integer(gsub(fireSenseUtils::yearTxt, "", out$fireYear))
  expectedDataYear <- max.col(outer(fireYearNum, dataYears, ">=") * 1, ties.method = "last")
  expect_identical(out$year, dataYears[expectedDataYear])
  expect_identical(out$fuel, out$year)
})
