## A hindcast (fireSense_dataPrepPredict, vegetationYear = 1985) overwrites sim$rstLCC_RTM and
## sim$flammableRTM with 1985 before prepare_SpreadFit() runs. Whether treed wetland is a spread
## covariate (minCovariateProp) and the fuel-class land-cover shares must still come from the most
## recent data year, or the spread formula gets a term the fit does not have.

lccRaster <- function(vals) terra::rast(nrows = 2, ncols = 5, vals = vals)

test_that("the land-cover shares are those of the latest year in the lists, not of sim$rstLCC_RTM", {
  lcc2010 <- lccRaster(c(rep(1L, 8), 2L, 2L))   # class 2 is 20%
  lcc2020 <- lccRaster(c(rep(1L, 9), 2L))       # class 2 is 10%
  flam <- lccRaster(rep(1L, 10))
  sim <- list(rstLCCs = list(`2010` = lcc2010, `2020` = lcc2020),
              flammableRTMs = list(`2010` = flam, `2020` = flam))
  expected <- c(`1` = 0.9, `2` = 0.1)
  expect_equal(latestLccFlammableShare(sim), expected)

  ## the shared objects a hindcast overwrites with another year change nothing
  sim$rstLCC_RTM <- lcc2010
  sim$flammableRTM <- flam
  expect_equal(latestLccFlammableShare(sim), expected)
})

test_that("neither dataPrepBuild nor prepare_SpreadFit reads the shared year's rasters for shares", {
  for (f in c("dataPrepBuild", "prepare_SpreadFit")) {
    txt <- paste(deparse(body(get(f))), collapse = " ")
    expect_true(grepl("latestLccFlammableShare(sim)", txt, fixed = TRUE), info = f)
    expect_false(grepl("lccFlammableShare(sim$rstLCC_RTM", txt, fixed = TRUE), info = f)
  }
})
