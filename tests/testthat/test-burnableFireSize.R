## The spread fit's landscape has burnable pixels only (prepare_SpreadFit() keeps the pixels of `landcoverDTs`,
## and fireSenseUtils' objective gives the others no spread), but each fire's target size was its whole polygon
## (POLY_HA). ELF 10.3.1_2, 2012: a 7,516 ha polygon with 397 ha on burnable pixels, a size no simulated fire
## could reach.

## Two fires of 2001 on 6.25 ha pixels. Fire 11 burned pixels 1-5, of which only 1 and 2 are burnable; fire 12
## burned 9-11, all burnable. Pixels 6-8 and 12 are their unburned buffers.
toyBuffers <- function() {
  list(year2001 = data.table::data.table(pixelID = 1:12,
                                         buffer = c(1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 1L, 1L, 0L),
                                         ids = c(rep(11L, 8), rep(12L, 4))))
}
toyBurnable <- function() list(year2001 = c(1L, 2L, 6:12))
toyStarts <- function() list(year2001 = data.frame(FIRE_ID = c(11L, 12L), POLY_HA = c(500, 20)))

test_that("a fire mostly on non-burnable pixels gets its burnable area as its spread-fit target", {
  out <- spreadFitBurnableFires(toyBuffers(), toyStarts(), toyBurnable(), pixSizeHa = 6.25, minHa = 0,
                                idCol = "FIRE_ID")
  pts <- out$spreadFirePoints$year2001
  expect_equal(pts[[fireSenseUtils::burnableHaTxt]], c(2, 3) * 6.25)
  expect_equal(pts$POLY_HA, c(500, 20))        # the polygon area is kept
  expect_identical(out$fireBufferedListDT, toyBuffers())
})

test_that("fires with less burnable area than minHa leave the spread fit's fires, with a message, but keep their start", {
  expect_message(
    out <- spreadFitBurnableFires(toyBuffers(), toyStarts(), toyBurnable(), pixSizeHa = 6.25, minHa = 15,
                                  idCol = "FIRE_ID"),
    "1 of 2 fires")
  fb <- out$fireBufferedListDT$year2001
  expect_identical(unique(fb$ids), 12L)
  ## the start points stay: the escape data place a polygon's ignition at its start point
  expect_identical(out$spreadFirePoints$year2001$FIRE_ID, c(11L, 12L))
  ## a year left with no fire is dropped
  out <- suppressMessages(spreadFitBurnableFires(toyBuffers(), toyStarts(), toyBurnable(), pixSizeHa = 6.25,
                                                 minHa = 100, idCol = "FIRE_ID"))
  expect_length(out$fireBufferedListDT, 0L)
})

test_that("each fire year's burnable pixels are its vegetation data year's landcoverDT pixels", {
  allYears <- list("1985" = c("year1986", "year1987"), "1990" = "year1991")
  lcc <- list(year1985 = data.table::data.table(pixelID = 1:3), year1990 = data.table::data.table(pixelID = 4:5))
  expect_identical(burnablePixelsByFireYear(lcc, allYears),
                   list(year1986 = 1:3, year1987 = 1:3, year1991 = 4:5))
})

test_that("minBurnableFireSizeHa defaults to escapeSizeHa, the smallest fire the spread objective scores", {
  expect_identical(minBurnableFireSizeHa(list(minBurnableFireSizeHa = NA_real_, escapeSizeHa = 50)), 50)
  expect_identical(minBurnableFireSizeHa(list(minBurnableFireSizeHa = 6.25, escapeSizeHa = 50)), 6.25)
})
