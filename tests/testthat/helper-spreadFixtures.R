## Fixtures shared by several test files: test files do not see each other's functions, helpers are sourced into all.

## the whole-landscape fixture (test-wholeLandscapeSpreadCovariates.R, test-spreadClimateHinge.R)
wholeFixture <- function() {
  r <- terra::rast(nrows = 4, ncols = 5)
  years <- 2001:2004
  clim <- list(CMD = terra::rast(lapply(years, function(y) { x <- r; terra::values(x) <- seq_len(20) + y; x })))
  names(clim$CMD) <- paste0("year", years)
  tsd <- function(v) { x <- r; terra::values(x) <- v; x }
  tsds <- list(year2000 = tsd(c(rep(100, 19), 3)), year2002 = tsd(rep(100, 20)))
  allYears <- list(`2000` = c("year2001", "year2002"), `2002` = c("year2003", "year2004"))
  burnable <- seq_len(16)                                        # pixels 17:20 cannot burn
  burnablePixels <- setNames(rep(list(burnable), 4), unlist(allYears))
  buffers <- list(year2001 = data.table::data.table(pixelID = c(2L, 3L, 5L)),
                  year2002 = data.table::data.table(pixelID = integer(0)),    # a year without fire
                  year2003 = data.table::data.table(pixelID = c(3L, 16L, 17L)),
                  year2004 = data.table::data.table(pixelID = c(7L, 8L)))
  list(clim = clim, tsds = tsds, allYears = allYears, burnablePixels = burnablePixels, buffers = buffers,
       fires = list(`2003` = 9L))
}

## a ledger with this ELF (4.2.2, CMD) and a neighbour (4.2.1, CMDsm) (test-spreadClimateOwnELF.R, test-spreadClimateHinge.R)
ownAndNeighbourLedger <- function() {
  own <- oneELFLedger()                     # 4.2.2, fitted with CMD
  nb <- oneELFLedger()
  nb[[fireSenseUtils::polygonIDTxt]] <- "4.2.1"
  nbPars <- own$params[[1]]
  names(nbPars)[names(nbPars) == "CMD"] <- "CMDsm"
  ledger <- rbind(nb, own)
  ledger$params <- list(nbPars, own$params[[1]])
  ledger
}
