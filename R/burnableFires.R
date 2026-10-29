## The spread fit's landscape is the burnable pixels only (`landcoverDTs`; fireSenseUtils' objective gives the
## others no spread), so each fire's target size is its burned area on those pixels, not its polygon area. On
## ELF 10.3.1_2 the median fit fire was 6% burnable (2012: a 7,516 ha polygon with 397 ha burnable).

#' Each fire year's burnable pixels
#'
#' A fire year's pixels that can burn are those of its vegetation data year's `landcoverDTs` (flammable pixels).
#'
#' @param landcoverDTs list named `year<dataYear>` of `data.table`s with `pixelID`.
#' @param allYears list named by data year of the fire year names (`year<year>`) that use it, as `mod$allYears`.
#' @return list named by fire year (`year<year>`) of integer `pixelID`s.
burnablePixelsByFireYear <- function(landcoverDTs, allYears) {
  byDataYear <- Map(fy = allYears, dy = names(allYears), function(fy, dy) {
    rep(list(landcoverDTs[[paste0(fireSenseUtils::yearTxt, dy)]]$pixelID), length(fy)) |> setNames(fy)
  })
  unlist(unname(byDataYear), recursive = FALSE)
}

#' The fire ID column of fire points or polygons
#'
#' @param fires `sf` or `SpatVector` of fires.
#' @return character, the column whose name has `ID`, other than `GID`.
fireIDColumn <- function(fires) {
  grep("ID", names(fires), value = TRUE) |> setdiff("GID")
}

#' The burnable area of each spread-fit fire, and the fires with enough of it
#'
#' Adds to each fire's start point its burned area on burnable pixels (column `fireSenseUtils::burnableHaTxt`),
#' which `fireSense_spreadFit` uses as the fire's target size, and leaves the fires with less than `minHa` of it
#' out of `fireBufferedListDT`, the fires the spread model is fitted to. Their start points stay: the escape
#' data place a polygon's ignition at its start point (`fireSenseUtils::polygonIgnitions()`).
#'
#' @param fireBufferedListDT list named `year<year>` of `data.table`s of `pixelID`, `buffer` (1 burned) and `ids`.
#' @param spreadFirePoints list named `year<year>` of the fires' start points (`sf`, `SpatVector` or `data.frame`).
#' @param burnablePixels from [burnablePixelsByFireYear()].
#' @param pixSizeHa area of one pixel (ha).
#' @param minHa fires with less burnable area (ha) leave the spread fit.
#' @param idCol the fire ID column of `spreadFirePoints`, matching `ids`.
#' @return list of `fireBufferedListDT` (years with no fire left are dropped) and `spreadFirePoints`.
spreadFitBurnableFires <- function(fireBufferedListDT, spreadFirePoints, burnablePixels, pixSizeHa, minHa,
                                   idCol) {
  burnableHa <- Map(fb = fireBufferedListDT, fy = names(fireBufferedListDT), function(fb, fy) {
    burned <- fb$buffer == 1L
    nBurnable <- tapply(burned & fb$pixelID %in% burnablePixels[[fy]], fb$ids, sum)
    nBurnable * pixSizeHa
  })
  pts <- Map(p = spreadFirePoints, fy = names(spreadFirePoints), function(p, fy) {
    p[[fireSenseUtils::burnableHaTxt]] <- as.vector(burnableHa[[fy]][as.character(as.data.frame(p)[[idCol]])])
    p
  })
  keep <- Map(fb = fireBufferedListDT, ha = burnableHa, function(fb, ha) {
    fb[as.character(fb$ids) %in% names(ha)[ha >= minHa], ]
  })
  nFires <- length(unlist(burnableHa))
  nDropped <- nFires - sum(vapply(keep, function(fb) length(unique(fb$ids)), integer(1)))
  if (nDropped > 0)
    message("fireSense_dataPrepFit: ", nDropped, " of ", nFires, " fires have less than ", minHa,
            " ha on burnable pixels and are left out of the spread fit")
  list(fireBufferedListDT = keep[vapply(keep, NROW, integer(1)) > 0L], spreadFirePoints = pts)
}

#' The smallest burnable area (ha) of a spread-fit fire
#'
#' @param Par the module's parameters.
#' @return `minBurnableFireSizeHa`, or `escapeSizeHa` when it is `NA`.
minBurnableFireSizeHa <- function(Par) {
  if (is.na(Par$minBurnableFireSizeHa)) Par$escapeSizeHa else Par$minBurnableFireSizeHa
}

#' The data-year groups that have a fire year with fires
#'
#' @param allYears list named by data year of fire year names, as `mod$allYears`.
#' @param yearsWithFire fire year names (`year<year>`) that have fires.
#' @return `allYears` without the groups that have none.
yearGroupsWithFire <- function(allYears, yearsWithFire) {
  allYears[vapply(allYears, function(ay) any(yearsWithFire %in% ay), logical(1))]
}
