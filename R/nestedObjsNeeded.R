#' Objects passed from the parent sim to the nested `Biomass_borealDataPrep` runs
#'
#' The ecoregion layers (`ecoregionLayer`, `ecoregionRst`) are deliberately not passed: the fit's
#' vegetation reconstruction uses `Biomass_borealDataPrep`'s own default ecoregions, not the local
#' regions (e.g. BEC zones from `localEcozones`) that other modules put in the parent sim.
#'
#' @param available character; names of the objects in the parent sim.
#' @return the subset of `available` to pass on.
#' @noRd
nestedObjsNeeded <- function(available) {
  intersect(available,
            c("firePerimeters",
              "rasterToMatch", "studyArea",
              "rstLCCs",
              "standAgeMaps",
              "studyArea_biomassParam", "rasterToMatch_biomassParam", #needed by BBDP
              "species", "speciesTable", "sppEquiv"))
}
