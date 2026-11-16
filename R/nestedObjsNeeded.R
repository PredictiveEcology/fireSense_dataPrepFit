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

#' The `.seed` parameter of the nested `Biomass_borealDataPrep` runs
#'
#' SpaDES.core calls `set.seed()` for an event named in `.seed`; the draws of the nested module
#' (age and biomass model subsets, replacement land cover) all happen in its `init` event.
#'
#' @param seed numeric or `NA`; the `vegInitSeed` parameter.
#' @param current the `.seed` the nested module already has (kept when `seed` is `NA`).
#' @return a named list of seeds per event.
#' @noRd
nestedVegSeed <- function(seed, current = NULL) {
  if (is.null(seed) || is.na(seed)) return(current)
  utils::modifyList(as.list(current), list(init = seed))
}
