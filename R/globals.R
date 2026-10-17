## Names that R CMD check cannot see: `Par` is defined in the module's environment by SpaDES.core, the
## rest are data.table columns used in `[` calls, or objects `list2env()` puts in the function environment.
## (Called without `utils::` because tests that read module code take the head of every top-level call as a name.)
globalVariables(c(".", "B", "FuelClass", "Par", "age", "assignedFuelClass", "auc", "burned", "cell",
                      "chosen", "i.FuelClass", "lcc", "mod", "newLCC", "nonForestedLCCGroups", "pixelID",
                      "pixelIndex", "rho", "species", "speciesCode", "sppEquiv", "totalBiomass", "var"))
