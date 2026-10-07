## The nested Biomass_borealDataPrep runs were given the parent's `ecoregionLayer` when one existed
## (e.g. built by localEcozones), so the fit used local regions depending on module load order.

test_that("the parent's ecoregion layers are not passed to the nested runs", {
  available <- c("ecoregionLayer", "ecoregionRst", "studyArea", "rasterToMatch", "sppEquiv", "other")
  out <- nestedObjsNeeded(available)
  expect_false(any(c("ecoregionLayer", "ecoregionRst") %in% out))
  expect_setequal(out, c("studyArea", "rasterToMatch", "sppEquiv"))
})
