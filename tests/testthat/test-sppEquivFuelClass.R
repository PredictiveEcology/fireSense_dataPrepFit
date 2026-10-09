## Estimating fuel classes dropped the species without cohorts from sppEquiv (a right join to the
## species that have cohorts) while sppNameVector and sppColorVect kept them, and Biomass_core's
## sppHarmonize later stopped on the length mismatch.

test_that("sppEquiv keeps every species when some have no cohorts; the others get their fuel class", {
  sppEquiv <- data.table::data.table(LandR = c("Abie_bal", "Pice_gla", "Pice_mar", "Pinu_ban"),
                                     EN_generic_short = c("Fir", "Spruce", "Spruce", "Pine"),
                                     FuelClass = c("old1", "old2", "old3", "old4"))
  sppNameVector <- sppEquiv$LandR
  assigned <- data.table::data.table(species = c("Pice_mar", "Pinu_ban"),
                                     assignedFuelClass = c("SpruceFuel", "PineFuel"))
  data.table::setnames(assigned, c("LandR", "FuelClass"))
  before <- data.table::copy(sppEquiv)

  out <- setSppEquivFuelClass(sppEquiv, assigned, "LandR", "FuelClass")

  expect_identical(NROW(out), NROW(before))
  expect_identical(out$LandR, sppNameVector)
  expect_identical(out$FuelClass, c(NA, NA, "SpruceFuel", "PineFuel"))
  expect_identical(colnames(out), colnames(before))
  expect_identical(out$EN_generic_short, before$EN_generic_short)
  expect_identical(sppEquiv$FuelClass, before$FuelClass) # input untouched
})
