## Root cause: fireSense_dataPrepFit.R:81 hard-coded `nonflammableLCC = c(0, 20, 31, 32, 33)`,
## the NTEMS non-flammable codes, missing SCANFI's combined rock/exposed code (30). Land cover
## built by `fireSenseUtils::makeFireSenseLCC()` from SCANFI carries that code, so rock entered
## fits as flammable non-forest (20% of ELF 14.4's "flammable" pixels in the 2026-09 land cover).
## The default now comes from `fireSenseUtils::fireSenseNonflammableLCC`, the single source of
## truth `makeFireSenseLCC()` itself also uses.

test_that("nonflammableLCC defaults to fireSenseUtils::fireSenseNonflammableLCC", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  default <- md$parameters$default[md$parameters$paramName == "nonflammableLCC"][[1]]
  expect_identical(default, fireSenseUtils::fireSenseNonflammableLCC)
})

test_that("the flammability step marks SCANFI rock/exposed (class 30) non-flammable", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  default <- md$parameters$default[md$parameters$paramName == "nonflammableLCC"][[1]]

  ## Toy land cover: coniferous forest (flammable), SCANFI rock/exposed, water.
  toyLCC <- terra::rast(nrows = 1, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 1,
                        crs = "EPSG:3978")
  terra::values(toyLCC) <- c(210L, 30L, 20L)

  flammable <- LandR::defineFlammable(toyLCC, nonFlammClasses = default)
  expect_equal(as.vector(terra::values(flammable, mat = FALSE)), c(1, 0, 0))
})
