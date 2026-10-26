## 2026-10-26: prepare_SpreadFitFire_Vector() kept a year's fire polygons only when its points' IDs were in the
## same order: `poly[[IDvar]] %in% point[[IDvar]]` compared two data.frames. 6.6.1 lost 64 of 127 polygons, and
## with them every escape of 9 years (escapes are counted from these polygons since #82).

test_that("withIDsIn keeps the rows whose ID is in the other vector, in any order", {
  sq <- function(x0) sprintf("POLYGON ((%f 0, %f 0, %f 100, %f 100, %f 0))", x0, x0 + 100, x0 + 100, x0, x0)
  polys <- terra::vect(sq(c(0, 200, 400)), crs = "EPSG:3978")
  polys$FIRE_ID <- c(11, 22, 33)
  pts <- terra::vect(cbind(c(450, 50, 250), 50), crs = "EPSG:3978")
  pts$FIRE_ID <- c(33, 11, 22) # same IDs, another order, as harmonizeFireData() returns them
  ## the old form: one value for the whole column
  expect_length(polys[["FIRE_ID"]] %in% pts[["FIRE_ID"]], 1L)
  expect_false(polys[["FIRE_ID"]] %in% pts[["FIRE_ID"]])
  expect_equal(nrow(withIDsIn(polys, pts, "FIRE_ID")), 3)
  ## a polygon without a point goes, the others stay
  kept <- withIDsIn(polys, pts[pts$FIRE_ID != 22, ], "FIRE_ID")
  expect_identical(kept$FIRE_ID, c(11, 33))
  expect_identical(withIDsIn(pts, polys[1, ], "FIRE_ID")$FIRE_ID, 11)
})
