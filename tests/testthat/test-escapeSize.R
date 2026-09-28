## A fire counts as escaped when it reached `escapeSizeHa` (default 50): the escape model's response and the
## fires the spread model is fitted to. Before, any fire larger than one pixel (about 6 ha) counted.

fires <- data.frame(SIZE_HA = c(0.1, 3, 6, 12, 49.9, 50, 51, 800), id = 1:8)

test_that("an escaped fire reached escapeSizeHa and is larger than one pixel", {
  expect_identical(escapedFires(fires, escapeSizeHa = 50, pixSizeHa = 5.76)$id, 6:8)
  ## a threshold below one pixel still needs the fire to be larger than a pixel (the old rule)
  expect_identical(escapedFires(fires, escapeSizeHa = 0, pixSizeHa = 5.76)$id, 3:8)
})

test_that("the sizes of the fires that did not escape are kept, in ha", {
  expect_equal(nonEscapedFireSizes(fires, 50), c(0.1, 3, 6, 12, 49.9))
  expect_equal(nonEscapedFireSizes(data.frame(SIZE_HA = c(NA, 0, 5)), 50), 5)
})

test_that("the rules work on terra points, as the module's fire points are", {
  skip_if_not_installed("terra")
  v <- terra::vect(cbind(1:8, 1:8), atts = fires)
  expect_identical(escapedFires(v, 50, 5.76)$id, 6:8)
  expect_equal(nonEscapedFireSizes(v, 50), c(0.1, 3, 6, 12, 49.9))
})

test_that("the module uses escapeSizeHa for the escape response, the spread fires and the small-fire sizes", {
  src <- readLines(file.path(modulePath, moduleName, paste0(moduleName, ".R")))
  expect_true(any(grepl("escapedFires(sim$ignitionFirePoints, Par$escapeSizeHa", src, fixed = TRUE)))
  expect_true(any(grepl("escapedFires(x, Par$escapeSizeHa, pixSizeHa", src, fixed = TRUE)))
  expect_true(any(grepl("nonEscapedFireSizes(sim$ignitionFirePoints, Par$escapeSizeHa)", src, fixed = TRUE)))
})

## PR #44 filtered spreadFirePoints with escapedFires() (bigger than one pixel AND at least escapeSizeHa)
## but left spreadFirePolys filtered only by pixel size (`x[[haColname]] > pixSizeHa`). A fire between
## one pixel and escapeSizeHa then survived in the polygons but not in the points, so the two lists no
## longer named the same fires per year and fireSenseUtils::harmonizeFireData()'s nfires_poly/nfires_point
## check stopped with "spread fire point and poly harmonization error in dataPrepFit". Real case: ELF
## 4.2.2, 1985-1989, polys 1,6,28,12,34 vs points 0,5,25,9,28.
test_that("spreadFirePolys are filtered with escapedFires too, so they keep the same fires as spreadFirePoints", {
  src <- readLines(file.path(modulePath, moduleName, paste0(moduleName, ".R")))
  polysFilterLine <- grep("sim$spreadFirePolys <- lapply(sim$spreadFirePolys", src, fixed = TRUE)
  expect_length(polysFilterLine, 1)
  polysFilterBody <- src[polysFilterLine:(polysFilterLine + 3)]
  expect_true(any(grepl("escapedFires(x, Par$escapeSizeHa, pixSizeHa, sizeCol = haColname)",
                        polysFilterBody, fixed = TRUE)))
  ## the old, pixel-size-only filter must be gone
  expect_false(any(grepl("x[[haColname]] > pixSizeHa", polysFilterBody, fixed = TRUE)))
})

test_that("escapedFires keeps the same fires whether applied to points or to polygons", {
  ## a fire of 20 ha: bigger than one pixel (5.76 ha) but below escapeSizeHa (50 ha), the size class
  ## that exposed the mismatch between the points and polys filters
  points <- list(
    "2000" = data.frame(FIRE_ID = 1:3, SIZE_HA = c(20, 60, 3)),
    "2001" = data.frame(FIRE_ID = 4:5, SIZE_HA = c(80, 0.5))
  )
  polys <- lapply(points, identity) # same fires, geometry does not matter here

  pointsFiltered <- lapply(points, escapedFires, escapeSizeHa = 50, pixSizeHa = 5.76, sizeCol = "SIZE_HA")
  polysFiltered <- lapply(polys, escapedFires, escapeSizeHa = 50, pixSizeHa = 5.76, sizeCol = "SIZE_HA")
  expect_identical(lapply(pointsFiltered, `[[`, "FIRE_ID"), lapply(polysFiltered, `[[`, "FIRE_ID"))

  ## the old poly filter (pixel size only) kept the 20 ha fire that escapedFires() drops for points
  polysOldFilter <- lapply(polys, function(x) x[x[["SIZE_HA"]] > 5.76, ])
  expect_false(identical(lapply(pointsFiltered, `[[`, "FIRE_ID"), lapply(polysOldFilter, `[[`, "FIRE_ID")))
})
