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
