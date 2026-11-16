## The nested Biomass_borealDataPrep runs drew random numbers with no seed, so the same ELF
## got different cohortData on different machines. `vegInitSeed` seeds them through `.seed`.

test_that("vegInitSeed becomes the nested module's .seed for init; NA leaves it as it was", {
  expect_identical(nestedVegSeed(1), list(init = 1))
  expect_identical(nestedVegSeed(7, list(init = 3, other = 2)), list(init = 7, other = 2))
  expect_null(nestedVegSeed(NA_real_))
  expect_identical(nestedVegSeed(NA, list(init = 3)), list(init = 3))
})

test_that("the module defines vegInitSeed with a fixed default seed", {
  pars <- moduleParamTable()
  default <- pars$default[[match("vegInitSeed", pars$paramName)]]
  expect_equal(default, 1)
  expect_true(is.numeric(default))
})

test_that("a seeded init event draws the same subset for the same seed, another for another; cached", {
  skip_if_not_installed("LandR")
  data.table::setDTthreads(1)
  modDir <- file.path(withr::local_tempdir(), "modules")
  dir.create(file.path(modDir, "toyDraw"), recursive = TRUE)
  writeLines(c(
    'defineModule(sim, list(name = "toyDraw", description = "toy", keywords = "toy", childModules = character(0),',
    '  authors = person("A", "B", email = "a@b.c", role = c("aut", "cre")), spatialExtent = terra::ext(rep(0, 4)),',
    '  loadOrder = list(after = NULL, before = NULL), version = list(toyDraw = "0.0.1"), timeframe = as.POSIXlt(c(NA, NA)),',
    '  timeunit = "year", citation = list(), documentation = list(), reqdPkgs = list("data.table", "LandR"),',
    '  parameters = bindrows(defineParameter(".useCache", "character", "init", NA, NA, "")),',
    '  inputObjects = bindrows(), outputObjects = bindrows(createsOutput("picked", "data.table", "")) ))',
    'doEvent.toyDraw <- function(sim, eventTime, eventType) {',
    '  if (eventType == "init") {',
    '    DT <- data.table::data.table(g = rep(1:2, each = 200), v = 1:400)',
    '    sim$picked <- LandR::subsetDT(DT, by = "g", doSubset = 5)',
    '  }',
    '  invisible(sim)',
    '}'), file.path(modDir, "toyDraw", "toyDraw.R"))
  run <- function(seed, cachePath) {
    pars <- if (is.na(seed)) list() else list(toyDraw = list(.seed = nestedVegSeed(seed)))
    SpaDES.core::simInitAndSpades(
      modules = "toyDraw", params = pars, times = list(start = 0, end = 0),
      paths = list(modulePath = modDir, cachePath = cachePath, inputPath = cachePath,
                   outputPath = cachePath))$picked$v
  }
  withr::local_options(list(reproducible.useCache = TRUE, spades.useCache = "all"))
  cp <- withr::local_tempdir()
  a <- run(1, cp); b <- run(1, cp); c2 <- run(2, cp)
  expect_identical(a, b)
  expect_false(identical(a, c2))
  ## the same seed in a fresh cache gives the same draw: it is the seed, not a cache replay
  expect_identical(run(1, withr::local_tempdir()), a)
  ## unseeded draws vary (a fresh cache each time)
  u <- replicate(3, run(NA, withr::local_tempdir()), simplify = FALSE)
  expect_gt(length(unique(u)), 1L)
})

test_that("the outer Cache() of the nested simInitAndSpades() keys on the seed in its params", {
  skip_if_not_installed("LandR")
  data.table::setDTthreads(1)
  modDir <- file.path(withr::local_tempdir(), "modules")
  dir.create(file.path(modDir, "toyDraw"), recursive = TRUE)
  writeLines(c(
    'defineModule(sim, list(name = "toyDraw", description = "toy", keywords = "toy", childModules = character(0),',
    '  authors = person("A", "B", email = "a@b.c", role = c("aut", "cre")), spatialExtent = terra::ext(rep(0, 4)),',
    '  loadOrder = list(after = NULL, before = NULL), version = list(toyDraw = "0.0.1"), timeframe = as.POSIXlt(c(NA, NA)),',
    '  timeunit = "year", citation = list(), documentation = list(), reqdPkgs = list("data.table"),',
    '  parameters = bindrows(defineParameter(".useCache", "character", "init", NA, NA, "")),',
    '  inputObjects = bindrows(), outputObjects = bindrows()))',
    'doEvent.toyDraw <- function(sim, eventTime, eventType) invisible(sim)'),
    file.path(modDir, "toyDraw", "toyDraw.R"))
  cp <- withr::local_tempdir()
  ## the same call shape as runBorealDP_forCohortData(): params in, paths omitted from the key
  nested <- function(seed) {
    pathsLocal <- list(modulePath = modDir, cachePath = cp, inputPath = cp, outputPath = withr::local_tempdir())
    parms <- list(toyDraw = list(.seed = nestedVegSeed(seed)))
    SpaDES.core::simInitAndSpades(paths = pathsLocal, params = parms, times = list(start = 0, end = 0),
                                  modules = "toyDraw") |>
      reproducible::Cache(.functionName = "toyNested", omitArgs = "paths", cachePath = cp)
    nrow(reproducible::showCache(cp, userTags = "toyNested", verbose = -2))
  }
  nested(1); nested(1)
  one <- length(unique(reproducible::showCache(cp, userTags = "toyNested", verbose = -2)$cacheId))
  nested(2)
  two <- length(unique(reproducible::showCache(cp, userTags = "toyNested", verbose = -2)$cacheId))
  expect_gt(one, 0L)
  expect_equal(two, 2L * one) # a second seed adds an entry; a repeat of the first adds none
})
