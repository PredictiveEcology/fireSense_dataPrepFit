## runBorealDP_forCohortData() keys its cached Biomass_borealDataPrep run on a digest of the modules'
## code. It digested `dir(..., pattern = ".R$")`, which returns bare file names; the digest of a name
## that does not resolve ignores the file's content, so an edited module still hit the old cached run.
test_that("moduleCodeDigest changes when a module's code changes", {
  modulePath <- withr::local_tempdir()
  modDir <- file.path(modulePath, "someModule")
  dir.create(file.path(modDir, "R"), recursive = TRUE)
  writeLines("x <- 1", file.path(modDir, "someModule.R"))
  writeLines("f <- function() 1", file.path(modDir, "R", "helpers.R"))
  withr::local_dir(withr::local_tempdir()) # the working directory is not the module path

  before <- moduleCodeDigest(modulePath, "someModule")
  expect_identical(moduleCodeDigest(modulePath, "someModule"), before)

  writeLines("f <- function() 2", file.path(modDir, "R", "helpers.R"))
  expect_false(identical(moduleCodeDigest(modulePath, "someModule"), before))
})
