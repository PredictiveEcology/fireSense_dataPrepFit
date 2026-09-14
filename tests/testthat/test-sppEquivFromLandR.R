## .inputObjects takes sppEquiv from LandR::speciesInStudyArea()$sppEquiv (LandR >= 1.2.0.9012),
## the same table fireSense_ELFs uses, instead of building its own slightly different one (it
## kept `_Spp` genus entries and did not merge Engelmann spruce). This checks the wiring: the
## statement, parsed from the module file, with the LandR call stubbed.

sppEquivStatement <- function() {
  exprs <- parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE)
  found <- list()
  walk <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], as.name("<-")) && identical(deparse(x[[2]]), "sim$sppEquiv") &&
          grepl("speciesInStudyArea", paste(deparse(x[[3]]), collapse = " ")))
        found[[length(found) + 1L]] <<- x
      for (i in seq_along(x)[-1L]) if (is.call(x[[i]])) walk(x[[i]])
    }
  }
  for (e in exprs) walk(e)
  found
}

test_that(".inputObjects takes sppEquiv from LandR::speciesInStudyArea()", {
  stmt <- sppEquivStatement()
  expect_length(stmt, 1L)
  skip_if(length(stmt) != 1L)

  sentinel <- data.frame(LandR = "Abie_ama")
  args <- NULL
  local_mocked_bindings(speciesInStudyArea = function(...) {
    args <<- list(...)
    list(speciesList = "ABIE_AMA", sppEquiv = sentinel)
  }, .package = "LandR")

  env <- new.env(parent = globalenv())
  env$sim <- new.env()
  env$sim$studyArea <- "the study area"
  env$Par <- list(sppEquivCol = "LandR")
  env$inputPath <- function(sim) tempdir()
  eval(stmt[[1]], env)

  expect_identical(env$sim$sppEquiv, sentinel)
  expect_identical(args$sppEquivCol, "LandR")
  expect_identical(args$studyArea, "the study area")
})
