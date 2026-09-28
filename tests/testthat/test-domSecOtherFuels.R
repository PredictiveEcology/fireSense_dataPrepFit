## fireSense_dataPrepFit.R ~690 (prepare_SpreadFit()): fireSenseCovariatesCreate() was never given
## `rstLCC`, so `treedWetland` never appeared, and there was no way to ask for the new
## dom/sec/other AGB fuel representation. `fuelCovariates` (default "domSecOther") now picks the
## representation, and the dominant/secondary classes are chosen once per ELF (not once per data
## year, or every prediction would build a different pair of columns) via
## fireSenseUtils::chooseDomSecFuelClasses(), stored in sim$fuelClassRoles.

## generic AST walk: apply `f` to every call node in `expr`, recursing into every part
walkCalls <- function(expr, f) {
  if (is.call(expr)) {
    f(expr)
    for (part in as.list(expr)) if (!missing(part) && !is.null(part)) walkCalls(part, f)
  }
  invisible(NULL)
}

test_that("fuelCovariates defaults to domSecOther, with species as the only other choice", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  def <- stats::setNames(md$parameters$default, md$parameters$paramName)
  expect_identical(def$fuelCovariates, c("domSecOther", "species"))
})

test_that("reqdPkgs floors fireSenseUtils at >= 0.2.3.9053 (chooseDomSecFuelClasses etc.)", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  fsu <- grep("fireSenseUtils", md$reqdPkgs, value = TRUE)
  expect_length(fsu, 1L)
  expect_match(fsu, "0\\.2\\.3\\.9053")
})

test_that("prepare_SpreadFit() picks fuelClassRoles once per ELF via chooseDomSecFuelClasses()", {
  body <- body(prepare_SpreadFit)
  found <- FALSE
  walkCalls(body, function(e) {
    if (identical(deparse(e[[1]]), "<-") && is.call(e[[2]]) &&
        identical(deparse(e[[2]]), "sim$fuelClassRoles") &&
        any(grepl("chooseDomSecFuelClasses", deparse(e[[3]]), fixed = TRUE)))
      found <<- TRUE
  })
  expect_true(found)
})

test_that("chooseDomSecFuelClasses() is asked about one representative data year, not every year", {
  ## "once per ELF": cohortData/pixelGroupMap/flammableRTM/landcoverDT are each a single tail()
  ## element, not the per-data-year lists (sim$cohortDatas etc.) fireSenseCovariatesCreate()'s
  ## Map() call uses.
  body <- body(prepare_SpreadFit)
  callTxt <- NULL
  walkCalls(body, function(e) {
    if (identical(deparse(e[[1]]), "<-") && is.call(e[[2]]) &&
        identical(deparse(e[[2]]), "sim$fuelClassRoles") &&
        any(grepl("chooseDomSecFuelClasses", deparse(e[[3]]), fixed = TRUE)))
      callTxt <<- gsub("\\s+", " ", paste(deparse(e[[3]]), collapse = " "))
  })
  expect_match(callTxt, "tail(sim$cohortDatas, 1)", fixed = TRUE)
  expect_false(grepl("cohortData = sim$cohortDatas,", callTxt, fixed = TRUE))
})

test_that("the fireSenseCovariatesCreate() Map() call passes rstLCC, fuelCovariates, domClass, secClass", {
  body <- body(prepare_SpreadFit)
  mapCall <- NULL
  walkCalls(body, function(e) {
    if (identical(deparse(e[[1]]), "Map") && !is.null(names(e)) && "f" %in% names(e) &&
        grepl("fireSenseCovariatesCreate", deparse(e[["f"]]), fixed = TRUE))
      mapCall <<- e
  })
  expect_false(is.null(mapCall))
  expect_true("rstLCC" %in% names(mapCall))
  expect_identical(deparse(mapCall[["rstLCC"]]), "sim$rstLCCs")

  moreArgsTxt <- gsub("\\s+", " ", paste(deparse(mapCall[["MoreArgs"]]), collapse = " "))
  expect_match(moreArgsTxt, "fuelCovariates = fuelCovariates", fixed = TRUE)
  expect_match(moreArgsTxt, "domClass = sim$fuelClassRoles$domClass", fixed = TRUE)
  expect_match(moreArgsTxt, "secClass = sim$fuelClassRoles$secClass", fixed = TRUE)
})
