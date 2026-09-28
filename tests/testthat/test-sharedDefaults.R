## Root cause: fireSense_dataPrepFit.R and fireSense_dataPrepPredict.R each hard-coded the same
## defaults (forestedLCC, cutoffForYoungAge, nonForestCanBeYoungAge, flammabilityThreshold,
## fuelClassCol, igAggFactor) independently, as nonflammableLCC used to (PR #49). A fit only
## matches its predictions if both modules use the same values, so the defaults now come from
## fireSenseUtils::fireSenseSharedDefaults (fireSenseUtils PR #105), the single source of truth
## fireSense_dataPrepPredict also uses.

test_that("shared-default parameters equal the fireSenseUtils constants", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  paramDefault <- function(name) md$parameters$default[md$parameters$paramName == name][[1]]

  expect_identical(paramDefault("forestedLCC"), fireSenseUtils::fireSenseForestedLCC)
  expect_identical(paramDefault("cutoffForYoungAge"), fireSenseUtils::fireSenseYoungAgeCutoff)
  expect_identical(paramDefault("nonForestCanBeYoungAge"), fireSenseUtils::fireSenseNonForestCanBeYoungAge)
  expect_identical(paramDefault("flammabilityThreshold"), fireSenseUtils::fireSenseFlammabilityThreshold)
  expect_identical(paramDefault("fuelClassCol"), fireSenseUtils::fireSenseFuelClassCol)
  expect_identical(paramDefault("igAggFactor"), fireSenseUtils::fireSenseIgAggFactor)
  expect_identical(paramDefault("scanfiVersion"), fireSenseUtils::fireSenseSCANFIVersion)
})

## the module's cached makeFireSenseLCC() call, as in test-lccCacheExtra.R
makeFireSenseLCCCallExpr <- function() {
  mainFiles <- list.files(moduleRoot, pattern = "\\.R$", full.names = TRUE)
  mainFiles <- mainFiles[vapply(mainFiles, function(f)
    any(grepl("defineModule\\(", readLines(f, warn = FALSE))), logical(1))]
  expect_length(mainFiles, 1L)
  found <- list()
  walk <- function(e) {
    if (is.call(e)) {
      if (identical(deparse(e[[1]]), "fireSenseUtils::makeFireSenseLCC"))
        found[[length(found) + 1L]] <<- e
      for (part in as.list(e)) if (!missing(part) && !is.null(part)) walk(part)
    }
    invisible(NULL)
  }
  for (ex in parse(mainFiles)) walk(ex)
  found
}

test_that("every makeFireSenseLCC() call passes scanfiVersion", {
  calls <- makeFireSenseLCCCallExpr()
  expect_length(calls, 1L)
  args <- names(calls[[1]])
  expect_true("scanfiVersion" %in% args)
  expect_identical(deparse(calls[[1]][["scanfiVersion"]]), "P(sim)$scanfiVersion")
})
