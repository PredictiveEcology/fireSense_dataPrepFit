## reproducible::Cache() digests only the called function's own code, so the cached makeFireSenseLCC()
## call has to name the functions IT calls in `.cacheExtra`. Which functions those are depends on
## `lccSource` -- a run-time option -- so a list written out here is wrong for the other source and goes
## stale when the default moves. It did: this module pinned LandR::prepInputs_NTEMS_LCC_FAO and kept it
## after fireSenseUtils #64 made SCANFI the default, so SCANFI and CWIM changes stopped invalidating the
## cache while NTEMS changes invalidated it for nothing. fireSenseUtils::makeFireSenseLCCDeps() owns the
## answer; this pins that we ask it rather than guessing.

## the .cacheExtra expression from the module's cached makeFireSenseLCC() call
lccCacheExtraExpr <- function() {
  mainFiles <- list.files(moduleRoot, pattern = "\\.R$", full.names = TRUE)
  mainFiles <- mainFiles[vapply(mainFiles, function(f)
    any(grepl("defineModule\\(", readLines(f, warn = FALSE))), logical(1))]
  expect_length(mainFiles, 1L)
  found <- list()
  walk <- function(e) {
    if (is.call(e)) {
      nms <- names(e)
      ## `x |> Cache(...)` parses to Cache(x, ...), so the piped call is the first argument
      if (identical(deparse(e[[1]]), "Cache") && !is.null(nms) && ".cacheExtra" %in% nms &&
          any(grepl("makeFireSenseLCC", deparse(e[[2]]))))
        found[[length(found) + 1L]] <<- e[[".cacheExtra"]]
      for (part in as.list(e)) if (!missing(part) && !is.null(part)) walk(part)
    }
    invisible(NULL)
  }
  for (ex in parse(mainFiles)) walk(ex)
  found
}

test_that("the cached makeFireSenseLCC() call asks fireSenseUtils which functions to digest", {
  exprs <- lccCacheExtraExpr()
  expect_length(exprs, 1L)
  txt <- deparse(exprs[[1]])
  expect_match(paste(txt, collapse = " "), "makeFireSenseLCCDeps")
  ## the hard-coded list is what went stale; it must not come back
  expect_false(any(grepl("prepInputs_NTEMS_LCC_FAO|prepInputs_SCANFI_LCC_FAO", txt)))
})

test_that("that .cacheExtra changes with lccSource, so a SCANFI run cannot reuse an NTEMS cache", {
  skip_if_not(requireNamespace("fireSenseUtils", quietly = TRUE) &&
                "makeFireSenseLCCDeps" %in% getNamespaceExports("fireSenseUtils"),
              "fireSenseUtils without makeFireSenseLCCDeps() (PredictiveEcology/fireSenseUtils#65)")
  skip_if_not("prepInputs_CWIM" %in% getNamespaceExports("LandR"),
              "LandR without prepInputs_CWIM (PredictiveEcology/LandR#228)")
  expr <- lccCacheExtraExpr()[[1]]
  scanfi <- withr::with_options(list(fireSense.lccSource = "SCANFI"), eval(expr))
  ntems  <- withr::with_options(list(fireSense.lccSource = "NTEMS"),  eval(expr))
  expect_true(all(vapply(scanfi, is.function, logical(1))))
  ## the user-visible consequence, and what the hard-coded list got wrong
  expect_false(identical(reproducible::.robustDigest(scanfi), reproducible::.robustDigest(ntems)))
})
