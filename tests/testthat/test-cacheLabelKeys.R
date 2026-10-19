## Phase-1 data preparation reran for every replicate although nothing in the data had changed, because
## two things that only label a run were part of the cache keys: the download folder `destinationPath`
## (inputPath(sim) holds the replicate folder) of makeFireSenseLCC() and prepInputsStandAgeMap(), and the
## parameter `fitFileSuffix` (it picks the spread-fit ledger in Init(), which is never cached; the ledger
## itself reaches the cached events through `.cacheExtra`).

mainFile <- function() {
  f <- list.files(moduleRoot, pattern = "\\.R$", full.names = TRUE)
  f[vapply(f, function(x) any(grepl("defineModule\\(", readLines(x, warn = FALSE))), logical(1))]
}

## the Cache() calls (`x |> Cache(...)` parses to Cache(x, ...)) whose piped call is to `fn`
cacheCallsOf <- function(fn) {
  found <- list()
  walk <- function(e) {
    if (is.call(e)) {
      if (identical(deparse(e[[1]]), "Cache") && length(e) > 1 && is.call(e[[2]]) &&
          identical(deparse(e[[2]][[1]]), fn))
        found[[length(found) + 1L]] <<- e
      for (part in as.list(e)) if (!missing(part) && !is.null(part)) walk(part)
    }
    invisible(NULL)
  }
  for (ex in parse(mainFile(), keep.source = FALSE)) walk(ex)
  found
}

test_that("makeFireSenseLCC() and prepInputsStandAgeMap() are cached without their download folder", {
  for (fn in c("fireSenseUtils::makeFireSenseLCC", "prepInputsStandAgeMap")) {
    calls <- cacheCallsOf(fn)
    expect_length(calls, 1L)
    expect_true("destinationPath" %in% names(calls[[1]][[2]]))
    expect_true("destinationPath" %in% eval(calls[[1]][["omitArgs"]]), info = fn)
  }
})

defaultCacheArgs <- function() {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  md$parameters$default[[which(md$parameters$paramName == ".useCacheArgs")]]
}

test_that("fitFileSuffix is left out of the cache key of every cached event", {
  args <- defaultCacheArgs()
  for (ev in c(".inputObjects", "dataPrepBuild", "prepIgnitionFitData", "prepEscapeFitData", "prepSpreadFitData"))
    expect_true("fitFileSuffix" %in% args[[ev]]$omitParams, info = ev)
})

test_that("only fitFileSuffix is omitted, it is a parameter of the module, and no cached code reads it", {
  args <- defaultCacheArgs()
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  omitted <- unique(unlist(lapply(args, `[[`, "omitParams")))
  expect_identical(omitted, "fitFileSuffix")
  expect_true(omitted %in% md$parameters$paramName)
  ## only Init() reads it, and Init() is never cached
  users <- Filter(function(d) "fitFileSuffix" %in% all.names(d) || any(grepl("fitFileSuffix", deparse(d))),
                  lapply(as.list(parse(mainFile(), keep.source = FALSE)), function(e)
                    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
                        !identical(deparse(e[[2]]), "defineModule")) e else NULL))
  users <- Filter(Negate(is.null), users)
  expect_identical(vapply(users, function(d) as.character(d[[2]]), ""), "Init")
  expect_true("init" %in% md$parameters$default[[which(md$parameters$paramName == ".neverCache")]])
})
