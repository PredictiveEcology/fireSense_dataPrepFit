## Every argument handed to reproducible::Cache() must be one Cache() knows. A stray one (the
## piped call after `Map(f = fuelClassPrep, ...)` passed `.omitArgs`, a typo for `omitArgs`) is
## silently ignored while caching is on -- the two rasters it meant to omit were digested every
## time -- and with caching off (`spades.useCache = "eventsOnly"`, i.e. reproducible.useCache =
## FALSE) it flips Cache()'s "uses dots" heuristic, whose bypass then calls the already-evaluated
## result as a function: 'could not find function "FUN"'. That killed an ELF on 2026-09-13.
test_that("no Cache() call in the module passes a misspelled dot-prefixed argument", {
  ## the main module file is the one holding defineModule(); it is named for the module, which
  ## need not be the name of the checkout directory (a git worktree, for instance)
  mainFiles <- list.files(moduleRoot, pattern = "\\.R$", full.names = TRUE)
  mainFiles <- mainFiles[vapply(mainFiles, function(f)
    any(grepl("defineModule\\(", readLines(f, warn = FALSE))), logical(1))]
  expect_length(mainFiles, 1L)
  files <- c(mainFiles, list.files(file.path(moduleRoot, "R"), pattern = "\\.R$", full.names = TRUE))
  code <- unlist(lapply(files, readLines, warn = FALSE))
  found <- unique(unlist(regmatches(code, gregexpr("\\.[A-Za-z]+(?= *= )", code, perl = TRUE))))
  stray <- intersect(found, c(".omitArgs", ".userTags", ".useCache", ".useCloud", ".notOlderThan",
                              ".showSimilar", ".quick", ".cacheId"))
  expect_identical(stray, character(0))
  expect_false(any(grepl("\\.omitArgs", code)))
})
