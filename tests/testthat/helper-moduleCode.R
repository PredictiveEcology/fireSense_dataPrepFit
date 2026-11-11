## The module's own functions, read from its source (the main file and R/), without running it.

## The definitions `name <- function(...)` of the module, as the `function` calls, by name.
moduleDefs <- function() {
  files <- c(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"),
             list.files(testthat::test_path("..", "..", "R"), "\\.R$", full.names = TRUE))
  isDef <- function(d) is.call(d) && as.character(d[[1]]) %in% c("<-", "=") && is.name(d[[2]]) &&
    is.call(d[[3]]) && identical(d[[3]][[1]], as.name("function"))
  defs <- Filter(isDef, unlist(lapply(files, function(f) as.list(parse(f, keep.source = FALSE))),
                               recursive = FALSE))
  stats::setNames(lapply(defs, function(d) d[[3]]), vapply(defs, function(d) as.character(d[[2]]), ""))
}

## The module functions reached from `start` (themselves included), through the module's own functions.
moduleFunsFrom <- function(start, defs = moduleDefs()) {
  seen <- character()
  todo <- start
  while (length(todo)) {
    f <- todo[1]
    todo <- todo[-1]
    if (f %in% seen) next
    seen <- c(seen, f)
    todo <- c(todo, setdiff(intersect(all.names(defs[[f]]), names(defs)), seen))
  }
  seen
}
