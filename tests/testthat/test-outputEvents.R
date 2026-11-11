## SpaDES.core (>= the version that reads `createsOutput(events = )`) restores, on a cache hit of an event,
## only the outputs declared for that event type. An output declared for the wrong event is restored by the
## wrong cache hit (a stale copy over a fresh one) or never restored. So `events` must be exactly the events
## whose code assigns the object. This reads the code: `SpaDES.core:::.findElement(type = "assign")` finds the objects each function
## assigns into `sim`, and an event's functions are those its branch of `doEvent` reaches (`moduleFunsFrom()`).
test_that("each output's events are exactly the events that assign it", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_true("events" %in% names(md$outputObjects))
  defs <- moduleDefs()
  ## .findElement() takes `sim$x$y <-` for a write of `y`, and not `=`, `names(sim$x) <-` or `sim$x[i] <-`: cut every
  ## assignment target back to its `sim$x` / `sim[["x"]]` root first
  toRoot <- function(x) {
    if (!is.call(x)) return(x)
    for (i in seq_along(x)) if (!identical(x[[i]], quote(expr = ))) x[i] <- list(toRoot(x[[i]]))  # `x[ , 1]` has an empty argument
    if (is.name(x[[1]]) && as.character(x[[1]]) %in% c("<-", "=") && length(x) == 3) {
      x[[1]] <- as.name("<-")
      while (is.call(x[[2]]) && !(as.character(x[[2]][[1]])[1] %in% c("$", "[[") && identical(x[[2]][[2]], as.name("sim"))))
        x[[2]] <- x[[2]][[2]]
    }
    x
  }
  written <- function(funs) unique(unlist(lapply(funs, function(f) SpaDES.core:::.findElement(toRoot(defs[[f]][[3]]), type = "assign"))))
  ## the event dispatcher's branches: init -> Init, prepSpreadFitData -> prepare_SpreadFit (which calls the
  ## spread fire-data builders), and so on; `.inputObjects` is its own function
  branches <- as.list(Find(function(x) is.call(x) && identical(x[[1]], as.name("switch")),
                           as.list(defs[[paste0("doEvent.", moduleName)]][[3]])))[-(1:2)]
  branches <- branches[nzchar(names(branches))]
  start <- lapply(branches, function(b) intersect(all.names(b), names(defs)))
  start[[".inputObjects"]] <- ".inputObjects"
  expect_true(all(lengths(start) > 0))
  writers <- lapply(start, function(f) written(moduleFunsFrom(f, defs)))
  ## not an assignment to `sim$x`: `.inputObjects` puts the land cover it builds into the simList with
  ## `list2env(<named list>, envir = envir(sim))` (`standAgeMaps`, `propFlammables` and `rstLCCs`; of these only
  ## `rstLCCs` is an output)
  writers[[".inputObjects"]] <- union(writers[[".inputObjects"]], "rstLCCs")
  outputs <- md$outputObjects
  for (i in seq_len(nrow(outputs))) {
    nm <- outputs$objectName[i]
    expected <- names(writers)[vapply(writers, function(w) nm %in% w, logical(1))]
    expect_setequal(outputs$events[[i]], expected)
  }
})
