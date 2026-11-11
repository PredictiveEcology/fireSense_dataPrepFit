## SpaDES.core (>= the version that reads `createsOutput(events = )`) restores, on a cache hit of an event,
## only the outputs declared for that event type. An output declared for the wrong event is restored by the
## wrong cache hit (a stale copy over a fresh one) or never restored. So `events` must be exactly the events
## whose code assigns the object. This reads the code: every `sim$x <-`, `sim[["x"]] <-`, `sim$x$y <-`,
## `names(sim$x) <-`, `attr(sim$x, ..) <-` is a write of `x`, and the event is the one that calls the function.
test_that("each output's events are exactly the events that assign it", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_true("events" %in% names(md$outputObjects))
  ## the module's functions are in its main file, or in R/ once convertToPackage() has moved them
  files <- c(file.path(moduleRoot, paste0(moduleName, ".R")),
             list.files(file.path(moduleRoot, "R"), pattern = "\\.R$", full.names = TRUE))
  exprs <- unlist(lapply(files, function(f) as.list(parse(f, keep.source = FALSE))), recursive = FALSE)
  funs <- list()
  for (e in exprs)
    if (is.call(e) && is.name(e[[1]]) && as.character(e[[1]]) %in% c("<-", "=") && is.name(e[[2]]) &&
        is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function")))
      funs[[as.character(e[[2]])]] <- e[[3]]
  rootOf <- function(lhs) {  # the `x` of the `sim$x` or `sim[["x"]]` at the root of an assignment target
    while (is.call(lhs)) {
      if (is.name(lhs[[1]]) && as.character(lhs[[1]]) %in% c("$", "[[") && identical(lhs[[2]], as.name("sim")))
        return(as.character(lhs[[3]]))
      lhs <- lhs[[2]]
    }
    NULL
  }
  written <- function(fn) {
    out <- character()
    walk <- function(x) {
      if (is.call(x)) {
        if (is.name(x[[1]]) && as.character(x[[1]]) %in% c("<-", "=") && length(x) == 3)
          out <<- c(out, rootOf(x[[2]]))
        for (i in seq_along(x)) if (!identical(x[[i]], quote(expr = ))) walk(x[[i]])  # `x[ , 1]` has an empty argument
      }
    }
    walk(fn[[3]])  # a `function` call: formals, body
    unique(out)
  }
  ## the event dispatcher calls: init -> Init, prepSpreadFitData -> prepare_SpreadFit (which calls the
  ## spread fire-data builders), and so on; `.inputObjects` is its own function
  eventFuns <- list(init = "Init", dataPrepBuild = "dataPrepBuild",
                    prepIgnitionFitData = "prepare_IgnitionFit", prepEscapeFitData = "prepare_EscapeFit",
                    prepSpreadFitData = c("prepare_SpreadFit", "prepare_SpreadFitFire_Vector"),
                    plotAndMessage = "plotAndMessage", cleanUp = "cleanUpMod",
                    .inputObjects = ".inputObjects")
  expect_true(all(unlist(eventFuns) %in% names(funs)))
  writers <- lapply(eventFuns, function(f) unique(unlist(lapply(funs[f], written))))
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
