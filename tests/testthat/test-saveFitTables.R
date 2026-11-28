## `saveFitTables` writes the ignition / escape fit tables to outputPath(sim) after they are built, so analyses
## outside the fit can use the exact table a fit was trained on. It is its own uncached event: the tables come
## from the simList, whether `prepIgnitionFitData` just ran or its cache restored them.

tableSim <- function(saveFitTables, outDir, ELF = "6.6.2") {
  sim <- suppressMessages(SpaDES.core::simInit(paths = list(outputPath = outDir)))
  sim@params <- list(fireSense_dataPrepFit = list(saveFitTables = saveFitTables))
  sim@current <- data.table::data.table(moduleName = "fireSense_dataPrepFit")
  sim$.ELFind <- ELF
  sim$fireSense_ignitionCovariates <- data.table::data.table(
    pixelID = c(3L, 5L, 3L), year = c(2001, 2001, 2002), ignitions = c(0L, 1L, 2L),
    CMDsm = c(1.5, 2.5, 3.5), coreFraction = c(1, 0.5, 1))
  sim$fireSense_escapeCovariates <- data.frame(sim$fireSense_ignitionCovariates, escapes = c(0L, 1L, 1L)) |> data.table::as.data.table()
  sim
}

## the module's functions get `P()` from the running simList; a test supplies it
runSave <- function(sim) {
  env <- new.env(parent = environment(saveFitTables))
  env$P <- function(sim) sim@params$fireSense_dataPrepFit
  f <- saveFitTables
  environment(f) <- env
  f(sim)
}

test_that("saveFitTables writes each requested table, identical to the one in the simList, and registers it", {
  out <- withr::local_tempdir()
  sim <- tableSim(c("ignition", "escape"), out)
  sim <- runSave(sim)
  files <- list.files(out, "\\.qs2$")
  expect_length(files, 2L)
  expect_true(all(grepl("6.6.2", files, fixed = TRUE)))
  ign <- qs2::qs_read(file.path(out, grep("ignition", files, value = TRUE)))
  esc <- qs2::qs_read(file.path(out, grep("escape", files, value = TRUE)))
  expect_equal(ign, sim$fireSense_ignitionCovariates)
  expect_equal(esc, sim$fireSense_escapeCovariates)
  expect_setequal(basename(SpaDES.core::outputs(sim)$file), files)
})

test_that("saveFitTables writes only the tables asked for, and nothing by default", {
  out <- withr::local_tempdir()
  runSave(tableSim("ignition", out))
  expect_length(list.files(out), 1L)
  out2 <- withr::local_tempdir()
  runSave(tableSim(character(0), out2))
  expect_length(list.files(out2), 0L)
  tbl <- moduleParamTable()
  expect_identical(tbl$default[[which(tbl$paramName == "saveFitTables")]], character(0))
  expect_false(tbl$cacheKey[[which(tbl$paramName == "saveFitTables")]])   # changing it does not re-run cached events
})

test_that("saveFitTables does not change the tables", {
  sim <- tableSim(c("ignition", "escape"), withr::local_tempdir())
  before <- data.table::copy(sim$fireSense_ignitionCovariates)
  sim <- runSave(sim)
  expect_identical(sim$fireSense_ignitionCovariates, before)
})

test_that("a table the run did not build is skipped with a message, not an error", {
  out <- withr::local_tempdir()
  sim <- tableSim("escape", out)
  sim$fireSense_escapeCovariates <- NULL
  expect_message(runSave(sim), "fireSense_escapeCovariates")
  expect_length(list.files(out), 0L)
})

test_that("an unknown table name is refused", {
  expect_error(runSave(tableSim("spread", withr::local_tempdir())), "saveFitTables")
})

test_that("the event is scheduled only when tables are requested, runs after the prep events, and is never cached", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  declared <- md$parameters$default[[which(md$parameters$paramName == ".neverCache")]]
  expect_true("saveFitTables" %in% declared)   # a cache hit would skip the write
  expect_identical(saveFitTablesEvents(character(0)), character(0))
  expect_identical(saveFitTablesEvents("ignition"), "saveFitTables")
})

test_that("init schedules the event after the prep events (priority 2, they are 1) and the dispatcher runs it", {
  defs <- moduleDefs()
  init <- paste(deparse(defs[["doEvent.fireSense_dataPrepFit"]]), collapse = "\n")
  expect_match(init, "saveFitTablesEvents\\(P\\(sim\\)\\$saveFitTables\\)")
  expect_match(init, "\"saveFitTables\", eventPriority = 2|ev, eventPriority = 2")
  expect_match(init, "saveFitTables = \\{\\s*sim <- saveFitTables\\(sim\\)")
})
