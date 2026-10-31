## A project that sets no `.useCache` gets the module's default. It must cache what the project script
## used to switch on for this module: `.inputObjects` and the four data-preparation events. (The
## project's `.globals$.useCache` also listed `init`, `initPlot`, `estimateThreshold`, `spreadFitPrepare` and
## `checkData`; this module's own value beat it, and `init` is never cached, see `.neverCache`.)

md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
useCacheDefault <- md$parameters$default[[which(md$parameters$paramName == ".useCache")]]
neverCache <- md$parameters$default[[which(md$parameters$paramName == ".neverCache")]]
code <- readLines(file.path(moduleRoot, paste0(moduleName, ".R")))

test_that("the default .useCache is .inputObjects and every data-preparation event of this module", {
  expect_setequal(useCacheDefault, c(".inputObjects", "dataPrepBuild", "prepEscapeFitData",
                                     "prepIgnitionFitData", "prepSpreadFitData"))
})

test_that("each name in the default .useCache is an event or .inputObjects of this module", {
  exprs <- parse(text = code, keep.source = FALSE)
  assigned <- vapply(exprs, function(e) if (is.call(e) && as.character(e[[1]])[1] %in% c("=", "<-")) as.character(e[[2]])[1] else "", "")
  doEvent <- exprs[[which(assigned == paste0("doEvent.", moduleName))]][[3]]
  switchCall <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("switch")), as.list(doEvent[[3]]))[[1]]
  events <- setdiff(names(switchCall), "")
  expect_true(".inputObjects" %in% assigned)
  expect_true(all(setdiff(useCacheDefault, ".inputObjects") %in% events))
})

test_that("the default .useCache does not name an event that .neverCache forbids", {
  expect_length(intersect(useCacheDefault, neverCache), 0L)
})
