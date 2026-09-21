## Metadata that the code never uses misleads a project into setting or supplying it. Six such
## entries were removed at 1.2.0.9007, with the non-xgboost ignition formula (it needed an object,
## `ignitionClimate`, that was never defined). `.useCache` and `.useCacheArgs` are read by SpaDES.core.
test_that("every declared parameter and input is read, and every declared output is assigned", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  files <- c(file.path(moduleRoot, paste0(moduleName, ".R")),
             list.files(file.path(moduleRoot, "R"), pattern = "\\.R$", full.names = TRUE))
  code <- unlist(lapply(files, readLines, warn = FALSE))
  code <- code[!grepl("defineParameter\\(|expectsInput\\(|createsOutput\\(|^\\s*#", code)]
  used <- function(patterns) vapply(patterns, function(p) any(grepl(p, code, perl = TRUE)), logical(1))

  params <- setdiff(md$parameters$paramName, c(".useCache", ".useCacheArgs"))
  notRead <- params[!used(paste0("(\\$|\")", gsub(".", "\\.", params, fixed = TRUE), "\\b"))]
  expect_identical(notRead, character(0))

  inputs <- md$inputObjects$objectName
  expect_identical(inputs[!used(paste0("\\b", inputs, "\\b"))], character(0))

  outputs <- md$outputObjects$objectName
  expect_identical(outputs[!used(paste0("sim\\$", outputs, "(\\[.*\\])? *<-"))], character(0))

  expect_false(any(grepl("modelAlgorithm|fireSense_ignitionFormula", code)))
})
