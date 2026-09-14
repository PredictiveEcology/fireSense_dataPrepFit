## A NULL `sim$sppEquiv` (left by fireSense_ELFs for ELFs without Engelmann spruce) used to
## reach Biomass_speciesData inside runBorealDP_forCohortData() and stop in LandR with
## "object 'LandR' not found". Two places let it through:
##   1. .inputObjects skipped building sppEquiv because suppliedElsewhere() was TRUE, even
##      though the value was NULL;
##   2. runBorealDP_forCohortData() passed the NULL on as `objects =`, and simInit applies
##      `objects` after .inputObjects, so it replaced the table Biomass_speciesData built.
##
## Both need a full simList to run, so these evaluate the relevant pieces of the module's
## own source (parsed, not sourced) with small stand-ins.

moduleExprs <- function() {
  parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE)
}

moduleFunctionBody <- function(name) {
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                  identical(x[[2]], as.name(name)), moduleExprs())
  stopifnot(length(def) == 1L)
  body(eval(def[[1]][[3]]))
}

test_that("runBorealDP_forCohortData does not pass NULL objects into the nested simInit", {
  bod <- as.list(moduleFunctionBody("runBorealDP_forCohortData"))
  ## the statements that build `objsNeeded` (and `ecoFile`, which they use), in order
  assigns <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                      is.name(x[[2]]) && as.character(x[[2]]) %in% c("ecoFile", "objsNeeded"), bod)
  expect_gte(length(assigns), 3L)

  simEnv <- new.env()
  assign("sppEquiv", NULL, envir = simEnv)
  assign("studyArea", "a study area", envir = simEnv)
  e <- new.env()
  e$sim <- simEnv
  e$envir <- function(x) x
  for (a in assigns) eval(a, e)

  expect_false("sppEquiv" %in% names(e$objsNeeded))
  expect_identical(e$objsNeeded[["studyArea"]], "a study area")
})

test_that(".inputObjects builds sppEquiv when the supplied value is NULL", {
  bod <- as.list(moduleFunctionBody(".inputObjects"))
  ifs <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("if")) &&
                  any(grepl('suppliedElsewhere\\("sppEquiv"', deparse(x[[2]]))), bod)
  expect_length(ifs, 1L)
  cond <- ifs[[1]][[2]]

  e <- new.env()
  e$suppliedElsewhere <- function(...) TRUE  # another module declares sppEquiv as an output
  e$sim <- list(sppEquiv = NULL)
  expect_true(eval(cond, e))                  # NULL: build it

  e$sim <- list(sppEquiv = data.frame(LandR = "Pice_mar"))
  expect_false(eval(cond, e))                 # a real table: leave it alone
})
