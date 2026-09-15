## .inputObjects builds rstLCCs and propFlammables only when nothing else supplies rstLCCs, but its per-year
## function returned `LCC$lcc` and `LCC$flammableProp` whether or not it had built them. Supplying rstLCCs
## therefore stopped every run with "object 'LCC' not found". A supplied rstLCCs must be left as it is.
## The function needs a full simList to run in place, so this evaluates it from the module's parsed source.

moduleExprs <- function() {
  parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE)
}

moduleFunctionBody <- function(name) {
  def <- Filter(function(x) is.call(x) && (identical(x[[1]], as.name("<-")) || identical(x[[1]], as.name("="))) &&
                  identical(x[[2]], as.name(name)), moduleExprs())
  stopifnot(length(def) == 1L)
  body(eval(def[[1]][[3]]))
}

test_that(".inputObjects keeps a supplied rstLCCs instead of failing on land cover it did not build", {
  bod <- as.list(moduleFunctionBody(".inputObjects"))
  outsAssign <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                         identical(x[[2]], as.name("outs")), bod)
  expect_length(outsAssign, 1L)
  perYearExpr <- Filter(function(a) is.call(a) && identical(a[[1]], as.name("function")),
                        as.list(outsAssign[[1]][[3]])[-1])
  expect_length(perYearExpr, 1L)

  e <- new.env()
  e$doRstLCCs <- FALSE      # another module, or the user, supplied rstLCCs
  e$doStandAgeMaps <- FALSE
  e$sim <- list(standAgeMaps = list(year2000 = "standAgeMap for 2000"))
  perYear <- eval(perYearExpr[[1]], envir = e)

  expect_identical(perYear(dyChar = "year2000", dy = 2000),
                   list(standAgeMaps = "standAgeMap for 2000"))
})
