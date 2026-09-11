## The fire-buffer steps used to fork one worker per fire year (parallel::mclapply in
## fireSenseUtils::rasterFireBufferDT, parallel::mcMap in fireSenseUtils::bufferToArea).
## In a fits job process every forked child hung forever, waiting on a lock with 0 s of CPU
## from the moment it was forked, so each job stalled at prepSpreadFitData. Both calls must
## run sequentially until forking is safe in that process.
test_that("the fire-buffer steps do not fork", {
  exprs <- parse(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"), keep.source = FALSE)
  fnBody <- function(name) {
    def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                    identical(x[[2]], as.name(name)), exprs)
    expect_length(def, 1L)
    def[[1]][[3]]
  }
  coresValues <- function(body) {
    out <- list()
    walk <- function(x) {
      if (is.call(x)) {
        if (identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name("nCores")))
          out[[length(out) + 1L]] <<- x[[3]]
        for (i in seq_along(x)[-1L]) if (is.call(x[[i]])) walk(x[[i]])
      }
    }
    walk(body)
    out
  }
  for (fn in c("prepare_SpreadFitFire_Raster", "prepare_SpreadFitFire_Vector")) {
    vals <- coresValues(fnBody(fn))
    expect_length(vals, 1L)
    expect_identical(vals[[1]], 1L, info = fn)
  }
})
