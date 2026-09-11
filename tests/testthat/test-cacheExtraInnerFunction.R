## reproducible::Cache digests only the called function's own code, so a cached call keeps returning its
## old result when a function it calls changes. The land-cover and stand-age calls therefore pass the
## LandR functions their results depend on in `.cacheExtra`. This pins the behaviour those calls rely on.
test_that("a function in .cacheExtra makes Cache re-run when that function's code changes", {
  withr::local_options(reproducible.cachePath = withr::local_tempdir(), reproducible.verbose = -2)
  inner <- function(x) x + 1
  outer <- function(x) inner(x)

  expect_identical(reproducible::Cache(outer(1)), 2)
  expect_identical(reproducible::Cache(outer(1), .cacheExtra = list(inner)), 2)

  inner <- function(x) x + 100
  ## without it, the old result comes back
  expect_identical(reproducible::Cache(outer(1)), 2)
  ## with it, the call runs again
  expect_identical(reproducible::Cache(outer(1), .cacheExtra = list(inner)), 101)
})
