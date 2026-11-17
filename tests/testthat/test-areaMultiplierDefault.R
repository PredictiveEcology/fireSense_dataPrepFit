## Defaults production already passes: areaMultiplier = 10 x fireSenseUtils::multiplier, as a call that is
## identical() to bquote(function(size) .(10) * fireSenseUtils::multiplier(size)) (so the cache digest is unchanged).


test_that("areaMultiplier defaults to 10 x fireSenseUtils::multiplier", {
  params <- moduleParamTable()
  i <- match("areaMultiplier", params$paramName)
  default <- params$default[[i]]
  expect_true("call" %in% params$paramClass[[i]])
  fn <- eval(default)
  for (size in c(1, 10, 100, 1e4, 1e6))
    expect_equal(fn(size), 10 * fireSenseUtils::multiplier(size))
})

test_that("the default areaMultiplier is identical() to the value production passes, with the same digest", {
  ## evaluated from text without srcrefs, as an Rscript run does (testthat parses test files with srcrefs,
  ## which a `function` call would keep as its 4th element)
  production <- eval(parse(text = "bquote(function(size) .(10) * fireSenseUtils::multiplier(size))",
                           keep.source = FALSE))
  params <- moduleParamTable()
  default <- params$default[[match("areaMultiplier", params$paramName)]]
  expect_identical(default, production)
  expect_identical(reproducible::.robustDigest(default), reproducible::.robustDigest(production))
})
