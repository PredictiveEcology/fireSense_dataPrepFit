## `spread = "auto"` chooses among `spreadClimateCandidates`, no longer among the ignition variables.
## Summer maximum temperature (Tmax_sm) and precipitation (PPT_sm) are offered; the ignition list is unchanged.

yearlyRaster <- function(values, years) {
  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  s <- terra::rast(lapply(values, function(v) terra::setValues(r, rep(v, 4))))
  names(s) <- paste0(fireSenseUtils::yearTxt, years)
  s
}

test_that("the default candidates are the default ignition variables plus Tmax_sm and PPT_sm; ignition is unchanged", {
  def <- fireSenseUtils::defaultClimateVariablesForFire
  expect_identical(def$ignition, c("CMD", "cumMDC", "CMD_sm", "CMD_sp"))
  expect_identical(withSpreadCandidates(def, c(def$ignition, "Tmax_sm", "PPT_sm"))$ignition, def$ignition)
  expect_setequal(withSpreadCandidates(def, c(def$ignition, "Tmax_sm", "PPT_sm"))$spread,
                  c(def$ignition, "Tmax_sm", "PPT_sm"))
})

test_that("the candidates are not added to the layers when spread is not auto", {
  cv <- list(ignition = c("CMD", "CMD_sm"), spread = "CMD")
  expect_identical(withSpreadCandidates(cv, c("CMD", "Tmax_sm")), cv)
})

test_that("Tmax_sm and PPT_sm are offered to selectSpreadClimateVariable, whose layers are prepared", {
  skip_if_not_installed("climateData")
  cv <- withSpreadCandidates(fireSenseUtils::defaultClimateVariablesForFire, c("CMD", "Tmax_sm", "PPT_sm"))
  cl <- fireSenseUtils::fireClimateLayers(cv, historicalYears = 1990:1995)
  expect_true(all(c("historical_CMD", "historical_Tmaxsm", "historical_PPTsm") %in% names(cl)))
  ## and what the ignition models take is the untouched list
  expect_false(any(c("Tmax_sm", "PPT_sm") %in% fireSenseUtils::defaultClimateVariablesForFire$ignition))
})

years <- 2001:2012
fireBuffers <- function(burned) {          # per year, `burned` pixels in buffer 1
  setNames(lapply(burned, function(n) data.frame(buffer = c(rep(1, n), rep(0, 5)))), paste0(fireSenseUtils::yearTxt, years))
}
burned <- c(1:9, 100, 200, 300)

test_that("auto offers Tmax_sm and PPT_sm, and chooses the one tracking bad years; PPT is scored as decreasing with dryness", {
  rasters <- list(CMD = yearlyRaster(c(9, 1, 8, 2, 7, 3, 6, 4, 5, 5, 5, 5), years),
                  Tmaxsm = yearlyRaster(c(9, 1, 8, 2, 7, 3, 6, 4, 5, 10, 11, 12), years),
                  PPTsm = yearlyRaster(-c(9, 1, 8, 2, 7, 3, 6, 4, 5, 10, 11, 12), years))
  sel <- autoSpreadClimate(rasters, fireBuffers(burned), years, c("CMD", "Tmax_sm", "PPT_sm"))
  expect_setequal(sel$var, c("CMD", "Tmaxsm", "PPTsm"))
  expect_equal(sel[var == "PPTsm", auc], 1)
  expect_equal(sel[var == "Tmaxsm", auc], 1)
  ## a precipitation layer that falls in the bad years is chosen over one that rises
  rasters$PPTsm <- yearlyRaster(c(9, 1, 8, 2, 7, 3, 6, 4, 5, 10, 11, 12), years)
  sel2 <- autoSpreadClimate(rasters[c("CMD", "PPTsm")], fireBuffers(burned), years, c("CMD", "PPT_sm"))
  expect_lt(sel2[var == "PPTsm", auc], 0.5)
  expect_identical(sel2[chosen == TRUE, var], "CMD")
})

test_that("a run with spreadClimateCandidates = 'CMD' chooses CMD", {
  rasters <- list(CMD = yearlyRaster(c(9, 1, 8, 2, 7, 3, 6, 4, 5, 10, 11, 12), years),
                  Tmaxsm = yearlyRaster(c(9, 1, 8, 2, 7, 3, 6, 4, 5, 10, 11, 12) + 100, years))
  expect_identical(autoSpreadClimate(rasters, fireBuffers(burned), years, "CMD")[chosen == TRUE, var], "CMD")
  expect_identical(autoSpreadClimate(rasters, fireBuffers(burned), years, "CMD")$var, "CMD")
})

test_that("spreadClimateCandidates is a declared parameter whose default is the ignition list plus Tmax_sm and PPT_sm", {
  txt <- paste(readLines(file.path(moduleRoot, "fireSense_dataPrepFit.R")), collapse = "\n")
  expect_match(txt, 'defineParameter\\("spreadClimateCandidates", "character"')
  expect_match(txt, 'fireSenseUtils::defaultClimateVariablesForFire\\$ignition, "Tmax_sm", "PPT_sm"', fixed = FALSE)
  expect_match(txt, "P\\(sim\\)\\$spreadClimateCandidates")
})
