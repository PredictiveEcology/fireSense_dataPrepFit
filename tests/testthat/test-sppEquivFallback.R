## A zero-row sppEquiv is fireSense_ELFs' statement that the ELF has no tree species. The
## .inputObjects fallback used to rebuild the species list whenever NROW was 0, which (for ELF
## 3.2.5, 2026-09-12) found species over this module's studyArea and sent the nested
## Biomass_borealDataPrep run down the with-species path. Only a missing table may fall back.
test_that(".inputObjects rebuilds sppEquiv only when it is missing, not when it has zero rows", {
  src <- readLines(testthat::test_path("..", "..", "fireSense_dataPrepFit.R"))
  cond <- grep("suppliedElsewhere(\"sppEquiv\", sim", src, value = TRUE, fixed = TRUE)
  expect_length(cond, 1L)
  expect_match(cond, "is.null(sim$sppEquiv)", fixed = TRUE)
  expect_false(grepl("NROW(sim$sppEquiv) == 0", cond, fixed = TRUE))
})
