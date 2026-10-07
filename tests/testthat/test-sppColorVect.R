## Init() read every ledger row that intersects the study area (the own ELF and its neighbours) and
## stopped with "sppColVect has unique colour values for a single species" when a species had a
## different colour in two rows. Colours only drive plots: keep the own ELF's colour and carry on.

ids <- c("4.2.2", "13.1")
test_that("a species with different colours in two ledger rows keeps the own ELF's colour", {
  rows <- list(c(Abie_bal = "#111111", Pice_mar = "#222222", Mixed = "#999999"),   # neighbour
               c(Abie_bal = "#AAAAAA", Pice_gla = "#333333", Mixed = "#999999"))   # own ELF
  expect_no_error(out <- suppressMessages(mergeSppColorVects(rows, "13.1", ids)))
  expect_identical(out[["Abie_bal"]], "#AAAAAA")
  expect_identical(sort(names(out)), sort(c("Abie_bal", "Pice_gla", "Pice_mar", "Mixed")))
  expect_identical(names(out)[length(out)], "Mixed")
  expect_false(anyDuplicated(names(out)) > 0)
  expect_message(mergeSppColorVects(rows, "13.1", ids), "Abie_bal")
})

test_that("without the own ELF in the ledger the first row's colour is kept, and agreeing rows are silent", {
  rows <- list(c(Abie_bal = "#111111", Mixed = "#999999"), c(Abie_bal = "#AAAAAA", Mixed = "#999999"))
  expect_identical(suppressMessages(mergeSppColorVects(rows, NULL, ids[1:2]))[["Abie_bal"]], "#111111")
  same <- list(c(Abie_bal = "#111111", Mixed = "#999999"), c(Abie_bal = "#111111", Mixed = "#999999"))
  expect_no_message(mergeSppColorVects(same, "A", c("A", "B")))
})
