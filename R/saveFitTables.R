## The fit tables `saveFitTables` can write: parameter value -> the simList object
fitTableObjects <- c(ignition = "fireSense_ignitionCovariates", escape = "fireSense_escapeCovariates")

#' The events to schedule for `saveFitTables`
#'
#' @param which character, `P(sim)$saveFitTables`.
#' @return `"saveFitTables"` if any table is asked for, otherwise `character(0)`.
saveFitTablesEvents <- function(which) {
  if (length(which)) "saveFitTables" else character(0)
}

#' Write the ignition and/or escape fit tables to `outputPath(sim)`
#'
#' Runs after the prep events, as its own uncached event, so it sees the tables whether `prepIgnitionFitData`
#' just built them or its cache restored them. It only reads them: no table and no cache key changes. Each
#' file is `<object>_<ELF>.qs2` (`sim$.ELFind`; no ELF suffix without one) and is listed in `outputs(sim)`.
#' A requested table the run did not build (`skipUnusedFitPrep`) is reported, not an error.
#'
#' @param sim a `simList`.
#' @return the `simList`, invisibly, with the files registered as outputs.
saveFitTables <- function(sim) {
  which <- P(sim)$saveFitTables
  if (!all(which %in% names(fitTableObjects)))
    stop("saveFitTables must be among ", paste(shQuote(names(fitTableObjects)), collapse = ", "),
         "; got: ", paste(shQuote(setdiff(which, names(fitTableObjects))), collapse = ", "))
  elf <- if (length(sim$.ELFind)) paste0("_", paste(sim$.ELFind, collapse = "_"))
  for (obj in fitTableObjects[unique(which)]) {
    if (is.null(sim[[obj]])) {
      message("saveFitTables: ", obj, " was not built in this run; nothing written")
      next
    }
    f <- file.path(SpaDES.core::outputPath(sim), paste0(obj, elf, ".qs2"))
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    qs2::qs_save(sim[[obj]], f)
    sim <- SpaDES.core::registerOutputs(f, sim)
  }
  invisible(sim)
}
