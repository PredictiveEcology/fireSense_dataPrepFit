#' The ledger row of the run's own ELF
#'
#' @param ownID character; the own ELF (`sim$.ELFind`; may be `NULL`).
#' @param rowIDs character; the ELF of each row, as `fireSenseUtils::polygonIDTxt` in the ledger.
#' @return the index of the own ELF's row, or `NA` when `ownID` is `NULL` or has no row.
#' @noRd
ownLedgerRow <- function(ownID, rowIDs) {
  if (is.null(ownID) || !length(ownID)) return(NA_integer_)
  match(as.character(ownID)[1], rowIDs)
}

#' One colour per species from the SpreadFit ledger rows
#'
#' Each ledger row (the own ELF and the neighbours whose polygons intersect the study area) carries a
#' species colour vector, and the same species can have a different colour in different rows. Colours
#' only drive plots, so the own ELF's row wins; otherwise the first row that has the species does.
#' "Mixed" is last.
#'
#' @param sppColorVects list of named colour vectors, one per ledger row.
#' @param ownID character; the own ELF (`sim$.ELFind`; may be `NULL`).
#' @param rowIDs character; the ELF of each row, as `fireSenseUtils::polygonIDTxt` in the ledger.
#' @return a named colour vector, one entry per species.
#' @noRd
mergeSppColorVects <- function(sppColorVects, ownID, rowIDs) {
  ownIndex <- ownLedgerRow(ownID, rowIDs)
  ord <- if (is.na(ownIndex)) seq_along(sppColorVects) else c(ownIndex, setdiff(seq_along(sppColorVects), ownIndex))
  scv <- unlist(unname(sppColorVects[ord]))
  differ <- names(which(tapply(scv, names(scv), function(x) length(unique(x))) > 1L))
  if (length(differ))
    message("fireSense_dataPrepFit: species colours differ among ledger rows; keeping the own ELF's (or the first row's) for: ",
            paste(differ, collapse = ", "))
  scv <- scv[!duplicated(names(scv))]
  scv <- scv[order(names(scv) == "Mixed", names(scv))]
  scv
}
