#' Select the ledger columns that exist
#'
#' A ledger written before a column was added (e.g. `covMinMax_spread`) lacks it, and
#' `ledger[colNames]` would error.
#'
#' @param ledger data.frame or `sf`, the SpreadFit ledger (`sim$spreadFitPreRun`).
#' @param colNames character, the columns wanted.
#' @return `ledger` with only the columns in `colNames` that it has (an `sf` keeps its geometry).
ledgerColumns <- function(ledger, colNames) {
  ledger[intersect(colNames, names(ledger))]
}
