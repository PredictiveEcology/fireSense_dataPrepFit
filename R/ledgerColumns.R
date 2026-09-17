## The ledger columns this module reads, keeping only those the ledger has: rows written before a column
## was added (e.g. covMinMax_spread, 2026-09-17) lack it, and `x[colNames]` would error.
ledgerColumns <- function(ledger, colNames) {
  ledger[intersect(colNames, names(ledger))]
}
