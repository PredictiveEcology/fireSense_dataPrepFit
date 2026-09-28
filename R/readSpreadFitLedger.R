#' Read the SpreadFit ledger rows for a study area
#'
#' With `spreadFitFilename = "latest"`, the rows come from `fireSenseUtils::latestSpreadFits()`: for
#' each polygon, the most recent ledger file that has it. `CacheGeo()` then picks this study area's
#' rows, as it does from a named file. The combined rows go to a local file named by their content
#' (`prepInputs()` would otherwise find its checksum changed and try to download it).
#'
#' @param spreadFitFilename The module's `spreadFitFilename`: a file in `cloudFolderID`, or `"latest"`.
#' @param cloudFolderID The Google Drive folder holding the ledger files.
#' @param domain The study area, an `sf` object.
#' @param destinationPath Where the ledger files are downloaded.
#' @return What `CacheGeo()` returns: the ledger rows within `domain`, or `NULL`.
readSpreadFitLedger <- function(spreadFitFilename, cloudFolderID, domain, destinationPath) {
  purge <- 7                                         # download the named file again: it changes on Drive
  if (identical(spreadFitFilename, "latest")) {
    latest <- fireSenseUtils::latestSpreadFits(cloudFolderID, destinationPath)
    if (is.null(latest))
      return(NULL)
    tf <- tempfile(fileext = ".rds")
    on.exit(unlink(tf), add = TRUE)
    saveRDS(latest, tf)
    spreadFitFilename <- paste0("fireSenseParams_latest_", substr(tools::md5sum(tf), 1, 12), ".rds")
    file.copy(tf, file.path(destinationPath, spreadFitFilename), overwrite = TRUE)
    cloudFolderID <- NULL
    purge <- FALSE                                   # local, with nothing to download it from
  }
  CacheGeo(cloudFolderID = cloudFolderID, targetFile = spreadFitFilename, purge = purge,
           domain = domain, action = "nothing", useCache = FALSE,
           destinationPath = destinationPath, bufferOK = TRUE)
}
