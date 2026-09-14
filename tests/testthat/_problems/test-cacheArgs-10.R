# Extracted from test-cacheArgs.R:10

# test -------------------------------------------------------------------------
files <- c(file.path(moduleRoot, paste0(moduleName, ".R")),
             list.files(file.path(moduleRoot, "R"), pattern = "\\.R$", full.names = TRUE))
code <- unlist(lapply(files, readLines, warn = FALSE))
