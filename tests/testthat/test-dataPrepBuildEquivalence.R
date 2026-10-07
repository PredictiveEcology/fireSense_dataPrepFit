## The split of init into init + dataPrepBuild must change no result. This runs the module as it was
## before the split and as it is now, from the same captured inputs, each first with a cold cache
## and then with a warm one (a warm run loads dataPrepBuild from the cache), and compares every
## object and the module's mod$ contents: pre and post cold, and pre and post warm.
##
## Opt-in, because it needs a real study area and Google Drive access (init asks Drive for an
## existing fit). Set
##   DPF_EQUIVALENCE_FIXTURE    a fixture written by dpfFixtureFromSim() (helper-dataPrepBuild.R)
##                              from a simList stopped just before this module's init
## and optionally
##   DPF_EQUIVALENCE_INPUTPATH  inputPath for the runs (default: the captured run's, reusing downloads)
##   DPF_EQUIVALENCE_TMP        where the four runs' caches and outputs go (default: tempdir())
##   DPF_EQUIVALENCE_PRE_REF    git ref of the pre-split module (default below)

preSplitRef <- "58cfb1a4" # development just before the split

test_that("dataPrepBuild gives the pre-split results, with a cold and with a warm cache", {
  fixtureFile <- Sys.getenv("DPF_EQUIVALENCE_FIXTURE")
  skip_if(!nzchar(fixtureFile), "DPF_EQUIVALENCE_FIXTURE is not set")
  skip_if(!nzchar(Sys.which("git")), "git is not available")
  withr::local_options(spades.useCache = "eventsOnly", reproducible.useMemoise = FALSE)

  root <- withr::local_tempdir(tmpdir = Sys.getenv("DPF_EQUIVALENCE_TMP", tempdir()))
  ## SpaDES.core finds a module in a directory named for it
  preRoot <- file.path(root, "pre")
  dir.create(file.path(preRoot, dpfModule), recursive = TRUE)
  tarFile <- file.path(root, "pre.tar")
  preRef <- Sys.getenv("DPF_EQUIVALENCE_PRE_REF", preSplitRef)
  expect_identical(system2("git", c("-C", shQuote(moduleRoot), "archive", "-o", shQuote(tarFile), preRef)), 0L)
  utils::untar(tarFile, exdir = file.path(preRoot, dpfModule))
  postRoot <- file.path(root, "post")
  dir.create(postRoot)
  expect_true(file.symlink(moduleRoot, file.path(postRoot, dpfModule)))

  fixture <- qs2::qs_read(fixtureFile)
  inputPath <- Sys.getenv("DPF_EQUIVALENCE_INPUTPATH", fixture$inputPath)
  run <- function(modulePath, cacheName) {
    res <- dpfRunFromFixture(fixture, modulePath = modulePath,
                             cachePath = file.path(root, cacheName),
                             outputPath = file.path(root, paste0(cacheName, "_outputs")),
                             inputPath = inputPath)
    list(loadedBuild = any(grepl("Loaded! Cached result from previous doEvent.fireSense_dataPrepFit::dataPrepBuild",
                                 res$messages, fixed = TRUE)),
         result = dpfResult(res$sim))
  }

  preCold <- run(preRoot, "cache_pre")
  preWarm <- run(preRoot, "cache_pre")
  postCold <- run(postRoot, "cache_post")
  postWarm <- run(postRoot, "cache_post")

  expect_false(postCold$loadedBuild) # a cold cache builds
  expect_true(postWarm$loadedBuild)  # a warm one loads dataPrepBuild
  ## init now records in mod$ whether the user supplied the fuel objects; the pre-split module had no such entry
  intended <- "mod$userSuppliedFuelObjs"
  dpfExpectNone(setdiff(dpfDifferences(postCold$result, preCold$result), intended),
                "differ from pre-split, cold cache")
  dpfExpectNone(setdiff(dpfDifferences(postWarm$result, preWarm$result), intended),
                "differ from pre-split, warm cache")
})
