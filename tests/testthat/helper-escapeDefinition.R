## fireSenseUtils::escapeDefinitionOf(sim) reads the parameters of a simList. A test that gives a prepare function a
## plain environment for `sim` gives it the definition of `Par` (its `escapeSizeHa` and `fireCauses`) instead.
mockEscapeDefinition <- function(Par, .env = parent.frame()) {
  testthat::local_mocked_bindings(
    escapeDefinitionOf = function(sim) fireSenseUtils::escapeDefinition(Par$escapeSizeHa, Par$fireCauses),
    .package = "fireSenseUtils", .env = .env)
}
