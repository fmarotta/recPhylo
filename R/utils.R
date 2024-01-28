#' Get the path to an example RecPhyloXML
#'
#' @param name The name of the example.
#'
#' @returns
#' If name is `NULL`, a vector of possible examples. Otherwise, the
#' path to the requested example.
#'
#' @export
recphylo_example <- function(name = NULL) {
  if (is.null(name)) {
    dir(system.file("extdata", package = "recPhyloParse"))
  } else {
    system.file("extdata", name, package = "recPhyloParse", mustWork = TRUE)
  }
}
