#' @export
recphylo_example <- function(path = NULL) {
  if (is.null(path)) {
    dir(system.file("extdata", package = "recPhyloParse"))
  } else {
    system.file("extdata", path, package = "recPhyloParse", mustWork = TRUE)
  }
}
