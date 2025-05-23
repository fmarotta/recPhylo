% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloXML.R
\name{traverse_clades}
\alias{traverse_clades}
\title{Recursively Traverse Clades in a Phylogenetic Tree}
\usage{
traverse_clades(clade, FUN, ..., .order = "pre")
}
\arguments{
\item{clade}{A list representing a clade of a phylogenetic tree. Each clade is expected to have a \code{clade} element which is a list of sub-clades.}

\item{FUN}{A function to be applied to each clade. It should accept the clade object as first argument and potentially additional arguments \code{...}.}

\item{...}{Additional arguments to be passed to \code{FUN}.}

\item{.order}{A character string indicating the order of traversal: "pre" for pre-order (default) and "post" for post-order.}
}
\value{
A list of results from applying \code{FUN} to each clade in the specified order.
}
\description{
This function traverses the clades of a phylogenetic tree, applying a specified function to each clade.
}
\note{
If FUN returns \code{NULL}, the value will not be included in the output. In other words,
}
\examples{
clade <- example_phyloXML_clade()
extract_clade_name <- function(clade) {
  clade$name
}
traverse_clades(clade, extract_clade_name)
traverse_clades(clade, extract_clade_name, .order = "post")

}
