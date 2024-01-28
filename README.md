# recPhyloParse

<!-- badges: start -->
  [![R-CMD-check](https://github.com/fmarotta/recPhyloParse/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fmarotta/recPhyloParse/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

recPhyloParse provides a way to import reconciliated phylogenies into R in a way that makes it easy to plot the trees with ggplot2. The starting point is an XML file following the [recPhyloXML format](http://phylariane.univ-lyon1.fr/recphyloxml/), which will be parsed into an [R6 object](https://r6.r-lib.org/index.html). During the parsing, the coordinates and attributes of all the plot elements will be appropriately extracted. The package provides methods for converting the trees into data.frames, so that standard ggplot2 geoms can be used.

## Installation

You can install the development version of recPhyloParse from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fmarotta/recPhyloParse")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(recPhyloParse)

recphylo <- RecPhylo$new("my_reconciliated_tree.xml")

spTree <- as.data.frame(recphylo$spTree)
recGeneTree <- as.data.frame(recphylo$recGeneTree)
spEdges <- get_spedges(recphylo$spTree)
recGeneEdges <- get_gedges(recphylo$recGeneTree)

ggplot() +
  geom_line(data = spEdges, aes(x, y, group = group)) +
  geom_line(data = recGeneEdges, aes(x, y, group = group))
```

## Bugs and missing features

* We support only one gene tree per species tree.
* BifurcationOut events are currently not supported.
* We don't have the ability to (easily) flip the children of a node (both gene and species tree).
* We make no attempt to minimise the edge crossings for lateral transfers.
* We do not add clade attributes and internal elements as columns in the resulting data.frame (we could add all these: http://www.phyloxml.org/documentation/version_1.20/phyloxml.html)
