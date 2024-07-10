
<!-- README.md is generated from README.Rmd. Please edit that file -->

# recPhylo

<!-- badges: start -->

[![R-CMD-check](https://github.com/fmarotta/recPhylo/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fmarotta/recPhylo/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/fmarotta/recPhylo/branch/main/graph/badge.svg)](https://app.codecov.io/gh/fmarotta/recPhylo?branch=main)
<!-- badges: end -->

`recPhylo` provides a way to import simple or reconciliated phylogenies
into R in a way that makes it easy to plot the trees with `ggplot2`. The
starting point is an XML file following either the
[phyloXML](http://www.phyloxml.org/) or the [recPhyloXML
format](http://phylariane.univ-lyon1.fr/recphyloxml/), which will be
parsed into specific R objects. Dedicated `layout_*()` functions can
convert these objects into data.frames ready to be used with ggplot.

Contrary to `ggtree` or `ggraph`, with this package you don’t get a plot
directly from the tree object, but rather a pair of data.frames
representing the computed layouts of nodes and edges. It is your
responsibility, then, to use the appropriate ggplot geoms (typically
geom_point() or geom_line()). The advantages of this approach are:

- We don’t need to duplicate existing geoms (e.g. we don’t need to wrap
  geom_point() inside geom_node_point(), we can simply reuse
  geom_point())
- We can make use of special geoms defined in other packages
  (e.g. geom_link2() from ggforce, without having to define a special
  geom for phylogeny objects)
- We have better control over the layout and annotation: the x and y
  coordinates of all the objects are transparently saved in the layout
  data.frames, so it’s straightforward to add custom ggplot elements
  corresponding to those coordinates (e.g. a tree and a multiple
  sequence alignment).
- We can make use of other `ggplot2` tricks, like coord_polar(), to
  create new layouts “for free”.
- Should you want to use a graphics package different from `ggplot2`,
  you can reuse the layout computed by recPhylo.

This package does not have as many features as `ggtree`. Use it if you
are comfortable with `ggplot2` and like to have full control over your
figures.

## Installation

You can install the development version of recPhylo from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fmarotta/recPhylo")
```

## Example

### recPhylo objects

This is a basic example which shows you how to plot a reconciliated
phylogenetic tree.

``` r

library(recPhylo)
library(ggplot2)

recphylo_file <- system.file("extdata", "recphylo_example_1.xml", package = "recPhylo")
recphylo_xml <- read_recPhyloXML(recphylo_file)
#> Warning: Phylogeny doesn't have a 'rooted' attribute. Inferred 'TRUE'.
#> Warning: Phylogeny doesn't have a 'rooted' attribute. Inferred 'TRUE'.
recphylo_layout <- layout_recphylo(recphylo_xml, branch_length_scale = 5)
#> Warning: No branch length was found in clade ABCD. Setting it to 1
#> automatically.

ggplot() +
  geom_line(data = recphylo_layout$spEdges, aes(x, y, group = group)) +
  geom_line(data = recphylo_layout$recGeneEdges, aes(x, y, group = group, color = lineage, linetype = leg_type), show.legend = F) +
  scale_linetype_manual(values = c("loss_vertical" = 2, "transferBack" = 3), na.value = 1) +
  theme_void() +
  NULL
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

### Simple phylogeny objects

Although the primary goal of the package was to plot reconciliated
trees, for which there was no alternative in R, it is also possible to
plot simple phylogenies.

``` r

library(recPhylo)
library(ggplot2)

phylo_file <- system.file("extdata", "phylo_example_1.xml", package = "recPhylo")
phylo_xml <- read_phyloXML(phylo_file)[[1]]
phylo_xml <- add_annot(phylo_xml, data.frame(
  name = find_descendants(phylo_xml$clade, "Bacteria"),
  is_bacterium = TRUE
))
phylo_layout <- layout_phylogeny(phylo_xml, branch_length_scale = 5)
#> Warning: No branch length was found in clade 'unnamed'. Setting it to 1
#> automatically.

ggplot() +
  geom_line(data = phylo_layout$edges, aes(x, y, group = group, color = is_bacterium)) +
  theme_void() +
  NULL
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

## Features

- Parse and draw phyloXML or recPhyloXML trees with ggplot2
- Three linear layouts: real branch lengths, uniform branch lengths,
  aligned leaves
- Corresponding circular layouts “for free” thanks to coord_polar()
- Manually adjust the space between clades
- Swap left and right children of a species node
- Swap children of a duplication event
- Import branch lengths from another species tree file
- Add custom annotation for species nodes, gene nodes, species edges,
  gene edges
- Convenience methods to extract data from phylogenies or decorate them
  with additional attributes

## Bugs and missing features

- bifurcationOut events are not supported (but planned)
- We don’t automatically scale the plot to make it look nicer (planned)
- We make no attempt to minimise the edge crossings for lateral
  transfers (though you can manually swap left/right children to make
  the tree look better)
