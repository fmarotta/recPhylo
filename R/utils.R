#' Get examples for recPhyloXML or phyloXML objects
#'
#' @returns An example of the requested type.
#'
#' @export
#' @name phyloXML_examples
example_recPhyloXML_file <- function(name = NULL) {
  system.file("extdata", "recphyloxml_syntax_example_1.recphyloxml", package = "recPhylo", mustWork = TRUE)
}

#' @export
#' @rdname phyloXML_examples
example_phyloXML_file <- function(name = NULL) {
  system.file("extdata", "phyloxml_syntax_example_1.phyloxml", package = "recPhylo", mustWork = TRUE)
}

#' @export
#' @rdname phyloXML_examples
example_recPhyloXML <- function() {
  read_recPhyloXML(example_recPhyloXML_file())
}

#' @export
#' @rdname phyloXML_examples
example_phyloXML <- function() {
  read_phyloXML(example_phyloXML_file())
}

#' @export
#' @rdname phyloXML_examples
example_phyloXML_phylogeny <- function() {
  phyloXML <- example_phyloXML()
  phyloXML[[1]]
}

#' @export
#' @rdname phyloXML_examples
example_phyloXML_clade <- function() {
  phylo <- example_phyloXML_phylogeny()
  phylo$clade
}

#' Convert phyloXML layout to data frame
#'
#' This function converts a phyloXML layout (a list structure) into a data frame.
#'
#' @param cl A list representing the phyloXML layout. The list sholud have multiple nested lists under the element named 'children'.
#'
#' @return A data frame containing the concatenated elements of the phyloXML layout.
#'
#' @examples
#' cl <- list(name = "root", value = 1, children = list(
#'               list(name = "child1", value = 2, children = list()),
#'               list(name = "child2", value = 3, children = list())))
#' as.data.frame.phyloXML_layout(cl)
#'
#' @export
as.data.frame.phyloXML_layout <- function(cl) {
  fields <- setdiff(names(cl), "children")
  rbind(
    as.data.frame(cl[fields]),
    Reduce(rbind, lapply(cl$children, as.data.frame))
  )
}

#' Merge Phylogeny Data with Layout Data
#'
#' This function merges a `phyloXML_phylogeny` list with a layout data frame based on a common column 'name'.
#' Additionally, it merges annotation data if present in the phylogeny.
#'
#' @param phylogeny A list of class `phyloXML_phylogeny, potentially annotated with an `annot` attribute.
#' @param layout A data frame representing the layout data.
#'
#' @return A data frame resulting from the merge of the phylogeny and layout data frames, with annotation data if available.
#'
#' @examples
#' phylogeny <- example_phyloXML_phylogeny()
#' layout <- data.frame(name = c("unnamed", "Octopus", "Bacteria", "E. coli", "B. subtilis"), position = c(100, 100, 200, 150, 250))
#' attr(phylogeny, "annot") <- data.frame(name = c("Octopus"), annotation = c("annot1"))
#' merge_layout(phylogeny, layout)
merge_layout <- function(phylogeny, layout) {
  # NOTE: phyolgeny lists are automatically casted to data.frame
  res <- merge(
    phylogeny,
    layout,
    by = "name"
  )
  annot <- attr(phylogeny, "annot", exact = TRUE)
  if (!is.null(annot)) {
    res <- merge(res, annot, all.x = TRUE)
  }
  res
}

# XXX: I'm trying to upstream this to the ape package! https://github.com/emmanuelparadis/ape/pull/120
#' Save a "phylo" object to a file in phyloXML format
#'
#' @param tree A phylogenetic tree (e.g. from ape::read.tree()).
#' @param path Path to the file where the tree will be written.
#'
#' @export
write_phyloXML <- function(tree, path) {
  phyloxml <- xml2::xml_new_root("phyloxml",
    xmlns = "http://www.phyloxml.org",
    `xmlns:xsi` = "http://www.w3.org/2001/XMLSchema-instance",
    `xsi:schemaLocation` = "http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd"
  )
  phylogeny <- xml2::xml_add_child(phyloxml, "phylogeny", rooted = "false")
  xml2::xml_add_child(phylogeny, "name", "my_tree")
  xml2::xml_add_child(phylogeny, "description", "beautiful")
  root_idx <- unique(tree$edge[! tree$edge[, 1] %in% tree$edge[, 2], 1])
  stopifnot(length(root_idx) == 1)
  clades <- phylo_to_xml_clades(root_idx, tree)
  if (!is.null(tree$root.edge)) {
    xml2::xml_set_attr(clades, "branch_length", tree$root.edge)
  }
  xml2::xml_add_child(phylogeny, clades)
  xml2::write_xml(phyloxml, path)
}

phylo_to_xml_clades <- function(parent_idx, tree) {
  n_tips <- tree$Nnodes + 1 + is.rooted(tree)
  parent <- read_xml("<clade></clade>")
  if (!is.null(tree$tip.label) && parent_idx <= length(tree$tip.label)) {
    xml2::xml_add_child(parent, "name", tree$tip.label[parent_idx])
  } else if (!is.null(tree$node.label) && parent_idx > n_tips) {
    xml2::xml_add_child(parent, "name", tree$node.label[parent_idx - n_tips])
  } else {
    xml2::xml_add_child(parent, "name", paste("node", parent_idx, sep = "_"))
  }
  which_children <- which(tree$edge[, 1] == parent_idx)
  if (length(which_children) == 0) {
    return(parent)
  }
  lapply(which_children, function(which_child) {
    child_idx <- tree$edge[which_child, 2]
    child <- phylo_to_xml_clades(child_idx, tree)
    if (!is.null(tree$edge.length)) {
      branch_length <- tree$edge.length[which_child]
      xml2::xml_set_attr(child, "branch_length", branch_length)
    }
    xml2::xml_add_child(parent, child)
  })
  return(parent)
}
