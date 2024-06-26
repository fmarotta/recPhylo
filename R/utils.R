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


#' Find the Most Recent Common Ancestor (MRCA) in a Phylogenetic List
#'
#' This function searches through a binary phylogenetic tree represented
#' as a list for the most recent common ancestor (MRCA) of a specified
#' set of nodes. The tree is expected to be provided as a recursive list
#' structure with `left_child` and `right_child` elements defining the
#' bifurcations, and `name` elements defining the node names. Nodes in
#' the tree should be uniquely named.
#'
#' @param phylo_list A list representing the phylogenetic tree. This
#'   list should have elements `name` for the node names, and
#'   `left_child` and `right_child` for the branches, recursively
#' @param nodes A character vector with the names of the nodes for which
#'   the MRCA is sought
#'
#' @return Returns the name of the MRCA node if found; otherwise, `NA`.
#' Additionally, if `NA` is returned, an attribute "which_nodes" is
#' added to the result. This attribute is a logical vector indicating
#' which of the nodes specified in the `nodes` argument are present in
#' the subtree where the search ended.
#'
#' @examples
#' # Construct a simple phylogenetic tree
#' tree <- list(name="root",
#'              left_child=list(name="A",
#'                              left_child=NULL,
#'                              right_child=NULL),
#'              right_child=list(name="BC",
#'                               left_child=list(name="B",
#'                                               left_child=NULL,
#'                                               right_child=NULL),
#'                               right_child=list(name="C",
#'                                                left_child=NULL,
#'                                                right_child=NULL)))
#' # Find the MRCA of nodes "A" and "B"
#' mrca(tree, c("A", "B"))
#'
#' @export
mrca <- function(phylo_list, nodes) {
  if (is.null(phylo_list)) {
    ret <- NA
    attr(ret, "which_nodes") <- rep(FALSE, length(nodes))
    return(ret)
  }
  mrca_left <- mrca(phylo_list$left_child, nodes)
  mrca_right <- mrca(phylo_list$right_child, nodes)
  # if there are no duplicated node names, _left and _right can't be both true at the same time.
  if (!is.na(mrca_left)) {
    return(mrca_left)
  } else if (!is.na(mrca_right)) {
    return(mrca_right)
  }
  set <- nodes %in% phylo_list$name | attr(mrca_left, "which_nodes") | attr(mrca_right, "which_nodes")
  if (all(set)) {
    return(phylo_list$name)
  } else {
    ret <- NA
    attr(ret, "which_nodes") <- set
    return(ret)
  }
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

nextperm <- function(vec) {
  # TODO: correct implementation!
  rev(vec)
}
