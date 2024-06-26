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
#'                                                right_child=NULL))
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

