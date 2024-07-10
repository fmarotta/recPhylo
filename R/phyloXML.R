# NOTE: gene and species names will be made unique
# NOTE: phylogenies will be given a unique idx field (their index)

#' Read a phyloXML file
#'
#' @description Read a tree from a file in phyloXML format and parse it
#' into a list of phylogeny objects.
#'
#' @param x A string, a connection, or a raw vector (see xml2::read_xml).
#' @param ... Additional arguments passed on to xml2::read_xml.
#'
#' @return A list of the phylogenies that were parsed from the file.
#'
#' @section Phylogeny objects:
#' Each phylogeny is a list of class `phyloXML_phylogeny` with one field
#' for each attribute or child of the `<phylogeny>` tag. The `clade` field
#' is a nested list structure of class `phyloXML_clade`, with one field
#' for each attribute of child of the `<clade>` tag.
#'
#' @note The `<name>` of the nodes in the tree will be made unique, if necessary, by adding
#' suffixes of the form `#1`, `#2`, and so on. Unnamed clades will be given
#' "unnamed" as their name. This is because the functions in this package
#' heavily rely on the clades being uniquely named.
#'
#' @examples
#' path_to_xml <- example_phyloXML_file()
#' read_phyloXML(path_to_xml)
#'
#' @export
read_phyloXML <- function(x, ...) {
  xml <- xml2::read_xml(x, ...)
  parse_phyloXML(xml)
}

#' Read a recPhyloXML file
#'
#' @inherit read_phyloXML
#'
#' @description Read a reconciliated tree from a file in recPhyloXML format and parse it
#' into a list of spTree and recGeneTrees objects.
#'
#' @return A list with two elements: spTree and recGeneTrees. spTree is a simple phylogeny, recGeneTrees is a list of phylogenies.
#'
#' @examples
#' path_to_xml <- example_phyloXML_file()
#' read_phyloXML(path_to_xml)
#'
#' @export
read_recPhyloXML <- function(x, ...) {
  xml <- xml2::read_xml(x, ...)
  parse_recPhyloXML(xml)
}

parse_phyloXML <- function(xml) {
  xml2::xml_ns_strip(xml)
  make_unique_names(xml)
  phylogenies <- xml2::xml_find_all(xml, "./phylogeny")
  l <- lapply(seq_along(phylogenies), function(phy_idx) {
    phy <- parse_phylogeny(phylogenies[[phy_idx]], type = "phyloXML")
    phy$idx <- phy_idx
    phy
  })
  class(l) <- c("phyloXML", class(l))
  l
}

parse_recPhyloXML <- function(xml) {
  xml2::xml_ns_strip(xml)
  make_unique_names(xml2::xml_find_first(xml, "./spTree"))
  make_unique_recGeneTree_names(xml2::xml_find_first(xml, "./recGeneTree"))
  spTree_xml <- xml2::xml_find_all(xml, "./spTree/phylogeny")
  recGeneTree_xml <- xml2::xml_find_all(xml, "./recGeneTree/phylogeny")
  stopifnot(length(spTree_xml) <= 1)
  stopifnot(length(recGeneTree_xml) >= 1)
  recGeneTrees <- lapply(seq_along(recGeneTree_xml), function(phy_idx) {
    phy <- parse_phylogeny(recGeneTree_xml[[phy_idx]], "recGeneTree")
    phy$idx <- phy_idx
    phy
  })
  stopifnot(all(sapply(recGeneTrees, `[[`, "rooted")))
  # TODO: give names to recGeneTrees
  spTree <- if (length(spTree_xml) == 0) {
    # TODO: infer species tree. (Is it even possible?)
    stop("Missing spTree.")
  } else if (length(spTree_xml) == 1) {
    parse_phylogeny(spTree_xml[[1]], "spTree")
  } else {
    stop("Missing spTree.")
  }
  l <- list(
    spTree = spTree,
    recGeneTrees = recGeneTrees
  )
  class(l) <- c("recPhyloXML", class(l))
  l
}

parse_phylogeny <- function(phylogeny_xml, type) {
  clade <- if (type == "phyloXML") {
    parse_phyloXML_clade(xml2::xml_find_first(phylogeny_xml, "./clade"))
  } else if (type == "spTree") {
    parse_spTree_clade(xml2::xml_find_first(phylogeny_xml, "./clade"))
  } else if (type == "recGeneTree") {
    parse_recGeneTree_clade(xml2::xml_find_first(phylogeny_xml, "./clade"))
  } else {
    stop("Unsupported phylogeny type.")
  }
  rooted <- as.logical(xml2::xml_attr(phylogeny_xml, "rooted"))
  if (is.na(rooted)) {
    rooted <- length(clade$clade) == 2
    warning("Phylogeny doesn't have a 'rooted' attribute. Inferred '", rooted, "'.", call. = FALSE)
  }
  phylogeny <- list(
    name = xml2::xml_text(xml2::xml_find_first(phylogeny_xml, "./name")),
    id = xml2::xml_text(xml2::xml_find_first(phylogeny_xml, "./id")),
    description = xml2::xml_text(xml2::xml_find_first(phylogeny_xml, "./description")),
    date = as.Date(xml2::xml_text(xml2::xml_find_first(phylogeny_xml, "./date"))),
    rooted = rooted,
    rerootable = as.logical(xml2::xml_attr(phylogeny_xml, "rerootable")),
    branch_length_unit = xml2::xml_attr(phylogeny_xml, "branch_length_unit"),
    type = xml2::xml_attr(phylogeny_xml, "type"),
    # confidence = xml2::xml_find_all(phylogeny_xml, "./confidence"),
    # clade_relation = xml2::xml_find_all(phylogeny_xml, "./clade_relation"),
    # sequence_relation = xml2::xml_find_all(phylogeny_xml, "./sequence_relation"),
    # property = xml2::xml_find_all(phylogeny_xml, "./property"),
    clade = clade
  )
  class(phylogeny) <- c("phyloXML_phylogeny", class(phylogeny))
  phylogeny
}

parse_clade <- function(clade_xml) {
  name <- xml2::xml_text(xml2::xml_find_first(clade_xml, "./name"))
  branch_length_attr <- as.numeric(xml2::xml_attr(clade_xml, "branch_length"))
  branch_length_tag <- xml2::xml_double(xml2::xml_find_first(clade_xml, "branch_length"))
  branch_length <- if (!is.na(branch_length_attr) && !is.na(branch_length_tag)) {
    warning("clade '", name, "' defines branch length both with an attribute and with a tag, which is ambiguous. I chose the attribute value, ", branch_length_attr, ".", call. = FALSE)
    branch_length_attr
  } else if (!is.na(branch_length_attr)) {
    branch_length_attr
  } else if (!is.na(branch_length_tag)) {
    branch_length_tag
  } else {
    NA_real_
  }
  clade <- list(
    branch_length = branch_length,
    name = name,
    # id_source = xml2::xml_double(xml2::xml_find_first(clade_xml, "./id_source")),
    collapse = as.logical(xml2::xml_attr(clade_xml, "collapse")),
    # confidence = xml2::xml_text(xml2::xml_find_all(clade_xml, "./confidence")),
    width = xml2::xml_double(xml2::xml_find_first(clade_xml, "./width")),
    # color = xml2::xml_text(xml2::xml_find_first(clade_xml, "./color")),
    # taxonomy = xml2::xml_text(xml2::xml_find_all(clade_xml, "./taxonomy")),
    # sequence = xml2::xml_text(xml2::xml_find_all(clade_xml, "./sequence")),
    # events = xml2::xml_text(xml2::xml_find_first(clade_xml, "./events")),
    # binary_characters = xml2::xml_text(xml2::xml_find_first(clade_xml, "./binary_characters")),
    # distribution = xml2::xml_text(xml2::xml_find_all(clade_xml, "./distribution")),
    date = as.Date(xml2::xml_text(xml2::xml_find_first(clade_xml, "./date")))
    # reference = xml2::xml_text(xml2::xml_find_all(clade_xml, "./reference")),
    # property = xml2::xml_text(xml2::xml_find_all(clade_xml, "./property")),
  )
  class(clade) <- c("phyloXML_clade", class(clade))
  clade
}

parse_phyloXML_clade <- function(clade_xml) {
  clade <- parse_clade(clade_xml)
  clade$clade <- lapply(xml2::xml_find_all(clade_xml, "./clade"), parse_phyloXML_clade)
  return(clade)
}

parse_spTree_clade <- function(clade_xml) {
  clade <- parse_clade(clade_xml)
  clade$clade <- lapply(xml2::xml_find_all(clade_xml, "./clade"), parse_spTree_clade)
  # clade$geography <- ... TODO
  return(clade)
}

parse_recGeneTree_clade <- function(clade_xml) {
  clade <- parse_clade(clade_xml)
  event_xml <- xml2::xml_find_first(clade_xml, "./eventsRec/*[self::leaf or self::duplication or self::loss or self::branchingOut or self::speciation]")
  clade$event_type <- xml2::xml_name(event_xml)
  clade$event_location <- xml2::xml_attr(event_xml, "speciesLocation")
  # clade$geography <- ... TODO
  clade$clade <- lapply(xml2::xml_find_all(clade_xml, "./clade"), parse_recGeneTree_clade)
  if (length(xml2::xml_find_first(clade_xml, "./eventsRec/*[self::transferBack]")) > 0) {
    parent <- clade
    parent$name <- paste(parent$name, "transferBack", sep = "@")
    parent$event_type <- "transferBack"
    parent$clade <- list(clade)
    return(parent)
  }
  return(clade)
}

# S3 interop

#' @export
as.data.frame.phyloXML <- function(phylo_xml) {
  Reduce(rbind, lapply(phylo_xml, as.data.frame))
}

#' @export
as.data.frame.phyloXML_phylogeny <- function(phy) {
  phy_fields <- setdiff(names(phy), "clade")
  phy_values <- phy[phy_fields]
  names(phy_values) <- paste("phylogeny", phy_fields, sep = "_")
  # TODO: consider removing fields that are NA in the phylogeny
  cbind(as.data.frame(phy$clade), phy_values)
}

#' @export
as.data.frame.phyloXML_clade <- function(cl) {
  fields <- setdiff(names(cl), "clade")
  Reduce(rbind, traverse_clades(cl, function(x) list(as.data.frame(x[fields]))))
}

#' @export
as.data.frame.recPhyloXML <- function(recphylo_xml) {
  sp_df <- as.data.frame(recphylo_xml$spTree)
  sp_df$tree_type <- "spTree"
  recGene_df <- as.data.frame(recphylo_xml$recGeneTrees)
  recGene_df$tree_type <- "recGeneTree"
  rbind(sp_df, recGene_df)
}

#' @export
print.phyloXML <- function(phylo_xml) {
  n_phylogenies <- length(phylo_xml)
  cat("phyloXML object with", n_phylogenies, if (n_phylogenies == 1) "phylogeny.\n" else "phylogenies.\n")
  lapply(seq_along(phylo_xml), function(idx) print(phylo_xml[[idx]], idx = idx))
  invisible(NULL)
}

#' @export
print.phyloXML_phylogeny <- function(phy, idx = NULL) {
  fields <- setdiff(names(phy), c("clade", "idx"))
  leaf_status <- traverse_clades(phy$clade, function(cl) {
    length(cl$clade) == 0
  })
  if (!is.null(idx)) {
    cat("[[", idx, "]] ", sep = "")
  }
  cat("phyloXML_phylogeny object with ", length(leaf_status), " nodes (", sum(leaf_status), " leaves).\n", sep = "")
  str(phy[fields])
  annot <- attr(phy, "annot", exact = TRUE)
  cat("Attributes:\n")
  if (!is.null(annot)) {
    str(annot)
  } else {
    cat("NULL\n")
  }
  invisible(NULL)
}

#' @export
print.phyloXML_clade <- function(cl, indent = "") {
  fields <- setdiff(names(cl), c("clade"))
  cat(indent, cl$name, "\n", sep = "")
  lapply(cl$clade, print, indent = paste0(indent, "  "))
  invisible(NULL)
}

#' @export
print.recPhyloXML <- function(recphylo_xml) {
  n_gene_trees <- length(recphylo_xml$recGeneTrees)
  cat("recPhyloXML object with 1 spTree and", n_gene_trees, if (n_gene_trees == 1) "recGeneTree.\n" else "recGeneTrees.\n")
  cat("spTree:\n")
  print(recphylo_xml$spTree)
  cat("recGeneTrees:\n")
  lapply(seq_along(recphylo_xml$recGeneTrees), function(idx) print(recphylo_xml$recGeneTrees[[idx]], idx = idx))
  invisible(NULL)
}

# Public utils

#' Decorate Clades with a Specific Attribute
#'
#' This function adds or modify an attribute of each node. The attribute is computed by a user-defined function which takes the node as input. decorate_clades() returns a new object of `phyloXML_clade` class.
#'
#' @param clade A list representing a phylogenetic clade. Each node in the clade is expected to have a `clade` element containing its child nodes.
#' @param attr A character string representing the name of the attribute to be added to each node.
#' @param FUN A function that computes the value of the attribute for each node. The function should accept the clade as its first argument and any additional parameters through `...`.
#' @param ... Additional parameters to be passed to `FUN`.
#' 
#' @return The input `clade` list, with each node augmented by the computed attribute.
#'
#' @examples
#' # Create a sample clade
#' clade <- example_phyloXML_clade()
#' # Define a function to compute an attribute
#' compute_attr <- function(clade, factor = 1) {
#'   length(clade$clade) * factor
#' }
#' # Decorate the clade with the new attribute
#' decorated_clade <- decorate_clades(clade, "node_count", compute_attr, factor = 2)
#'
#' @export
decorate_clades <- function(clade, attr, FUN, ...) {
  clade[[attr]] <- FUN(clade, ...)
  clade$clade <- lapply(clade$clade, decorate_clades, attr, FUN, ...)
  clade
}

#' Recursively Traverse Clades in a Phylogenetic Tree
#'
#' @description This function traverses the clades of a phylogenetic tree, applying a specified function to each clade.
#'
#' @param clade A list representing a clade of a phylogenetic tree. Each clade is expected to have a `clade` element which is a list of sub-clades.
#' @param FUN A function to be applied to each clade. It should accept the clade object as first argument and potentially additional arguments `...`.
#' @param ... Additional arguments to be passed to `FUN`.
#' @param .order A character string indicating the order of traversal: "pre" for pre-order (default) and "post" for post-order.
#'
#' @return A list of results from applying `FUN` to each clade in the specified order.
#'
#' @note If FUN returns `NULL`, the value will not be included in the output. In other words, 
#'
#' @examples
#' clade <- example_phyloXML_clade()
#' extract_clade_name <- function(clade) {
#'   clade$name
#' }
#' traverse_clades(clade, extract_clade_name)
#' traverse_clades(clade, extract_clade_name, .order = "post")
#'
#' @export
traverse_clades <- function(clade, FUN, ..., .order = "pre") {
  item <- FUN(clade, ...)
  children <- lapply(clade$clade, traverse_clades, FUN, .order = .order, ...)
  if (.order == "pre") {
    c(item, unlist(children, recursive = FALSE))
  } else if (.order == "post") {
    c(unlist(children, recursive = FALSE), item)
  } else {
    stop("Unsupported `.order`. Please choose either 'pre' or 'post'.")
  }
}

#' Find the Most Recent Common Ancestor (MRCA)
#'
#' This function identifies the MRCA of a set of nodes. It traverses through the clades and returns the MRCA, if possible.
#'
#' @param clade A list of class `phyloXML_clade` representing a
#' phylogenetic tree.
#' @param node_names A character vector containing the names of the
#' nodes for which the MRCA is to be found.
#'
#' @return The function returns the name of the MRCA of the specified
#' nodes, if found. If the MRCA is not found, the function returns `NA`
#' with an attribute "which_node_names" indicating the nodes that were
#' not found in the tree.
#'
#' @examples
#' clade <- example_phyloXML_clade()
#' find_mrca(clade, c("B. subtilis", "E. coli"))
#' find_mrca(clade, c("B. subtilis", "Bogus"))
#'
#' @export
find_mrca <- function(clade, node_names) {
  mrca_clade <- find_mrca_clade(clade, node_names)
  if (length(mrca_clade) == 1 && is.na(mrca_clade)) {
    return(mrca_clade)
  }
  mrca_clade$name
}

find_mrca_clade <- function(clade, node_names) {
  mrcas <- lapply(clade$clade, find_mrca_clade, node_names)
  if (sum(!is.na(mrcas)) == 1) {
    return(mrcas[[which(!is.na(mrcas))]])
  }
  set <- node_names %in% clade$name
  if (length(mrcas) > 0) {
   set <- set | Reduce(`|`, lapply(mrcas, attr, "which_node_names"))
  }
  if (all(set)) {
    return(clade)
  } else {
    ret <- NA
    attr(ret, "which_node_names") <- set
    return(ret)
  }
}

#' @export
find_descendants <- function(clade, ancestor_name) {
  dec_clade <- decorate_clades(clade, "clade", function(cl) {
    if (cl$name == ancestor_name || !is.null(cl$descendant)) {
      lapply(cl$clade, function(child) {
        child$descendant <- TRUE
        child
      })
    } else {
      cl$clade
    }
  })
  traverse_clades(dec_clade, function(cl) {
    if (cl$name == ancestor_name || !is.null(cl$descendant)) {
      return(cl$name)
    }
  })
}

#' @export
select_phylogeny <- function(phyloxml, idx) {
  # To be used with pipes, e.g. phylo |> select_phylogeny(1) |> import_branch_lengths(ape) |> bla...
  phyloxml[[idx]]
}

# Public branch length utils

# XXX: importing branch lengths makes sense only for individual phylogenies, not for phyloXML objects (which are lists of phylogenies). Why would all the phylogenies in the list have the same branch lengths? Anyway, one can easily use lapply in that case.
# But maybe we should provide a helper function that returns the whole phyloXML, so that one doesn't have to repeat the indexing... e.g. we could do
# import_branch_lengths <- function(phyloXML, phylo, which.phylogeny) {
#   phyloXML[[which.phylogeny]] <- import_branch_lengths(phyloXML[[which.phylogeny]], phylo)
#   return(phyloXML)
# }

#' Import Branch Lengths from a `phylo` Object
#'
#' @description This function takes branch lengths from an `ape::phylo`
#' object and imports them into the corresponding nodes of a
#' `phyloXML_phylogeny` or `recPhyloXML` object. For a `recPhyloXML`,
#' the branch lengths are imported only to the species tree.
#'
#' @param x Either a `phyloXML_phylogeny` or a `recPhyloXML` object to which the branch lengths are to be imported.
#' @param phylo A `phylo` object containing the branch lengths.
#' 
#' @return The modified phylogeny with the branch lengths from the `phylo` object.
#'
#' @details The `edge.length` and `node.label` fields of the `phylo` object must exist.
#'
#' @examples
#' \dontrun{
#' phyloxml <- example_phyloXML()
#' ape_tree <- ape::read.tree("path/to/tree.nwk")
#' phyloxml[[1]] <- import_branch_lengths(phyloxml[[1]], ape_tree)
#' }
#'
#' \dontrun{
#' recphyloxml <- example_recPhyloXML()
#' ape_tree <- ape::read.tree("path/to/tree.nwk")
#' recphyloxml <- import_branch_lengths(recphyloxml, ape_tree)
#' }
#'
#' @export
import_branch_lengths <- function(x, phylo) {
  UseMethod("import_branch_lengths")
}

#' @export
import_branch_lengths.phyloXML_phylogeny <- function(x, phylo) {
  if (is.null(phylo$edge.length)) {
    stop("`edge.length` not found in phylo object.")
  }
  if (is.null(phylo$node.label)) {
    # TODO: use the tree structure to infer them.
    stop("`node.label` not found in phylo object.")
  }
  species_labels <- c(phylo$tip.label, phylo$node.label)
  species_table <- data.frame(
    parent = species_labels[phylo$edge[, 1]],
    child = species_labels[phylo$edge[, 2]],
    branchLength = phylo$edge.length
  )
  x$clade <- decorate_clades(x$clade, "branch_length", function(cl) {
    branch_length <- species_table[species_table$child == cl$name, ]$branchLength
    if (length(branch_length) == 1) {
      branch_length
    } else {
      NA
    }
  })
  # Set the root branch length
  x$clade$branch_length <- phylo$root.edge
  x
}

#' @export
import_branch_lengths.recPhyloXML <- function(x, phylo) {
  x$spTree <- import_branch_lengths(x$spTree, phylo)
  x
}

# Public flipping utils

#' Flip Children of a Specified Clade
#'
#' @description This function permutes the order of child clades. Every time it is called on the same parent clade, the order of its children becomes the next permutation in lexycographic order, coming back to the initial values when all possibilities have been exhausted. For example, for a node with three children, the initial order will be c(1, 2, 3); calling flip_children() once will give c(1, 3, 2), calling it again c(2, 1, 3), then c(2, 3, 1), c(3, 1, 2), c(3, 2, 1), and finally it will come back to c(1, 2, 3). The optional argument `perm` allows to specify a target permutation directly.
#'
#' @param x A `phyloXML` or `phyloXML_phylogeny` object containing the clades.
#' @param name A character string specifying the name of the clade whose children are to be flipped.
#' @param perm An optional numeric vector specifying the new permutation order of the children. If `NULL`, the next permutation in lexycographic order is used.
#' 
#' @return A new object with the child clades of the specified clade reordered.
#'
#' @details When applied to a `phyloXML` object, all of the phylogenies will be affected (if they have a clade with the specified name).
#
#' @note Flipping the order of children does NOT change the topology of the tree. It matters only for visualization purposes.
#'
#' @examples
#' phylogeny <- example_phyloXML_phylogeny()
#' flipped_phylogeny <- flip_children(phylogeny, "Bacteria")
flip_children <- function(x, name, perm = NULL) {
  UseMethod("flip_children")
}

#' @export
flip_children.phyloXML <- function(x, name, perm = NULL) {
  # TODO: consider whether this method makes sense.
  # Consider consistency when it's applied to recGeneTrees.
  lapply(x, flip_children, name = name, perm = perm)
}

#' @export
flip_children.phyloXML_phylogeny <- function(x, name, perm = NULL) {
  x$clade <- flip_children(x$clade, name = name, perm = perm)
  x
}

flip_children.phyloXML_clade <- function(x, name, perm = NULL) {
  decorate_clades(x, "clade", function(cl) {
    if (cl$name != name || length(cl$clade) <= 1) {
      # Not our target, don't change anything
      return(cl$clade)
    }
    current_order <- if (is.null(attr(cl, "children_permutation", exact = TRUE))) {
      seq_along(cl$clade)
    } else {
      attr(cl, "children_permutation")
    }
    new_order <- if (is.null(perm)) {
      nextperm(current_order)
    } else {
      stopifnot(is.numeric(perm))
      stopifont(length(perm) == length(cl$clade))
      perm
    }
    attr(cl, "children_permutation") <- new_order
    if (is.null(new_order)) {
      new_order <- seq_along(cl$clade)
    }
    cl$clade[order(current_order)][new_order]
  })
}

#' Flip the children of a clade in a spTree
#'
#' @inherit flip_children
#'
#' @param x A `recPhyloXML` object.
#'
#' @return The `recPhyloXML` object with the modified spTree.
#'
#' @details Only the spTree attribute will be affected.
#'
#' @export
flip_sp_children <- function(x, name, perm = NULL) {
  x$spTree <- flip_children(x$spTree, name = name, perm = perm)
  # TODO: Consider flipping the genes as well.
  # On second thought it seems too complex... I'd have to flip all the descendants genes as well because they are symmetrically flipped to
  # another side... Also, what to do for branchingOuts?
  x
}

#' Flip the children of a clade in recGeneTrees
#'
#' @inherit flip_children
#'
#' @param x A `recPhyloXML` object.
#'
#' @return The `recPhyloXML` object with the modified recGeneTrees.
#'
#' @details All of the recGeneTrees phylogenies will be affected, if they have a clade of the specified name.
#'
#' @export
flip_recGene_children <- function(x, name, perm = NULL) {
  # TODO: if name is anything but a duplication, flipping will not have any
  # effect. If it's a speciation, we should advise user to flip the
  # event_location species instead (or do it ourselves). If it's a
  # branchingOut, flip the mrca of the event_location species.
  # TODO: why couldn't we do x$recGeneTrees[[1]] <- flip_children(x$recGeneTrees[[1]]
  x$recGeneTrees <- lapply(x$recGeneTrees, flip_children, name = name, perm = perm)
  x
}

# Public annot utils

#' Add annotation to phyloXML objects
#'
#' @description Provide a data.frame with a `name` column matching the names of the clades in the phylogeny, and the new columns will be available in the layout for plotting.
#'
#' @param x A `phyloXML` or `phyloXML_phylogeny` object.
#' @param df A data.frame with a `name` column, to be matched with the nodes in the phylogeny.
#'
#' @return The modified `phyloXML` or `phyloXML_phylogeny` object.
#'
#' @details If x is a `phyloXML` object, the annotation is added to all of its phylogenies.
#'
#' @note If some annotation was already present, the new data.frame is
#' merged with the old one by the `name` column, keeping the existing
#' annotation. Use clear_annot() to remove all annotations.
#'
#' @export
add_annot <- function(x, df) {
  UseMethod("add_annot")
}

#' @export
add_annot.phyloXML <- function(x, df) {
  lapply(x, add_annot, df = df)
}

#' @export
add_annot.phyloXML_phylogeny <- function(x, df) {
  stopifnot(inherits(df, "data.frame"))
  if (is.null(attr(x, "annot", exact = TRUE))) {
    attr(x, "annot") <- df
  } else {
    attr(x, "annot") <- merge(attr(x, "annot", exact = TRUE), df, by = "name", all = TRUE)
  }
  x
}

#' Add annotation to recPhyloXML objects
#'
#' @inheritParams add_annot
#' @inherit add_annot description
#'
#' @param x An object of `recPhyloXML` class.
#'
#' @details add_sp_annot() adds the annotation to the spTree while add_recGene_annot() to the recGeneTrees.
#'
#' @note If some annotation was already present, the new data.frame is
#' merged with the old one by the `name` column, keeping the existing
#' annotation. Use clear_sp_annot() or clear_recGene_annot() to remove
#' all annotations.
#'
#' @export
add_sp_annot <- function(x, df) {
  # attr(x, "sp_annot") <- df  # duplicated attrs are not good
  x$spTree <- add_annot(x$spTree, df)
  x
}

#' @export
#' @rdname add_sp_annot
add_recGene_annot <- function(x, df) {
  # attr(x, "recGene_annot") <- df
  x$recGeneTrees <- add_annot(x$recGeneTrees, df)
  x
}

#' @export
clear_annot <- function(x) {
  UseMethod("clear_annot")
}

#' @export
clear_annot.phyloXML <- function(x) {
  lapply(x, clear_annot, df = df)
}

#' @export
clear_annot.phyloXML_phylogeny <- function(x) {
  attr(x, "annot") <- NULL
  x
}

#' @export
clear_sp_annot <- function(x) {
  # attr(x, "sp_annot") <- df  # duplicated attrs are not good
  x$spTree <- clear_annot(x$spTree)
  x
}

#' @export
clear_recGene_annot <- function(x) {
  # attr(x, "recGene_annot") <- df
  x$recGeneTrees <- clear_annot(x$recGeneTrees)
  x
}

# Private utils

make_unique_names <- function(root_xml) {
  # Give a name to unnamed clades
  unnamed <- xml2::xml_find_all(root_xml, "//clade[not(name)]")
  xml2::xml_add_child(unnamed, xml2::read_xml("<name>unnamed</name>"))
  # We also need to make sure that each gene parent has a unique name.
  # In many cases, the name is just "NULL".
  # Sometimes an internal node in the gene tree has the same name as a leaf node.
  # This is done by generax for who knows what reason.
  # We need to make sure that this doesn't happen.
  # Also, we make sure that leaves get the name without the # suffix and
  # the internal nodes get the suffix.
  # For phyloXML, the names must be unique across phylogenies.
  old_names <- xml2::xml_find_all(root_xml, "//name")
  new_names <- make.unique(xml2::xml_text(old_names), sep = "#")
  lapply(seq_along(old_names), function(i) {
    xml2::xml_text(old_names[i]) <- new_names[i]
  })
  invisible(NULL)
}

make_unique_recGeneTree_names <- function(root_xml) {
  # Give a name to unnamed clades
  unnamed <- xml2::xml_find_all(root_xml, "//clade[not(name)]")
  xml2::xml_add_child(unnamed, xml2::read_xml("<name>unnamed</name>"))
  # We also need to make sure that each gene parent has a unique name.
  # In many cases, the name is just "NULL".
  # Sometimes an internal node in the gene tree has the same name as a leaf node.
  # This is done by generax for who knows what reason.
  # We need to make sure that this doesn't happen.
  # Also, we make sure that leaves get the name without the # suffix and
  # the internal nodes get the suffix.
  # For recGeneTrees, the names must be unique across phylogenies.
  gene_names <- xml2::xml_find_all(root_xml, "//name")
  gene_events <- sapply(gene_names, function(g) {
    xml2::xml_name(xml2::xml_find_first(g, "../eventsRec/*[self::leaf or self::duplication or self::loss or self::branchingOut or self::speciation]"))
  })
  gene_events <- factor(gene_events, levels = c("leaf", "loss", "speciation", "duplication", "branchingOut"), ordered = TRUE)
  gene_names <- gene_names[order(gene_events)]
  new_names <- make.unique(xml2::xml_text(gene_names), sep = "#")
  lapply(seq_along(gene_names), function(i) {
    xml2::xml_text(gene_names[i]) <- new_names[i]
  })
  invisible(NULL)
}

nextperm <- function(vec) {
  stopifnot(length(vec) > 0)
  pivots <- which(vec[2:length(vec)] > vec[1:(length(vec)-1)])
  if (length(pivots) == 0) {
    return(NULL)
  }
  pivot <- max(pivots)
  tmp <- vec[pivot]
  successor <- which.min(vec[seq_along(vec) > pivot & vec > vec[pivot]]) + pivot
  vec[pivot] <- vec[successor]
  vec[successor] <- tmp
  vec[(pivot+1):length(vec)] <- rev(vec[(pivot+1):length(vec)])
  vec
}
