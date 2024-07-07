# NOTE: gene and species names will be made unique
# NOTE: phylogenies will be given a unique idx field (their index)

read_phyloXML <- function(x, ...) {
  xml <- read_xml(x, ...)
  parse_phyloXML(xml)
}

read_recPhyloXML <- function(x, ...) {
  xml <- read_xml(x, ...)
  parse_recPhyloXML(xml)
}

parse_phyloXML <- function(xml) {
  make_unique_names(xml)
  phylogenies <- xml_find_all(xml, "./phylogeny")
  l <- lapply(seq_along(phylogenies), function(phy_idx) {
    phy <- parse_phylogeny(phylogenies[[phy_idx]], type = "phyloXML")
    phy$idx <- phy_idx
    phy
  })
  class(l) <- c("phyloXML", class(l))
  l
}

parse_recPhyloXML <- function(xml) {
  make_unique_names(xml_find_first(xml, "./spTree"))
  make_unique_recGeneTree_names(xml_find_first(xml, "./recGeneTree"))
  spTree_xml <- xml_find_all(xml, "./spTree/phylogeny")
  recGeneTree_xml <- xml_find_all(xml, "./recGeneTree/phylogeny")
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
    parse_phyloXML_clade(xml_find_first(phylogeny_xml, "./clade"))
  } else if (type == "spTree") {
    parse_spTree_clade(xml_find_first(phylogeny_xml, "./clade"))
  } else if (type == "recGeneTree") {
    parse_recGeneTree_clade(xml_find_first(phylogeny_xml, "./clade"))
  } else {
    stop("Unsupported phylogeny type.")
  }
  rooted <- as.logical(xml_attr(phylogeny_xml, "rooted"))
  if (is.na(rooted)) {
    rooted <- length(clade$clade) == 2
    warning("Phylogeny doesn't have a 'rooted' attribute. Inferred '", rooted, "'.")
  }
  phylogeny <- list(
    name = xml_text(xml_find_first(phylogeny_xml, "./name")),
    id = xml_text(xml_find_first(phylogeny_xml, "./id")),
    description = xml_text(xml_find_first(phylogeny_xml, "./description")),
    date = as.Date(xml_text(xml_find_first(phylogeny_xml, "./date"))),
    rooted = rooted,
    rerootable = as.logical(xml_attr(phylogeny_xml, "rerootable")),
    branch_length_unit = xml_attr(phylogeny_xml, "branch_length_unit"),
    type = xml_attr(phylogeny_xml, "type"),
    # confidence = xml_find_all(phylogeny_xml, "./confidence"),
    # clade_relation = xml_find_all(phylogeny_xml, "./clade_relation"),
    # sequence_relation = xml_find_all(phylogeny_xml, "./sequence_relation"),
    # property = xml_find_all(phylogeny_xml, "./property"),
    clade = clade
  )
  class(phylogeny) <- c("phyloXML_phylogeny", class(phylogeny))
  phylogeny
}

parse_clade <- function(clade_xml) {
  name <- xml_text(xml_find_first(clade_xml, "./name"))
  branch_length_attr <- as.numeric(xml_attr(clade_xml, "branch_length"))
  branch_length_tag <- xml_double(xml_find_first(clade_xml, "branch_length"))
  branch_length <- if (!is.na(branch_length_attr) && !is.na(branch_length_tag)) {
    warning("clade '", name, "' defines branch length both with an attribute and with a tag, which is ambiguous. I chose the attribute value, ", branch_length_attr, ".")
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
    # id_source = xml_double(xml_find_first(clade_xml, "./id_source")),
    collapse = as.logical(xml_attr(clade_xml, "collapse")),
    # confidence = xml_text(xml_find_all(clade_xml, "./confidence")),
    width = xml_double(xml_find_first(clade_xml, "./width")),
    # color = xml_text(xml_find_first(clade_xml, "./color")),
    # taxonomy = xml_text(xml_find_all(clade_xml, "./taxonomy")),
    # sequence = xml_text(xml_find_all(clade_xml, "./sequence")),
    # events = xml_text(xml_find_first(clade_xml, "./events")),
    # binary_characters = xml_text(xml_find_first(clade_xml, "./binary_characters")),
    # distribution = xml_text(xml_find_all(clade_xml, "./distribution")),
    date = as.Date(xml_text(xml_find_first(clade_xml, "./date")))
    # reference = xml_text(xml_find_all(clade_xml, "./reference")),
    # property = xml_text(xml_find_all(clade_xml, "./property")),
  )
  class(clade) <- c("phyloXML_clade", class(clade))
  clade
}

parse_phyloXML_clade <- function(clade_xml) {
  clade <- parse_clade(clade_xml)
  clade$clade <- lapply(xml_find_all(clade_xml, "./clade"), parse_phyloXML_clade)
  return(clade)
}

parse_spTree_clade <- function(clade_xml) {
  clade <- parse_clade(clade_xml)
  clade$clade <- lapply(xml_find_all(clade_xml, "./clade"), parse_spTree_clade)
  # clade$geography <- ... TODO
  return(clade)
}

parse_recGeneTree_clade <- function(clade_xml) {
  clade <- parse_clade(clade_xml)
  event_xml <- xml_find_first(clade_xml, "./eventsRec/*[self::leaf or self::duplication or self::loss or self::branchingOut or self::speciation]")
  clade$event_type <- xml_name(event_xml)
  clade$event_location <- xml_attr(event_xml, "speciesLocation")
  # clade$geography <- ... TODO
  clade$clade <- lapply(xml_find_all(clade_xml, "./clade"), parse_recGeneTree_clade)
  if (length(xml_find_first(clade_xml, "./eventsRec/*[self::transferBack]")) > 0) {
    parent <- clade
    parent$name <- paste(parent$name, "transferBack", sep = "@")
    parent$event_type <- "transferBack"
    parent$clade <- list(clade)
    return(parent)
  }
  return(clade)
}

# data.frame interop

as.data.frame.phyloXML <- function(phylo_xml) {
  lapply(seq_along(phylo_xml), function(phy_idx) {
    as.data.frame(phylo_xml[[phy_idx]])
  })
}

as.data.frame.phyloXML_phylogeny <- function(phy) {
  phy_fields <- setdiff(names(phy), "clade")
  phy_values <- phy[phy_fields]
  names(phy_values) <- paste("phylogeny", phy_fields, sep = "_")
  # TODO: consider removing fields that are NA in the phylogeny
  cbind(as.data.frame(phy$clade), phy_values)
}

as.data.frame.phyloXML_clade <- function(cl) {
  fields <- setdiff(names(cl), "clade")
  rbind(
    as.data.frame(cl[fields]),
    Reduce(rbind, lapply(cl$clade, as.data.frame))
  )
}

as.data.frame.recPhyloXML <- function(recphylo_xml) {
  sp_df <- as.data.frame(recphylo_xml$spTree)
  sp_df$tree_type <- "spTree"
  recGene_df <- as.data.frame(recphylo_xml$recGeneTrees)
  recGene_df$tree_type <- "recGeneTree"
  rbind(sp_df, recGene_df)
}

# Public utils

decorate_clades <- function(clade, attr, FUN, ...) {
  clade[[attr]] <- FUN(clade, ...)
  clade$clade <- lapply(clade$clade, decorate_clades, attr, FUN, ...)
  clade
}

traverse_clades <- function(clade, FUN, ..., .order = "pre") {
  item <- FUN(clade, ...)
  children <- lapply(clade$clade, traverse_clades, FUN, .order = .order, ...)
  if (.order == "pre") {
    c(item, unlist(children))
  } else if (.order == "post") {
    c(unlist(children), item)
  } else {
    stop("Unsupported `.order`. Please choose either 'pre' or 'post'.")
  }
}

# Public branch length utils

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

import_branch_lengths.recPhyloXML <- function(x, phylo) {
  x$spTree <- import_branch_lengths(x$spTree, phylo)
  x
}

# Public flipping utils

flip_children.phyloXML <- function(x, name, perm = NULL) {
  # TODO: consider whether this method makes sense.
  # Consider consistency when it's applied to recGeneTrees.
  lapply(x, flip_children, name = name, perm = perm)
}

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
      seq_along(cl$children)
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
      new_order <- seq_along(cl$children)
    }
    cl$clade[order(current_order)][new_order]
  })
}

flip_sp_children <- function(x, name, perm = NULL) {
  x$spTree <- flip_children(x$spTree, name = name, perm = perm)
  # TODO: Consider flipping the genes as well.
  # On second thought it seems too complex... I'd have to flip all the descendants genes as well because they are symmetrically flipped to
  # another side... Also, what to do for branchingOuts?
  x
}

flip_recGene_children <- function(x, name, perm = NULL) {
  # TODO: if name is anything but a duplication, flipping will not have any
  # effect. If it's a speciation, we should advise user to flip the
  # event_location species instead (or do it ourselves). If it's a
  # branchingOut, flip the mrca of the event_location species.
  # TODO: why couldn't we do x$recGeneTrees[[1]] <- flip_children(x$recGeneTrees[[1]]
  x$recGeneTrees <- lapply(x$recGeneTrees, flip_children, name = name, perm = perm)
}

# Public annot utils

add_annot <- function(x, df) {
  UseMethod("add_annot")
}

add_annot.phyloXML <- function(x, df) {
  invisible(lapply(x, add_annot, df = df))
}

add_annot.phyloXML_phylogeny <- function(x, df) {
  attr(x, "annot") <- df
  invisible(x)
}

add_sp_annot <- function(x, df) {
  # attr(x, "sp_annot") <- df  # duplicated attrs are not good
  add_annot(x$spTree, df)
  invisible(x)
}

add_recGene_annot <- function(x, df) {
  # attr(x, "recGene_annot") <- df
  add_annot(x$recGeneTrees, df)
  invisible(x)
}

# Private utils

make_unique_names <- function(root_xml) {
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
