#' R6 Class Phylo
#'
#' @description
#' Class used to import and manipulate phylogenetic trees in [PhyloXML format](https://en.wikipedia.org/wiki/PhyloXML).
#'
#'
#' @param use_branch_length What to use as branch length
#' @param x_padding Distance between species in the plot
#' @param branch_length_scale Multiplier for the branch length
#' @param use_y_shift Whether to use the y_shift correction
#'
#' @export
Phylo <- R6::R6Class("Phylo",
  active = list(
    #' @field spList Species tree in list form.
    spList = function(value) {
      if (missing(value)) private$.spList
      else stop("Can't assign value to spList")
    },
    #' @field spNodes data.frame with the nodes of the species tree.
    spNodes = function(value) {
      if (missing(value)) {
        Reduce(\(x, y) merge(x, y, all.x = T), private$.spNodesAnnot, private$.spNodes)
      } else {
        stop("Can't assign value to spNodes")
      }
    },
    #' @field spEdges data.frame with the edges of the species tree.
    spEdges = function(value) {
      if (missing(value)) {
        Reduce(\(x, y) merge(x, y, all.x = T), private$.spEdgesAnnot, private$.spEdges)
      } else {
        stop("Can't assign value to spEdges")
      }
    }
  ),
  public = list(
    #' @description Create a new Phylo object from a PhyloXML file.
    #'
    #' @param xml_file Path to the PhyloXML file.
    #'
    #' @returns A new Phylo object.
    #'
    #' @examples
    #' xml_file <- system.file("extdata", "example_1.phyloxml", package = "recPhyloParse")
    #' p <- Phylo$new(xml_file)
    initialize = function(xml_file, use_branch_length = "branch_length", x_padding = 1, branch_length_scale = 1, use_y_shift = TRUE) {
      if ("xml_document" %in% class(xml_file)) {
        private$phylo_xml <- xml2::xml_unserialize(xml2::xml_serialize(xml_file, NULL))
      } else if (is.character(xml_file)) {
        private$phylo_xml <- xml2::read_xml(xml_file)
      } else {
        stop("`xml_file` must be either an xml_document or the path to an xml file.")
      }
      xml2::xml_ns_strip(private$phylo_xml)
      if (! (is.character(use_branch_length) || is.numeric(use_branch_length) || isFALSE(use_branch_length))) {
        stop("`use_branch_length` can only be a character vector, numeric vector, or FALSE.")
      }
      private$config <- list(
        use_y_shift = use_y_shift,
        use_branch_length = use_branch_length,
        x_padding = x_padding,
        branch_length_scale = branch_length_scale
      )
      all_species_names <- xml2::xml_text(xml2::xml_find_all(private$phylo_xml, "//phylogeny//clade/name"))

      self$redraw()
    },
    #' @description Redraw the reconciliated tree, changing parameters.
    #'
    #' @examples
    #' xml_file <- system.file("extdata", "example_1.recphyloxml", package = "recPhyloParse")
    #' rp <- RecPhylo$new(xml_file)
    #' rp$redraw(x_padding = 5)
    redraw = function(x_padding, use_branch_length, branch_length_scale, use_y_shift) {
      for (param in names(as.list(match.call())[-1])) {
        value <- get(param)
        # print(paste(param, "is", value))
        if (value != private$config[[param]]) {
          private$config[[param]] <- value
        }
      }
      spRoot <- xml2::xml_find_first(private$phylo_xml, "phylogeny/clade")
      private$max_x <- private$max_y <- 0
      private$.spList <- private$parse_sptree(spRoot, y_start = -private$config$branch_length_scale)
      private$.spNodes <- as.data.frame(private$.spList)
      private$.spEdges <- get_spedges(private$.spList)
      invisible(self)
    },
    #' @description Import branch lengths from another species tree.
    #'
    #' @param phylo Object of class "phylo" (e.g. from ape's read.tree())
    import_branch_lengths = function(phylo) {
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
      lapply(xml2::xml_find_all(private$phylo_xml, "phylogeny//clade"), function(cl) {
        child_name <- xml2::xml_text(xml2::xml_find_first(cl, "name"))
        branch_length <- species_table[species_table$child == child_name, ]$branchLength
        if (! length(branch_length) || is.null(branch_length)) {
          branch_length <- 1
        }
        xml2::xml_set_attr(cl, "branch_length", branch_length)
      })
      self$redraw()
    },
    #' @description Add species nodes attributes
    #'
    #' @param df data.frame to be merged with spNodes
    add_species_nodes_annotation = function(df) {
      private$.spNodesAnnot <- append(private$.spNodesAnnot, list(df))
    },
    #' @description Add species edges attributes
    #'
    #' @param df data.frame to be merged with spEdges
    add_species_edges_annotation = function(df) {
      private$.spEdgesAnnot <- append(private$.spEdgesAnnot, list(df))
    },
    #' @description Swap the children of a species node
    #'
    #' @details Flipping is idempotent: doing it twice restores the original state.
    #'
    #' @param spnames names of the species whose children should be flipped
    flip_species_children = function(spnames) {
      # For efficiency, the first time we should just reflect, then for subsequent redraws we should reverse the children during the parsing.
      # For now we just reverse the children and parse again.
      if (any(! spnames %in% names(private$.spNodesFlips))) {
        warning(setdiff(spnames, names(private$.spNodesFlips)), " are not species names in this tree and they will be ignored.")
        spnames <- intersect(spnames, names(private$.spNodesFlips))
      }
      for (sp in spnames) {
        private$.spNodesFlips[[sp]] <- nextperm(private$.spNodesFlips[[sp]])
      }
      self$redraw()
    },
    #' @description Print config, spNodes, and recGeneNodes.
    print = function() {
      cat("<config>\n")
      str(private$config)
      cat("\n<spNodes>\n")
      str(self$spNodes)
    },
    #' @description Print summary information about the tree.
    summarize = function() {
      cat("Phylo object with ", nrow(self$spNodes), " species (", sum(self$spNodes$is_leaf), " leaves)\n", sep = "")
    },
    #' @description Save the tree in PhyloXML format.
    #'
    #' @param file Path to the file that will be written
    #' @param ... Additional options passed to xml2::write_xml()
    write = function(file, ...) {
      xml2::write_xml(private$phylo_xml, file)
    },
    #' @description Create a test plot of the tree.
    #'
    #' @details
    #' This method is used for testing purposes, just to check that all the
    #' elements are OK. Feel free to copy the code and edit it to make your
    #' own pretty plots.
    plot = function() {
      if (requireNamespace("ggplot2", quietly = T)) {
        auxpoints <- data.frame(
          x = c(private$max_x, self$spNodes[is.na(self$spNodes$parent), "x"]),
          y = c(private$max_y + private$config$branch_length_scale, self$spNodes[is.na(self$spNodes$parent), "y"] - private$config$branch_length_scale)
        )
        ggplot2::ggplot() +
          ggplot2::geom_line(data = self$spEdges, ggplot2::aes(x, y, group = group), lineend = "round") +
          ggplot2::geom_point(data = self$spNodes, ggplot2::aes(x, y)) +
          ggplot2::geom_text(data = self$spNodes, ggplot2::aes(x, y, label = name)) +
          ggplot2::geom_point(data = auxpoints, ggplot2::aes(x, y), alpha = 0) +
          # coord_polar() +
          # coord_flip() +
          # theme_void() +
          NULL
      } else {
        stop("Please install `ggplot2` before using the plot() method.")
      }
    }
  ),
  private = list(
    .spList = NULL,
    .spNodes = NULL,
    .spEdges = NULL,
    .spNodesAnnot = list(),
    .spEdgesAnnot = list(),
    .spNodesFlips = list(),
    phylo_xml = NULL,
    config = list(),
    warnings = list(negative_branch_height = T, missing_branch_length = T),
    max_x = 0,
    max_y = 0,
    side = list(),
    parse_sptree = function(spnode, parent = NA, side = "root", y_start = 0, warn = T) {
      name <- xml2::xml_text(xml2::xml_find_first(spnode, "name"))
      # Find all gene events that occur in this species
      half_x_thickness <- 0
      half_y_thickness <- 0
      # Branch length
      if (is.numeric(private$config$use_branch_length)) {
        branch_length <- private$config$use_branch_length * private$config$branch_length_scale
        y <- y_start + branch_length
      } else if (is.character(private$config$use_branch_length)) {
        branch_length <- as.numeric(xml2::xml_attr(spnode, private$config$use_branch_length))
        if (is.na(branch_length)) {
          branch_length <- xml2::xml_double(xml2::xml_find_first(spnode, private$config$use_branch_length))
        }
        if (is.na(branch_length)) {
          branch_length <- 1
          if (isTRUE(private$warnings$missing_branch_length)) {
            private$warnings$missing_branch_length <- FALSE
            warning("'", private$config$use_branch_length, "' not found in clade ", name, ". Setting branch length to 1 automatically.", call. = FALSE)
          }
        }
        branch_length <- branch_length * private$config$branch_length_scale
        y <- y_start + branch_length
      } else if (isFALSE(private$config$use_branch_length)) {
        branch_length <- private$config$branch_length_scale
        y <- 0
      } else {
        stop("Unsupported type for `use_branch_length`.")
      }
      y_shift <- 0
      children <- xml2::xml_find_all(spnode, "./clade")
      if (length(children) == 0) {
        # It's an extant so we can have gene leaves but not speciations
        # The number of events is the x_thickness
        left_child <- right_child <- NULL
        x <- private$max_x + half_x_thickness
        private$max_x <- x + half_x_thickness + private$config$x_padding
        is_leaf <- TRUE
      } else if (length(children) == 2) {
        if (name %in% private$.spNodesFlips) {
          children <- rev(children)
        }
        left_child <- private$parse_sptree(children[[1]], parent = name, side = "left", y_start = y, warn = warn)
        right_child <- private$parse_sptree(children[[2]], parent = name, side = "right", y_start = y, warn = warn)
        x <- (left_child$x + left_child$half_x_thickness + right_child$x - right_child$half_x_thickness) / 2
        is_leaf <- FALSE
        if (private$config$use_branch_length == FALSE) {
          y <- min(left_child$y, right_child$y) - private$config$branch_length_scale
        }
        # left_child_min_branch_height <- sum(xml2::xml_name(private$internal_events[[left_child$name]]) %in% c("duplication", "loss", "branchingOut", "transferBack"))
        # right_child_min_branch_height <- sum(xml2::xml_name(private$internal_events[[right_child$name]]) %in% c("duplication", "loss", "branchingOut", "transferBack"))
        left_child_min_branch_height <- 0
        right_child_min_branch_height <- 0
        if (private$config$use_y_shift == T) {
          y_shift <- min(
            left_child$y - half_y_thickness - y - left_child$half_y_thickness + left_child$y_shift - left_child_min_branch_height,
            right_child$y - half_y_thickness - y - right_child$half_y_thickness + right_child$y_shift - right_child_min_branch_height,
            0
          )
        }
      } else {
        stop("This species tree is not binary.")
      }
      if (y > private$max_y) {
        private$max_y <- y
      }
      # Return the list
      splist <- list(
        name = name,
        parent = parent,
        side = side,
        is_leaf = is_leaf,
        x = x,
        y = y,
        half_x_thickness = half_x_thickness,
        half_y_thickness = half_y_thickness,
        y_shift = y_shift,
        left_child = left_child,
        right_child = right_child
      )
      class(splist) <- c("recphylo_spTree", class(splist))
      return(splist)
    }
  )
)



