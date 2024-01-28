as.data.frame.recphylo_spTree <- function(l) {
  if (is.null(l)) {
    return(NULL)
  }
  fields <- c("name", "parent", "side", "is_leaf", "x", "y", "half_x_thickness", "half_y_thickness", "y_shift")
  rbind(
    as.data.frame(l[fields]),
    as.data.frame(l$left_child),
    as.data.frame(l$right_child)
  )
}

as.data.frame.recphylo_recGeneTree <- function(l) {
  if (is.null(l)) {
    return(NULL)
  }
  fields <- c("id", "name", "event_type", "event_location", "parent", "side", "lineage", "x", "y")
  rbind(
    as.data.frame(l[fields]),
    as.data.frame(l$left_child),
    as.data.frame(l$right_child)
  )
}

# TODO: https://arxiv.org/pdf/2008.08960.pdf

# Pass aes(y = -y) everywhere and then coord_polar() to make an inverted radial plot. or use scale_y_reverse() and coord_polar().

#' @export
RecPhylo <- R6::R6Class("RecPhylo",
  public = list(
    spList = NULL,
    spNodes = NULL,
    spEdges = NULL,
    gList = NULL,
    recGeneNodes = NULL,
    recGeneEdges = NULL,
    initialize = function(xml_file, use_branch_length = "branch_length", use_y_shift = TRUE, x_padding = 1, branch_length_scale = 1) {
      if ("xml_document" %in% class(xml_file)) {
        private$recphylo_xml <- xml2::xml_unserialize(xml2::xml_serialize(xml_file, NULL))
      } else if (is.character(xml_file)) {
        private$recphylo_xml <- xml2::read_xml(xml_file)
      } else {
        stop("`xml_file` must be either an xml_document or the path to an xml file.")
      }
      xml2::xml_ns_strip(private$recphylo_xml)
      if (! (is.character(use_branch_length) || is.numeric(use_branch_length) || isFALSE(use_branch_length))) {
        stop("`use_branch_length` can only be a character vector, numeric vector, or FALSE.")
      }
      private$config <- list(
        use_y_shift = use_y_shift,
        use_branch_length = use_branch_length,
        x_padding = x_padding,
        branch_length_scale = branch_length_scale
      )
      species_names <- xml2::xml_text(xml2::xml_find_all(private$recphylo_xml, "//spTree/phylogeny//clade/name"))
      private$internal_events <- sapply(species_names, function(sp) {
          xml2::xml_find_all(private$recphylo_xml, paste0("recGeneTree//*[@speciesLocation='", sp, "'] | recGeneTree//transferBack[following-sibling::*[@speciesLocation='", sp, "']]"))
      }, simplify = FALSE)
      self$redraw()
    },
    redraw = function(x_padding, use_branch_length, branch_length_scale, use_y_shift) {
      for (param in names(as.list(match.call())[-1])) {
        value <- get(param)
        # print(paste(param, "is", value))
        if (value != private$config[[param]]) {
          private$config[[param]] <- value
        }
      }
      spRoot <- xml2::xml_find_first(private$recphylo_xml, "spTree//clade")
      gRoot <- xml2::xml_find_first(private$recphylo_xml, "recGeneTree//clade")
      private$max_x <- private$max_y <- 0
      private$warnings$negative_branch_height <- T
      self$spList <- private$parse_sptree(spRoot, y_start = -private$config$branch_length_scale)
      private$calc_branch_heights(self$spList)
      self$gList <- private$parse_gtree(gRoot)
      self$spNodes <- as.data.frame(self$spList)
      self$spEdges <- get_spedges(self$spList)
      self$recGeneNodes <- as.data.frame(self$gList)
      self$recGeneEdges <- get_gedges(self$gList)
      invisible(self)
    },
    flip_children_species = function(sp) {
      # TODO. should be relatively easy, just reflect across the node's x all
      # the x coords of downstream genes and species. but we need to take that
      # into account when we redraw, so we should save it somewhere...
    },
    flip_children_gene = function(g) {
      # TODO. more tricky, may entail flipping species as well.
    },
    print = function() {
      cat("<config>\n")
      str(private$config)
      cat("\n<spNodes>\n")
      str(self$spNodes)
      cat("\n<recGeneNodes>\n")
      str(self$recGeneNodes)
    },
    summary = function() {
      cat("RecPhylo object with ", nrow(self$spNodes), " species (", sum(self$spNodes$is_leaf), " leaves) and ", sum(self$recGeneNodes$event_type == "leaf"), " genes (", nrow(self$recGeneNodes), " events)\n", sep = "")
    },
    plot = function() {
      if (requireNamespace("ggplot2", quietly = T)) {
        auxpoints <- data.frame(
          x = c(private$max_x, self$spNodes[is.na(self$spNodes$parent), "x"]),
          y = c(private$max_y + private$config$branch_length_scale, self$spNodes[is.na(self$spNodes$parent), "y"] - private$config$branch_length_scale)
        )
        ggplot2::ggplot() +
          ggplot2::geom_line(data = self$spEdges, ggplot2::aes(x, y, group = group)) +
          ggplot2::geom_point(data = self$spNodes, ggplot2::aes(x, y)) +
          ggplot2::geom_text(data = self$spNodes, ggplot2::aes(x, y, label = name)) +
          ggplot2::geom_point(data = self$recGeneNodes, ggplot2::aes(x, y)) +
          ggplot2::geom_line(data = self$recGeneEdges, ggplot2::aes(x, y, group = group, color = gsub("l+$", "", lineage), linetype = event_type), show.legend = F) +
          ggplot2::geom_point(data = auxpoints, ggplot2::aes(x, y), alpha = 0) +
          ggplot2::scale_linetype_manual(values = c("loss_v" = 4, "transferBack" = 3), na.value = 1) +
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
    recphylo_xml = NULL,
    config = list(),
    warnings = list(negative_branch_height = T, missing_branch_length = T),
    max_x = 0,
    max_y = 0,
    internal_events = list(),
    intra_species_x = list(),
    intra_species_y = list(),
    intra_species_h = list(),
    side = list(),
    branch_height_fraction = list(),
    parse_sptree = function(spnode, parent = NA, side = "root", y_start = 0, warn = T) {
      name <- xml2::xml_text(xml2::xml_find_first(spnode, "name"))
      # Find all gene events that occur in this species
      events <- private$internal_events[[name]]
      half_x_thickness <- max(sum(xml2::xml_name(events) %in% c("speciation", "loss", "leaf")) + 1, 2) / 2
      half_y_thickness <- max(sum(xml2::xml_name(events) == "speciation") + 1, 2) / 2
      # Branch length
      if (is.numeric(private$config$use_branch_length)) {
        branch_length <- private$config$use_branch_length * private$config$branch_length_scale
        y <- y_start + branch_length
      } else if (is.character(private$config$use_branch_length)) {
        branch_length <- xml2::xml_double(xml2::xml_find_first(spnode, private$config$use_branch_length))
        if (is.na(branch_length)) {
          branch_length <- 1
          if (isTRUE(private$warnings$missing_branch_length)) {
            private$warnings$missing_branch_length <- FALSE
            warning("'", private$config$use_branch_length, "' not found in clade ", name, ". Setting branch length to 1 automatically.", call. = FALSE)
          }
        }
        branch_length <- branch_length * private$config$branch_length_scale
        y <- y_start + branch_length
      } else {
        branch_length <- private$config$branch_length_scale
        y <- 0
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
        left_child <- private$parse_sptree(children[[1]], parent = name, side = "left", y_start = y, warn = warn)
        right_child <- private$parse_sptree(children[[2]], parent = name, side = "right", y_start = y, warn = warn)
        x <- (left_child$x + left_child$half_x_thickness + right_child$x - right_child$half_x_thickness) / 2
        is_leaf <- FALSE
        if (private$config$use_branch_length == FALSE) {
          y <- min(left_child$y, right_child$y) - private$config$branch_length_scale
        }
        left_child_min_branch_height <- sum(xml2::xml_name(private$internal_events[[left_child$name]]) %in% c("duplication", "loss", "branchingOut", "transferBack"))
        right_child_min_branch_height <- sum(xml2::xml_name(private$internal_events[[right_child$name]]) %in% c("duplication", "loss", "branchingOut", "transferBack"))
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
    },
    calc_branch_heights = function(splist, parent_y = 0, parent_half_y_thickness = 0, parent_y_shift = 0, warn = TRUE) {
      if (is.null(splist)) {
        return(invisible(self))
      }
      private$side[[splist$name]] <- splist$side
      branch_height <- splist$y + splist$y_shift - parent_y - parent_y_shift - splist$half_y_thickness - parent_half_y_thickness
      min_branch_height <- sum(xml2::xml_name(private$internal_events[[splist$name]]) %in% c("duplication", "loss", "branchingOut", "transferBack"))
      if (isTRUE(private$warnings$negative_branch_height) && min_branch_height + parent_y_shift < 0) {
        private$warnings$negative_branch_height = FALSE
        warning("Negative branch heights: the tree will be gibberish. Please increase `branch_length_scale`.", call. = FALSE)
      }
      private$branch_height_fraction[[splist$name]] <- branch_height / (min_branch_height + 1)
      private$intra_species_h[[splist$name]] <- parent_y + parent_y_shift + parent_half_y_thickness + private$branch_height_fraction[[splist$name]]
      if (splist$side == "left") {
        private$intra_species_x[[splist$name]] <- splist$x - splist$half_x_thickness + 1
      } else {
        private$intra_species_x[[splist$name]] <- splist$x + splist$half_x_thickness - 1
      }
      if (splist$is_leaf == TRUE) {
        private$intra_species_y[[splist$name]] <- splist$y
      } else {
        private$intra_species_y[[splist$name]] <- splist$y - splist$half_y_thickness + splist$y_shift + 1
      }
      private$calc_branch_heights(splist$left_child, splist$y, splist$half_y_thickness, splist$y_shift, warn)
      private$calc_branch_heights(splist$right_child, splist$y, splist$half_y_thickness, splist$y_shift, warn)
      return(invisible(self))
    },
    parse_gtree = function(gnode, parent = NA, side = "root", lineage = "r") {
      # * Depth first is the only way to preserve the order. If a branch was at the top in one species node, it must be at the top in the descendant species as well. So it's also pre-order: we visit the node itself first, then the left child, then the right child.
      # * If it's a speciation, we put it in the first (most external) available x-y coord IN THE INTERNAL LEVEL (not in the branch height level). Then we increment x and y for the internal level, and we increment x for the branch height level.
      # * If it's a duplication, we put it in the first available x-y coord in the branch-height level, then we increment x and y for the branch height level, we don't touch the internal level, so the next event will fall exactly in place.
      # * If it's a leaf, we put it in the first availalbe x coord (y is same as the species leaf), and increment x-coord. we also increment x-coord in the branch-height level.
      # * so every time we increment the x-coord for the internal level, we do so for the branch too.
      # * If it's a transfer and we are in the same species, same as duplication.
      # * If it's a loss, we increment the x-coord of both internal and branch levels, we increment the y of the branch level.
      name <- xml2::xml_text(xml2::xml_find_first(gnode, "name"))
      event_xml <- xml2::xml_find_first(gnode, "./eventsRec/*[self::leaf or self::duplication or self::loss or self::branchingOut or self::speciation]")
      event_type <- xml2::xml_name(event_xml)
      event_location <- xml2::xml_attr(event_xml, "speciesLocation")
      id <- paste(name, parent, sep = "@")
      x <- private$intra_species_x[[event_location]]
      if (event_type == "speciation") {
        y <- private$intra_species_y[[event_location]]
        private$increment_x(event_location)
        private$increment_y(event_location)
      } else if (event_type %in% c("duplication", "branchingOut")) {
        y <- private$intra_species_h[[event_location]]
        private$increment_h(event_location)
      } else if (event_type == "leaf") {
        y <- private$intra_species_y[[event_location]]
        private$increment_x(event_location)
      } else if (event_type == "loss") {
        y <- private$intra_species_h[[event_location]]
        private$increment_x(event_location)
        private$increment_h(event_location)
      }
      children <- xml2::xml_find_all(gnode, "./clade")
      if (length(children) == 0) {
        left_child <- right_child <- NULL
      } else if (length(children) == 2) {
        if (event_type == "duplication") {
          lineage1 <- paste0(lineage, "l")
          lineage2 <- paste0(lineage, "r")
        } else {
          lineage1 <- lineage2 <- lineage
        }
        left_child <- private$parse_gtree(
          children[[1]], parent = id, side = "left", lineage = lineage1)
        right_child <- private$parse_gtree(
          children[[2]], parent = id, side = "right", lineage = lineage2)
      } else {
        stop("This gene tree is not binary.")
      }
      glist <- list(
        id = id,
        name = name,
        parent = parent,
        event_type = event_type,
        event_location = event_location,
        side = side,
        x = x,
        y = y,
        lineage = lineage,
        left_child = left_child,
        right_child = right_child
      )
      class(glist) <- c("recphylo_recGeneTree", class(glist))
      # If there was a transferBack, add this node as a left child of the transfer
      if (length(xml2::xml_find_first(gnode, "./eventsRec/*[self::transferBack]")) > 0) {
        new_id <- paste(id, "transferBack", sep = "@")
        glist$parent <- new_id
        glist <- list(
          id = new_id,
          name = paste(name, "transferBack", sep = "@"),
          parent = parent,
          event_type = "transferBack",
          event_location = event_location,
          side = side,
          x = x,
          y = private$intra_species_h[[event_location]],
          lineage = lineage,
          left_child = glist,
          right_child = NULL
        )
        class(glist) <- c("recphylo_recGeneTree", class(glist))
        private$increment_h(event_location)
      }
      return(glist)
    },
    increment_x = function(spname) {
      private$intra_species_x[[spname]] <- private$intra_species_x[[spname]] + 2*(private$side[[spname]] == "left")-1
    },
    increment_y = function(spname) {
      private$intra_species_y[[spname]] <- private$intra_species_y[[spname]] + 1
    },
    increment_h = function(spname) {
      private$intra_species_h[[spname]] <- private$intra_species_h[[spname]] + private$branch_height_fraction[[spname]]
    }
  )
)
