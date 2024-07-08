# NOTE: Pass aes(y = -y) everywhere and then coord_polar() to make an inverted radial plot. or use scale_y_reverse() and coord_polar().

# @export
RecPhyloLayout <- R6::R6Class("RecPhyloLayout",
  public = list(
    initialize = function(recPhylo, use_branch_length = TRUE, x_padding = 1, branch_length_scale = 1, use_y_shift = TRUE) {
      stopifnot("recPhyloXML" %in% class(recPhylo))
      private$.recPhylo <- recPhylo
      private$set_config(use_branch_length, x_padding, branch_length_scale, use_y_shift)
      species_names <- traverse_clades(recPhylo$spTree$clade, function(x) x$name)
      private$internal_events <- sapply(species_names, function(sp) {
        lapply(recPhylo$recGeneTrees, function(rgt) {
          traverse_clades(rgt$clade, function(x) {
            if (x$event_location == sp) {
              x$event_type
            }
          })
        })
      }, simplify = FALSE)
      self$redraw()
    },
    redraw = function() {
      private$max_x <- 0
      private$warnings$negative_branch_height <- TRUE
      private$layout_spTree <- private$species_pipes_layout(private$.recPhylo$spTree$clade)
      private$inner_species_pipes_layout(private$layout_spTree)
      private$layout_recGeneTrees <- lapply(private$.recPhylo$recGeneTrees, function(rgt) {
        private$gene_tree_layout(rgt$clade)
      })
      invisible(self)
    },
    plot = function() {
      if (!requireNamespace("ggplot2", quietly = T)) {
        stop("Please install `ggplot2` before using the plot() method.")
      }
      spNodes <- self$spLayout$nodes
      spEdges <- self$spLayout$edges
      geneNodes <- self$recGeneLayout$nodes
      geneEdges <- self$recGeneLayout$edges
      auxpoints <- data.frame(
        x = c(private$max_x, private$layout_spTree$x),
        y = c(max(spNodes$y) + private$config$branch_length_scale, private$layout_spTree$y - private$config$branch_length_scale)
      )
      ggplot2::ggplot() +
        ggplot2::geom_line(data = spEdges, ggplot2::aes(x, y, group = group), lineend = "round") +
        ggplot2::geom_point(data = spNodes, ggplot2::aes(x, y)) +
        ggplot2::geom_text(data = spNodes, ggplot2::aes(x, y, label = name)) +
        ggplot2::geom_point(data = geneNodes, ggplot2::aes(x, y)) +
        ggplot2::geom_line(data = geneEdges, ggplot2::aes(x, y, group = group, color = lineage, linetype = paste(event_type, leg, sep = "_")), lineend = "round", show.legend = F) +
        ggplot2::geom_point(data = auxpoints, ggplot2::aes(x, y), alpha = 0) +
        ggplot2::scale_linetype_manual(values = c("loss_vertical" = 2, "transferBack_transfer" = 3), na.value = 1) +
        # coord_polar() +
        # coord_flip() +
        # theme_void() +
        NULL
    }
  ),
  active = list(
    recPhylo = function(value) {
      if (!missing(value)) {
        stop("Cannot assign to recPhylo. Initialize a new object instead.")
      }
      private$.recPhylo
    },
    spLayout = function(value) {
      if (!missing(value)) {
        stop("Can't assign to spLayout.")
      }
      list(
        nodes = merge_layout(private$.recPhylo$spTree, as.data.frame(private$layout_spTree)),
        edges = merge_layout(private$.recPhylo$spTree, get_species_pipes_edges(private$layout_spTree))
      )
    },
    recGeneLayout = function(value) {
      if (!missing(value)) {
        stop("Can't assign to recGeneLayout.")
      }
      list(
        nodes = Reduce(rbind, lapply(seq_along(private$.recPhylo$recGeneTrees), function(phylogeny_idx) {
          merge_layout(
            private$.recPhylo$recGeneTrees[[phylogeny_idx]],
            as.data.frame(private$layout_recGeneTrees[[phylogeny_idx]])
          )
        })),
        edges = Reduce(rbind, lapply(seq_along(private$.recPhylo$recGeneTrees), function(phylogeny_idx) {
          merge_layout(
            private$.recPhylo$recGeneTrees[[phylogeny_idx]],
            get_gene_tree_edges(private$layout_recGeneTrees[[phylogeny_idx]])
          )
        }))
      )
    }
  ),
  private = list(
    .recPhylo = list(),
    config = list(),
    warnings = list(missing_branch_length = TRUE, negative_branch_height = TRUE),
    max_x = 0,
    internal_events = list(),
    layout_spTree = list(),
    layout_recGeneTrees = list(),
    inner_layout = list(),
    set_config = function(use_branch_length, x_padding, branch_length_scale, use_y_shift) {
      if (! (is.numeric(use_branch_length) || is.logical(use_branch_length))) {
        stop("`use_branch_length` can only be a numeric or logical value.")
      } else if (length(use_branch_length) != 1) {
        stop("`use_branch_length` must have length one.")
      }
      private$config <- list(
        use_branch_length = use_branch_length,
        x_padding = x_padding,
        branch_length_scale = branch_length_scale,
        use_y_shift = use_y_shift
      )
    },
    species_pipes_layout = function(clade, child_idx = 0, y_start = 0, parent_half_y_thickness = 0) {
      # We need x, y, half_x_thickness, half_y_thickness, branch_height, is_leaf
      events <- unlist(private$internal_events[[clade$name]])
      half_x_thickness <- max(sum(events %in% c("speciation", "loss", "leaf")) + 1, 2) / 2
      half_y_thickness <- max(sum(events == "speciation") + 1, 2) / 2
      min_branch_height <- sum(events %in% c("duplication", "loss", "branchingOut", "transferBack"))
      branch_length <- if (is.numeric(private$config$use_branch_length)) {
        private$config$use_branch_length * private$config$branch_length_scale
      } else if (isTRUE(private$config$use_branch_length)) {
        if (!is.na(clade$branch_length)) {
          clade$branch_length * private$config$branch_length_scale
        } else {
          if (isTRUE(private$warnings$missing_branch_length)) {
            warning("No branch length was found in clade ", clade$name, ". Setting branch length to 1 automatically.", call. = FALSE)
            private$warnings$missing_branch_length <- FALSE
          }
          private$config$branch_length_scale
        }
      } else if (isFALSE(private$config$use_branch_length)) {
        private$config$branch_length_scale
      } else {
        stop("Unsupported type for `use_branch_length`.")
      }
      y <- if (isFALSE(private$config$use_branch_length)) {
        0
      } else {
        y_start + branch_length
      }
      children <- lapply(seq_along(clade$clade), function(child_idx) {
        private$species_pipes_layout(clade$clade[[child_idx]], child_idx, y, half_y_thickness)
      })
      y_shift <- 0
      if (length(children) == 0) {
        is_leaf <- TRUE
        x <- private$max_x + half_x_thickness
        private$max_x <<- x + half_x_thickness + private$config$x_padding
      } else {
        is_leaf <- FALSE
        x <- (children[[1]]$x + children[[1]]$half_x_thickness + children[[length(children)]]$x - children[[length(children)]]$half_x_thickness) / 2
        if (isFALSE(private$config$use_branch_length)) {
          y <- min(sapply(children, `[[`, "y")) - branch_length
        }
        if (private$config$use_y_shift) {
          y_shifts <- lapply(children, function(child) {
            child$y - child$half_y_thickness - child$branch_height - y - half_y_thickness
          })
          y_shift <- min(c(0, unlist(y_shifts)))
          # y-shifting a parent changes the branch heights of its children
          # when use_branch_length isFALSE, children's height must also be recomputed
        }
        children <- lapply(children, function(child) {
          child$side <- if (child$x <= x) "left" else "right"
          child$branch_height <- child$y - child$half_y_thickness - (y + half_y_thickness + y_shift)
          child
        })
      }
      branch_height <- if (private$config$use_y_shift) {
        max(min_branch_height, branch_length - half_y_thickness - parent_half_y_thickness)
      } else {
        branch_length - half_y_thickness - parent_half_y_thickness
      }
      l <- list(
        name = clade$name,
        child_idx = child_idx,
        side = "root",
        is_leaf = is_leaf,
        x = x,
        y = y,
        half_x_thickness = half_x_thickness,
        half_y_thickness = half_y_thickness,
        branch_height = branch_height,
        y_shift = y_shift,
        children = children
      )
      class(l) <- c("phyloXML_layout", class(l))
      l
    },
    inner_species_pipes_layout = function(splist) {
      if (isTRUE(private$warnings$negative_branch_height) && splist$branch_height + splist$y_shift < 0) {
        warning("Negative branch heights: the tree will be gibberish. Please increase `branch_length_scale`.", call. = FALSE)
        private$warnings$negative_branch_height <- FALSE
      }
      n_branch_events <- sum(unlist(private$internal_events[[splist$name]]) %in% c("duplication", "loss", "branchingOut", "transferBack"))
      branch_height_fraction <- splist$branch_height / (n_branch_events + 1)
      inner_h <- splist$y - splist$half_y_thickness + splist$y_shift - splist$branch_height + branch_height_fraction
      inner_x <- if (splist$side == "left") {
        splist$x - splist$half_x_thickness + 1
      } else {
        splist$x + splist$half_x_thickness - 1
      }
      inner_y <- if (isTRUE(splist$is_leaf)) {
        splist$y
      } else {
        splist$y - splist$half_y_thickness + splist$y_shift + 1
      }
      private$inner_layout[[splist$name]] = list(
        side = splist$side,
        branch_height_fraction = branch_height_fraction,
        inner_h = inner_h,
        inner_x = inner_x,
        inner_y = inner_y
      )
      lapply(splist$children, private$inner_species_pipes_layout)
      invisible(NULL)
    },
    gene_tree_layout = function(clade, child_idx = 0, lineage = "0") {
      # * Depth first is the only way to preserve the order. If a branch was at the top in one species node, it must be at the top in the descendant species as well. So it's also pre-order: we visit the node itself first, then the left child, then the right child.
      # * If it's a speciation, we put it in the first (most external) available x-y coord IN THE INTERNAL LEVEL (not in the branch height level). Then we increment x and y for the internal level, and we increment x for the branch height level.
      # * If it's a duplication, we put it in the first available x-y coord in the branch-height level, then we increment x and y for the branch height level, we don't touch the internal level, so the next event will fall exactly in place.
      # * If it's a leaf, we put it in the first availalbe x coord (y is same as the species leaf), and increment x-coord. we also increment x-coord in the branch-height level.
      # * so every time we increment the x-coord for the internal level, we do so for the branch too.
      # * If it's a transfer and we are in the same species, same as duplication.
      # * If it's a loss, we increment the x-coord of both internal and branch levels, we increment the y of the branch level.
      inner_layout <- private$inner_layout[[clade$event_location]]
      x <- inner_layout$inner_x
      y <- if (clade$event_type %in% c("speciation", "leaf")) {
        y <- inner_layout$inner_y
      } else if (clade$event_type %in% c("duplication", "branchingOut", "transferBack", "loss")) {
        y <- inner_layout$inner_h
      } else {
        stop("Unsupported event type.")
      }
      if (clade$event_type %in% c("speciation", "leaf", "loss")) {
        private$inner_layout[[clade$event_location]]$inner_x <- inner_layout$inner_x + 2*(inner_layout$side == "left") - 1
      }
      if (clade$event_type == "speciation") {
        private$inner_layout[[clade$event_location]]$inner_y <- inner_layout$inner_y + 1
      }
      if (clade$event_type %in% c("duplication", "branchingOut", "transferBack", "loss")) {
        private$inner_layout[[clade$event_location]]$inner_h <- inner_layout$inner_h + inner_layout$branch_height_fraction
      }
      children_permutation <- if (is.null(private$permutations_recGeneTree[[clade$name]])) {
        seq_along(clade$clade)
      } else {
        private$permutations_recGeneTree[[clade$name]]
      }
      children <- lapply(children_permutation, function(child_idx) {
        new_lineage <- if (clade$event_type != "duplication") {
          lineage
        } else if (child_idx == 1) {
          lineage
        } else {
          paste(lineage, clade$name, sep = ":-:")
        }
        private$gene_tree_layout(clade$clade[[child_idx]], child_idx, new_lineage)
      })
      l <- list(
        name = clade$name,
        child_idx = child_idx,
        lineage = lineage,
        is_leaf = clade$event_type == "leaf",
        x = x,
        y = y,
        layout_event = clade$event_type,
        children = children
      )
      class(l) <- c("phyloXML_layout", class(l))
      l
    }
  )
)
