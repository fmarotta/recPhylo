#' @export
PhylogenyLayout <- R6::R6Class("PhylogenyLayout",
  public = list(
    initialize = function(phylogeny, use_branch_length = TRUE, x_padding = 1, branch_length_scale = 1) {
      stopifnot("phyloXML_phylogeny" %in% class(phylogeny))
      private$.phylogeny <- phylogeny
      private$set_config(use_branch_length, x_padding, branch_length_scale)
      self$redraw()
    },
    redraw = function() {
      private$max_x <- 0
      private$layout_phylogeny <- private$simple_phylo_layout(private$.phylogeny$clade)
      invisible(self)
    },
    testplot = function() {
      if (!requireNamespace("ggplot2", quietly = T)) {
        stop("Please install `ggplot2` before using the plot() method.")
      }
      nodes <- self$layout$nodes
      edges <- self$layout$edges
      auxpoints <- data.frame(
        x = c(private$max_x, private$layout_phylogeny$x),
        y = c(max(nodes$y) + private$config$branch_length_scale, private$layout_phylogeny$y - private$config$branch_length_scale)
      )
      ggplot2::ggplot() +
        ggplot2::geom_line(data = edges, ggplot2::aes(x, y, group = group), lineend = "round") +
        ggplot2::geom_point(data = nodes, ggplot2::aes(x, y)) +
        ggplot2::geom_text(data = nodes, ggplot2::aes(x, y, label = name)) +
        # coord_polar() +
        # coord_flip() +
        # theme_void() +
        NULL
    }
  ),
  active = list(
    phylogeny = function(value) {
      if (!missing(value)) {
        stop("Cannot assign to phylogeny. Initialize a new object instead.")
      }
      private$.phylogeny
    },
    layout = function(value) {
      if (!missing(value)) {
        stop("Cannot assign to layout.")
      }
      list(
        nodes = merge_layout(private$.phylogeny, as.data.frame(private$layout_phylogeny)),
        edges = merge_layout(private$.phylogeny, get_phylogeny_edges_link(private$layout_phylogeny))
      )
    }
  ),
  private = list(
    .phylogeny = list(),
    config = list(),
    warnings = list(missing_branch_length = TRUE),
    max_x = 0,
    layout_phylogeny = list(),
    set_config = function(use_branch_length, x_padding, branch_length_scale) {
      if (! (is.numeric(use_branch_length) || is.logical(use_branch_length))) {
        stop("`use_branch_length` can only be a numeric or logical value.")
      } else if (length(use_branch_length) != 1) {
        stop("`use_branch_length` must have length one.")
      }
      private$config <- list(
        use_branch_length = use_branch_length,
        x_padding = x_padding,
        branch_length_scale = branch_length_scale
      )
    },
    simple_phylo_layout = function(clade, child_idx = 0, y_start = 0) {
      # Branch length
      if (is.numeric(private$config$use_branch_length)) {
        branch_length <- private$config$use_branch_length * private$config$branch_length_scale
        y <- y_start + branch_length
      } else if (isTRUE(private$config$use_branch_length)) {
        branch_length <- clade$branch_length
        if (is.na(branch_length)) {
          branch_length <- 1
          if (isTRUE(private$warnings$missing_branch_length)) {
            private$warnings$missing_branch_length <- FALSE
            warning("No branch length was found in clade '", clade$name, "'. Setting it to 1 automatically.", call. = FALSE)
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
      children <- lapply(seq_along(clade$clade), function(child_idx) {
        private$simple_phylo_layout(clade$clade[[child_idx]], child_idx, y)
      })
      if (length(children) == 0) {
        # It's an extant so we can have gene leaves but not speciations
        # The number of events is the x_thickness
        x <- private$max_x
        private$max_x <- x + private$config$x_padding
        is_leaf <- TRUE
      } else {
        x <- (children[[1]]$x + children[[length(children)]]$x) / 2
        is_leaf <- FALSE
        if (isFALSE(private$config$use_branch_length)) {
          y <- min(sapply(children, `[[`, "y")) - branch_length
        }
      }
      # Return the list
      l <- list(
        name = clade$name,
        child_idx = child_idx,
        is_leaf = is_leaf,
        x = x,
        y = y,
        branch_length = branch_length,
        children = children
      )
      class(l) <- c("phyloXML_layout", class(l))
      return(l)
    }
  )
)

#' @export
layout_phylogeny <- function(phylogeny, ...) {
  PhylogenyLayout$new(phylogeny, ...)$layout
}
