RecPhylo <- R6::R6Class("RecPhylo",
  public = list(
    spTree = list(),
    recGeneTree = list(),
    initialize = function(xml_file, x_padding = 1, branch_length_scale = 1) {
      recphylo_xml <- read_xml(xml_file)
      xml_ns_strip(recphylo_xml)
      sproot <- xml_find_first(recphylo_xml, "spTree//clade")
      self$spTree <- private$parse_sptree(sproot, recphylo_xml, x_padding = x_padding, branch_length_scale = branch_length_scale)
      groot <- xml_find_first(recphylo_xml, "recGeneTree//clade")
      self$recGeneTree <- private$parse_gtree(groot)
    }
  ),
  private = list(
    side = list(),
    n_internal_events = list(),
    branch_height = list(),
    intra_species_x = list(),
    intra_species_y = list(),
    intra_branch_y = list(),
    max_x = 0,
    parse_sptree = function(spnode, recphylo_xml,
                             x_padding = 1, branch_length_scale = 1,
                             parent = NA, side = "root", y_start = 0) {
      name <- xml_text(xml_find_first(spnode, "name"))
      branch_length <- as.numeric(xml_text(xml_find_first(spnode, "branch_length")))
      if (is.na(branch_length)) {
        branch_length <- 1
      }
      branch_length <- branch_length * branch_length_scale
      y_coord <- y_start + branch_length
      # Find all gene events that occur in this species
      events <- xml_find_all(recphylo_xml, paste0("recGeneTree//*[@speciesLocation='", name, "']"))
      x_thickness <- max(sum(xml_name(events) %in% c("speciation", "loss", "leaf")) + 1, 2)
      y_thickness <- max(sum(xml_name(events) == "speciation") + 1, 2)
      n_internal_events <- sum(xml_name(events) %in% c("duplication", "loss", "branchingOut"))
      children <- xml_find_all(spnode, "./clade")
      if (length(children) == 0) {
        # It's an extant so we can have gene leaves but not speciations
        # The number of events is the x_thickness
        left_child <- right_child <- NULL
        x_coord <- private$max_x + x_thickness / 2
        private$max_x <- x_coord + x_thickness / 2 + x_padding
        y_shift <- 0  # Leaves cannot be y-shifted
        is_leaf <- TRUE
      } else if (length(children) == 2) {
        left_child <- private$parse_sptree(
          children[[1]], recphylo_xml,
          parent = name, y_start = y_coord, side = "left",
          x_padding = x_padding, branch_length_scale = branch_length_scale)
        right_child <- private$parse_sptree(
          children[[2]], recphylo_xml,
          parent = name, y_start = y_coord, side = "right",
          x_padding = x_padding, branch_length_scale = branch_length_scale)
        x_coord <- (left_child$x + left_child$half_x_thickness + right_child$x - right_child$half_x_thickness) / 2
        y_shift <- min(
          (left_child$y - left_child$half_y_thickness - left_child$n_internal_events + left_child$y_shift) - (y_coord + y_thickness / 2),
          (right_child$y - right_child$half_y_thickness - right_child$n_internal_events + right_child$y_shift) - (y_coord + y_thickness / 2)
        )
        if (y_shift > 0) {
          y_shift <- 0
        }
        # We need to know the parent's y_thickness and y_shift before we can know the branch_height of the children, but the y_shift of the parent depends on that of the children
        left_child$branch_height <- left_child$y - y_coord - left_child$half_y_thickness - y_thickness / 2 - y_shift
        right_child$branch_height <- right_child$y - y_coord - right_child$half_y_thickness - y_thickness / 2 - y_shift
        private$branch_height[[left_child$name]] <- left_child$branch_height
        private$branch_height[[right_child$name]] <- right_child$branch_height
        private$intra_branch_y[[left_child$name]] <- left_child$y - left_child$half_y_thickness + left_child$branch_height * (1 / (left_child$n_internal_events + 1) - 1)
        private$intra_branch_y[[right_child$name]] <- right_child$y - right_child$half_y_thickness + right_child$branch_height * (1 / (right_child$n_internal_events + 1) - 1)
        is_leaf <- FALSE
      } else {
        stop("This species tree is not binary.")
      }
      private$side[[name]] <- side
      private$n_internal_events[[name]] <- n_internal_events
      private$branch_height[[name]] <- 0
      if (side == "left") {
        private$intra_species_x[[name]] <- x_coord - x_thickness / 2 + 1
      } else {
        private$intra_species_x[[name]] <- x_coord + x_thickness / 2 - 1
      }
      if (is_leaf == TRUE) {
        private$intra_species_y[[name]] <- y_coord
      } else {
        private$intra_species_y[[name]] <- y_coord - y_thickness / 2 + y_shift + 1
      }
      # Return the list
      splist <- list(
        name = name,
        parent = parent,
        side = side,
        is_leaf = is_leaf,
        x = x_coord,
        y = y_coord,
        half_x_thickness = x_thickness / 2,
        half_y_thickness = y_thickness / 2,
        n_internal_events = n_internal_events,
        branch_height = 0,
        y_shift = y_shift,
        left_child = left_child,
        right_child = right_child
      )
      class(splist) <- c("recphylo_spTree", class(splist))
      return(splist)
    },
    parse_gtree = function(gnode, parent = NA, side = "root", lineage = "r") {
      # * No. Depth first is the only way to preserve the order. If a branch was at the top in one species node, it must be at the top in the descendant species as well. So it's also pre-order: we visit the node itself first, then the left child, then the right child.
      # * If it's a speciation, we put it in the first (most external) available x-y coord IN THE INTERNAL LEVEL (not in the branch height level). Then we increment x and y for the internal level, and we increment x for the branch height level.
      # * If it's a duplication, we put it in the first available x-y coord in the branch-height level, then we increment x and y for the branch height level, we don't touch the internal level, so the next event will fall exactly in place.
      # * If it's a leaf, we put it in the first availalbe x coord (y is same as the species leaf), and increment x-coord. we also increment x-coord in the branch-height level.
      # * so every time we increment the x-coord for the internal level, we do so for the branch too.
      # * If it's a transfer and we are in the same species, same as duplication.
      # * If it's a loss, we increment the x-coord of both internal and branch levels, we increment the y of the branch level.
      name <- xml_text(xml_find_first(gnode, "name"))
      event_xml <- xml_find_first(gnode, "./eventsRec/*[self::leaf or self::duplication or self::loss or self::branchingOut or self::speciation]")
      event_type <- xml_name(event_xml)
      event_location <- xml_attr(event_xml, "speciesLocation")
      id <- paste(name, parent, sep = "@")
      if (event_type == "speciation") {
        x <- private$intra_species_x[[event_location]]
        y <- private$intra_species_y[[event_location]]
        private$increment_x(event_location)
        private$increment_splevel_y(event_location)
      } else if (event_type %in% c("duplication", "branchingOut")) {
        x <- private$intra_species_x[[event_location]]
        y <- private$intra_branch_y[[event_location]]
        private$increment_brlevel_y(event_location)
      } else if (event_type == "leaf") {
        x <- private$intra_species_x[[event_location]]
        y <- private$intra_species_y[[event_location]]
        private$increment_x(event_location)
      } else if (event_type == "loss") {
        x <- private$intra_species_x[[event_location]]
        y <- private$intra_branch_y[[event_location]]
        private$increment_x(event_location)
        private$increment_brlevel_y(event_location)
      }
      children <- xml_find_all(gnode, "./clade")
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
      return(glist)
    },
    increment_x = function(spname) {
      private$intra_species_x[[spname]] <- private$intra_species_x[[spname]] + 2*(private$side[[spname]] == "left")-1
      invisible(self)
    },
    increment_splevel_y = function(spname) {
      private$intra_species_y[[spname]] <- private$intra_species_y[[spname]] + 1
      invisible(self)
    },
    increment_brlevel_y = function(spname) {
      private$intra_branch_y[[spname]] <- private$intra_branch_y[[spname]] + private$branch_height[[spname]] / (private$n_internal_events[[spname]] + 1)
      invisible(self)
    }
  )
)

as.data.frame.recphylo_spTree <- function(l) {
  if (is.null(l)) {
    return(NULL)
  }
  fields <- c("name", "parent", "side", "is_leaf", "x", "y", "half_x_thickness", "half_y_thickness", "y_shift", "branch_height", "n_internal_events")
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

get_spedges = function(splist) {
  if (splist$is_leaf == TRUE) {
    data.frame(  # The pipe closed at one end
      group = paste(splist$name, "close"),
      x = c(splist$x - splist$half_x_thickness, splist$x - splist$half_x_thickness, splist$x + splist$half_x_thickness, splist$x + splist$half_x_thickness),
      y = c(splist$y - splist$half_y_thickness, splist$y, splist$y, splist$y - splist$half_y_thickness)
    )
  } else {
    rbind(
      data.frame(  # The inner 'U'
        group = paste(splist$name, "in"),
        x = c(splist$left_child$x + splist$left_child$half_x_thickness, splist$left_child$x + splist$left_child$half_x_thickness, splist$right_child$x - splist$right_child$half_x_thickness, splist$right_child$x - splist$right_child$half_x_thickness),
        y = c(splist$left_child$y - splist$left_child$half_y_thickness, splist$y + splist$half_y_thickness + splist$y_shift, splist$y + splist$half_y_thickness + splist$y_shift, splist$right_child$y - splist$right_child$half_y_thickness)
      ),
      data.frame(  # The left shoulder
        group = paste(splist$name, "outl"),
        x = c(splist$left_child$x - splist$left_child$half_x_thickness, splist$left_child$x - splist$left_child$half_x_thickness, splist$x - splist$half_x_thickness),
        y = c(splist$left_child$y - splist$left_child$half_y_thickness, splist$y - splist$half_y_thickness + splist$y_shift, splist$y - splist$half_y_thickness + splist$y_shift)
      ),
      data.frame(  # The right shoulder
        group = paste(splist$name, "outr"),
        x = c(splist$x + splist$half_x_thickness, splist$right_child$x + splist$right_child$half_x_thickness, splist$right_child$x + splist$right_child$half_x_thickness),
        y = c(splist$y - splist$half_y_thickness + splist$y_shift, splist$y - splist$half_y_thickness + splist$y_shift, splist$right_child$y - splist$right_child$half_y_thickness)
      ),
      get_spedges(splist$left_child),
      get_spedges(splist$right_child)
    )
  }
}

get_gedges <- function(glist, parent = NULL) {
  # each node is responsible to connect with its parent.
  if (is.null(glist)) {
    return(NULL)
  }
  if (is.null(parent)) {
    item <- data.frame(
      x = c(glist$x, glist$x),
      y = c(glist$y, glist$y - 2),
      group = glist$id,
      event_type = glist$event_type,
      side = glist$side,
      lineage = glist$lineage
    )
  } else {
    if (glist$event_type == "loss") {
      item <- rbind(
        data.frame(
          x = c(glist$x, glist$x),
          y = c(glist$y, parent$y),
          group = paste(glist$id, "v", sep = "_"),
          event_type = paste(glist$event_type, "v", sep = "_"),
          side = glist$side,
          lineage = glist$lineage
        ),
        data.frame(
          x = c(glist$x, parent$x),
          y = c(parent$y, parent$y),
          group = paste(glist$id, "h", sep = "_"),
          event_type = paste(glist$event_type, "h", sep = "_"),
          side = glist$side,
          lineage = glist$lineage
        )
      )
    } else if (parent$x > glist$x) {
      item <- data.frame(
        x = c(glist$x, glist$x, parent$x),
        y = c(glist$y, parent$y, parent$y),
        group = glist$id,
        event_type = glist$event_type,
        side = glist$side,
        lineage = glist$lineage
      )
    } else {
      item <- data.frame(
        x = c(parent$x, glist$x, glist$x),
        y = c(parent$y, parent$y, glist$y),
        group = glist$id,
        event_type = glist$event_type,
        side = glist$side,
        lineage = glist$lineage
      )
    }
  }
  rbind(
    item,
    get_gedges(glist$left_child, glist),
    get_gedges(glist$right_child, glist)
  )
}



ex1 <- RecPhylo$new("inst/extdata/example_1.recphyloxml", branch_length_scale = 10)
# ex1 <- RecPhylo$new("/g/scb/bork/marotta/prj/mycogenes_survey/notebooks/2023-11-05_rasmus_figure/recphylo/generax/selected_sequences/reconciliations/selected_sequences_reconciliated_fixed.xml", branch_length_scale = 60, x_padding = 10)
sptree <- as.data.frame(ex1$spTree)
gtree <- as.data.frame(ex1$recGeneTree)
spedges <- get_spedges(ex1$spTree)
gedges <- get_gedges(ex1$recGeneTree)
# ggplot() +
#   geom_line(data = spedges, aes(x, y, group = group)) +
#   # geom_point(data = sptree, aes(x, y, shape = event_type), color = "red", alpha = 0.4) +
#   geom_line(data = gedges, aes(x, y, group = group, color = gsub("l+$", "", lineage), linetype = event_type == "loss_v"), show.legend = F) +
#   # geom_text(data = gtree, aes(x, y, label = name), angle = 50, hjust = 0, nudge_y = 0.3) +
#   # geom_point(data = data.frame(x = 17, y = 10), aes(x, y)) +
#   # coord_polar() +
#   # theme_void() +
#   NULL

