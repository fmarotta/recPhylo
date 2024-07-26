# TODO: have these functions return a list of phyloXML_layout class

get_species_pipes_edges <- function(splist) {
  item <- if (splist$is_leaf) {
    data.frame(
      name = splist$name,
      group = paste(splist$name, "close"),
      x = c(splist$x - splist$half_x_thickness, splist$x - splist$half_x_thickness, splist$x + splist$half_x_thickness, splist$x + splist$half_x_thickness),
      y = c(splist$y - splist$half_y_thickness - splist$branch_height, splist$y, splist$y, splist$y - splist$half_y_thickness - splist$branch_height)
    )
  } else {
    rbind(
      data.frame(
        name = splist$name,
        group = paste(splist$name, "inl"),
        x = c(splist$children[[1]]$x + splist$children[[1]]$half_x_thickness, splist$x),
        y = splist$y + splist$half_y_thickness + splist$y_shift
      ),
      data.frame(
        name = splist$name,
        group = paste(splist$name, "inr"),
        x = c(splist$x, splist$children[[length(splist$children)]]$x - splist$children[[length(splist$children)]]$half_x_thickness),
        y = splist$y + splist$half_y_thickness + splist$y_shift
      ),
      data.frame(
        name = splist$name,
        group = paste(splist$name, "outl"),
        x = c(splist$children[[1]]$x - splist$children[[1]]$half_x_thickness, splist$children[[1]]$x - splist$children[[1]]$half_x_thickness, splist$x - splist$half_x_thickness),
        y = c(splist$y + splist$half_y_thickness + splist$y_shift, splist$y - splist$half_y_thickness + splist$y_shift, splist$y - splist$half_y_thickness + splist$y_shift)
      ),
      data.frame(
        name = splist$name,
        group = paste(splist$name, "outr"),
        x = c(splist$x + splist$half_x_thickness, splist$children[[length(splist$children)]]$x + splist$children[[length(splist$children)]]$half_x_thickness, splist$children[[length(splist$children)]]$x + splist$children[[length(splist$children)]]$half_x_thickness),
        y = c(splist$y - splist$half_y_thickness + splist$y_shift, splist$y - splist$half_y_thickness + splist$y_shift, splist$y + splist$half_y_thickness + splist$y_shift)
      ),
      data.frame(
        name = splist$name,
        group = paste(splist$name, "steml"),
        x = c(splist$x - splist$half_x_thickness, splist$x - splist$half_x_thickness),
        y = c(splist$y - splist$half_y_thickness - splist$branch_height, splist$y - splist$half_y_thickness)
      ),
      data.frame(
        name = splist$name,
        group = paste(splist$name, "stemr"),
        x = c(splist$x + splist$half_x_thickness, splist$x + splist$half_x_thickness),
        y = c(splist$y - splist$half_y_thickness, splist$y - splist$half_y_thickness - splist$branch_height)
      )
    )
  }
  Reduce(rbind, c(list(item), lapply(splist$children, get_species_pipes_edges)))
}

get_gene_tree_edges <- function(glist, parent = NULL) {
  # each node is responsible to connect with its parent.
  if (is.null(parent)) {
    item <- data.frame(
      name = glist$name,
      group = glist$name,
      leg_type = "vertical",
      lineage = glist$lineage,
      x = c(glist$x, glist$x),
      y = c(glist$y, glist$y - 2)
    )
  } else if (glist$layout_event == "loss") {
    item <- rbind(
      data.frame(
        name = glist$name,
        group = paste(glist$name, "v", sep = "_"),
        leg_type = "loss_vertical",
        lineage = glist$lineage,
        x = c(glist$x, glist$x),
        y = c(glist$y, parent$y)
      ),
      data.frame(
        name = glist$name,
        group = paste(glist$name, "h", sep = "_"),
        leg_type = "loss_horizontal",
        lineage = glist$lineage,
        x = c(glist$x, parent$x),
        y = c(parent$y, parent$y)
      )
    )
  } else if (glist$layout_event == "transferBack") {
    item <- data.frame(
      name = glist$name,
      group = glist$name,
      leg_type = "transferBack",
      lineage = glist$lineage,
      x = c(glist$x, parent$x),
      y = c(glist$y, parent$y)
    )
  } else if (parent$x >= glist$x) {
    item <- data.frame(
      name = glist$name,
      group = glist$name,
      leg_type = "elbow",
      lineage = glist$lineage,
      x = c(glist$x, glist$x, parent$x),
      y = c(glist$y, parent$y, parent$y)
    )
  } else if (parent$x < glist$x) {
    item <- data.frame(
      name = glist$name,
      group = glist$name,
      leg_type = "elbow",
      lineage = glist$lineage,
      x = c(parent$x, glist$x, glist$x),
      y = c(parent$y, parent$y, glist$y)
    )
  } else {
    stop("Unexpected node type: ", glist)
  }
  Reduce(rbind, c(list(item), lapply(glist$children, get_gene_tree_edges, glist)))
}

get_phylogeny_edges_comb <- function(cl, parent = NULL) {
  item <- data.frame(
    name = cl$name,
    group = cl$name,
    leg = "vertical",
    x = cl$x,
    xend = cl$x,
    y = cl$y - cl$branch_length,
    yend = cl$y
  )
  if (length(cl$children) == 0) {
    return(item)
  }
  children_range <- range(sapply(cl$children, `[[`, "x"))
  rbind(
    item,
    data.frame(
      name = cl$name,
      group = paste(cl$name, "brace", sep = "@"),
      leg = "horizontal",
      x = children_range[1],
      xend = children_range[2],
      y = cl$y,
      yend = cl$y
    ),
    Reduce(rbind, lapply(cl$children, get_phylogeny_edges_comb, cl))
  )
}

get_phylogeny_edges_link <- function(cl, parent = NULL) {
  item <- if (is.null(parent)) {
    data.frame(
      name = cl$name,
      group = cl$name,
      x = cl$x,
      xend = cl$x,
      y = cl$y - cl$branch_length,
      yend = cl$y
    )
  } else {
    data.frame(
      name = cl$name,
      group = cl$name,
      x = parent$x,
      xend = cl$x,
      y = parent$y,
      yend = cl$y
    )
  }
  if (length(cl$children) == 0) {
    return(item)
  }
  rbind(
    item,
    Reduce(rbind, lapply(cl$children, get_phylogeny_edges_link, cl))
  )
}
