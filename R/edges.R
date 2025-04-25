# TODO: have these functions return a list of phyloXML_layout class

get_spTree_edges_elbow <- function(splist, parent = NULL) {
  item <- if (splist$is_leaf) {
    # For leaf nodes, we need to convert the 4 points into 3 line segments
    data.frame(
      name = rep(splist$name, 3),
      group = rep(paste(splist$name, "close"), 3),
      x = c(
        splist$x - splist$half_x_thickness,                  # Point 1 to 2
        splist$x - splist$half_x_thickness,                  # Point 2 to 3
        splist$x + splist$half_x_thickness                   # Point 3 to 4
      ),
      y = c(
        splist$y - splist$half_y_thickness - splist$branch_height,  # Point 1 to 2
        splist$y,                                            # Point 2 to 3
        splist$y                                             # Point 3 to 4
      ),
      xend = c(
        splist$x - splist$half_x_thickness,                  # Point 1 to 2
        splist$x + splist$half_x_thickness,                  # Point 2 to 3
        splist$x + splist$half_x_thickness                   # Point 3 to 4
      ),
      yend = c(
        splist$y,                                            # Point 1 to 2
        splist$y,                                            # Point 2 to 3
        splist$y - splist$half_y_thickness - splist$branch_height # Point 3 to 4
      )
    )
  } else {
    rbind(
      # Inner left segment (2 points = 1 segment)
      data.frame(
        name = splist$name,
        group = paste(splist$name, "inl"),
        x = splist$children[[1]]$x + splist$children[[1]]$half_x_thickness,
        y = splist$y + splist$half_y_thickness + splist$y_shift,
        xend = splist$x,
        yend = splist$y + splist$half_y_thickness + splist$y_shift
      ),
      # Inner right segment (2 points = 1 segment)
      data.frame(
        name = splist$name,
        group = paste(splist$name, "inr"),
        x = splist$x,
        y = splist$y + splist$half_y_thickness + splist$y_shift,
        xend = splist$children[[length(splist$children)]]$x - splist$children[[length(splist$children)]]$half_x_thickness,
        yend = splist$y + splist$half_y_thickness + splist$y_shift
      ),
      # Outer left segment (3 points = 2 segments)
      data.frame(
        name = rep(splist$name, 2),
        group = rep(paste(splist$name, "outl"), 2),
        x = c(
          splist$children[[1]]$x - splist$children[[1]]$half_x_thickness,  # Point 1 to 2
          splist$children[[1]]$x - splist$children[[1]]$half_x_thickness   # Point 2 to 3
        ),
        y = c(
          splist$y + splist$half_y_thickness + splist$y_shift,             # Point 1 to 2
          splist$y - splist$half_y_thickness + splist$y_shift              # Point 2 to 3
        ),
        xend = c(
          splist$children[[1]]$x - splist$children[[1]]$half_x_thickness,  # Point 1 to 2
          splist$x - splist$half_x_thickness                               # Point 2 to 3
        ),
        yend = c(
          splist$y - splist$half_y_thickness + splist$y_shift,             # Point 1 to 2
          splist$y - splist$half_y_thickness + splist$y_shift              # Point 2 to 3
        )
      ),
      # Outer right segment (3 points = 2 segments)
      data.frame(
        name = rep(splist$name, 2),
        group = rep(paste(splist$name, "outr"), 2),
        x = c(
          splist$x + splist$half_x_thickness,                              # Point 1 to 2
          splist$children[[length(splist$children)]]$x + splist$children[[length(splist$children)]]$half_x_thickness  # Point 2 to 3
        ),
        y = c(
          splist$y - splist$half_y_thickness + splist$y_shift,             # Point 1 to 2
          splist$y - splist$half_y_thickness + splist$y_shift              # Point 2 to 3
        ),
        xend = c(
          splist$children[[length(splist$children)]]$x + splist$children[[length(splist$children)]]$half_x_thickness, # Point 1 to 2
          splist$children[[length(splist$children)]]$x + splist$children[[length(splist$children)]]$half_x_thickness  # Point 2 to 3
        ),
        yend = c(
          splist$y - splist$half_y_thickness + splist$y_shift,             # Point 1 to 2
          splist$y + splist$half_y_thickness + splist$y_shift              # Point 2 to 3
        )
      ),
      # Left stem segment (2 points = 1 segment)
      data.frame(
        name = splist$name,
        group = paste(splist$name, "steml"),
        x = splist$x - splist$half_x_thickness,
        y = splist$y - splist$half_y_thickness - splist$branch_height,
        xend = splist$x - splist$half_x_thickness,
        yend = splist$y - splist$half_y_thickness
      ),
      # Right stem segment (2 points = 1 segment)
      data.frame(
        name = splist$name,
        group = paste(splist$name, "stemr"),
        x = splist$x + splist$half_x_thickness,
        y = splist$y - splist$half_y_thickness,
        xend = splist$x + splist$half_x_thickness,
        yend = splist$y - splist$half_y_thickness - splist$branch_height
      )
    )
  }
  # Recursively process child nodes and combine results
  Reduce(rbind, c(list(item), lapply(splist$children, get_spTree_edges_elbow)))
}

get_recGene_edge_root <- function(glist, branch_length = 2) {
  data.frame(
    name = glist$name,
    group = glist$name,
    leg_type = "vertical",
    lineage = glist$lineage,
    x = glist$x,
    xend = glist$x,
    y = glist$y - branch_length,
    yend = glist$y
  )
}

get_recGene_edge_elbow <- function(glist, parent) {
  rbind(
    data.frame(
      name = glist$name,
      group = paste(glist$name, "h", sep = "_"),
      leg_type = "horizontal",
      lineage = glist$lineage,
      x =  parent$x,
      xend = glist$x,
      y = parent$y,
      yend = parent$y
    ),
    data.frame(
      name = glist$name,
      group = paste(glist$name, "v", sep = "_"),
      leg_type = "vertical",
      lineage = glist$lineage,
      x = glist$x,
      xend = glist$x,
      y = parent$y,
      yend = glist$y
    )
  )
}

get_recGene_edge_link <- function(glist, parent) {
  data.frame(
    name = glist$name,
    group = paste(glist$name, "l", sep = "_"),
    leg_type = "lateral",
    lineage = glist$lineage,
    x =  parent$x,
    xend = glist$x,
    y = parent$y,
    yend = glist$y
  )
}

get_recGene_edges <- function(glist, parent = NULL) {
  # each node is responsible to connect with its parent.
  item <- if (is.null(parent)) {
    get_recGene_edge_root(glist)
  } else if (glist$layout_event == "transferBack") {
    get_recGene_edge_link(glist, parent)
  } else if (glist$layout_event == "bifurcationOut") {
    get_recGene_edge_link(glist, parent)
  } else {
    get_recGene_edge_elbow(glist, parent)
  }
  Reduce(rbind, c(
    list(item),
    lapply(glist$children, get_recGene_edges, glist)
  ))
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
