# TODO: have these functions return a list of phyloXML_layout class

get_spTree_edges_link <- function(splist, parent = NULL) {
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
  Reduce(rbind, c(list(item), lapply(splist$children, get_spTree_edges_link)))
}

get_spTree_edges <- function(splist, parent = NULL) {
  if (splist$is_leaf) {
    return(
      data.frame(
        name = splist$name,
        group = splist$name,
        x = c(splist$x - splist$half_x_thickness, splist$x - splist$half_x_thickness, splist$x + splist$half_x_thickness),
        y = c(splist$y, splist$y, splist$y),
        xend = c(splist$x - splist$half_x_thickness, splist$x + splist$half_x_thickness, splist$x + splist$half_x_thickness),
        yend = c(splist$y - splist$half_y_thickness, splist$y, splist$y - splist$half_y_thickness)
      )
    )
  }
  items <- lapply(splist$children, function(child) {
    if (child$side == "left") {
      data.frame(
        name = child$name,
        group = paste(splist$name, child$name),
        x = c(splist$x - splist$half_x_thickness, splist$x),
        y = c(splist$y - splist$half_y_thickness, splist$y + splist$half_y_thickness),
        xend = c(child$x - child$half_x_thickness, child$x + child$half_x_thickness),
        yend = c(child$y - child$half_y_thickness, child$y - child$half_y_thickness)
      )
    } else if (child$side == "right") {
      data.frame(
        name = child$name,
        group = paste(splist$name, child$name),
        x = c(splist$x + splist$half_x_thickness, splist$x),
        y = c(splist$y - splist$half_y_thickness, splist$y + splist$half_y_thickness),
        xend = c(child$x + child$half_x_thickness, child$x - child$half_x_thickness),
        yend = c(child$y - child$half_y_thickness, child$y - child$half_y_thickness)
      )
    }
  })
  if (splist$side == "root") {
    items <- c(
      items,
      list(
        data.frame(
          name = splist$name,
          group = splist$name,
          x = c(splist$x - splist$half_x_thickness, splist$x + splist$half_x_thickness),
          y = c(splist$y - splist$half_y_thickness, splist$y - splist$half_y_thickness),
          xend = c(splist$x - splist$half_x_thickness, splist$x + splist$half_x_thickness),
          yend = c(splist$y - splist$branch_height, splist$y - splist$branch_height)
        )
      )
    )
  }
  Reduce(rbind, c(items, lapply(splist$children, get_spTree_edges)))
}
# ggplot(get_spTree_edges(splist), aes(x, y, xend = xend, yend = yend)) + geom_diagonal(orientation = "y")

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
    x = parent$x %||% glist$x,
    xend = glist$x,
    y = parent$y %||% (glist$y - 1),
    yend = glist$y
  )
}

get_recGene_edges <- function(glist, parent = NULL) {
  item <- get_recGene_edge_link(glist, parent)
  return(
    Reduce(rbind, c(
      list(item),
      lapply(glist$children, get_recGene_edges, glist)
    ))
  )
  # lol trollolol
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
