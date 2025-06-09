# TODO: have these functions return a list of phyloXML_layout class

get_spTree_edges <- function(splist, parent = NULL) {
  item <- if (splist$side == "left") {
    data.frame(
      name = splist$name,
      group = paste(parent$name, splist$name),
      side = splist$side,
      x = c(parent$x - parent$half_x_thickness, parent$x),
      y = c(parent$y - parent$half_y_thickness + parent$y_shift, parent$y + parent$half_y_thickness + parent$y_shift),
      xend = c(splist$x - splist$half_x_thickness, splist$x + splist$half_x_thickness),
      yend = c(splist$y - splist$half_y_thickness, splist$y - splist$half_y_thickness)
    )
  } else if (splist$side == "right") {
    data.frame(
      name = splist$name,
      group = paste(parent$name, splist$name),
      side = splist$side,
      x = c(parent$x + parent$half_x_thickness, parent$x),
      y = c(parent$y - parent$half_y_thickness + parent$y_shift, parent$y + parent$half_y_thickness + parent$y_shift),
      xend = c(splist$x + splist$half_x_thickness, splist$x - splist$half_x_thickness),
      yend = c(splist$y - splist$half_y_thickness, splist$y - splist$half_y_thickness)
    )
  } else if (splist$side == "root") {
    data.frame(
      name = splist$name,
      group = splist$name,
      side = splist$side,
      x = c(splist$x - splist$half_x_thickness, splist$x + splist$half_x_thickness),
      y = c(splist$y - splist$half_y_thickness + splist$y_shift, splist$y - splist$half_y_thickness + splist$y_shift),
      xend = c(splist$x - splist$half_x_thickness, splist$x + splist$half_x_thickness),
      yend = c(splist$y - splist$branch_height + splist$y_shift, splist$y - splist$branch_height + splist$y_shift)
    )
  }
  if (splist$is_leaf) {
    item$yend <- item$yend + splist$half_y_thickness
    item <- rbind(
      item,
      data.frame(
        name = splist$name,
        group = splist$name,
        side = "cap",
        x = splist$x - splist$half_x_thickness,
        y = splist$y,
        xend = splist$x + splist$half_x_thickness,
        yend = splist$y
      )
    )
  }
  Reduce(rbind, c(list(item), lapply(splist$children, get_spTree_edges, splist)))
}

get_recGene_edges <- function(glist, parent = NULL) {
  item <- data.frame(
    name = glist$name,
    group = paste(glist$name, "l", sep = "_"),
    side = if (is.null(parent)) "root" else if (glist$x < parent$x) "left" else "right",
    leg_type = "lateral",
    lineage = glist$lineage,
    x = parent$x %||% glist$x,
    xend = glist$x,
    y = parent$y %||% (glist$y - 1),
    yend = glist$y
  )
  return(
    Reduce(rbind, c(
      list(item),
      lapply(glist$children, get_recGene_edges, glist)
    ))
  )
}

get_phylogeny_edges_link <- function(cl, parent = NULL) {
  item <- if (is.null(parent)) {
    data.frame(
      name = cl$name,
      group = cl$name,
      side = cl$side,
      x = cl$x,
      xend = cl$x,
      y = cl$y - cl$branch_length,
      yend = cl$y
    )
  } else {
    data.frame(
      name = cl$name,
      group = cl$name,
      side = cl$side,
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

get_phylogeny_edges_comb <- function(cl, parent = NULL) {
  item <- data.frame(
    name = cl$name,
    group = cl$name,
    side = cl$side,
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
      side = cl$side,
      leg = "horizontal",
      x = children_range[1],
      xend = children_range[2],
      y = cl$y,
      yend = cl$y
    ),
    Reduce(rbind, lapply(cl$children, get_phylogeny_edges_comb, cl))
  )
}
