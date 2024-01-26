
get_spedges <- function(splist) {
  if (splist$is_leaf == TRUE) {
    data.frame(  # The pipe closed at one end
      group = paste(splist$name, "close"),
      x = c(splist$x - splist$half_x_thickness, splist$x - splist$half_x_thickness, splist$x + splist$half_x_thickness, splist$x + splist$half_x_thickness),
      y = c(splist$y - splist$half_y_thickness, splist$y, splist$y, splist$y - splist$half_y_thickness)
    )
  } else {
    rbind(
      if (is.na(splist$parent)) {
        data.frame(  # The root's pipe
          group = rep(c(paste(splist$name, "rootl"), paste(splist$name, "rootr")), each = 2),
          x = c(splist$x - splist$half_x_thickness, splist$x - splist$half_x_thickness, splist$x + splist$half_x_thickness, splist$x + splist$half_x_thickness),
          y = c(splist$y - splist$half_y_thickness + splist$y_shift, splist$y - splist$half_y_thickness + splist$y_shift - 2, splist$y - splist$half_y_thickness + splist$y_shift, splist$y - splist$half_y_thickness + splist$y_shift - 2)
        )
      } else {
        NULL
      },
      data.frame(  # The inner 'U' (left part)
        group = paste(splist$name, "inl"),
        x = c(splist$left_child$x + splist$left_child$half_x_thickness, splist$left_child$x + splist$left_child$half_x_thickness, splist$x),
        y = c(splist$left_child$y - splist$left_child$half_y_thickness, splist$y + splist$half_y_thickness + splist$y_shift, splist$y + splist$half_y_thickness + splist$y_shift)
      ),
      data.frame(  # The inner 'U' (right part)
        group = paste(splist$name, "inr"),
        x = c(splist$x, splist$right_child$x - splist$right_child$half_x_thickness, splist$right_child$x - splist$right_child$half_x_thickness),
        y = c(splist$y + splist$half_y_thickness + splist$y_shift, splist$y + splist$half_y_thickness + splist$y_shift, splist$right_child$y - splist$right_child$half_y_thickness)
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
