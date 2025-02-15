GeomElbow <- ggproto("GeomElbow", GeomSegment,
  draw_panel = function(data, panel_params, coord) {
    # Create a data frame for each part
    upper_segment <- data
    upper_segment$yend <- data$y

    lower_segment <- data
    lower_segment$x <- data$xend

    # Use the draw_panel method of GeomSegment
    grid::gList(
      GeomSegment$draw_panel(upper_segment, panel_params, coord),
      GeomSegment$draw_panel(lower_segment, panel_params, coord)
    )
  }
)

geom_elbow <- function(mapping = NULL, data = NULL, stat = "identity",
                       position = "identity", ..., arrow = NULL,
                       arrow.fill = NULL, lineend = "butt", linejoin = "round",
                       na.rm = FALSE, show.legend = NA,
                       inherit.aes = TRUE) {
  layer(
    geom = GeomElbow, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(arrow = arrow, arrow.fill = arrow.fill, lineend = lineend, linejoin = linejoin, na.rm = na.rm, ...)
  )
}
