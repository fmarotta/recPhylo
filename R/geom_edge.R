GeomElbow <- ggplot2::ggproto("GeomElbow", ggplot2::GeomSegment,
  draw_panel = function(data, panel_params, coord, direction = "hv",
                        arrow = NULL, arrow.fill = NULL, lineend = "butt",
                        linejoin = "round", na.rm = FALSE) {
    # Create a data frame for each part
    first_segment <- data
    second_segment <- data
    if (direction == "hv") {
      first_segment$yend <- data$y
      second_segment$x <- data$xend
    } else if (direction == "vh") {
      first_segment$xend <- data$x
      second_segment$y <- data$yend
    } else {
      stop("`direction` must be either 'hv' or 'vh', not '", direction, "'.")
    }

    # Use the draw_panel method of GeomSegment
    grid::gList(
      ggplot2::GeomSegment$draw_panel(first_segment, panel_params, coord, arrow = NULL, lineend = lineend, linejoin = linejoin, na.rm = na.rm),
      ggplot2::GeomSegment$draw_panel(second_segment, panel_params, coord, arrow = arrow, arrow.fill = arrow.fill, lineend = lineend, linejoin = linejoin, na.rm = na.rm)
    )
  }
)

#' geom_elbow: Draw Elbow Segments
#'
#' `geom_elbow()` creates elbow-shaped segments connecting two points by
#' introducing a new intermediate point at the bend of the elbow.
#'
#' @inherit ggplot2::geom_segment
#'
#' @param direction Direction of stairs: 'vh' for vertical then horizontal,
#'     'hv' for horizontal then vertical.
#'
#' @examples
#' library(ggplot2)
#' df <- data.frame(x1 = 2.62, x2 = 3.57, y1 = 21.0, y2 = 15.0)
#' ggplot(df) +
#'   geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment")) +
#'   geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2, colour = "curve")) +
#'   geom_elbow(aes(x = x1, y = y1, xend = x2, yend = y2, colour = "elbow"))
#'
#' ggplot(df) +
#'   geom_elbow(
#'     aes(x = x1, y = y1, xend = x2, yend = y2),
#'     direction = "vh",
#'     arrow = arrow(type = "closed"),
#'     arrow.fill = "blue"
#'   )
#'
#' @import ggplot2
#'
#' @export
geom_elbow <- function(mapping = NULL, data = NULL, stat = "identity",
                       position = "identity", ..., direction = "hv", arrow = NULL,
                       arrow.fill = NULL, lineend = "butt", linejoin = "round",
                       na.rm = FALSE, show.legend = NA,
                       inherit.aes = TRUE) {
  ggplot2::layer(
    geom = GeomElbow, mapping = mapping, data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(direction = direction, arrow = arrow, arrow.fill = arrow.fill, lineend = lineend, linejoin = linejoin, na.rm = na.rm, ...)
  )
}
