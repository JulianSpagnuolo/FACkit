starLegend <- function (labels, colors = grDevices::rainbow(length(labels)), main = "")
{
  graphics::plot(1, type = "n", xlab = "", ylab = "", xlim = c(-10, 10), ylim = c(-3, 3),
                 asp = 1, bty = "n", xaxt = "n", yaxt = "n", main = main)
  graphics::stars(matrix(c(1:(2 * length(labels))), nrow = 2), col.segments = colors,
                  locations = c(0, 0), draw.segments = TRUE, add = TRUE, inches = FALSE)
  n <- length(labels)
  angle <- 2 * pi/n
  angles <- seq(angle/2, 2 * pi, by = angle)
  left <- (angles > (pi/2) & angles < (3 * pi/2))
  x <- c(2, -2)[left + 1]
  y_tmp <- c(seq(-2, 2, by = 4/(sum(!left) + 1))[-c(1, sum(!left) + 2)],
             seq(2, -2, by = -4/(sum(left) + 1))[-c(1, sum(left) + 2)])
  y <- shiftFunction(y_tmp, max((cummax(y_tmp) < 0) * seq_along(y_tmp)))
  for (i in seq_along(labels))
    {
    graphics::text(x = x[i], y = y[i], labels = labels[i], adj = c(as.numeric(left)[i], 0.5), cex = 0.5)
    graphics::lines(x = c(x[i] + c(-0.2, 0.2)[left[i] + 1], c(1.5, -1.5)[left[i] + 1], cos(angles[i])),
                    y = c(y[i], y[i], sin(angles[i])), col = colors[i], lwd = 2)
    }
}
