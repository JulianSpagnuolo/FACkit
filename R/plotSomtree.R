plot.somtree <- function(somtree, type=c("stars","pies"), var.data, col.pal, layout)
{
  if(somtree$algorithm == "GrowSOM")
  {
    if(type =="stars")
    {
      add.vertex.shape("star", clip = igraph:::igraph.shape.noclip, plot = mystars,
                       parameters = list(vertex.data = NULL, vertex.cP = col.pal, vertex.scale = T,
                                         vertex.size=NULL))
      oldpar <- graphics::par(no.readonly = TRUE)
      graphics::par(mar = c(1, 1, 1, 1))
      graphics::layout(matrix(c(1, 2), 1, 2, byrow = TRUE), widths = c(1, 2), heights = c(1))
      starLegend(colnames(var.data), col.pal)
      plot.igraph(somtree$mst,vertex.shape = "star", vertex.label = NA,
                  vertex.size = somtree$map$nodes$size, vertex.data = var.data, vertex.cP = col.pal,
                  vertex.scale = T, layout=layout, edge.lty = 1)
      graphics::par(oldpar)
      graphics::layout(1)
    }

    if(type == "pies")
    {
      t <- table(factor(somtree$map$mapped$bmn, levels = seq_along(somtree$map$nodes$size)),
                 factor(var.data, levels = c(levels(var.data), "empty")))
      t[rowSums(t) == 0, "empty"] <- 1
      data <- unlist(apply(t, 1, list), recursive = FALSE)
      colors <- list(c(col.pal, "#000000"))
      oldpar <- graphics::par(no.readonly = TRUE)
      graphics::par(mar = c(1, 1, 1, 1))
      graphics::layout(matrix(c(1, 2), 1, 2, byrow = TRUE), widths = c(1, 2), heights = c(1))
      graphics::plot.new()
      legend("center", legend = levels(var.data), fill = colors[[1]], cex = 0.7, ncol = 1, bty = "n")
      igraph::plot.igraph(somtree$mst, vertex.shape = "pie", vertex.label = NA, vertex.size = somtree$map$nodes$size, vertex.pie = data, vertex.pie.color = colors, layout = layout, edge.lty = 1)
      graphics::par(oldpar)
      graphics::layout(1)
    }
  }
  if(somtree$algorithm == "kohonen")
  {
    if(type =="stars")
    {
      add.vertex.shape("star", clip = igraph:::igraph.shape.noclip, plot = mystars,
                       parameters = list(vertex.data = NULL, vertex.cP = col.pal, vertex.scale = T,
                                         vertex.size=NULL))
      oldpar <- graphics::par(no.readonly = TRUE)
      graphics::par(mar = c(1, 1, 1, 1))
      graphics::layout(matrix(c(1, 2), 1, 2, byrow = TRUE), widths = c(1, 2), heights = c(1))
      starLegend(colnames(var.data), col.pal)
      plot.igraph(somtree$mst,vertex.shape = "star", vertex.label = NA,
                  vertex.size = somtree$map$nodes$size, vertex.data = var.data, vertex.cP = col.pal,
                  vertex.scale = T, layout=layout, edge.lty = 1)
      graphics::par(oldpar)
      graphics::layout(1)
    }

    if(type == "pies")
    {
      t <- table(factor(somtree$map$unit.classif, levels = seq_along(somtree$map$nodes$size)),
                 factor(var.data, levels = c(levels(var.data), "empty")))
      t[rowSums(t) == 0, "empty"] <- 1
      data <- unlist(apply(t, 1, list), recursive = FALSE)
      colors <- list(c(col.pal, "#000000"))
      oldpar <- graphics::par(no.readonly = TRUE)
      graphics::par(mar = c(1, 1, 1, 1))
      graphics::layout(matrix(c(1, 2), 1, 2, byrow = TRUE), widths = c(1, 2), heights = c(1))
      graphics::plot.new()
      legend("center", legend = levels(var.data), fill = colors[[1]], cex = 0.7, ncol = 1, bty = "n")
      igraph::plot.igraph(somtree$mst, vertex.shape = "pie", vertex.label = NA, vertex.size = 5, vertex.pie = data, vertex.pie.color = colors, layout = layout, edge.lty = 1)
      graphics::par(oldpar)
      graphics::layout(1)
    }
  }
}
