mystars <- function(coords, params)
{
  data <- params("vertex","data")
  cP <- params("vertex","cP")
  scale <- params("vertex","scale")
  size <- params("vertex","size")
  graphics::stars(data, locations = coords, labels=NULL, scale = scale, len = size, col.segments = cP, draw.segments = TRUE,
                  mar=c(0,0,0,0), add=TRUE, inches=FALSE)
  graphics::symbols(x = coords[,1], y=coords[,2], circles = size, inches=FALSE, bg="transparent", bty="n", add=TRUE)
}
