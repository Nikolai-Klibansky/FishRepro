#' Compute diameter of a circle from its area
#'
#' Compute diameter of a circle from its area
#' @param A area of a circle
#' @author Nikolai Klibansky
#' @export
#' @examples
#' # Example
#' x <- rnorm(10)
#' y <- rnorm(10)
#' n <- sample(1:100,size=10,replace=TRUE)
#' cex <- 3
#' cex.scaled <- a2d_circ(cex*(n/max(n))) # Scale cex (diameter) of points so area is proportional to n
#' plot(x,y,cex=cex)
#' plot(x,y,cex=cex.scaled)
#'

a2d_circ <- function(A){2*sqrt(A/pi)}
