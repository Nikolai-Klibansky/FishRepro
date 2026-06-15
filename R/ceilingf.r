#' Compute the greatest multiple of a factor f which is greater than or equal to a value x
#'
#' Compute the greatest multiple of a factor f which is greater than or equal to a value x
#' @param x numeric vector
#' @param f factor
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Example
#'}
#'
ceilingf <- function(x,f) {ceiling(x/f)*f}
