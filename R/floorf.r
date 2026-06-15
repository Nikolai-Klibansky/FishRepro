#' Compute the greatest multiple of a factor f which is less than or equal to a value x
#'
#' Compute the greatest multiple of a factor f which is less than or equal to a value x
#' @param x numeric vector
#' @param f factor
#' @author Nikolai Klibansky
#' @export
#' @examples
#' \dontrun{
#' # Example
#'}
#'
floorf <- function(x,f) {floor(x/f)*f}
