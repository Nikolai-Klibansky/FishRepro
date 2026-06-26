#' Plateau model with four parameters
#'
#' Plateau model with four parameters (\eqn{y = a/(1+exp(-(x-(m-0.5))/s))-a/(1+exp(-(x-((m-0.5)+(d+1)))/s))}). A plateau model is type of double logistic model where the slope of the first (left) logistic function is positive and the slope of the second (right) logistic function is negative, resulting in a domed shape. If the slopes are steep enough and the points of inflection are far enough apart, then the model will have a flat-topped section between the points of inflection.
#' @param x numeric vector
#' @param m point of inflection for left logistic function
#' @param s scale parameter shared by both logistic functions
#' @param d difference between left and right points of inflection. When fit to binary event data, average duration of events.
#' @param a asymptote of each logistic function
#' @author Nikolai Klibansky
#' @export
#' @examples
#' # Example

pla_msda <- function(x,m,s,d,a) {
  a/(1+exp(-(x-(m-0.5))/s))-a/(1+exp(-(x-((m-0.5)+(d+1)))/s))}
