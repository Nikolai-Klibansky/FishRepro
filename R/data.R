#' Red snapper batch fecundity-length data
#'
#' Batch fecundity, length, and related biological data, provided to the 2017 (SEDAR 41) stock assessment of the Southeast US Atlantic stock of red snapper.
#'
#' @format ## `data_fb_L`
#' \describe{
#'   \item{Source}{Agency which collected the data (FWRI = Florida Fish and Wildlife Research Institute, SCDNR = South Carolina Department of Natural Resources)}
#'   \item{Date}{Date (MM/DD/YYYY) of data collection}
#'   \item{Month}{Month of data collection}
#'   \item{TL}{Fish total length (mm)}
#'   \item{TW}{Fish total weight (g)}
#'   \item{SW}{Fish somatic weight (g)}
#'   \item{Age}{Fish age (yr)}
#'   \item{fb}{Fish batch fecundity (eggs)}
#' }
#' @source <https://sedarweb.org/documents/s41dw49-estimates-of-reproductive-activity-in-red-snapper-by-size-season-and-time-of-day-with-nonlinear-models/>
"data_fb_L"

#' Red snapper spawning indicator-size data
#'
#' Red snapper spawning indicator-size data, from an earlier draft of: Porch, C. E., Fitzhugh, G. R., Lang, E. T., Lyon, H. M., & Linton, B. C. (2015). Estimating the dependence of spawning frequency on size and age in Gulf of Mexico Red Snapper. Marine and Coastal Fisheries, 7(1), 233-245. Similar to Table 1 of the published manuscript.
#'
#' @format ## `data_SI_L`
#' \describe{
#'   \item{Age}{Age of females}
#'   \item{N.f}{Number of females sampled}
#'   \item{P.SI}{Proportion of females with spawing indicators}
#'   \item{P.SI.adj}{Proportion of females with spawing indicators, adjusted by multiplying by a correction factor (CF) based on a spawning indicator duration of 34 hours \eqn{CF = (10+24)/24} }
#'   \item{N.f.SI}{Number of females sampled with spawning indicators}
#'   \item{L}{Fish length (mm) estimated from ages. \eqn{L = 848.4797*(1-exp(-0.2188*(Age-(-0.0611))))]}. SEDAR 31 2.9.1 pg 32 Eq 4. converted to mm and g by B0'=B0*(0.0393701^B1)/0.00220462}
#' }
#' @source <https://academic.oup.com/mcf/article/7/1/233/7827111?guestAccessKey=>
"data_SI_L"

#' Variance-covariance matrix for parameters in four parameter logistic model fit to red snapper spawning indicator-size data.
#'
#' Covariance matrix for parameters in four parameter logistic model fit \eqn{y=f+(a-f)/(1+exp(-(x-m)/s))} to red snapper spawning indicator-size data stored in data_SI_L. The logistic model was fit with with \code{\link[stats]{optim}}, minimizing the negative log-likelihood. Running \code{\link[base]{solve}} on the hessian matrix returns the variance-covariance matrix among the model parameters.
#'
#' @format ## `covmat_SI_L`
#' \describe{
#'   \item{m}{Parameter m. point of inflection}
#'   \item{s}{Parameter s. scale parameter}
#'   \item{a}{Parameter a. upper asymptote}
#'   \item{f}{Parameter f. lower asymptote}
#' }
#' @source <https://academic.oup.com/mcf/article/7/1/233/7827111?guestAccessKey=>
"covmat_SI_L"

#' Seasonal patterns in red snapper spawning activity
#'
#' Data shows the proportions of female red snapper with spawning indicators by time of year, from Table 1 of: Fitzhugh, G.R., E.T. Lang, and H. Lyon. 2012. Expanded Annual Stock Assessment Survey 2011: Red Snapper Reproduction. SEDAR31-DW07. SEDAR, North Charleston, SC. 33 pp.
#'
#' @format ## `data_SI_d`
#' \describe{
#'   \item{Period}{Numeric ID of sampling period, in sequence.}
#'   \item{Start.date}{Start date of time period}
#'   \item{N.dates}{Number of sampling dates in each time period}
#'   \item{N.f}{Number of females sampled}
#'   \item{N.f.SI}{Number of females sampled with spawning indicators}
#'   \item{P.SI}{Proportion of females with spawing indicators}
#'   \item{Age.mn}{Mean age of females}
#'   \item{Age.rng}{Age range of females}
#' }
#' @source <https://sedarweb.org/documents/sedar31-dw07-expanded-annual-stock-assessment-survey-2011-red-snapper-reproduction/>
"data_SI_d"

#' Seasonal patterns in red snapper spawning activity (raw version)
#'
#' Data converted from `data_SI_d` so that each row represents data for one fish and SI is binary
#'
#' @format ## `data_SI_d_bin`
#' \describe{
#'   \item{Start.date}{Start date of time period}
#'   \item{SI}{Presence (1) or absence (0) of spawing indicators}
#' }
#' @source <https://sedarweb.org/documents/sedar31-dw07-expanded-annual-stock-assessment-survey-2011-red-snapper-reproduction/>
"data_SI_d_bin"

#' Diel patterns in red snapper spawning activity
#'
#' Data shows the proportions of female red snapper with spawning indicators by time of day, from data captured from Figure 1b of: Jackson, M. W., Nieland, D. L., & Cowan Jr, J. H. (2006). Diel spawning periodicity of red snapper Lutjanus campechanus in the northern Gulf of Mexico. Journal of Fish Biology, 68(3), 695-706.
#'
#' @format ## `data_SI_h`
#' \describe{
#'   \item{CruiseNum}{Numeric ID of sampling cruise, in sequence.}
#'   \item{Cruise}{Alphanumeric cruise ID}
#'   \item{Hour}{Hour of sampling}
#'   \item{Fmat}{Number of mature females sampled}
#'   \item{HO}{Number of females with hydrated oocytes}
#'   \item{POF}{Number of females with post-ovulatory follicles}
#' }
#' @source <https://onlinelibrary.wiley.com/doi/abs/10.1111/j.0022-1112.2006.00951.x>
"data_SI_h"

#' Diel patterns in red snapper spawning activity (raw version)
#'
#' Data converted from `data_SI_h` so that each row represents data for one fish and SI is binary
#'
#' @format ## `data_SI_h_bin`
#' \describe{
#'   \item{CruiseNum}{Numeric ID of sampling cruise, in sequence.}
#'   \item{Cruise}{Alphanumeric cruise ID}
#'   \item{Hour}{Hour of sampling}
#'   \item{Fmat}{Number of mature females sampled}
#'   \item{HO}{Presence (1) or absence (0) of hydrated oocytes}
#'   \item{POF}{Presence (1) or absence (0) of post-ovulatory follicles}
#'   \item{Hour.shf}{Shifted hour of sampling, used for plotting distribution of POF. Hours 0-8 are shifted by +24 hours.}
#' }
#' @source <https://onlinelibrary.wiley.com/doi/abs/10.1111/j.0022-1112.2006.00951.x>
"data_SI_h_bin"
