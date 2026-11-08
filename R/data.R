#' Data for borders elimination simulation in Anderson, Larch an Yotov (2018)
#'
#' @description
#' This data is used to replicate the results of Anderson, Larch an Yotov (2018) for
#' border elimination.
#'
#' The data has been directly ported from the original Stata dataset. The original data was obtained by contacting the
#' authors diretly.
#'
#' @format A data.frame of 10 variables and 1681 observations for the year 2002.
#' \itemize{
#'   \item exporter ISO-3 codes for the exporters
#'   \item importer ISO-3 codes for the importers
#'   \item trade trade in United States dollars
#'   \item ln_dist natural logarithm of the distance between countries in kilometers (between capitols)
#'   \item cntg dummy variable for land continuity (1 means contiguous countries)
#'   \item brdr dummy variable for shared border (1 means shared border)
#'   \item output COMPLETE
#'   \item expndr COMPLETE
#'   \item expndr_deu COMPLETE
#'   \item pair exporter and importer pair (for robust standard errors)
#' }
#'
#' @references
#' Anderson, James E., Mario Larch, and Yoto V. Yotov. "GEPPML: General equilibrium analysis with PPML."
#' \emph{The World Economy} 41, no. 10 (2018): 2750-2782.
"aly2018"
