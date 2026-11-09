#' @title Price change
#' @detail OMR and IMR are power transforms of the price indexes, the paper reports the change in
#' the price indexes
#' @param new price vector from current fit
#' @param old price vector from previous fit
#' @export
price_change <- function(new, old) {
  ((new / old)^(1 / (1 - sigma)) - 1) * 100
}

rgdp <- function(output, imr) output / imr^(1 / (1 - sigma))

exports <- function(output, expenditure, t, imr, omr) {
  x <- outer(output, expenditure) * t / outer(omr, imr)
  rowSums(x) - diag(x)
}
