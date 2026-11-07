#' @title Add estimated fixed effects
#' @details Adds the exponentiated fixed effects to the model data. Capybara returns normalized
#' fixed effects so the importer fixed effect of the reference country is zero.
#' @param data input data
#' @param fit model to recover the fixed effects from
#' @param exp_name exporters' names
#' @param imp_name importers' names
#' @export
add_fes <- function(data, fit, exp_name, imp_name) {
  fe_imp <- fit$fixed_effects$importer
  fe_exp <- fit$fixed_effects$exporter + fe_imp[ref_country]
  fe_imp <- fe_imp - fe_imp[ref_country]

  data[data.table(exporter = names(fe_exp), fe = exp(fe_exp)), on = "exporter", (exp_name) := i.fe]
  data[data.table(importer = names(fe_imp), fe = exp(fe_imp)), on = "importer", (imp_name) := i.fe]
}

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
