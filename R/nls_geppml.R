#' General equilibrium PPML (Anderson, Larch and Yotov, 2018) with a non-linear solver
#'
#' @description
#' Estimates the gravity equation with PPML, builds the baseline and
#' counterfactual (no borders) trade costs and solves the multilateral
#' resistances and factory-gate prices as a system of non-linear equations,
#' as in the original MATLAB code (ge_ppml_indexes.m). The results should match
#' those from the regression-based approach.
#'
#' @param data a data.frame with the columns exporter, importer, trade, output,
#' expndr and the regressors (e.g., \code{aly18}). Output is the sum of trade by
#' exporter and expenditure the sum of trade by importer.
#' @param formula the gravity specification, with the regressors being columns of
#' \code{data}. The default is \code{trade ~ ln_dist + cntg + brdr | exporter + importer}.
#' @param border name of the border variable that is set to zero in the
#' counterfactual (internal trade costs are kept as in the baseline).
#' @param removed optional logical vector (one value per row of \code{data}) with the
#' observations in which the border is removed in the counterfactual, e.g., to
#' remove it only in one direction. By default it is removed in all observations.
#' @param sigma the elasticity of substitution.
#' @param ref_country the reference country, its inward multilateral resistance
#' is normalised to one.
#' @param tol tolerance for the norm of the residuals of the non-linear systems.
#' @param maxiter maximum number of iterations of each system.
#' @param estibrated a boolean indicating whether the gravity error (observed over
#' predicted trade) is treated as a component of the trade costs, so the baseline
#' matches observed trade (estibrated trade costs). Zero flows have zero trade
#' costs, so they must be kept in \code{data} to have a balanced panel.
#' @param verbose a boolean indicating whether convergence information should be printed.
#'
#' @return A data.table with the changes (in percentage) in exports, real GDP,
#' inward and outward multilateral resistances, and factory-gate prices for the
#' conditional and full endowment scenarios.
#'
#' @examples
#' \dontrun{
#' nls_geppml(aly18, ref_country = "0-DEU")
#' }
#'
#' @importFrom capybara fepoisson
#' @export
nls_geppml <- function(data,
                       formula = trade ~ ln_dist + cntg + brdr | exporter + importer,
                       border = "brdr",
                       removed = TRUE,
                       sigma = 7,
                       ref_country = "0-DEU",
                       tol = 1e-12,
                       maxiter = 1400,
                       estibrated = FALSE,
                       verbose = FALSE) {
  dt <- copy(as.data.table(data))

  if (!ref_country %in% dt$exporter) {
    stop("'ref_country' not found in 'exporter'", call. = FALSE)
  }

  fit <- fepoisson(
    formula,
    data = dt,
    control = list(dev_tol = 1e-12, center_tol = 1e-12)
  )
  b <- coef(fit)
  if (!border %in% names(b)) {
    stop("'border' must be a regressor of 'formula'", call. = FALSE)
  }

  dt[, removed_ij := as.numeric(removed)]
  if (estibrated) {
    dt[, error_ij := trade / predict(fit)]
  }

  # the reference country goes last, as in the MATLAB code
  countries <- sort(unique(dt$exporter))
  countries <- c(setdiff(countries, ref_country), ref_country)
  n <- length(countries)

  dt[, exporter := factor(exporter, levels = countries)]
  dt[, importer := factor(importer, levels = countries)]
  setorderv(dt, c("exporter", "importer"))
  if (nrow(dt) != n * n) {
    stop("'data' must be a balanced exporter-importer panel", call. = FALSE)
  }

  xb <- as.vector(as.matrix(dt[, names(b), with = FALSE]) %*% b)
  internal <- dt$exporter == dt$importer
  t_bln <- exp(xb)
  t_ctrf <- ifelse(internal, t_bln, exp(xb - b[border] * dt[[border]] * dt$removed_ij))
  if (estibrated) {
    t_bln <- t_bln * dt$error_ij
    t_ctrf <- t_ctrf * dt$error_ij
  }
  t_bln <- matrix(t_bln, n, n, byrow = TRUE)
  t_ctrf <- matrix(t_ctrf, n, n, byrow = TRUE)

  q <- dt[, list(q = output[1]), by = exporter]$q
  e <- dt[, list(e = expndr[1]), by = importer]$e
  phi <- e / q

  sol <- nls_geppml_(t_bln, t_ctrf, q, phi, sigma, tol, maxiter, verbose)

  # the solver works with the inverse of the multilateral resistances
  total_bln <- sum(q)
  total_full <- sum(sol$pstar * q)
  imr_bln <- 1 / sol$imr_bln
  omr_bln <- total_bln / sol$omr_bln
  imr_cdl <- 1 / sol$imr_cdl
  omr_cdl <- total_bln / sol$omr_cdl
  imr_full <- 1 / sol$imr_full
  omr_full <- total_full / sol$omr_full

  exp_bln <- exports(q, e, t_bln, imr_bln, omr_bln)
  exp_cdl <- exports(q, e, t_ctrf, imr_cdl, omr_cdl)
  exp_full <- exports(sol$pstar * q, phi * sol$pstar * q, t_ctrf, imr_full, omr_full)

  out <- data.table(
    country = countries,
    exports_cdl = (exp_cdl / exp_bln - 1) * 100,
    exports_full = (exp_full / exp_bln - 1) * 100,
    rGDP_cdl = (rgdp(q, imr_cdl) / rgdp(q, imr_bln) - 1) * 100,
    rGDP_full = (rgdp(q * sol$pstar, imr_full) / rgdp(q, imr_bln) - 1) * 100,
    imr_cdl = price_change(imr_cdl, imr_bln),
    imr_full = price_change(imr_full, imr_bln),
    omr_cdl = price_change(omr_cdl, omr_bln),
    omr_full = price_change(omr_full, omr_bln),
    prices_full = (sol$pstar - 1) * 100
  )
  out[, country := ifelse(country == ref_country, "DEU", country)]
  setorderv(out, "country")
  attr(out, "convergence") <- setNames(sol$norm, c("baseline", "conditional", "full"))
  out
}
